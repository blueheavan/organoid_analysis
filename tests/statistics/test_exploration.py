from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from organoid_analysis.statistics.exploration import (
    assess_data_quality,
    categorize_prolate_oblate,
    check_equal_variance,
    check_normality,
    cluster_characterization,
    cluster_feature_importance,
    cohens_d,
    compare_two_groups,
    create_dataset_inventory,
    detect_feature_columns,
    detect_outliers_iqr,
    effect_size_label,
    field_coefficient_of_variation,
    find_grouping_column,
    perform_kmeans_clustering,
    prepare_binary_data,
    simplify_feature_name,
    train_binary_classifiers,
    train_multiclass_classifier,
)


def _make_morphology_df(n=600) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    control = pd.DataFrame(
        {
            "Condition": "Isotonic",
            "Well": "W1",
            "Field": np.repeat([1, 2, 3], n // 3),
            "Value of Default_Volume_50.0_μm": rng.normal(100, 20, n),
            "Value of Default_Roundness_50.0_μm": rng.normal(0.8, 0.1, n),
        }
    )
    stressed = pd.DataFrame(
        {
            "Condition": "Hypertonic",
            "Well": "W2",
            "Field": np.repeat([4, 5, 6], n // 3),
            "Value of Default_Volume_50.0_μm": rng.normal(60, 20, n),
            "Value of Default_Roundness_50.0_μm": rng.normal(0.5, 0.1, n),
        }
    )
    return pd.concat([control, stressed], ignore_index=True)


def _make_topology_df() -> pd.DataFrame:
    rng = np.random.default_rng(1)
    rows = []
    for well in ("A01", "A07"):
        for _ in range(120):
            base = 1.0 if well == "A01" else 2.0
            rows.append(
                {
                    "Well": well,
                    "Value of Default_n_nuclei_neighbors_50.0_μm": rng.normal(6 * base, 1),
                    "Value of Default_ratio_minor_by_major_50.0_μm": rng.normal(0.5, 0.1),
                }
            )
    return pd.DataFrame(rows)


class AnalysisHelpersTests(unittest.TestCase):
    def test_detect_feature_columns(self) -> None:
        df = pd.DataFrame(
            {
                "Well": ["A01"],
                "Value of Default_Volume_50.0_μm": [1.0],
                "other": [2.0],
            }
        )
        self.assertEqual(detect_feature_columns(df), ["Value of Default_Volume_50.0_μm"])

    def test_simplify_feature_name(self) -> None:
        self.assertEqual(
            simplify_feature_name("Value of Default_n_nuclei_neighbors_50.0_μm"),
            "n_nuclei_neighbors",
        )

    def test_find_grouping_column_prefers_condition(self) -> None:
        df = pd.DataFrame(
            {"Well": ["a", "b"], "Condition": ["x", "y"], "Value of Default_v": [1, 2]}
        )
        self.assertEqual(find_grouping_column(df), "Condition")

    def test_cohens_d_separates_clear_groups(self) -> None:
        d = cohens_d(pd.Series([10.0, 9.0, 11.0, 8.0, 12.0] * 6), pd.Series([0.0, -1.0, 1.0, -2.0, 2.0] * 6))
        self.assertGreater(abs(d), 1.0)

    def test_effect_size_label(self) -> None:
        self.assertEqual(effect_size_label(0.1), "Negligible")
        self.assertEqual(effect_size_label(0.3), "Small")
        self.assertEqual(effect_size_label(0.6), "Medium")
        self.assertEqual(effect_size_label(1.5), "Large")


class MorphologyStatsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.df = _make_morphology_df()
        self.features = detect_feature_columns(self.df)

    def test_compare_two_groups_detects_volume_effect(self) -> None:
        result = compare_two_groups(self.df, self.features, "Condition")
        self.assertIn("Cohens_d", result.columns)
        self.assertIn("Significant_Bonferroni", result.columns)
        volume = result[result["Feature"] == "Volume"].iloc[0]
        self.assertGreater(abs(volume["Cohens_d"]), 1.0)
        self.assertTrue(bool(volume["Significant_Bonferroni"]))

    def test_normality_and_variance(self) -> None:
        norm = check_normality(self.df, self.features, "Condition")
        var = check_equal_variance(self.df, self.features, "Condition")
        self.assertEqual(len(norm), len(self.features) * 2)
        self.assertEqual(len(var), len(self.features))

    def test_field_coefficient_of_variation_shape(self) -> None:
        cv = field_coefficient_of_variation(self.df, self.features, "Condition", "Field")
        self.assertTrue(set(cv["Condition"]).issubset({"Isotonic", "Hypertonic"}))
        self.assertIn("Field_CV_%", cv.columns)


class MLClassificationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.df = _make_topology_df()
        self.features = detect_feature_columns(self.df)

    def test_prepare_binary_data(self) -> None:
        X, y, label_map = prepare_binary_data(self.df, self.features, "Well", "A01", "A07")
        self.assertEqual(X.shape[1], len(self.features))
        self.assertEqual(set(np.unique(y)), {0, 1})
        self.assertEqual(label_map[1], "A01")

    def test_binary_classifiers_recover_separable_groups(self) -> None:
        X, y, label_map = prepare_binary_data(self.df, self.features, "Well", "A01", "A07")
        summary, fitted, _, _, _ = train_binary_classifiers(X, y, label_map)
        self.assertEqual(len(summary), 3)
        best_auc = summary["Test_ROC_AUC"].max()
        self.assertGreater(best_auc, 0.8)

    def test_feature_importance_present(self) -> None:
        X, y, label_map = prepare_binary_data(self.df, self.features, "Well", "A01", "A07")
        _, fitted, _, _, _ = train_binary_classifiers(X, y, label_map)
        from organoid_analysis.statistics.exploration import (
            feature_importance,
        )

        importance = feature_importance(fitted, self.features)
        self.assertEqual(len(importance), len(self.features))
        self.assertIn("Consensus", importance.columns)

    def test_multiclass(self) -> None:
        results = train_multiclass_classifier(self.df, self.features, "Well")
        self.assertGreaterEqual(results["accuracy"], 0.0)
        self.assertEqual(len(results["class_names"]), 2)


class ClusteringTests(unittest.TestCase):
    def setUp(self) -> None:
        rng = np.random.default_rng(2)
        n = 400
        self.df = pd.DataFrame(
            {
                "Value of Default_Prolate_ratio": np.r_[
                    rng.uniform(0.6, 1.0, n // 2), rng.uniform(0.0, 0.2, n // 2)
                ],
                "Value of Default_Oblate_ratio": np.r_[
                    rng.uniform(0.0, 0.2, n // 2), rng.uniform(0.6, 1.0, n // 2)
                ],
                "Value of Default_Mean_intensity": rng.normal(10, 2, n),
            }
        )
        self.features = detect_feature_columns(self.df)

    def test_prolate_oblate_categories(self) -> None:
        categorized, p_col, o_col = categorize_prolate_oblate(self.df, self.features)
        self.assertIsNotNone(p_col)
        self.assertIn("Shape_Category", categorized.columns)
        self.assertIn("Rod", set(categorized["Shape_Category"]))

    def test_kmeans_and_characterization(self) -> None:
        clustered, _, _, _ = perform_kmeans_clustering(self.df, self.features, 2)
        self.assertEqual(clustered["Cluster"].nunique(), 2)
        characterization = cluster_characterization(clustered, self.features, 2)
        self.assertTrue((characterization["eta_squared"] >= 0).all())

    def test_cluster_feature_importance(self) -> None:
        clustered, _, _, _ = perform_kmeans_clustering(self.df, self.features, 2)
        importance = cluster_feature_importance(clustered, self.features)
        self.assertEqual(len(importance), len(self.features))
        self.assertGreaterEqual(importance.attrs["test_accuracy"], 0.0)


class ExplorationTests(unittest.TestCase):
    def test_inventory(self) -> None:
        df = pd.DataFrame({"Well": ["a", "b"], "Value of Default_v": [1, np.nan]})
        inventory = create_dataset_inventory({"x": df})
        self.assertEqual(inventory["Numeric_Features"].iloc[0], 1)

    def test_outliers_flagged(self) -> None:
        df = pd.DataFrame({"A": list(range(10)) + [1000]})
        outliers = detect_outliers_iqr(df)
        self.assertEqual(len(outliers), 1)

    def test_quality_zero_variance(self) -> None:
        df = pd.DataFrame({"A": [1.0, np.nan, 1.0], "B": [1.0, np.nan, 3.0]})
        quality = assess_data_quality(df)
        self.assertIn("A", quality["Column"].tolist())


if __name__ == "__main__":
    unittest.main()
