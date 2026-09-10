"""Reusable statistical analysis for the 3D organoid pipeline.

This module centralises the analysis workflows that live in Tutorials 2-5 so
they can be reused from the Streamlit UI or standalone scripts.  Functions are
pure (data in -> DataFrame/model out) and deliberately free of printing so they
are easy to unit test.

Dataset layout used across the tutorials:
  * Feature columns start with ``Value of Default_`` (or ``Value of``).
  * Metadata columns describe Batch / Sample / Well / Field / Cell.
  * Condition / grouping columns vary by dataset (e.g. ``Condition`` in the
    Fig-3 osmotic-stress data, ``Well`` in the Fig-4 topology data).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from organoid_analysis.segmentation.paths import DATA_DIR

FEATURE_PREFIXES = ("Value of Default_", "Value of")


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------


def find_dataset_file(name: str, data_directory: Path = DATA_DIR) -> Path | None:
    """Return the path to ``name`` in the data directory, if present.

    Returns ``None`` (instead of raising) when the file is missing so callers
    can show a friendly message instead of crashing.
    """
    path = Path(data_directory) / name
    return path if path.exists() else None


def detect_feature_columns(df: pd.DataFrame) -> list[str]:
    """Return columns that hold cell features (start with ``Value of``)."""
    return [col for col in df.columns if col.startswith(FEATURE_PREFIXES)]


def find_grouping_column(df: pd.DataFrame) -> str | None:
    """Pick the most likely grouping column (Condition > Sample > Well)."""
    lowered = {col: col.lower() for col in df.columns}
    for keyword in ("condition", "treatment", "sample", "well"):
        for col, low in lowered.items():
            if keyword in low and df[col].nunique() >= 2:
                return col
    return None


def simplify_feature_name(feature: str) -> str:
    """Strip the ``Value of Default_`` prefix and any radius suffix."""
    name = feature
    for prefix in FEATURE_PREFIXES:
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    for suffix in ("_in_50.0_μm", "_50.0_μm", "_in_20.0_μm", "_20.0_μm"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


# ---------------------------------------------------------------------------
# Generic descriptive helpers (Tutorial 2)
# ---------------------------------------------------------------------------


def create_dataset_inventory(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Build a summary table describing every loaded sheet.

    Args:
        frames: mapping of ``"file - sheet"`` -> DataFrame.

    Returns:
        DataFrame with one row per sheet.
    """
    rows = []
    for label, df in frames.items():
        numeric = df.select_dtypes(include=[np.number]).columns
        well_col = next((c for c in df.columns if "well" in c.lower()), None)
        rows.append(
            {
                "Dataset": label,
                "Total_Cells": df.shape[0],
                "Total_Features": df.shape[1],
                "Numeric_Features": len(numeric),
                "N_Wells": df[well_col].nunique() if well_col else "N/A",
                "Missing_Data_%": round(
                    float(df.isnull().sum().sum() / max(df.size, 1) * 100), 2
                ),
            }
        )
    return pd.DataFrame(rows)


def assess_data_quality(df: pd.DataFrame) -> pd.DataFrame:
    """Return per-column missing-value counts alongside '% missing'."""
    missing = df.isnull().sum()
    missing_pct = (missing / len(df)) * 100
    quality = pd.DataFrame(
        {
            "Column": df.columns,
            "Missing_Count": missing.values,
            "Missing_Percent": np.round(missing_pct.values, 2),
            "Zero_Variance": [
                col in df.columns
                and pd.api.types.is_numeric_dtype(df[col])
                and df[col].std() == 0
                for col in df.columns
            ],
        }
    )
    return quality[quality["Missing_Count"] > 0].sort_values(
        "Missing_Percent", ascending=False
    )


def detect_outliers_iqr(
    df: pd.DataFrame, n_features_to_show: int = 10
) -> pd.DataFrame:
    """Flag outliers per numeric feature using the 1.5 x IQR rule."""
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    rows = []
    for col in numeric_cols:
        q1 = df[col].quantile(0.25)
        q3 = df[col].quantile(0.75)
        iqr = q3 - q1
        lower = q1 - 1.5 * iqr
        upper = q3 + 1.5 * iqr
        n_outliers = int(
            ((df[col] < lower) | (df[col] > upper)).sum()
        )
        rows.append(
            {
                "Feature": col,
                "N_Outliers": n_outliers,
                "Percent_Outliers": round(n_outliers / len(df) * 100, 2),
                "Lower_Bound": round(float(lower), 4),
                "Upper_Bound": round(float(upper), 4),
            }
        )
    outlier_df = pd.DataFrame(rows)
    outlier_df = outlier_df[outlier_df["N_Outliers"] > 0].sort_values(
        "Percent_Outliers", ascending=False
    )
    return outlier_df.head(n_features_to_show)


# ---------------------------------------------------------------------------
# Tutorial 3: morphological statistics (two-group comparison)
# ---------------------------------------------------------------------------


def _safe_shapiro(data: pd.Series, sample_size: int = 5000) -> float:
    """Shapiro p-value on a size-capped sample (avoids over-rejection)."""
    values = data.dropna()
    if len(values) > sample_size:
        values = values.sample(sample_size, random_state=42)
    if len(values) < 3 or not np.isfinite(values).all() or values.nunique() < 2:
        return np.nan
    return float(stats.shapiro(values).pvalue)


def check_normality(
    df: pd.DataFrame, feature_cols: list[str], condition_col: str
) -> pd.DataFrame:
    """Per-condition Shapiro-Wilk normality results."""
    rows = []
    for col in feature_cols:
        for condition in df[condition_col].unique():
            p = _safe_shapiro(df.loc[df[condition_col] == condition, col])
            rows.append(
                {
                    "Feature": simplify_feature_name(col),
                    "Condition": condition,
                    "P_value": p,
                    "Is_Normal": p > 0.05 if np.isfinite(p) else pd.NA,
                }
            )
    return pd.DataFrame(rows)


def check_equal_variance(
    df: pd.DataFrame, feature_cols: list[str], condition_col: str
) -> pd.DataFrame:
    """Levene's test for equal variance between conditions."""
    conditions = df[condition_col].unique()
    if len(conditions) < 2:
        raise ValueError("At least two conditions are required.")
    rows = []
    for col in feature_cols:
        groups = [df.loc[df[condition_col] == c, col].dropna() for c in conditions]
        if any(len(g) < 2 for g in groups):
            p_value = np.nan
        else:
            _, p_value = stats.levene(*groups)
        rows.append(
            {
                "Feature": simplify_feature_name(col),
                "P_value": float(p_value),
                "Equal_Variance": p_value > 0.05 if np.isfinite(p_value) else pd.NA,
            }
        )
    return pd.DataFrame(rows)


def cohens_d(group1: pd.Series, group2: pd.Series) -> float:
    """Pooled Cohen's d effect size between two groups."""
    g1 = group1.dropna()
    g2 = group2.dropna()
    n1, n2 = len(g1), len(g2)
    if n1 < 2 or n2 < 2:
        return np.nan
    var1, var2 = g1.var(ddof=1), g2.var(ddof=1)
    pooled = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled == 0:
        return np.nan
    return float((g1.mean() - g2.mean()) / pooled)


def effect_size_label(d: float) -> str:
    """Classify Cohen's d magnitude into a labelled band."""
    if not np.isfinite(d):
        return "Not estimable"
    magnitude = abs(d)
    if magnitude < 0.2:
        return "Negligible"
    if magnitude < 0.5:
        return "Small"
    if magnitude < 0.8:
        return "Medium"
    return "Large"


def compare_two_groups(
    df: pd.DataFrame,
    feature_cols: list[str],
    group_col: str,
    group_order: tuple[str, str] | None = None,
) -> pd.DataFrame:
    """Two-group statistical comparison with Bonferroni adjustment.

    Uses a t-test when normality + equal-variance assumptions hold, otherwise a
    Mann-Whitney U test.  Also reports Cohen's d and Bonferroni significance.
    """
    conditions = sorted(df[group_col].unique())
    if len(conditions) < 2:
        raise ValueError("At least two groups are required.")
    group1, group2 = (group_order if group_order is not None else tuple(conditions))[:2]
    if group1 == group2:
        raise ValueError("Select two distinct groups")
    if group1 not in df[group_col].values or group2 not in df[group_col].values:
        raise ValueError(f"Groups {group1!r} / {group2!r} not found in data.")

    rows = []
    for col in feature_cols:
        g1 = df.loc[df[group_col] == group1, col].dropna()
        g2 = df.loc[df[group_col] == group2, col].dropna()
        if not np.isfinite(g1).all() or not np.isfinite(g2).all():
            raise ValueError(f"{col} contains nonfinite measurements")
        p_norm1 = _safe_shapiro(g1)
        p_norm2 = _safe_shapiro(g2)
        p_var = (
            float(stats.levene(g1, g2).pvalue)
            if len(g1) >= 2 and len(g2) >= 2
            else np.nan
        )
        assumptions_met = p_norm1 > 0.05 and p_norm2 > 0.05 and p_var > 0.05
        if len(g1) < 3 or len(g2) < 3:
            # The existing selection rule requires Shapiro-Wilk, which has
            # no result below three observations. Do not fabricate p=1 or
            # silently select an alternative statistical procedure.
            p_value = np.nan
            test_used = "NOT ASSESSED"
        elif assumptions_met:
            _, p_value = stats.ttest_ind(g1, g2, equal_var=True)
            test_used = "t-test"
        else:
            _, p_value = stats.mannwhitneyu(g1, g2, alternative="two-sided")
            test_used = "Mann-Whitney U"
        d = cohens_d(g1, g2)
        rows.append(
            {
                "Feature": simplify_feature_name(col),
                f"{group1}_Mean": round(float(g1.mean()), 4),
                f"{group1}_SD": round(float(g1.std()), 4),
                f"{group1}_N": int(len(g1)),
                f"{group2}_Mean": round(float(g2.mean()), 4),
                f"{group2}_SD": round(float(g2.std()), 4),
                f"{group2}_N": int(len(g2)),
                "Mean_Difference": round(float(g2.mean() - g1.mean()), 4),
                "Test_Used": test_used,
                "P_value": float(p_value),
                "Cohens_d": round(d, 4),
                "Effect_Size": effect_size_label(d),
            }
        )

    result = pd.DataFrame(rows)
    threshold = 0.05 / len(feature_cols) if feature_cols else 1.0
    result["Significant_Bonferroni"] = (result["P_value"] < threshold).astype("boolean")
    result.loc[~np.isfinite(result["P_value"]), "Significant_Bonferroni"] = pd.NA
    return result


def field_coefficient_of_variation(
    df: pd.DataFrame, feature_cols: list[str], condition_col: str, field_col: str
) -> pd.DataFrame:
    """Coefficient of variation (%) of per-field means, per condition."""
    field_means = df.groupby([condition_col, field_col])[feature_cols].mean()
    rows = []
    for col in feature_cols:
        for condition in df[condition_col].unique():
            values = field_means.loc[condition, col]
            mean = values.mean()
            std = values.std()
            rows.append(
                {
                    "Feature": simplify_feature_name(col),
                    "Condition": condition,
                    "Field_CV_%": round((std / mean) * 100, 2) if mean != 0 else np.nan,
                    "Field_Mean": round(float(mean), 4),
                    "Field_SD": round(float(std), 4),
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Tutorial 4: machine-learning classification of spatial topology
# ---------------------------------------------------------------------------


def prepare_binary_data(
    df: pd.DataFrame, feature_cols: list[str], group_col: str, well1: str, well2: str
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Filter two groups and build a binary feature matrix + labels."""
    binary = df[df[group_col].isin([well1, well2])].copy()
    binary["Label"] = (binary[group_col] == well1).astype(int)
    label_map = {1: well1, 0: well2}
    # Preserve missing values until after splitting: a full-data mean leaks
    # holdout information into training observations.
    X = binary[feature_cols].to_numpy(dtype=float)
    y = binary["Label"].values
    return X, y, label_map


def _check_training_features(X: np.ndarray) -> None:
    """A training mean cannot be estimated for an entirely missing feature."""
    if np.isnan(X).all(axis=0).any():
        raise ValueError("Training data contain an entirely missing feature; mean imputation is undefined")


def train_binary_classifiers(
    X: np.ndarray,
    y: np.ndarray,
    label_map: dict,
):
    """Train Logistic Regression, Random Forest, and XGBoost classifiers.

    Returns a compact summary DataFrame of per-model metrics plus the fitted
    models, scaler, and test predictions.  The heavy ML backends (scikit-learn,
    xgboost) are imported lazily so the module stays importable without them.

    Model hyperparameters below (n_estimators/max_depth/learning_rate) are
    untuned engineering defaults, not cross-validated for this pipeline's
    data; see docs/PARAMETERS.md ("Exploratory statistics/ML parameters").
    Exploratory tutorial workflow output, not a validated classifier.
    """
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, roc_auc_score
    from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    _check_training_features(X_train)
    scaler = make_pipeline(SimpleImputer(keep_empty_features=True), StandardScaler())
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    ratio = (y_train == 0).sum() / max((y_train == 1).sum(), 1)
    models = {
        "Logistic Regression": LogisticRegression(
            max_iter=1000, random_state=42, class_weight="balanced"
        ),
        "Random Forest": RandomForestClassifier(
            n_estimators=100, random_state=42, class_weight="balanced", max_depth=10
        ),
    }
    import logging

    try:
        from xgboost import XGBClassifier

        models["XGBoost"] = XGBClassifier(
            n_estimators=100,
            random_state=42,
            max_depth=5,
            learning_rate=0.1,
            scale_pos_weight=ratio,
        )
    except Exception as error:  # noqa: BLE001 - e.g. missing libomp on macOS
        logging.getLogger(__name__).warning(
            "XGBoost unavailable (%s); continuing without it.", error
        )
    n_splits = min(5, int(np.bincount(y_train).min()))
    if n_splits < 2:
        raise ValueError("At least two training observations per class are required for cross-validation")
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    for train_indices, _ in cv.split(X_train, y_train):
        _check_training_features(X_train[train_indices])

    summary_rows = []
    fitted = {}
    predictions = {}
    for name, model in models.items():
        model.fit(X_train_s, y_train)
        cv_scores = cross_val_score(
            make_pipeline(SimpleImputer(keep_empty_features=True), StandardScaler(), model),
            X_train, y_train, cv=cv, scoring="accuracy", error_score="raise",
        )
        y_pred = model.predict(X_test_s)
        acc = accuracy_score(y_test, y_pred)
        y_prob = model.predict_proba(X_test_s)[:, -1]
        auc = roc_auc_score(y_test, y_prob) if len(np.unique(y)) > 1 else 0.0
        fitted[name] = model
        predictions[name] = {"y_true": y_test, "y_pred": y_pred, "y_prob": y_prob}
        summary_rows.append(
            {
                "Model": name,
                "CV_Accuracy": round(float(cv_scores.mean()), 4),
                "CV_Std": round(float(cv_scores.std()), 4),
                "Test_Accuracy": round(float(acc), 4),
                "Test_ROC_AUC": round(float(auc), 4),
            }
        )

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        best = summary.loc[summary["CV_Accuracy"].idxmax(), "Model"]
        summary.attrs["best_model"] = best
        summary.attrs["selection_metric"] = "training CV accuracy"
        summary.attrs["evaluation_scope"] = "object split; biological generalization NOT ASSESSED"
    return summary, fitted, scaler, predictions, label_map


def feature_importance(models: dict, feature_names: list[str]) -> pd.DataFrame:
    """Combine importance from Logistic Regression, Random Forest, and XGBoost.

    Importance is normalised to [0, 1] per model and averaged into a consensus
    score for cross-model ranking.
    """
    columns = {}
    labels = [simplify_feature_name(f) for f in feature_names]
    if "Logistic Regression" in models:
        columns["Logistic_Regression"] = np.abs(models["Logistic Regression"].coef_[0])
    if "Random Forest" in models:
        columns["Random_Forest"] = models["Random Forest"].feature_importances_
    if "XGBoost" in models:
        columns["XGBoost"] = models["XGBoost"].feature_importances_

    if not columns:
        raise ValueError("No trained models with feature importance found.")

    importance = pd.DataFrame(columns, index=labels)
    normalized = importance.div(importance.max(axis=0), axis=1).fillna(0)
    importance["Consensus"] = normalized.mean(axis=1)
    return importance.sort_values("Consensus", ascending=False)


def train_multiclass_classifier(
    df: pd.DataFrame, feature_cols: list[str], label_col: str
) -> dict:
    """Train a multi-class Random Forest over all groups (e.g. wells).

    `n_estimators`/`max_depth` are untuned engineering defaults; see
    docs/PARAMETERS.md ("Exploratory statistics/ML parameters").
    """
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
    from sklearn.model_selection import train_test_split
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import LabelEncoder, StandardScaler

    le = LabelEncoder()
    y = le.fit_transform(df[label_col].astype(str).values)
    X = df[feature_cols].to_numpy(dtype=float)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    _check_training_features(X_train)
    scaler = make_pipeline(SimpleImputer(keep_empty_features=True), StandardScaler())
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    model = RandomForestClassifier(
        n_estimators=200, max_depth=15, random_state=42, class_weight="balanced", n_jobs=-1
    )
    model.fit(X_train_s, y_train)
    y_pred = model.predict(X_test_s)
    accuracy = accuracy_score(y_test, y_pred)

    return {
        "model": model,
        "scaler": scaler,
        "label_encoder": le,
        "accuracy": float(accuracy),
        "confusion_matrix": confusion_matrix(y_test, y_pred),
        "class_names": list(le.classes_),
        "y_test": y_test,
        "y_pred": y_pred,
        "report": classification_report(
            y_test, y_pred, target_names=list(le.classes_), digits=3, zero_division=0
        ),
    }


# ---------------------------------------------------------------------------
# Tutorial 5: heterogeneity & unsupervised clustering
# ---------------------------------------------------------------------------


def determine_optimal_clusters(X_scaled: np.ndarray, max_k: int = 8) -> tuple[int, list, list]:
    """Choose k by maximising the silhouette score (elbow + silhouette)."""
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    k_values = list(range(2, max_k + 1))
    inertias: list[float] = []
    silhouettes: list[float] = []
    for k in k_values:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10).fit(X_scaled)
        inertias.append(float(kmeans.inertia_))
        silhouettes.append(float(silhouette_score(X_scaled, kmeans.labels_)))
    optimal_k = k_values[int(np.argmax(silhouettes))]
    return optimal_k, inertias, silhouettes


def perform_kmeans_clustering(
    df: pd.DataFrame, feature_cols: list[str], n_clusters: int
) -> tuple[pd.DataFrame, np.ndarray, object, object]:
    """Standardise features and return the labelled DataFrame, scaled X, scaler, model."""
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler

    X = df[feature_cols].fillna(df[feature_cols].mean()).values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10).fit(X_scaled)
    clustered = df.copy()
    clustered["Cluster"] = kmeans.labels_
    return clustered, X_scaled, scaler, kmeans


def silhouette_for(df: pd.DataFrame, feature_cols: list[str]) -> float:
    """Silhouette score of the stored cluster labels (best in a testable form)."""
    from sklearn.metrics import silhouette_score

    X = df[feature_cols].fillna(df[feature_cols].mean()).values
    return float(silhouette_score(X, df["Cluster"]))


def categorize_prolate_oblate(
    df: pd.DataFrame, feature_cols: list[str]
) -> tuple[pd.DataFrame, str | None, str | None]:
    """Gate Prolate/Oblate ratios into Rod / Disk / Sphere categories.

    Thresholds are the 75th/25th percentile of the *currently loaded*
    dataset's own prolate/oblate distribution (a relative, descriptive split,
    not a biologically validated shape boundary) -- the same object can
    receive a different label depending on what else is loaded. See
    docs/PARAMETERS.md ("Exploratory statistics/ML parameters").
    """
    prolate = next((c for c in feature_cols if "prolate" in c.lower()), None)
    oblate = next((c for c in feature_cols if "oblate" in c.lower()), None)
    result = df.copy()
    if prolate is None or oblate is None:
        return result, None, None

    p_high = np.percentile(df[prolate].dropna(), 75)
    p_low = np.percentile(df[prolate].dropna(), 25)
    o_high = np.percentile(df[oblate].dropna(), 75)
    o_low = np.percentile(df[oblate].dropna(), 25)

    def category(row):
        p, o = row[prolate], row[oblate]
        if pd.isna(p) or pd.isna(o):
            return "Unknown"
        if p > p_high and o < o_low:
            return "Rod"
        if p < p_low and o > o_high:
            return "Disk"
        return "Sphere"

    result["Shape_Category"] = result.apply(category, axis=1)
    return result, prolate, oblate


def cluster_characterization(
    df: pd.DataFrame, feature_cols: list[str], n_clusters: int
) -> pd.DataFrame:
    """ANOVA + eta-squared ranking of features that separate clusters."""
    X = df[feature_cols].fillna(df[feature_cols].mean()).values
    clusters = df["Cluster"].values
    rows = []
    for idx, feature in enumerate(feature_cols):
        data = X[:, idx]
        groups = [data[clusters == c] for c in range(n_clusters)]
        if all(len(g) >= 2 for g in groups):
            f_stat, p_value = stats.f_oneway(*groups)
            ss_between = sum(len(g) * (np.mean(g) - np.mean(data)) ** 2 for g in groups)
            ss_total = np.sum((data - np.mean(data)) ** 2)
            eta = ss_between / ss_total if ss_total > 0 else np.nan
        else:
            f_stat, p_value, eta = np.nan, np.nan, np.nan
        rows.append(
            {
                "Feature": simplify_feature_name(feature),
                "F_statistic": round(float(f_stat), 4),
                "p_value": float(p_value),
                "eta_squared": round(float(eta), 4),
            }
        )
    comparison = pd.DataFrame(rows).sort_values("eta_squared", ascending=False)
    return comparison


def cluster_feature_importance(df: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    """Random Forest classification importances predicting cluster membership.

    `n_estimators`/`max_depth` are untuned engineering defaults; see
    docs/PARAMETERS.md ("Exploratory statistics/ML parameters").
    """
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.model_selection import train_test_split
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    X = df[feature_cols].to_numpy(dtype=float)
    y = df["Cluster"].values
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    _check_training_features(X_train)
    scaler = make_pipeline(SimpleImputer(keep_empty_features=True), StandardScaler())
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    rf = RandomForestClassifier(n_estimators=200, max_depth=15, random_state=42, n_jobs=-1)
    rf.fit(X_train, y_train)
    importance = pd.DataFrame(
        {
            "Feature": [simplify_feature_name(f) for f in feature_cols],
            "Importance": rf.feature_importances_,
        }
    ).sort_values("Importance", ascending=False)
    importance.attrs["test_accuracy"] = float(rf.score(X_test, y_test))
    importance.attrs["train_accuracy"] = float(rf.score(X_train, y_train))
    return importance
