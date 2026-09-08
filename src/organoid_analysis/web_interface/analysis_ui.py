"""Streamlit UI for the statistical analysis workflows (Tutorials 2-5)."""

from __future__ import annotations

import hashlib
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st

from organoid_analysis.experimental_statistics import phenotype_exploration as analysis
from organoid_analysis.segmentation.paths import DATA_DIR, DATA_DIR_EXISTS

EXPECTED_FILES = {
    "Fig-3 osmotic stress": "Fig-3-41592_2025_2685_MOESM10_ESM.xlsx",
    "Fig-4 spatial topology": "Fig-4-41592_2025_2685_MOESM11_ESM.xlsx",
    "Fig-5 tumor heterogeneity": "Fig-5-41592_2025_2685_MOESM12_ESM.xlsx",
}


def _available_datasets() -> dict[str, Path]:
    found: dict[str, Path] = {}
    for label, filename in EXPECTED_FILES.items():
        path = analysis.find_dataset_file(filename, DATA_DIR)
        if path is not None:
            found[label] = path
    return found


def _list_uploaded_excels() -> dict[str, Path]:
    available = _available_datasets()
    if available:
        return available
    if not DATA_DIR_EXISTS:
        return {}
    found: dict[str, Path] = {}
    for path in sorted(DATA_DIR.glob("*.xlsx")) + sorted(DATA_DIR.glob("*.xls")):
        found[path.stem] = path
    return found


_MAX_EXCEL_UPLOAD_BYTES = 100 * 1024 * 1024

EXAMPLE_TABLE = pd.DataFrame({
    "Sample": ["S1", "S1", "S2", "S2"],
    "Well": ["A1", "A1", "B1", "B1"],
    "Condition": ["Control", "Control", "Treatment", "Treatment"],
    "Value of Default_Volume": [120.4, 131.2, 98.7, 250.6],
    "Value of Default_Sphericity": [0.91, 0.88, 0.85, 0.72],
})


def _save_uploaded_excels(uploads) -> dict[str, Path]:
    """Persist ``st.file_uploader`` results to disk so the rest of this module
    can treat them exactly like files found in ``data/`` (Excel readers need a
    real path/seekable file, not the upload's in-memory buffer alone)."""
    saved: dict[str, Path] = {}
    if not uploads:
        return saved
    upload_dir = _session_upload_dir()
    for upload in uploads:
        if upload.size > _MAX_EXCEL_UPLOAD_BYTES:
            st.warning(f"Skipped {Path(upload.name).name}: file exceeds the 100 MB upload limit.")
            continue
        data = upload.getbuffer()
        digest = hashlib.sha256(data).hexdigest()
        safe_name = Path(upload.name).name
        suffix = Path(safe_name).suffix.lower()
        path = upload_dir / f"{digest}{suffix}"
        if not path.exists():
            temporary = upload_dir / f".{digest}.tmp"
            temporary.write_bytes(data)
            temporary.replace(path)
        saved[Path(safe_name).stem] = path
    return saved


def _session_upload_dir() -> Path:
    """Return a private server-generated directory for this Streamlit session."""
    key = "stats_upload_directory"
    if key not in st.session_state:
        st.session_state[key] = tempfile.mkdtemp(prefix="organoid-stats-")
    path = Path(st.session_state[key])
    path.mkdir(mode=0o700, parents=True, exist_ok=True)
    return path


def _object_features_dataset() -> dict[str, Path]:
    """Expose the current segmentation's per-object feature table (populated in
    the 'Segmentation results' tab) as one more selectable dataset here, so
    Tutorials 2-5 can run on it without an external Excel file. Its numeric
    columns are renamed with the ``Value of`` prefix ``detect_feature_columns``
    looks for; it carries no Condition/Well metadata (a single segmentation run
    has no such grouping), so group-comparison tools will ask the user to pick
    a column manually, same as they do for any dataset without an obvious one.
    """
    features = st.session_state.get("features")
    if features is None or features.empty:
        return {}
    renamed = features.reset_index().rename(
        columns=lambda c: c if c == "Label" else f"Value of {c}"
    )
    upload_dir = _session_upload_dir()
    path = upload_dir / "current_segmentation_object_features.xlsx"
    export_key = st.session_state.get("feature_key", id(features))
    if st.session_state.get("stats_feature_export_key") != export_key or not path.exists():
        temporary = upload_dir / ".current_segmentation_object_features.tmp.xlsx"
        renamed.to_excel(temporary, index=False, sheet_name="Sheet1")
        temporary.replace(path)
        st.session_state["stats_feature_export_key"] = export_key
    return {"Current segmentation - Object features": path}


def _resolve_dataframe(path: Path, sheet: str | None) -> pd.DataFrame:
    # ``sheet_name=None`` reads every sheet into a dict, not a single
    # DataFrame, so a named sheet must be probed first before falling back to
    # scanning every sheet name for one with feature columns.
    if sheet is not None:
        if analysis.detect_feature_columns(pd.read_excel(path, sheet_name=sheet, nrows=0)):
            return pd.read_excel(path, sheet_name=sheet)
    for candidate in pd.ExcelFile(path).sheet_names:
        probe = pd.read_excel(path, sheet_name=candidate, nrows=0)
        if analysis.detect_feature_columns(probe):
            return pd.read_excel(path, sheet_name=candidate)
    return pd.read_excel(path, sheet_name=sheet if sheet is not None else 0)


def render_analysis() -> None:
    st.title("Statistical Analysis")
    st.caption("Tutorials 2-5: exploration, morphology stats, ML classification, clustering.")
    with st.expander("Before running an analysis", expanded=False):
        st.markdown(
            "These tools analyze the selected Excel rows and are intended for exploration of the tutorial "
            "datasets. They do not infer the experimental unit from filenames. Record which rows belong to "
            "the same image, well, biological replicate, donor and acquisition batch before interpreting a "
            "result. Cell- or organoid-level rows are correlated within a field, so a large row count does "
            "not by itself provide a large biological sample size.\n\n"
            "**Explore** reports missing values and IQR outliers; an outlier is a review cue, not an automatic "
            "exclusion rule. **Morphology statistics** compares two selected groups and adjusts across features, "
            "but does not model nested wells or batches. **ML classification** measures separation of the selected "
            "labels, not causal biology or clinical performance; split data by independent source when possible. "
            "**Clustering** finds mathematical groups whose biological identity must be confirmed independently."
        )

    available = _list_uploaded_excels()

    uploads = st.file_uploader(
        "Upload Excel dataset(s)", type=["xlsx", "xls"], accept_multiple_files=True,
        help="Feature columns must start with `Value of Default_` (or `Value of`); "
        "see the format example below.",
    )
    available.update(_save_uploaded_excels(uploads))
    available.update(_object_features_dataset())

    with st.expander("Expected Excel format", expanded=False):
        st.caption(
            "Feature columns start with `Value of Default_` (or `Value of`). Metadata "
            "columns describe Batch/Sample/Well/Field/Cell. A Condition/Well-style column "
            "(2-4 distinct values) provides the grouping used for comparisons, e.g.:"
        )
        st.dataframe(EXAMPLE_TABLE, use_container_width=True)

    if not available:
        st.info(
            "No dataset available yet. Upload an Excel file above, place one in a "
            "`data/` directory and restart, or run a segmentation (its Object features "
            "table is picked up here automatically once computed)."
        )
        return

    with st.expander("Detected datasets", expanded=False):
        st.write({label: str(path) for label, path in available.items()})

    tab_explore, tab_morph, tab_ml, tab_cluster = st.tabs(
        ["Explore", "Morphology Stats", "ML Classification", "Clustering"]
    )

    with tab_explore:
        _render_explore(available)
    with tab_morph:
        _render_morphology(available)
    with tab_ml:
        _render_ml(available)
    with tab_cluster:
        _render_clustering(available)


# ---------------------------------------------------------------------------
# Tutorial 2: Exploration
# ---------------------------------------------------------------------------


def _render_explore(available: dict[str, Path]) -> None:
    st.subheader("Dataset exploration (Tutorial 2)")
    st.caption("Inspect completeness and distribution problems before filtering, testing or modeling.")
    selected = st.multiselect(
        "Files to include", list(available), default=list(available)[:1]
    )
    if not selected:
        st.warning("Select at least one dataset.")
        return

    frames: dict[str, pd.DataFrame] = {}
    for label in selected:
        path = available[label]
        for sheet in pd.ExcelFile(path).sheet_names:
            df = pd.read_excel(path, sheet_name=sheet)
            if df.shape[0] == 0:
                continue
            frames[f"{path.stem} - {sheet}"] = df

    st.markdown("**Inventory**")
    st.dataframe(analysis.create_dataset_inventory(frames))

    dataset_choice = st.selectbox("Dataset to inspect", list(frames), key="explore_ds")
    df = frames[dataset_choice]

    if st.button("Run quality & outlier checks", type="primary"):
        quality = analysis.assess_data_quality(df)
        st.markdown("**Missing values / zero variance**")
        if quality.empty:
            st.success("No missing values detected.")
        else:
            st.dataframe(quality)

        outliers = analysis.detect_outliers_iqr(df)
        st.markdown("**Outliers (1.5 x IQR)**")
        if outliers.empty:
            st.success("No outliers detected.")
        else:
            st.dataframe(outliers)


# ---------------------------------------------------------------------------
# Tutorial 3: Morphology statistics
# ---------------------------------------------------------------------------


def _render_morphology(available: dict[str, Path]) -> None:
    st.subheader("Morphology statistics (Tutorial 3)")
    st.caption(
        "Two-group, row-level exploratory comparison. It chooses a t-test only when its displayed "
        "normality and equal-variance checks pass; otherwise it uses Mann-Whitney U. Aggregate to the "
        "independent experimental unit before confirmatory inference."
    )
    options = {k: v for k, v in available.items()}
    label = st.selectbox("Dataset", list(options), key="morph_ds")
    path = options[label]
    try:
        df = _resolve_dataframe(path, None)
    except Exception as error:  # noqa: BLE001
        st.error(f"Could not load dataset: {error}")
        return

    feature_cols = analysis.detect_feature_columns(df)
    if not feature_cols:
        st.info("No feature columns ('Value of ...') found in this dataset.")
        return

    group_col = analysis.find_grouping_column(df)
    group_options = [c for c in df.columns if df[c].nunique() in (2, 3, 4)][:10]
    group_col = st.selectbox(
        "Grouping column (condition)",
        group_options if group_options else list(df.columns),
        index=group_options.index(group_col) if group_col in group_options else 0,
        key="morph_group",
    )
    groups = sorted(df[group_col].unique())
    g1 = st.selectbox("Group 1", groups, key="morph_g1", index=0)
    g2 = st.selectbox("Group 2", groups, key="morph_g2", index=min(1, len(groups) - 1))

    if st.button("Run statistical comparison", type="primary"):
        try:
            result = analysis.compare_two_groups(df, feature_cols, group_col, (g1, g2))
        except ValueError as error:
            st.error(str(error))
            return
        st.dataframe(result)
        _barplot_cohens_d(result)
        st.caption("Bonferroni threshold applied across all features.")


def _barplot_cohens_d(result: pd.DataFrame) -> None:
    plot = result.sort_values("Cohens_d").copy()
    colors = [
        "#e74c3c" if s else "#95a5a6" for s in plot["Significant_Bonferroni"]
    ]
    fig, ax = plt.subplots(figsize=(9, max(4, 0.5 * len(plot) + 2)))
    ax.barh(plot["Feature"], plot["Cohens_d"], color=colors, edgecolor="black")
    ax.axvline(0, color="black")
    for threshold in (0.2, 0.5, 0.8):
        ax.axvline(threshold, color="gray", linestyle="--", alpha=0.5)
        ax.axvline(-threshold, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Cohen's d (effect size)")
    ax.set_title("Effect sizes (red = significant after Bonferroni)")
    ax.grid(True, alpha=0.3, axis="x")
    st.pyplot(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Tutorial 4: ML classification
# ---------------------------------------------------------------------------


def _render_ml(available: dict[str, Path]) -> None:
    st.subheader("Machine-learning classification (Tutorial 4)")
    st.caption(
        "Models predict the selected grouping label from the selected features. Performance can be inflated "
        "when rows from one image, well or donor occur in both train and test sets; use source-grouped "
        "validation for research claims and an independent cohort for generalization."
    )
    label = st.selectbox("Dataset", list(available), key="ml_ds")
    path = available[label]
    try:
        df = _resolve_dataframe(path, None)
    except Exception as error:  # noqa: BLE001
        st.error(f"Could not load dataset: {error}")
        return
    feature_cols = analysis.detect_feature_columns(df)
    if not feature_cols:
        st.info("No feature columns found.")
        return

    group_col = st.selectbox(
        "Grouping column (e.g. Well)",
        [c for c in df.columns],
        key="ml_group",
    )
    groups = sorted(df[group_col].unique())

    mode = st.radio("Analysis mode", ["Binary comparison", "Multi-class"], key="ml_mode")
    if mode == "Binary comparison":
        if len(groups) < 2:
            st.warning("Need at least two groups.")
            return
        g1 = st.selectbox("Group 1", groups, key="ml_g1")
        g2 = st.selectbox("Group 2", groups, key="ml_g2", index=1)
        if st.button("Train binary classifiers", type="primary"):
            X, y, label_map = analysis.prepare_binary_data(
                df, feature_cols, group_col, g1, g2
            )
            summary, fitted, _, _, _ = analysis.train_binary_classifiers(X, y, label_map)
            st.success(f"Best model: {summary.attrs.get('best_model')}")
            st.dataframe(summary)
            importance = analysis.feature_importance(fitted, feature_cols)
            st.markdown("**Feature importance (consensus across models)**")
            st.dataframe(importance.head(15))
    else:
        if st.button("Train multi-class classifier", type="primary"):
            results = analysis.train_multiclass_classifier(df, feature_cols, group_col)
            st.metric("Accuracy", f"{results['accuracy']:.3f}")
            st.text(results["report"])
            _confusion_heatmap(results)


def _confusion_heatmap(results: dict) -> None:
    cm = np.asarray(results["confusion_matrix"])
    class_names = results["class_names"]
    fig, ax = plt.subplots(figsize=(max(6, len(class_names) * 0.8), max(5, len(class_names) * 0.6)))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=class_names,
        yticklabels=class_names,
        ax=ax,
    )
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title("Confusion matrix (multi-class)")
    st.pyplot(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Tutorial 5: Clustering
# ---------------------------------------------------------------------------


def _render_clustering(available: dict[str, Path]) -> None:
    st.subheader("Heterogeneity & clustering (Tutorial 5)")
    st.caption(
        "K-means always returns clusters. The silhouette score describes separation in the chosen feature "
        "space, not biological validity; compare clusters with imaging QC, known markers and independent replicates."
    )
    label = st.selectbox("Dataset", list(available), key="clus_ds")
    path = available[label]
    try:
        df = _resolve_dataframe(path, "5b (grey)") \
            if path.stem.startswith("Fig-5") else _resolve_dataframe(path, None)
    except Exception as error:  # noqa: BLE001
        st.error(f"Could not load dataset: {error}")
        return
    feature_cols = analysis.detect_feature_columns(df)
    if not feature_cols:
        st.info("No feature columns found.")
        return

    k_choice = st.selectbox("Number of clusters (k)", [2, 3, 4, 5, 6, 7, 8], index=2, key="clus_k")

    if st.button("Run clustering & characterization", type="primary"):
        optimal_k, inertias, silhouettes = analysis.determine_optimal_clusters(
            df[feature_cols].fillna(df[feature_cols].mean()).values, max_k=8
        )
        st.markdown(f"**Recommended k (highest silhouette): {optimal_k}**")

        clustered, X_scaled, _, _ = analysis.perform_kmeans_clustering(
            df, feature_cols, int(k_choice)
        )
        sil = analysis.silhouette_for(clustered, feature_cols)
        st.metric("Silhouette score", f"{sil:.3f}")
        st.write(clustered["Cluster"].value_counts().sort_index().rename("Cells per cluster"))

        characterization = analysis.cluster_characterization(clustered, feature_cols, int(k_choice))
        st.markdown("**Features separating clusters (ANOVA, eta-squared)**")
        st.dataframe(characterization.head(15))

        shaped, _, _ = analysis.categorize_prolate_oblate(clustered, feature_cols)
        if "Shape_Category" in shaped.columns:
            st.markdown("**Prolate/oblate shape categories**")
            st.write(shaped["Shape_Category"].value_counts())
