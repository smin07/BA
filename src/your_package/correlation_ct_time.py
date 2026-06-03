"""Helpers extracted from the RNA_ATAC_correlation_ct_time notebook."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


DEFAULT_WINDOW_LABELS = ("10kb", "20kb", "50kb", "100kb")
DEFAULT_CATEGORY_ORDER = ("sig. negative", "non-significant", "sig. positive")


def _import_numpy():
    import numpy as np

    return np


def _import_pandas():
    import pandas as pd

    return pd


def _parse_assigned_peaks(assigned_peaks_raw: Any) -> list[str]:
    import ast

    if isinstance(assigned_peaks_raw, str):
        parsed_peaks = ast.literal_eval(assigned_peaks_raw)
    else:
        parsed_peaks = list(assigned_peaks_raw)

    return list(parsed_peaks)


def classify_peak_gene_pair(row: Mapping[str, Any], padj_threshold: float = 0.05) -> str:
    """Classify a peak-gene pair using adjusted p-value and correlation sign."""

    padj = row["padj"]
    correlation = row["correlation"]

    if padj <= padj_threshold and correlation < 0:
        return "sig. negative"
    if padj <= padj_threshold and correlation > 0:
        return "sig. positive"
    return "non-significant"


def compute_window_peak_correlations(
    pb_atac_ct_time: Any,
    pb_rna_ct_time: Any,
    gene_peaks_df: Any,
    selected_peaks: Sequence[str] | None = None,
) -> dict[str, dict[str, dict[str, float]]]:
    """Compute peak-wise Pearson correlations for one window assignment table."""

    np = _import_numpy()
    from scipy.stats import pearsonr

    gene_peak_results: dict[str, dict[str, dict[str, float]]] = {}

    gene_peaks_indexed = gene_peaks_df.set_index("gene_id")

    # prepare fast membership test for any provided selected_peaks
    selected_set = set(selected_peaks) if selected_peaks is not None else None

    for gene_id in gene_peaks_indexed.index:
        assigned_peaks_raw = gene_peaks_indexed.loc[gene_id, "assigned_peaks"]
        assigned_peaks = [
            peak_id
            for peak_id in _parse_assigned_peaks(assigned_peaks_raw)
            if peak_id in pb_atac_ct_time.var_names
        ]

        # if caller provided a subset of peaks to consider, restrict to their intersection
        if selected_set is not None:
            assigned_peaks = [p for p in assigned_peaks if p in selected_set]

        peak_corrs: dict[str, dict[str, float]] = {}
        if len(assigned_peaks) == 0:
            gene_peak_results[gene_id] = peak_corrs
            continue

        gene_expression = np.asarray(pb_rna_ct_time[:, gene_id].X).ravel()

        for peak_id in assigned_peaks:
            peak_accessibility = np.asarray(pb_atac_ct_time[:, peak_id].X).ravel()

            if np.std(gene_expression) == 0 or np.std(peak_accessibility) == 0:
                peak_corrs[peak_id] = {"correlation": np.nan, "pval": np.nan, "padj": np.nan}
                continue

            correlation, pval = pearsonr(gene_expression, peak_accessibility)
            peak_corrs[peak_id] = {
                "correlation": float(correlation),
                "pval": float(pval),
                "padj": np.nan,
            }

        gene_peak_results[gene_id] = peak_corrs

    return gene_peak_results


def adjust_window_peak_correlations(
    gene_peak_results: dict[str, dict[str, dict[str, float]]],
    method: str = "fdr_bh",
) -> dict[str, dict[str, dict[str, float]]]:
    """Apply multiple-testing correction to a nested correlation results dictionary."""

    np = _import_numpy()
    from statsmodels.stats.multitest import multipletests

    pair_index: list[tuple[str, str]] = []
    pvals_all: list[float] = []

    for gene_id, peak_corrs in gene_peak_results.items():
        for peak_id, stats in peak_corrs.items():
            pval = stats.get("pval")
            if pval is not None and not np.isnan(pval):
                pair_index.append((gene_id, peak_id))
                pvals_all.append(float(pval))

    if pvals_all:
        _, padj_all, _, _ = multipletests(pvals_all, method=method)
        for (gene_id, peak_id), padj in zip(pair_index, padj_all):
            gene_peak_results[gene_id][peak_id]["padj"] = float(padj)

    return gene_peak_results


def compute_all_window_peak_correlations(
    pb_atac_ct_time: Any,
    pb_rna_ct_time: Any,
    window_assignments: Mapping[str, Any],
    method: str = "fdr_bh",
    selected_peaks: Sequence[str] | None = None,
) -> dict[str, dict[str, dict[str, dict[str, float]]]]:
    """Compute correlations for every window assignment table."""

    all_results: dict[str, dict[str, dict[str, dict[str, float]]]] = {}

    for window_label, gene_peaks_df in window_assignments.items():
        window_results = compute_window_peak_correlations(
            pb_atac_ct_time=pb_atac_ct_time,
            pb_rna_ct_time=pb_rna_ct_time,
            gene_peaks_df=gene_peaks_df,
            selected_peaks=selected_peaks,
        )
        all_results[window_label] = adjust_window_peak_correlations(window_results, method=method)

    return all_results


def build_correlation_dataframe(
    window_assignments: Mapping[str, dict[str, dict[str, dict[str, float]]]],
    window_labels: Sequence[str] = DEFAULT_WINDOW_LABELS,
) -> Any:
    """Convert nested correlation results into a tidy DataFrame."""

    np = _import_numpy()
    pd = _import_pandas()

    all_rows: list[dict[str, Any]] = []

    for window_label in window_labels:
        cor_results = window_assignments.get(window_label, {})
        for gene_id, peaks in cor_results.items():
            for peak_id, stats in peaks.items():
                correlation = stats.get("correlation", np.nan)
                pval = stats.get("pval", np.nan)
                padj = stats.get("padj", np.nan)

                if np.isfinite(correlation) and np.isfinite(pval) and 0 < pval <= 1:
                    all_rows.append(
                        {
                            "window": window_label,
                            "gene": gene_id,
                            "peak": peak_id,
                            "correlation": correlation,
                            "pvalue": pval,
                            "neglog10_pvalue": -np.log10(pval),
                            "padj": padj,
                            "neglog10_padj": -np.log10(padj) if np.isfinite(padj) and padj > 0 else np.nan,
                        }
                    )

    cor_res_df = pd.DataFrame(all_rows)
    if not cor_res_df.empty:
        cor_res_df["window"] = pd.Categorical(cor_res_df["window"], categories=list(window_labels), ordered=True)

    return cor_res_df


def add_correlation_categories(cor_res_df: Any, padj_threshold: float = 0.05) -> Any:
    """Add the category column used throughout the notebook plots."""

    pd = _import_pandas()

    categorized = cor_res_df.copy()
    categorized["category"] = categorized.apply(
        lambda row: classify_peak_gene_pair(row, padj_threshold=padj_threshold),
        axis=1,
    )
    categorized["category"] = pd.Categorical(
        categorized["category"],
        categories=list(DEFAULT_CATEGORY_ORDER),
        ordered=True,
    )
    return categorized


def aggregate_correlation_categories(cor_res_df: Any) -> Any:
    """Count peak-gene pairs per window, gene, and category."""

    return (
        cor_res_df.dropna(subset=["correlation", "padj"])
        .groupby(["window", "gene", "category"])
        .size()
        .reset_index(name="count")
    )


def compute_proportion_summary(
    cor_res_df: Any,
    window_labels: Sequence[str] = DEFAULT_WINDOW_LABELS,
    cutoffs: Sequence[float] = (0.001, 0.01, 0.05, 0.1, 0.2),
) -> tuple[Any, Any]:
    """Return the full cumulative proportion table and the cutoff summary table."""

    np = _import_numpy()
    pd = _import_pandas()

    prop_records: list[dict[str, Any]] = []
    summary_records: list[dict[str, Any]] = []

    thresholds = np.linspace(0, 1, 500)

    for window_label in window_labels:
        subset = cor_res_df[cor_res_df["window"] == window_label]["padj"].dropna()
        n_total = len(subset)
        if n_total == 0:
            continue

        for threshold in thresholds:
            prop_records.append(
                {
                    "window": window_label,
                    "threshold": threshold,
                    "proportion": (subset <= threshold).sum() / n_total,
                    "n_total": n_total,
                }
            )

        for cutoff in cutoffs:
            n_sig = (subset <= cutoff).sum()
            summary_records.append(
                {
                    "window": window_label,
                    "threshold": cutoff,
                    "n_total": n_total,
                    "n_sig": n_sig,
                    "proportion": n_sig / n_total,
                }
            )

    return pd.DataFrame(prop_records), pd.DataFrame(summary_records)


def plot_proportion_significant(
    cor_res_df: Any,
    output_path: str | Path | None = None,
    window_labels: Sequence[str] = DEFAULT_WINDOW_LABELS,
):
    """Plot the cumulative proportion of significant peak-gene pairs by window."""

    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import seaborn as sns

    prop_df, summary_df = compute_proportion_summary(cor_res_df, window_labels=window_labels)
    palette = sns.color_palette("colorblind", n_colors=len(window_labels))
    color_map = dict(zip(window_labels, palette))

    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    sns.set_style("whitegrid")

    ax = axes[0]
    for window_label in window_labels:
        wdf = prop_df[prop_df["window"] == window_label]
        if wdf.empty:
            continue
        ax.plot(wdf["threshold"], wdf["proportion"], label=window_label, color=color_map[window_label], linewidth=2.0)

    ax.axvline(0.05, color="crimson", linestyle="--", linewidth=1.2, alpha=0.85, label="padj = 0.05")
    ax.set_xlabel("padj threshold", fontsize=11)
    ax.set_ylabel("Proportion of pairs <= threshold", fontsize=11)
    ax.set_title("Cumulative proportion of significant\npeak-gene pairs", fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(title="Window", fontsize=9, title_fontsize=9)

    ax = axes[1]
    bar_df = summary_df[summary_df["threshold"] == 0.05].copy()
    bars = ax.bar(
        bar_df["window"],
        bar_df["proportion"],
        color=[color_map[w] for w in bar_df["window"]],
        edgecolor="white",
        linewidth=0.8,
        width=0.5,
    )

    for bar, (_, row) in zip(bars, bar_df.iterrows()):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.005,
            f"n={int(row['n_sig']):,}\n/{int(row['n_total']):,}",
            ha="center",
            va="bottom",
            fontsize=8.5,
        )

    ax.axhline(0.05, color="crimson", linestyle="--", linewidth=1.2, alpha=0.85, label="proportion = 0.05")
    ax.set_xlabel("Window size", fontsize=11)
    ax.set_ylabel("Proportion of pairs with padj <= 0.05", fontsize=11)
    ax.set_title("Proportion significant at padj <= 0.05\nper window size", fontsize=12)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=1))
    ax.set_ylim(0, bar_df["proportion"].max() * 1.25 if not bar_df.empty else 1)
    ax.legend(fontsize=9)

    fig.suptitle("Peak-gene pair significance across window sizes", y=1.02, fontsize=14)
    plt.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")

    return fig, axes, prop_df, summary_df


def plot_category_boxplot(
    cor_res_df: Any,
    output_path: str | Path | None = None,
    window_labels: Sequence[str] = DEFAULT_WINDOW_LABELS,
    category_order: Sequence[str] = DEFAULT_CATEGORY_ORDER,
):
    """Recreate the category boxplot from the notebook."""

    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import seaborn as sns

    np = _import_numpy()

    category_colors = {
        "sig. negative": "#D32F2F",
        "non-significant": "#90A4AE",
        "sig. positive": "#388E3C",
    }

    fig, axes = plt.subplots(1, len(window_labels), figsize=(22, 6), sharey=True)
    if len(window_labels) == 1:
        axes = [axes]

    for idx, (ax, window_label) in enumerate(zip(axes, window_labels)):
        subset = cor_res_df[cor_res_df["window"] == window_label].dropna(subset=["correlation", "padj"])
        if subset.empty:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center")
            ax.set_title(f"Window = {window_label}", fontsize=12)
            continue

        sns.boxplot(
            data=subset,
            x="category",
            y="correlation",
            order=list(category_order),
            palette=category_colors,
            width=0.5,
            linewidth=1.2,
            flierprops=dict(marker=".", markersize=1.5, alpha=0.2, markeredgewidth=0),
            ax=ax,
        )

        for cat in category_order:
            cat_data = subset[subset["category"] == cat]
            if len(cat_data) > 1000:
                cat_data = cat_data.sample(1000, random_state=42)
            x_pos = category_order.index(cat)
            ax.scatter(
                np.random.normal(x_pos, 0.08, size=len(cat_data)),
                cat_data["correlation"],
                color=category_colors[cat],
                alpha=0.25,
                s=4,
                linewidths=0,
                zorder=2,
            )

        ax.axhline(0.0, color="black", linestyle="--", linewidth=1.0, alpha=0.6)
        ax.set_title(f"Window = {window_label}", fontsize=12)
        ax.set_xlabel("")
        ax.set_ylabel("Pearson r" if idx == 0 else "", fontsize=11)
        ax.set_xticklabels(list(category_order), fontsize=9, rotation=15, ha="right")

    legend_patches = [mpatches.Patch(color=category_colors[cat], label=cat) for cat in category_order]
    fig.legend(handles=legend_patches, loc="lower center", ncol=3, fontsize=9, frameon=True, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("Pearson r distribution by significance category across window sizes", y=1.02, fontsize=14)
    plt.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")

    return fig, axes


def plot_peak_count_distribution(
    agg_cor_df: Any,
    output_path: str | Path | None = None,
    window_labels: Sequence[str] = DEFAULT_WINDOW_LABELS,
    category_order: Sequence[str] = DEFAULT_CATEGORY_ORDER,
):
    """Plot the grouped peak-count distribution across windows and categories."""

    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import seaborn as sns

    np = _import_numpy()
    pd = _import_pandas()

    category_colors = {
        "sig. negative": "#D32F2F",
        "non-significant": "#90A4AE",
        "sig. positive": "#388E3C",
    }

    plot_df = agg_cor_df.copy()
    plot_df["window"] = pd.Categorical(plot_df["window"], categories=list(window_labels), ordered=True)

    fig, ax = plt.subplots(figsize=(14, 6))

    sns.boxplot(
        data=plot_df.dropna(subset=["count"]),
        x="window",
        y="count",
        hue="category",
        hue_order=list(category_order),
        order=list(window_labels),
        palette=category_colors,
        width=0.6,
        linewidth=1.2,
        flierprops=dict(marker=".", markersize=1.5, alpha=0.2, markeredgewidth=0),
        ax=ax,
    )

    n_groups = len(category_order)
    total_width = 0.6
    box_width = total_width / n_groups

    for w_idx, window_label in enumerate(window_labels):
        w_subset = plot_df[plot_df["window"] == window_label]
        for c_idx, category in enumerate(category_order):
            vals = w_subset[w_subset["category"] == category]["count"]
            if vals.empty:
                continue

            x_pos = w_idx - total_width / 2 + box_width / 2 + c_idx * box_width
            y_med = vals.median()
            y_mean = vals.mean()
            ax.text(
                x=x_pos,
                y=ax.get_ylim()[0] + 0.00001 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                s=f"median ={int(y_med):,}\nmean ={int(y_mean):,}",
                ha="center",
                va="bottom",
                fontsize=6,
                fontweight="bold",
                color=category_colors[category],
            )

    ax.set_xlabel("Window size", fontsize=12)
    ax.set_ylabel("Number of peaks per gene", fontsize=12)
    ax.set_title("Distribution of peak counts per gene across window sizes (grouped by correlation category)", fontsize=13)

    legend_patches = [mpatches.Patch(color=category_colors[cat], label=cat) for cat in category_order]
    ax.legend(
        handles=legend_patches,
        title="Correlation Category\n(based on padj <= 0.05)",
        fontsize=9,
        title_fontsize=9,
        frameon=True,
        loc="upper left",
    )

    plt.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")

    return fig, ax


def load_ct_time_correlation_inputs(base_dir: str | Path) -> dict[str, Any]:
    """Load the notebook inputs from the repository layout used in this project."""

    pd = _import_pandas()
    import anndata

    base_path = Path(base_dir)
    results_dir = base_path / "BA_data" / "Rmd_n_notebook"
    pseudobulk_dir = base_path / "BA_data" / "Pseudobulks"

    return {
        "pb_atac_ct_time": anndata.read_h5ad(
            pseudobulk_dir / "ATAC" / "celltypes_times" / "agg_atac_ct_time.h5ad"
        ),
        "pb_rna_ct_time": anndata.read_h5ad(
            pseudobulk_dir / "RNA" / "celltypes_times" / "agg_rna_ct_time.h5ad"
        ),
        "gene_peaks_10kb": pd.read_csv(results_dir / "gene_peak_assignments_10kb.csv"),
        "gene_peaks_20kb": pd.read_csv(results_dir / "gene_peak_assignments_20kb.csv"),
        "gene_peaks_50kb": pd.read_csv(results_dir / "gene_peak_assignments_50kb.csv"),
        "gene_peaks_100kb": pd.read_csv(results_dir / "gene_peak_assignments_100kb.csv"),
    }