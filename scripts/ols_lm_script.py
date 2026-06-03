# OLS regression of gene expression on peak accessibility per gene
# Two-stage approach:
# 1) fit_regularized on all assigned peaks to rank/select features
# 2) refit standard OLS on the selected subset only
import numpy as np
import pandas as pd
import ast
from statsmodels.regression.linear_model import OLS
from statsmodels.tools import add_constant
from statsmodels.stats.multitest import multipletests

window_labels = ["10kb", "20kb", "50kb", "100kb"]

# ── per-window OLS result stores ──────────────────────────────────────────────
gene_peak_10kb_ols_results  = {}
gene_peak_20kb_ols_results  = {}
gene_peak_50kb_ols_results  = {}
gene_peak_100kb_ols_results = {}

ols_result_stores = {
    "10kb":  gene_peak_10kb_ols_results,
    "20kb":  gene_peak_20kb_ols_results,
    "50kb":  gene_peak_50kb_ols_results,
    "100kb": gene_peak_100kb_ols_results,
}

window_assignments_raw = {
    "10kb":  gene_peaks_10kb,
    "20kb":  gene_peaks_20kb,
    "50kb":  gene_peaks_50kb,
    "100kb": gene_peaks_100kb,
}

# ── 1. OLS loop per window → per gene ────────────────────────────────────────
for window_label, gene_peaks_df in window_assignments_raw.items():
    ols_results = ols_result_stores[window_label]

    for gene_id in gene_peaks_df["gene_id"]:
        # ── retrieve & parse assigned peaks ───────────────────────────────────
        assigned_peaks_raw = gene_peaks_df.loc[
            gene_peaks_df["gene_id"] == gene_id, "assigned_peaks"
        ].iloc[0]

        if isinstance(assigned_peaks_raw, str):
            assigned_peaks = ast.literal_eval(assigned_peaks_raw)
        else:
            assigned_peaks = assigned_peaks_raw

        # keep only peaks present in ATAC object
        assigned_peaks = [p for p in assigned_peaks if p in pb_atac_ct.var_names]

        peak_ols = {}
        if len(assigned_peaks) == 0:
            ols_results[gene_id] = peak_ols
            continue

        # ── response variable: gene expression (y) ────────────────────────────
        y = np.asarray(pb_rna_ct[:, gene_id].X).ravel().astype(float)

        # skip gene if expression is constant or sample size is too small
        if np.std(y) == 0 or y.shape[0] < 2:
            ols_results[gene_id] = peak_ols
            continue

        # ── predictor matrix: accessibility for all assigned peaks ────────────
        X_raw = np.asarray(pb_atac_ct[:, assigned_peaks].X).astype(float)
        n_samples = X_raw.shape[0]
        n_peaks_total = X_raw.shape[1]

        # select at most n_samples - 1 peaks, as requested
        n_select = max(1, min(n_peaks_total, n_samples - 1))

        # stage 1: regularized fit for feature ranking
        X_full = add_constant(X_raw, prepend=True)
        try:
            reg_result = OLS(y, X_full).fit_regularized(alpha=1.0, L1_wt=1.0)
            reg_params = np.asarray(reg_result.params)
        except Exception:
            # fallback: rank by simple correlation with y if regularized fit fails
            reg_params = np.zeros(n_peaks_total + 1, dtype=float)
            for idx in range(n_peaks_total):
                x_col = X_raw[:, idx]
                if np.std(x_col) == 0:
                    reg_params[idx + 1] = 0.0
                else:
                    corr = np.corrcoef(x_col, y)[0, 1]
                    reg_params[idx + 1] = 0.0 if np.isnan(corr) else corr

        # rank by the regularized coefficients themselves (excluding intercept)
        peak_scores = np.abs(reg_params[1:])
        ranked_peak_indices = np.argsort(peak_scores)[::-1]
        selected_peak_indices = ranked_peak_indices[:n_select]
        selected_peak_indices = np.sort(selected_peak_indices)
        selected_peaks = [assigned_peaks[idx] for idx in selected_peak_indices]

        # stage 2: refit standard OLS on the selected peaks only
        X_selected = X_raw[:, selected_peak_indices]
        X = add_constant(X_selected, prepend=True)
        model = OLS(y, X)
        result = model.fit()
        conf_int = result.conf_int()

        peak_ols["selection"] = {
            "n_samples": n_samples,
            "n_peaks_total": n_peaks_total,
            "n_selected": len(selected_peaks),
            "selected_peaks": selected_peaks,
            "regularized_params": reg_params.tolist(),
        }

        peak_ols["model"] = {
            "intercept": result.params[0],
            "r_squared": result.rsquared,
            "adj_r_squared": result.rsquared_adj,
            "f_statistic": result.fvalue,
            "f_pval": result.f_pvalue,
            "aic": result.aic,
            "bic": result.bic,
            "padj": np.nan,
        }

        for peak_idx, peak_id in enumerate(selected_peaks):
            param_idx = peak_idx + 1  # +1 because index 0 is intercept
            peak_ols[peak_id] = {
                "coef": result.params[param_idx],
                "std_err": result.bse[param_idx],
                "t_stat": result.tvalues[param_idx],
                "pval": result.pvalues[param_idx],
                "conf_int_lower": conf_int[param_idx, 0],
                "conf_int_upper": conf_int[param_idx, 1],
                "padj": np.nan,
            }

        ols_results[gene_id] = peak_ols

    print(f"[{window_label}] Two-stage OLS complete — {len(ols_results)} genes processed")



# ── 3. Build tidy summary DataFrame ──────────────────────────────────────────
all_ols_results = []

for window_label, ols_results in ols_result_stores.items():
    for gene_id, peak_ols in ols_results.items():
        for peak_id, stats in peak_ols.items():
            if peak_id in {"model", "selection"}:
                continue
            if not isinstance(stats, dict):
                continue

            pval = stats.get("pval", np.nan)
            padj = stats.get("padj", np.nan)

            if np.isfinite(pval) and 0 < pval <= 1:
                all_ols_results.append({
                    "window": window_label,
                    "gene": gene_id,
                    "peak": peak_id,
                    "coef": stats.get("coef", np.nan),
                    "std_err": stats.get("std_err", np.nan),
                    "t_stat": stats.get("t_stat", np.nan),
                    "conf_int_lower": stats.get("conf_int_lower", np.nan),
                    "conf_int_upper": stats.get("conf_int_upper", np.nan),
                    "pval": pval,
                    "neglog10_pval": -np.log10(pval),
                    "padj": padj,
                    "neglog10_padj": -np.log10(padj) if np.isfinite(padj) and padj > 0 else np.nan,
                })

ols_res_df = pd.DataFrame(all_ols_results)

ols_res_df["window"] = pd.Categorical(
    ols_res_df["window"],
    categories=window_labels,
    ordered=True,
)

print(f"\nOLS summary DataFrame: {ols_res_df.shape[0]:,} peak–gene pairs")
print(ols_res_df.head())