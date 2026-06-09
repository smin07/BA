import ast
import numpy as np
import pandas as pd

from statsmodels.regression.linear_model import OLS
from statsmodels.tools import add_constant


def compute_peak_gene_ols(
    rna_data,
    atac_data,
    gene_peaks_df: pd.DataFrame,
    window_label: str,
) -> tuple[dict, pd.DataFrame]:
    """
    Compute peak-gene OLS associations for a given genomic window.

    Parameters
    ----------
    rna_data
        AnnData object containing pseudobulk RNA expression.
        Genes must be present in `.var_names`.

    atac_data
        AnnData object containing pseudobulk ATAC accessibility.
        Peaks must be present in `.var_names`.

    gene_peaks_df : pd.DataFrame
        DataFrame containing:
            - gene_id
            - assigned_peaks

        `assigned_peaks` may be a Python list or a string
        representation of a list.

    window_label : str
        Label describing the genomic window
        (e.g. "10kb", "50kb", "100kb").

    Returns
    -------
    ols_results : dict
        Nested dictionary:

        {
            gene_id: {
                peak_id: {
                    "coef": float,
                    "pval": float
                }
            }
        }

    ols_df : pd.DataFrame
        Tidy DataFrame with columns:

            window
            gene
            peak
            coef
            pval
    """

    ols_results = {}

    gene_peaks_indexed = gene_peaks_df.set_index("gene_id")

    for gene_id in gene_peaks_indexed.index:

        assigned_peaks_raw = gene_peaks_indexed.loc[
            gene_id,
            "assigned_peaks"
        ]

        if isinstance(assigned_peaks_raw, str):
            assigned_peaks = ast.literal_eval(assigned_peaks_raw)
        else:
            assigned_peaks = list(assigned_peaks_raw)

        assigned_peaks = [
            p for p in assigned_peaks
            if p in atac_data.var_names
        ]

        if len(assigned_peaks) == 0:
            ols_results[gene_id] = {}
            continue

        y = (
            np.asarray(rna_data[:, gene_id].X)
            .ravel()
            .astype(float)
        )

        if np.std(y) == 0 or y.shape[0] < 2:
            ols_results[gene_id] = {}
            continue

        X_raw = (
            np.asarray(atac_data[:, assigned_peaks].X)
            .astype(float)
        )

        # Remove constant peaks
        peak_stds = X_raw.std(axis=0)

        constant_mask = peak_stds == 0

        if constant_mask.any():

            keep_mask = ~constant_mask

            assigned_peaks = [
                peak
                for peak, keep in zip(assigned_peaks, keep_mask)
                if keep
            ]

            X_raw = X_raw[:, keep_mask]
            peak_stds = peak_stds[keep_mask]

            if X_raw.shape[1] == 0:
                ols_results[gene_id] = {}
                continue

        # Z-score normalize ATAC predictors
        X_raw = (
            X_raw - X_raw.mean(axis=0)
        ) / peak_stds

        try:
            result = OLS(
                y,
                add_constant(X_raw, prepend=True)
            ).fit()

            pvalues = result.pvalues[1:]
            coefs = result.params[1:]

        except Exception:
            ols_results[gene_id] = {}
            continue

        ols_results[gene_id] = {
            peak_id: {
                "coef": float(coefs[i]),
                "pval": float(pvalues[i]),
            }
            for i, peak_id in enumerate(assigned_peaks)
        }

    print(
        f"[{window_label}] OLS complete — "
        f"{len(ols_results)} genes processed"
    )

    ols_df = pd.DataFrame(
        [
            {
                "window": window_label,
                "gene": gene_id,
                "peak": peak_id,
                "coef": stats["coef"],
                "pval": stats["pval"],
            }
            for gene_id, peak_data in ols_results.items()
            if peak_data
            for peak_id, stats in peak_data.items()
        ]
    )

    print(
        f"OLS DataFrame: "
        f"{ols_df.shape[0]:,} peak–gene pairs"
    )

    return ols_results, ols_df
