#loading filtered rna and atac data
import anndata
from anndata import AnnData
import pandas as pd

pb_atac_ct = anndata.read_h5ad("/home/fgsasse_lrs_1/Downloads/BA/BA_data/Pseudobulks/ATAC/celltypes/agg_atac_ct.h5ad")
pb_rna_ct = anndata.read_h5ad("/home/fgsasse_lrs_1/Downloads/BA/BA_data/Pseudobulks/RNA/celltypes/agg_rna_ct.h5ad")       

#Reindex the dataframes to ensure they are in the same order
pb_atac_ct = pb_atac_ct[pb_rna_ct.obs_names, :]
pb_rna_ct = pb_rna_ct[pb_atac_ct.obs_names, :]

#Calculate correlation matrix using pearsonr: genes (x) x peaks (y)
import numpy as np
from scipy.stats import pearsonr

# Extract data matrices - transpose so each row is a gene/peak
genes_data = pb_rna_ct.X.T  # genes × samples
peaks_data = pb_atac_ct.X.T  # peaks × samples

n_genes = genes_data.shape[0]
n_peaks = peaks_data.shape[0]

# Create correlation and p-value matrices
correlation_matrix = np.zeros((n_genes, n_peaks))
pvalue_matrix = np.zeros((n_genes, n_peaks))

# Calculate pearsonr for each gene-peak pair
for i in range(n_genes):
    for j in range(n_peaks):
        correlation_matrix[i, j], pvalue_matrix[i, j] = pearsonr(genes_data[i, :], peaks_data[j, :])

# Convert to DataFrame for visualization
correlation_df = pd.DataFrame(
    correlation_matrix,
    index=pb_rna_ct.var_names,
    columns=pb_atac_ct.var_names
)

print(f"Correlation matrix shape (genes × peaks): {correlation_df.shape}")
print("\nFirst 5 genes × first 5 peaks:")
print(correlation_df.iloc[:5, :5])

#Calculate correlation matrix using alternative methods: genes (x) x peaks (y) 
from scipy.stats import rankdata
import numpy as np

genes_data = pb_rna_ct.X.T  # genes × samples
peaks_data = pb_atac_ct.X.T  # peaks × samples


def pearson_corrmat(genes_data, peaks_data):
    gene_c = genes_data - genes_data.mean(axis=1, keepdims=True)
    peak_c = peaks_data - peaks_data.mean(axis=1, keepdims=True)
    num = gene_c @ peak_c.T
    gene_norm = np.linalg.norm(gene_c, axis=1, keepdims=True)
    peak_norm = np.linalg.norm(peak_c, axis=1, keepdims=True)
    return num / (gene_norm * peak_norm.T)

def spearman_corrmat(genes_data, peaks_data):
    # Rank each row
    genes_ranked = np.apply_along_axis(rankdata, 1, genes_data).astype(float)
    peaks_ranked = np.apply_along_axis(rankdata, 1, peaks_data).astype(float)
    # Then Pearson on ranks = Spearman
    return pearson_corrmat(genes_ranked, peaks_ranked)

r_pearson  = pearson_corrmat(genes_data, peaks_data)   # (n_genes, n_peaks)
r_spearman = spearman_corrmat(genes_data, peaks_data)  # (n_genes, n_peaks)
print(f"Pearson correlation matrix shape (genes × peaks): {r_pearson.shape}")
print(f"Spearman correlation matrix shape (genes × peaks): {r_spearman.shape}")