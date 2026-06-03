
import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
from datetime import datetime

# Set output directory
output_dir = "../../SEACells/outputs/RNA"
os.makedirs(output_dir, exist_ok=True)

# Load data
data_path = "../../BA_data/RNA/zf_multiome_atlas_full_RNA_v1_release.h5ad"
filtered_rna = sc.read(data_path)
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Loaded data: {filtered_rna.shape}")

gene_peaks_10kb = pd.read_csv("gene_peak_assignments_10kb.csv")
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Gene-peak assignments: {gene_peaks_10kb.shape}")

# Filter genes
genes_in_peaks = gene_peaks_10kb['gene_id'].unique()
filtered_rna = filtered_rna[:, filtered_rna.var_names.isin(genes_in_peaks)].copy()
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] After filtering, {filtered_rna.n_vars} genes remain for SEACells fitting")

## Subset to 10,000 cells (random sample)
#if filtered_rna.n_obs > 10000:
#    np.random.seed(42)
#    selected_cells = np.random.choice(filtered_rna.obs_names, 10000, replace=False)
#    filtered_rna = filtered_rna[selected_cells, :].copy()
#    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Subsetted to 10,000 cells: {filtered_rna.shape}")

# Copy counts to .raw
raw_filtered_rna = sc.AnnData(filtered_rna.layers['counts'].copy())
raw_filtered_rna.obs_names, raw_filtered_rna.var_names = filtered_rna.obs_names, filtered_rna.var_names
filtered_rna.X = raw_filtered_rna.X # to make sure any function using X by default works on raw counts
filtered_rna.raw = raw_filtered_rna # SEACells will use .raw for fitting, so we set it here
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Copied counts to .raw for SEACells fitting")

# Preprocessing
sc.pp.normalize_per_cell(filtered_rna)
sc.pp.log1p(filtered_rna)
sc.pp.highly_variable_genes(filtered_rna, n_top_genes=2000)
sc.tl.pca(filtered_rna, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(filtered_rna, n_neighbors=15, n_pcs=20)
sc.tl.umap(filtered_rna)

# Save UMAP plot
plt.figure()
sc.pl.umap(filtered_rna, show=False)
plt.savefig(os.path.join(output_dir, "umap_rna.png"))
plt.close()

# Save elbow plot
plt.figure()
sc.pl.pca_variance_ratio(filtered_rna, log=True, show=False)
plt.savefig(os.path.join(output_dir, "elbow_plot_rna.png"))
plt.close()

# SEACells parameters
n_SEACells = filtered_rna.n_obs // 75
build_kernel_on = 'X_pca'
n_waypoint_eigs = 10

# Fit SEACells
seacell_rna = filtered_rna.copy()
model = SEACells.core.SEACells(
    seacell_rna,
    build_kernel_on=build_kernel_on,
    n_SEACells=n_SEACells,
    n_waypoint_eigs=n_waypoint_eigs,
    convergence_epsilon=1e-5
)
model.construct_kernel_matrix()
M = model.kernel_matrix
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Kernel matrix shape: {M.shape}")

# Save the kernel matrix to a file for later use
np.save(os.path.join(output_dir, "kernel_matrix.npy"), M)
fn = os.path.join(output_dir, "kernel_matrix_cells.txt")
cells = seacell_rna.obs_names.values
with open(fn, 'w') as f:
    for cell in cells:
        f.write(f"{cell}\n")
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Saved kernel matrix for later use.")

# Fit the model
model.initialize_archetypes()
model.fit(min_iter=5, max_iter=100)
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] SEACells fitting complete.")

# Save model
#print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Saving model to {output_dir}...")
#model.save_model(output_dir)

# Save the SEACell assignments and model
seacells_fn = os.path.join(output_dir, "zf_multiome_atlas_full_RNA_v1_SEACells.h5ad")
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Saving SEACell assignments to {seacells_fn}...")
seacell_rna.write(seacells_fn)

# Save initialization plots
plt.figure()
SEACells.plot.plot_initialization(seacell_rna, model, plot_basis='X_pca', show=False)
plt.savefig(os.path.join(output_dir, "init_pca_rna.png"))
plt.close()

plt.figure()
SEACells.plot.plot_initialization(seacell_rna, model, plot_basis='X_umap', show=False)
plt.savefig(os.path.join(output_dir, "init_umap_rna.png"))
plt.close()

# Save convergence plot
plt.figure()
model.plot_convergence(show=False)
plt.savefig(os.path.join(output_dir, "convergence_rna.png"))
plt.close()

# Save assignment distribution plot
plt.figure(figsize=(3,2))
sns.histplot((model.A_.T > 0.1).sum(axis=1), bins=20)
plt.title('Non-trivial (> 0.1) assignments per cell')
plt.xlabel('# Non-trivial SEACell Assignments')
plt.ylabel('# Cells')
plt.savefig(os.path.join(output_dir, "assignment_distribution_rna.png"))
plt.close()

# Save top 5 assignment strengths heatmap
plt.figure(figsize=(3,2))
b = np.partition(model.A_.T, -5)
sns.heatmap(np.sort(b[:,-5:])[:, ::-1], cmap='viridis', vmin=0)
plt.title('Strength of top 5 strongest assignments')
plt.xlabel('$n^{th}$ strongest assignment')
plt.savefig(os.path.join(output_dir, "top5_assignments_rna.png"))
plt.close()

# Save SEACell sizes plot
plt.figure()
SEACells.plot.plot_SEACell_sizes(seacell_rna, bins=50, show=False)
plt.savefig(os.path.join(output_dir, "seacell_sizes_rna.png"))
plt.close()

# Save the SEACell assignments and model
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Saving SEACell assignments to {seacells_fn}...")
seacell_rna.write(seacells_fn)

# Summarize by SEACell
SEACell_ad = SEACells.core.summarize_by_SEACell(seacell_rna, SEACells_label='SEACell', summarize_layer='counts')

# Normalize and log1p summarized data
sc.pp.normalize_total(SEACell_ad)
sc.pp.log1p(SEACell_ad)

# Save summarized data
SEACell_ad.write(os.path.join(output_dir, "SEACell_summarized_RNA.h5ad"))

print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Analysis complete. All outputs and plots are saved in: {output_dir}")
