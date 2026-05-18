
import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
import SEACells
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
from datetime import datetime

# Set output directory
output_dir = "../../SEACells/outputs/ATAC"
os.makedirs(output_dir, exist_ok=True)

# Load data
data_path = "../ATAC/zf_multiome_atlas_full_ATAC_v1_release.h5ad"
filtered_atac = sc.read(data_path)
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Loaded data: {filtered_atac.shape}")

## Subset to 5,000 cells (random sample)
#if filtered_atac.n_obs > 5000:
#    np.random.seed(42)
#    selected_cells = np.random.choice(filtered_atac.obs_names, 5000, replace=False)
#    filtered_atac = filtered_atac[selected_cells, :].copy()
#    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Subsetted to 10,000 cells: {filtered_atac.shape}")

# Copy counts to .raw
raw_filtered_atac = sc.AnnData(filtered_atac.layers['counts'].copy())
raw_filtered_atac.obs_names, raw_filtered_atac.var_names = filtered_atac.obs_names, filtered_atac.var_names
filtered_atac.X = raw_filtered_atac.X # to make sure any function using X by default works on raw counts
#filtered_atac.X = filtered_atac.X.astype(float)  # ensure float type for normalization
filtered_atac.raw = raw_filtered_atac # SEACells will use .raw for fitting, so we set it here
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Copied counts to .raw for SEACells fitting")

# Preprocessing
mu.atac.pp.tfidf(filtered_atac, scale_factor=1e4)
mu.atac.tl.lsi(filtered_atac, n_comps=50)
sc.pp.neighbors(filtered_atac, use_rep="X_lsi", n_neighbors=15)
sc.tl.umap(filtered_atac)

# Save UMAP plot
plt.figure()
sc.pl.umap(filtered_atac, show=False)
plt.savefig(os.path.join(output_dir, "umap_atac.png"))
plt.close()

# SEACells parameters
n_SEACells = filtered_atac.n_obs // 75
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] n_SEACells={n_SEACells}")
build_kernel_on = 'X_lsi'
n_waypoint_eigs = 10

# Fit SEACells
seacell_atac = filtered_atac.copy()
model = SEACells.core.SEACells(
    seacell_atac,
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
cells = seacell_atac.obs_names.values
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
seacells_fn = os.path.join(output_dir, "zf_multiome_atlas_full_ATAC_v1_SEACells.h5ad")
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Saving SEACell assignments to {seacells_fn}...")
seacell_atac.write(seacells_fn)

# Save initialization plots
plt.figure()
SEACells.plot.plot_initialization(seacell_atac, model, plot_basis='X_lsi', show=False)
plt.savefig(os.path.join(output_dir, "init_lsi_atac.png"))
plt.close()

plt.figure()
SEACells.plot.plot_initialization(seacell_atac, model, plot_basis='X_umap', show=False)
plt.savefig(os.path.join(output_dir, "init_umap_atac.png"))
plt.close()

# Save convergence plot
plt.figure()
model.plot_convergence(show=False)
plt.savefig(os.path.join(output_dir, "convergence_atac.png"))
plt.close()

# Save assignment distribution plot
plt.figure(figsize=(3,2))
sns.histplot((model.A_.T > 0.1).sum(axis=1), bins=20)
plt.title('Non-trivial (> 0.1) assignments per cell')
plt.xlabel('# Non-trivial SEACell Assignments')
plt.ylabel('# Cells')
plt.savefig(os.path.join(output_dir, "assignment_distribution_atac.png"))
plt.close()

# Save top 5 assignment strengths heatmap
plt.figure(figsize=(3,2))
b = np.partition(model.A_.T, -5)
sns.heatmap(np.sort(b[:,-5:])[:, ::-1], cmap='viridis', vmin=0)
plt.title('Strength of top 5 strongest assignments')
plt.xlabel('$n^{th}$ strongest assignment')
plt.savefig(os.path.join(output_dir, "top5_assignments_atac.png"))
plt.close()

# Save SEACell sizes plot
plt.figure()
SEACells.plot.plot_SEACell_sizes(seacell_atac, bins=50, show=False)
plt.savefig(os.path.join(output_dir, "seacell_sizes_atac.png"))
plt.close()

# Save the SEACell assignments and model
print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Saving SEACell assignments to {seacells_fn}...")
seacell_atac.write(seacells_fn)

# Summarize by SEACell
SEACell_ad = SEACells.core.summarize_by_SEACell(seacell_atac, SEACells_label='SEACell', summarize_layer='counts')

# Normalize and log1p summarized data
sc.pp.normalize_total(SEACell_ad)
sc.pp.log1p(SEACell_ad)

# Save summarized data
SEACell_ad.write(os.path.join(output_dir, "SEACell_summarized_ATAC.h5ad"))

print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Analysis complete. All outputs and plots are saved in: {output_dir}")
