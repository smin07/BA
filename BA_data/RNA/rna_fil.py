import anndata
from anndata import AnnData

RNAdata = anndata.read_h5ad("/home/fgsasse_lrs_1/Downloads/BA/BA_data/RNA/zf_multiome_atlas_full_RNA_v1_release.h5ad")
RNAdata.obs

cells = []
with open("/home/fgsasse_lrs_1/Downloads/BA/BA_data/ATAC/cell_id.txt", "r") as f:
    for line in f:
        cells.append(line.strip())

cells


filtered_Rdata = RNAdata[cells, :]
filtered_Rdata.write("filtered_Rdata.h5ad")