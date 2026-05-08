import anndata
from anndata import AnnData

ATACdata = anndata.read_h5ad("/home/fgsasse_lrs_1/Downloads/BA/BA_data/ATAC/zf_multiome_atlas_full_ATAC_v1_release.h5ad", backed="r")
cells = []
with open("/home/fgsasse_lrs_1/Downloads/BA/BA_data/ATAC/cell_id.txt", "r") as f:
    for line in f:
        cells.append(line.strip())

cells

filtered_ATdata = ATACdata[cells, :].to_memory()
filtered_ATdata.write_h5ad("filtered_ATdata.h5ad")