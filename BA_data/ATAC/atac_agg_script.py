#loading filtered rna and atac data
import anndata
from anndata import AnnData
import pandas as pd

filtered_ATdata = anndata.read_h5ad("/home/fgsasse_lrs_1/Downloads/BA/BA_data/ATAC/zf_multiome_atlas_full_ATAC_v1_subset.h5ad")
filtered_ATdata.obs

#loading the columns of the rna data that should be added to the atac data
filtered_Rcol = pd.read_csv("/home/fgsasse_lrs_1/Downloads/BA/BA_data/ATAC/rna_fil.csv")
filtered_Rcol.head()


# before adding the column to the ATAC data, we need to make sure that the order of the cells in the ATAC data matches the order of the cells in the RNA data. We can do this by sorting both datasets by their cell IDs.
filtered_Rcol["cell_id"] == filtered_ATdata.obs.index

# Add RNA annotation to ATAC obs by matching on cell_id
rna_annotation = (
    filtered_Rcol.set_index("cell_id")["zebrafish_anatomy_ontology_class"]
    .reindex(filtered_ATdata.obs.index)
 )
filtered_ATdata.obs["annotation_ML"] = rna_annotation.values
filtered_ATdata.obs.head()

# creating an aggregation function to sum counts for each cell type and normalize by total counts
def agg_fun(mat, meta, fun="sum", peaks=None):
    import pandas as pd
    import numpy as np

    meta = pd.Series(meta)

    if fun != "sum":
        raise ValueError("Only fun='sum' is supported.")

    def agg_f(mask):
        mask = np.asarray(mask)  # safe for sparse indexing
        sub = mat[mask, :]
        col_sums = np.asarray(sub.sum(axis=0)).ravel()  # 1D array
        total = col_sums.sum()
        if total == 0:
            return np.zeros_like(col_sums, dtype=float)
        return col_sums / total

    groups = meta.unique()
    agg_result = pd.DataFrame(
        [agg_f(meta == group) for group in groups],
        index=groups,
        columns= peaks if peaks is not None else None,
    )
    return agg_result

# Keep matrix sparse to avoid memory errors & extract metadata
atac_ct_matrix = filtered_ATdata.X
cell_types = filtered_ATdata.obs["annotation_ML"].astype(str)
dev_times = filtered_ATdata.obs["dev_stage"].astype(str)
peaks = filtered_ATdata.var_names

#Aggregate by cell type and developmental stages
meta_ct_time = cell_types + "__" + dev_times #Build modified meta vector: cell type + developmental time
agg_atac_ct_time = agg_fun(atac_ct_matrix, meta_ct_time, fun="sum", peaks = peaks)
agg_atac_ct_time.head()

#Comparing the number of combinations of cell types and developmental stages ("meta_ct_time") with the index (row) of the aggregated data frame to find out which groups are (not) present in the aggregated data frame
import numpy as np

meta_groups = set(meta_ct_time.unique())
agg_groups = set(agg_atac_ct_time.index)

print("Missing in ct + time dataframe:", meta_groups - agg_groups)
print("Number of combinations of cell types and developmental stages:", len(meta_groups))
print("Number of groups in the ct + time dataframe:", len(agg_groups))

# Save the aggregated data frames to a new AnnData object
agg_atac_ct_time = AnnData(X=agg_atac_ct_time.values, obs=pd.DataFrame(index=agg_atac_ct_time.index), var=pd.DataFrame(index=agg_atac_ct_time.columns))
agg_atac_ct_time.write_h5ad("/home/fgsasse_lrs_1/Downloads/BA/BA_data/Pseudobulks/ATAC/celltypes_times/agg_atac_ct_time.h5ad")

#Subsetting some columns to get the mean values grouped by cell type and developmental stage for the metadata columns "nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"
metadata_lg = filtered_ATdata.obs[["nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "annotation_ML", "dev_stage"]].copy()
()

#group by cell type and developmental stage to get the mean values for the metadata columns "nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"
metadata_lg_ct_time = metadata_lg.groupby(["annotation_ML", "dev_stage"]).mean()


