scRNA-seq data preprocessing
================

## loading libraries

``` r
library(Seurat)
```

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 'SeuratObject' was built under R 4.5.2 but the current version is
    ## 4.5.3; it is recomended that you reinstall 'SeuratObject' as the ABI
    ## for R may have changed

    ## 'SeuratObject' was built with package 'Matrix' 1.7.4 but the current
    ## version is 1.7.5; it is recomended that you reinstall 'SeuratObject' as
    ## the ABI for 'Matrix' may have changed

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, t

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.2.1     ✔ readr     2.2.0
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ ggplot2   4.0.2     ✔ tibble    3.3.1
    ## ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
    ## ✔ purrr     1.2.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggplot2)
library(dplyr)
library(RColorBrewer)
```

## read in scRNA-seq data

``` r
#load RNA data with ReadMtx function from Seurat package
counts_matrix <- ReadMtx(
  "/home/fgsasse_lrs_1/Downloads/BA/BA_data/RNA/matrix.mtx",
  "/home/fgsasse_lrs_1/Downloads/BA/BA_data/RNA/cells.tsv",
  "/home/fgsasse_lrs_1/Downloads/BA/BA_data/RNA/features.tsv",
  cell.column = 1,
  feature.column = 1,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = TRUE,
  unique.features = TRUE,
  strip.suffix = FALSE
)

seu.obj <- CreateSeuratObject(counts = counts_matrix)
#saveRDS(seu.obj, file = "seurat_obj.rds")
```

# adding metadata to seurat object

``` r
#adding metadata to seurat object
metadata_old <- seu.obj@meta.data
metadata_old$cell_id <- rownames(metadata_old)

metadata_new <- read.delim("/home/fgsasse_lrs_1/Downloads/BA/BA_data/RNA/cells.tsv", header = FALSE)

colnames <- c("cell_id",'developmental_stage', 'dataset', 'zebrafish_anatomy_ontology_class',
       'zebrafish_anatomy_ontology_class_coarse', 'timepoint',
       'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt',
       'total_counts_nc', 'pct_counts_nc')
colnames(metadata_new) <- colnames
head(metadata_new)
```

    ##                cell_id developmental_stage dataset
    ## 1 AAACAGCCACCTAAGC-1_1           15somites  TDR118
    ## 2 AAACAGCCAGGGAGGA-1_1           15somites  TDR118
    ## 3 AAACAGCCATAGACCC-1_1           15somites  TDR118
    ## 4 AAACATGCAAACTCAT-1_1           15somites  TDR118
    ## 5 AAACATGCAAGGACCA-1_1           15somites  TDR118
    ## 6 AAACATGCAAGGATTA-1_1           15somites  TDR118
    ##   zebrafish_anatomy_ontology_class zebrafish_anatomy_ontology_class_coarse
    ## 1                        epidermis                               epidermis
    ## 2                       pronephros                              pronephros
    ## 3                        hindbrain                               hindbrain
    ## 4                      spinal_cord                             spinal_cord
    ## 5                    neural_optic2                            neural_optic
    ## 6               neural_floor_plate                      neural_floor_plate
    ##   timepoint n_genes_by_counts total_counts total_counts_mt pct_counts_mt
    ## 1     16hpf              2317         6522             160      2.453235
    ## 2     16hpf              2319         6100             245      4.016393
    ## 3     16hpf              3467        12581             779      6.191877
    ## 4     16hpf              2145         5642             265      4.696916
    ## 5     16hpf               838         2691             181      6.726124
    ## 6     16hpf              1703         4233             207      4.890149
    ##   total_counts_nc pct_counts_nc
    ## 1             726      11.13155
    ## 2            1051      17.22951
    ## 3            2542      20.20507
    ## 4            1075      19.05353
    ## 5             732      27.20178
    ## 6            1017      24.02551

``` r
metadata <- merge.data.frame(metadata_old, metadata_new, by = "cell_id", sort=FALSE)
rownames(metadata) <- metadata$cell_id

all(rownames(metadata) == rownames(seu.obj@meta.data)) # check if rownames of metadata match colnames of seurat object
```

    ## [1] TRUE

``` r
seu.obj@meta.data <- metadata
```

## QC and filtering

# comparing the number of cells in each cell type (zebrafish anatomy ontology class) with and without the developmental stages

``` r
#generating df of original data set
so_df <- seu.obj@meta.data

unique(so_df$timepoint)
```

    ## [1] "16hpf" "24hpf" "19hpf" "10hpf" "12hpf" "14hpf"

``` r
unique(so_df$zebrafish_anatomy_ontology_class )
```

    ##  [1] "epidermis"                    "pronephros"                  
    ##  [3] "hindbrain"                    "spinal_cord"                 
    ##  [5] "neural_optic2"                "neural_floor_plate"          
    ##  [7] "neural_crest2"                "PSM"                         
    ##  [9] "optic_cup"                    "lateral_plate_mesoderm"      
    ## [11] "midbrain_hindbrain_boundary2" "neural_telencephalon"        
    ## [13] "differentiating_neurons"      "muscle"                      
    ## [15] "fast_muscle"                  "heart_myocardium"            
    ## [17] "somites"                      "NMPs"                        
    ## [19] "epidermis2"                   "pharyngeal_arches"           
    ## [21] "floor_plate2"                 "hemangioblasts"              
    ## [23] "neural_posterior"             "floor_plate"                 
    ## [25] "tail_bud"                     "endoderm"                    
    ## [27] "midbrain_hindbrain_boundary"  "neural_crest"                
    ## [29] "neural_optic"                 "hematopoietic_vasculature"   
    ## [31] "endocrine_pancreas"           "hatching_gland"              
    ## [33] "neurons"                      "notochord"                   
    ## [35] "pronephros2"                  "enteric_neurons"             
    ## [37] "epidermis3"                   "epidermis4"                  
    ## [39] "neural"                       "primordial_germ_cells"

``` r
unique(so_df$dataset) #fish samples
```

    ## [1] "TDR118" "TDR119" "TDR124" "TDR125" "TDR126" "TDR127" "TDR128"

``` r
#numbers of cells in each cell type (= zebrafish anatomy ontology class)
cells_celltype1 <- so_df %>% group_by(zebrafish_anatomy_ontology_class) %>% summarise(n_cells = n())
cells_celltype1
```

    ## # A tibble: 40 × 2
    ##    zebrafish_anatomy_ontology_class n_cells
    ##    <chr>                              <int>
    ##  1 NMPs                                1904
    ##  2 PSM                                 4674
    ##  3 differentiating_neurons             2242
    ##  4 endocrine_pancreas                  3075
    ##  5 endoderm                            1881
    ##  6 enteric_neurons                      733
    ##  7 epidermis                           3327
    ##  8 epidermis2                          2561
    ##  9 epidermis3                          1840
    ## 10 epidermis4                          1768
    ## # ℹ 30 more rows

``` r
ct_n <- ggplot(cells_celltype1, aes(x=zebrafish_anatomy_ontology_class, y= n_cells, fill = zebrafish_anatomy_ontology_class))+
   geom_col(alpha = 0.7) +
   geom_text(aes(label = ifelse(n_cells == max(n_cells) | n_cells == min(n_cells), n_cells, "")),
             vjust = -0.4, size = 2.5, col = "darkblue", fontface = "bold") +
   theme_bw()+ 
   theme(axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "none") +
   labs(title = "Number of cells in each cell type before filtering",
                     x = "cell type",
                     y = "number of cells")
ct_n
```

![](scRNA_preprocessing_files/figure-gfm/QC%20filtering-1.png)<!-- -->

``` r
#numbers of cells of every density()#numbers of cells of every developmental stages in each cell type (n_cells#numbers of cells of every developmental stages in each cell type (= zebrafish anatomy ontology class)
cells_time1 <- so_df %>% group_by(zebrafish_anatomy_ontology_class, developmental_stage) %>% summarise(n_cells = n()) 
```

    ## `summarise()` has regrouped the output.
    ## ℹ Summaries were computed grouped by zebrafish_anatomy_ontology_class and
    ##   developmental_stage.
    ## ℹ Output is grouped by zebrafish_anatomy_ontology_class.
    ## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
    ## ℹ Use `summarise(.by = c(zebrafish_anatomy_ontology_class,
    ##   developmental_stage))` for per-operation grouping (`?dplyr::dplyr_by`)
    ##   instead.

``` r
cells_time1 
```

    ## # A tibble: 236 × 3
    ## # Groups:   zebrafish_anatomy_ontology_class [40]
    ##    zebrafish_anatomy_ontology_class developmental_stage n_cells
    ##    <chr>                            <chr>                 <int>
    ##  1 NMPs                             0somites                630
    ##  2 NMPs                             10somites               254
    ##  3 NMPs                             15somites               340
    ##  4 NMPs                             20somites                99
    ##  5 NMPs                             30somites                 6
    ##  6 NMPs                             5somites                575
    ##  7 PSM                              0somites               1316
    ##  8 PSM                              10somites               560
    ##  9 PSM                              15somites              1043
    ## 10 PSM                              20somites               368
    ## # ℹ 226 more rows

``` r
ct_t_n <- ggplot(cells_time1, aes(x=zebrafish_anatomy_ontology_class, y= n_cells, fill = developmental_stage))+
   geom_col(alpha = 0.7) +    
  geom_text(aes(label = ifelse(n_cells < 100, n_cells, "")),
             vjust = -0.4, size = 2.5, col = "darkblue", fontface = "bold") +
   facet_wrap(~ developmental_stage, ncol = 1)+ theme_bw(base_size = 6)+
   theme(axis.text.x = element_text(angle = 70, hjust = 1), axis.text = element_text(size = 8)) +
   labs(title = "Number of cells in each developmental stages in each cell type before filtering",
                     x = "cell type",
                     y = "number of cells")
ct_t_n
```

![](scRNA_preprocessing_files/figure-gfm/QC%20filtering-2.png)<!-- -->

``` r
#number of cells in each developmental stages of each cell type
c <- cells_time1 %>% group_by(zebrafish_anatomy_ontology_class, developmental_stage) %>% mutate(n_below100c = sum(n_cells < 100))
c %>% group_by(n_below100c) %>% count()
```

    ## # A tibble: 2 × 2
    ## # Groups:   n_below100c [2]
    ##   n_below100c     n
    ##         <int> <int>
    ## 1           0   196
    ## 2           1    40

Regarding the number of cells in each cell type: - only 87 cells of
“primordial germ” cells (min) - 5110 of “tail bud” cells (max)

Regarding the number of cells in each developmental stages of each cell
type: - in “15 somites stage” only 2 cell types below 100 cells (min) -
in “30 somites stage” 17 cell types below 100 cells (max)

\#Filtering the data by removing every cell type in each developmental
stages below 100 cells

``` r
#filtering out the celltype & timepoints below 100 cells
filtered <- so_df %>% group_by(zebrafish_anatomy_ontology_class, timepoint) %>% summarise(n_cells = n())%>% mutate(rm = n_cells <100)
```

    ## `summarise()` has regrouped the output.
    ## ℹ Summaries were computed grouped by zebrafish_anatomy_ontology_class and
    ##   timepoint.
    ## ℹ Output is grouped by zebrafish_anatomy_ontology_class.
    ## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
    ## ℹ Use `summarise(.by = c(zebrafish_anatomy_ontology_class, timepoint))` for
    ##   per-operation grouping (`?dplyr::dplyr_by`) instead.

``` r
filtered <- filtered %>% mutate(keep_cell = ifelse(!rm, paste(zebrafish_anatomy_ontology_class, timepoint), ""))

fil_cells <- filtered$keep_cell #cells we are keeping 

so_df$cell_name <- paste(so_df$zebrafish_anatomy_ontology_class, so_df$timepoint)

filtered_df <- so_df %>% filter(cell_name %in% fil_cells)

#comparing the filtered data with unfiltered data by plotting
cells_celltype2 <- filtered_df %>% group_by(zebrafish_anatomy_ontology_class) %>% summarise(n_cells = n())
cells_celltype2
```

    ## # A tibble: 39 × 2
    ##    zebrafish_anatomy_ontology_class n_cells
    ##    <chr>                              <int>
    ##  1 NMPs                                1799
    ##  2 PSM                                 4650
    ##  3 differentiating_neurons             2050
    ##  4 endocrine_pancreas                  3031
    ##  5 endoderm                            1853
    ##  6 enteric_neurons                      636
    ##  7 epidermis                           3327
    ##  8 epidermis2                          2560
    ##  9 epidermis3                          1840
    ## 10 epidermis4                          1763
    ## # ℹ 29 more rows

``` r
ct_n_fil <- ggplot(cells_celltype2, aes(x=zebrafish_anatomy_ontology_class, y= n_cells, fill = zebrafish_anatomy_ontology_class))+
   geom_col(alpha = 0.7) +
   geom_text(aes(label = ifelse(n_cells == max(n_cells) | n_cells == min(n_cells), n_cells, "")),
             vjust = -0.4, size = 2.5, col = "darkblue", fontface = "bold") +
   theme_bw()+ 
   theme(axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "none") +
   labs(title = "Number of cells in each cell type after filtering",
                     x = "cell type",
                     y = "number of cells")
ct_n| ct_n_fil 
```

![](scRNA_preprocessing_files/figure-gfm/plots%20of%20un-/filtered%20data-1.png)<!-- -->

``` r
#numbers of cells of every density()#numbers of cells of every developmental stages in each cell type (n_cells#numbers of cells of every developmental stages in each cell type (= zebrafish anatomy ontology class)
cells_time2 <- filtered_df %>% group_by(zebrafish_anatomy_ontology_class, developmental_stage) %>% summarise(n_cells = n()) 
```

    ## `summarise()` has regrouped the output.
    ## ℹ Summaries were computed grouped by zebrafish_anatomy_ontology_class and
    ##   developmental_stage.
    ## ℹ Output is grouped by zebrafish_anatomy_ontology_class.
    ## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
    ## ℹ Use `summarise(.by = c(zebrafish_anatomy_ontology_class,
    ##   developmental_stage))` for per-operation grouping (`?dplyr::dplyr_by`)
    ##   instead.

``` r
cells_time2 
```

    ## # A tibble: 196 × 3
    ## # Groups:   zebrafish_anatomy_ontology_class [39]
    ##    zebrafish_anatomy_ontology_class developmental_stage n_cells
    ##    <chr>                            <chr>                 <int>
    ##  1 NMPs                             0somites                630
    ##  2 NMPs                             10somites               254
    ##  3 NMPs                             15somites               340
    ##  4 NMPs                             5somites                575
    ##  5 PSM                              0somites               1316
    ##  6 PSM                              10somites               560
    ##  7 PSM                              15somites              1043
    ##  8 PSM                              20somites               368
    ##  9 PSM                              5somites               1363
    ## 10 differentiating_neurons          15somites               494
    ## # ℹ 186 more rows

``` r
ct_t_fil_n <- ggplot(cells_time2, aes(x=zebrafish_anatomy_ontology_class, y= n_cells, fill = developmental_stage))+
   geom_col(alpha = 0.7) +    
  geom_text(aes(label = ifelse(n_cells < 100, n_cells, "")),
             vjust = -0.4, size = 2.5, col = "darkblue", fontface = "bold") +
   facet_wrap(~ developmental_stage, ncol = 1)+ theme_bw(base_size = 6)+
   theme(axis.text.x = element_text(angle = 70, hjust = 1), axis.text = element_text(size = 8)) +
   labs(title = "Number of cells in each developmental stages in each cell type after filtering",
                     x = "cell type",
                     y = "number of cells")
ct_t_n | ct_t_fil_n 
```

![](scRNA_preprocessing_files/figure-gfm/plots%20of%20un-/filtered%20data-2.png)<!-- -->

``` r
#creating new Seurat object by filtered cell_ids
cell_id_fil <- filtered_df$cell_id

fil.seu.obj <- subset(seu.obj, cell_id %in% cell_id_fil)
```

\#plots of original and filtered data set regarding mt percentage and
total counts plot_mt \<- ggplot(so_df, aes(x = pct_counts_mt,
after_stat(scaled)))+ geom_density() + geom_density(aes(color =
fil_cell), alpha = 0.3)+ theme(axis.text.x = element_text(angle = 90,
hjust = 1)) + labs(title = “Distribution of mitochondrial gene
expression (%)”, x = “Percentage of mitochondrial gene expression”, y =
“Number of cells (%)”) plot_mt

plot_tc \<- ggplot(so_df, aes(x = total_counts, after_stat(scaled)))+
geom_density() + geom_density(aes(color = fil_cell)) + scale_x_log10() +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title =
“Distribution of total counts”, x = “Total counts”, y = “Number of cells
(%)”) plot_tc

plot_tc_mt \<- ggplot(so_df, aes(x = total_counts, y = pct_counts_mt))+
geom_point(aes(col = fil_cell), alpha = 0.01, size = 0.5)+
geom_smooth(method = “lm”, col = “black”, linewidth = 0.5, se = FALSE)+
scale_x_log10() + theme(axis.text.x = element_text(angle = 90, hjust =
1)) + labs(title = “Total counts vs percentage of mitochondrial gene
expression”, x = “Total counts”, y = “Percentage of mitochondrial gene
expression”) plot_tc_mt

The number of cells in each cell type decreases insanely after filtering
out as shown in the plots. Next step: See if marker genes of one cell
type can be still detected after filtering - barcode ids of cells before
and after filtering - gene marker analysis

## Seurat’s standard workflow with filtered data

\#Before continuing with the following analysis, we’ll visualize the
original data first

``` r
mycol <- c("#d9ed92", "#b5e48c", "#99d98c", "#76c893", "#52b69a", "#34a0a4", "#168aad", "#1a759f", "#1e6091", "#184e77" )
mycol2 <- c("#54478C", "#2C699A", "#048BA8", "#0DB39E", "#16DB93", "#83E377", "#B9E769", "#EFEA5A", "#F1C453", "#F29E4C")


#Visualizing the QC metrics regarding.. 
#whole dataset
VlnPlot(
  fil.seu.obj,
  features = c("nFeature_RNA", "nCount_RNA", "pct_counts_mt"),
  ncol = 3,
  alpha = 0.1,
  pt.size = 0,
  cols = mycol
) 
```

    ## Warning: Default search for "data" layer in "RNA" assay yielded no results;
    ## utilizing "counts" layer instead.

![](scRNA_preprocessing_files/figure-gfm/Plots%20regarding%20developmental%20stages%20and%20cell%20types-1.png)<!-- -->

``` r
#developmental stages
VlnPlot(
  fil.seu.obj,
  features = c("nFeature_RNA", "nCount_RNA", "pct_counts_mt"),
  ncol = 3,
  group.by = "timepoint",
  alpha = 0.1,
  pt.size = 0,
  cols = mycol
) + theme(legend.position = "right") + 
  scale_fill_discrete(name= "developmental stages") &
  theme( axis.text.x = element_blank()) 
```

    ## Warning: Default search for "data" layer in "RNA" assay yielded no results;
    ## utilizing "counts" layer instead.

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](scRNA_preprocessing_files/figure-gfm/Plots%20regarding%20developmental%20stages%20and%20cell%20types-2.png)<!-- -->

``` r
#cell types
VlnPlot(
  fil.seu.obj,
  features = c("nFeature_RNA", "nCount_RNA", "pct_counts_mt"),
  ncol = 3,
  group.by = "zebrafish_anatomy_ontology_class",
  alpha = 0.1,
  pt.size = 0
) +  
  theme(legend.position = "right") + 
  scale_fill_discrete(name= "cell type") &
  theme( axis.text.x = element_blank()) 
```

    ## Warning: Default search for "data" layer in "RNA" assay yielded no results;
    ## utilizing "counts" layer instead.

![](scRNA_preprocessing_files/figure-gfm/Plots%20regarding%20developmental%20stages%20and%20cell%20types-3.png)<!-- -->

``` r
#Scatterplot show good regression between genes (unique genes detected) and counts (total UMI counts) across timepoints
FeatureScatter(fil.seu.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.00001, group.by = "developmental_stage") + 
  geom_point(aes(color = "developmental_stage"), size = 0.1, alpha = 0.01) +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_log10() 
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](scRNA_preprocessing_files/figure-gfm/Plots%20regarding%20developmental%20stages%20and%20cell%20types-4.png)<!-- -->

``` r
  #the higher nFeature the nCounts should be, cause high nCounts represents higher sequencing depth

#Scatterplot 
FeatureScatter(fil.seu.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.00001, group.by = "zebrafish_anatomy_ontology_class") + 
  geom_point(aes(color = "zebrafish_anatomy_ontology_class"), size = 0.1, alpha = 0.01) +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_log10() 
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](scRNA_preprocessing_files/figure-gfm/Plots%20regarding%20developmental%20stages%20and%20cell%20types-5.png)<!-- -->

``` r
#ggplot(so_df, aes(x = nCount_RNA, y = nFeature_RNA)) +
  #geom_point(aes( color = developmental_stage), alpha = 0.01, size = 1) +
  #geom_smooth(method = "lm")+
  #theme(axis.text.x = element_text(angle = 70, hjust = 1))+
 # scale_color_manual(values = mycol) + theme_bw()

#ggplot(so_df, aes(x = nCount_RNA, y = nFeature_RNA, color = developmental_stage)) +
  #geom_point(alpha = 0.01, size = 0.5) +
  #facet_wrap(~ developmental_stage, nrow = 2, shrink = FALSE) +
  #geom_smooth(method = "lm", col = "black", linewidth = 0.5,  se = FALSE)+ theme_bw() +
  #theme(axis.text.x = element_text(angle = 70, hjust = 1))+
```

# Seurat workflow

1.  Normalization The standard normalization method in Seurat is a
    global-scaling normalization method (NormalizeData(), default
    LogNormalize), which transforms the feature counts for each cell by
    dividing by the total counts of the cell and multiplying by a scale
    factor (default = 10000) and then applying a natural log with a
    pseudocount of 1. This should help reduce variation from sequencing
    depth and somewhat prevent high abundance genes from dominating
    downstream analyses.

- we scaled the scale.factor = 1e6

2.  Finding Variable Features (marker genes) A subset of cells with high
    cell to cell variation can be selected using the variance
    stabilization transformation method of FindVariableFeatures(). This
    uses a mean-variance relationship to return 2,000 variable features
    to prioritize for downstream analyses.

3.  Scaling the data As a standard preprocessing step prior PCA, a
    linear transformation (‘scaling’) is applied with ScaleData(). The
    ScaleData() function:

- Shifts the expression of each gene, so that the mean expression across
  cells is 0
- Scales the expression of each gene, so that the variance across cells
  is 1 =\> This step gives equal weight in downstream analyses, so that
  highly-expressed genes do not dominate

4.  Performing linear Dimension Reduction via PCA - determining
    resolution & dimension

- maximizes variance (Global structure)

5.  Clustering based on the PCA rotations

6.  Performing non-linear dimension reduction via UMAP

- preserves local relationship

``` r
#Normalization of data: ((gene expression of each cell/total expression)*scale factor (10000)) -> log-transformation; log()
fil.seu.obj <- NormalizeData(fil.seu.obj, normalization.method = "LogNormalize", scale.factor = 1e6)
```

    ## Normalizing layer: counts

``` r
#Identification of high variable features (= genes) with "FindVariableFeatures()", which shows a high cell-to-cell variation
fil.seu.obj <- FindVariableFeatures(fil.seu.obj, selection.method = "vst")
```

    ## Finding variable features for layer counts

``` r
top10genes <- head(VariableFeatures(fil.seu.obj),10)

plot_vf <- VariableFeaturePlot(fil.seu.obj)
LabelPoints(plot = plot_vf, points = top10genes, repel = TRUE)
```

    ## When using repel, set xnudge and ynudge to 0 for optimal results

    ## Warning in scale_x_log10(): log-10 transformation introduced infinite values.

![](scRNA_preprocessing_files/figure-gfm/Seurat%20workflow-1.png)<!-- -->

``` r
#Identification of clusters (= cell population) with its differential expression between the clusters by "FindAllMarkers()"

#Scaling of the data to mitigate "irrelevant" sources of variation
#all_genes <- rownames(seu.obj.fil)
fil.seu.obj <- ScaleData(fil.seu.obj)
```

    ## Centering and scaling data matrix

``` r
#Linear dimension reduction: PCA & determining dimensionality of the data
fil.seu.obj <- RunPCA(fil.seu.obj, features = VariableFeatures(fil.seu.obj))
```

    ## PC_ 1 
    ## Positive:  hmgb2a, hspb1, si:ch211-152c2.3, apoc1, si:ch211-222l21.1, ptmab, ved, vox, szl, apoeb 
    ##     apela, cdx4, cdx1a, hmga1a, hes6, tagln2, zgc:174855, capn8, crabp2b, shisa2a 
    ##     tbx16, msgn1, zgc:174153, six7, her7, ncl, bmp7a, BX248122.1, dlx3b, chrd 
    ## Negative:  mdka, CT027638.1, bahcc1b, zgc:158463, foxp4, plxna1a, FO704772.3, agrn, arvcfb, gse1 
    ##     hdac4, pleca, fabp3, bnc2, sox6, fgfr2, rn7sk, zgc:110425, mef2d, zbtb16a 
    ##     col5a1, dacha, hist1h4l.6, col11a1a, traf4a, emp2, ttn.2, celf2, hbbe3, si:ch211-286o17.1 
    ## PC_ 2 
    ## Positive:  hbbe1.3, hbae1.3, hbbe1.2, hbbe3, hbae3, hbae1.1, hbae1.3.1, hbbe2, fabp7a, fabp3 
    ##     mdka, ncam1a, si:ch73-21g5.7, actc1b, CR318588.4, plp1a, nova2, fabp11a, hbae5, cox6b1 
    ##     tuba1a, tmsb, mylpfa, ncam1b, tuba1c, elavl3, brsk2b, myl1, tmem14ca, prdx2 
    ## Negative:  apoeb, arid3c, fn1a, NC-002333.17, hspb1, kif26ab, serpinh1b, cdh1, lama5, apoc1 
    ##     si:ch211-152c2.3, il17rd, si:dkey-261h17.1, fbn2b, nrp2b, tiam1b, tp63, cdx4, XKR4, itga5 
    ##     prickle1b, pls3, rarga, phc2a, fermt1, alcama, wnt5b, fn1b, tbx16, tfap2c 
    ## PC_ 3 
    ## Positive:  musk, cdh15, rc3h2, txlnba, ttn.1, ryr1b, si:ch1073-268j14.1, rbfox1l, neb, flncb 
    ##     chrng, nexn, ttn.2, smyd1b, si:dkeyp-69b9.3, akap6, obscnb, mybphb, ldb3a, rbm24a 
    ##     slc6a8, BX323458.2, zgc:92429, efemp2a, thrap3a, tns2b, tspan9a, tmem38a, mymk, atp1a2a 
    ## Negative:  zbtb16b, pcdh19, cntfr, ncam1a, gli3, fgfr2, zbtb16a, negr1, NC-002333.17, bahcc1b 
    ##     BX571942.1, foxp4, nova2, cadm1a, rfx4, pvrl2l, zfhx4, znrf3, sox13, col18a1a 
    ##     dacha, jag2b, ptprnb, hs6st1a, ptprn2, epha4b, chl1a, ephb2b, tmem108, efna5a 
    ## PC_ 4 
    ## Positive:  rfx4, ncam1a, lrig1, tenm4, robo1, bmpr1ba, NC-002333.17, enah, hs3st3b1b, cadm1a 
    ##     adgrv1, ptprnb, jag2b, negr1, nova2, notch1b, plekhh1, crb1, epha4b, sox13 
    ##     epha4a, nav2a, ncam1b, pkd1b, BX649294.1, FRMD5, sox5, BX571942.1, efnb3b, sema5ba 
    ## Negative:  adgra2, esama, hbae1.3, hbbe1.3, egfl7, hbbe1.2, cdh5, tie1, flt4, yrk 
    ##     hbae3, klhl4, hbae1.1, fli1a, hbae1.3.1, rasip1, fgd5a, fli1b, hbbe3, hbbe2 
    ##     runx1t1, arhgef9b, clec14a, zfpm1, tal1, mrc1a, si:ch211-248e11.2, tfr1a, si:ch211-227m13.1, stab2 
    ## PC_ 5 
    ## Positive:  lama5, pcdh19, rpe65b, tp63, fgfr2, col7a1, fermt1, tfap2c, ppap2d, itga3b 
    ##     znf385a, cldni, camk2d1, col14a1a, pcdh7b, palm3, dock5, bcam, map3k15, col18a1a 
    ##     slitrk6, plekha6, hs6st1a, actn1, arhgef25b, dag1, si:ch211-264f5.6, ecrg4b, frem3, myh9a 
    ## Negative:  XKR4, si:dkey-178k16.1, hoxc3a, scn8aa, stx1b, tbx16, draxin, samd10b, cadpsb, fn1b 
    ##     her1, tbx16l, elavl4, hoxc6b, atp1a3a, CACNA2D1, adam22, BX649294.1, stxbp1a, phc2a 
    ##     mllt3, mapk10, pcdh8, kif26ab, aplp1, rbm38, enox1, wnt5b, dlc, snap25a

``` r
ElbowPlot(fil.seu.obj)
```

![](scRNA_preprocessing_files/figure-gfm/Seurat%20workflow-2.png)<!-- -->

``` r
#Clustering 
fil.seu.obj <- FindNeighbors(fil.seu.obj, dims = 1:18)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
fil.seu.obj <- FindClusters(fil.seu.obj, resolution = c(0.25, 0.5, 0.75, 1))
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 93063
    ## Number of edges: 2984999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9446
    ## Number of communities: 16
    ## Elapsed time: 22 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 93063
    ## Number of edges: 2984999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9260
    ## Number of communities: 24
    ## Elapsed time: 21 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 93063
    ## Number of edges: 2984999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9136
    ## Number of communities: 30
    ## Elapsed time: 20 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 93063
    ## Number of edges: 2984999
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9039
    ## Number of communities: 35
    ## Elapsed time: 20 seconds

``` r
DimPlot(fil.seu.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
```

![](scRNA_preprocessing_files/figure-gfm/Seurat%20workflow-3.png)<!-- -->

``` r
#Non-linear dimension reduction: UMAP 
fil.seu.obj <- RunUMAP(fil.seu.obj, dims = 1:18)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 11:14:25 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 11:14:25 Read 93063 rows and found 18 numeric columns

    ## 11:14:25 Using Annoy for neighbor search, n_neighbors = 30

    ## 11:14:25 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 11:14:29 Writing NN index file to temp file /tmp/RtmpcxyIAi/file697039c28b7f
    ## 11:14:29 Searching Annoy index using 1 thread, search_k = 3000
    ## 11:14:49 Annoy recall = 100%
    ## 11:14:49 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 11:14:50 Initializing from normalized Laplacian + noise (using RSpectra)
    ## 11:14:53 Commencing optimization for 200 epochs, with 4144500 positive edges
    ## 11:14:53 Using rng type: pcg
    ## 11:15:16 Optimization finished

``` r
#saveRDS(fil.seu.obj, file = "/home/fgsasse_lrs_1/Downloads/BA/so_preprocessed")
#markers <- FindAllMarkers(seu.obj.fil) #ups dauert laenger als gedacht..
ct_plot <- DimPlot(fil.seu.obj, reduction = "umap", group.by = "zebrafish_anatomy_ontology_class", label = TRUE)
timepoint_plot <- DimPlot(fil.seu.obj, reduction = "umap", group.by = "timepoint", label = TRUE)

ct_plot|timepoint_plot
```

![](scRNA_preprocessing_files/figure-gfm/Seurat%20workflow-4.png)<!-- -->

## Pseudobulking scRNA-Seq data

After looking in to the single cell RNA data, we create a pseudobulk RNA
count matrix with raw counts grouped by cell type & developmental
stages. DE analysis can be performed by the function FindMarkers() but
these tests treat each cell as an independent replicate and ignore
inherent correlations between cells originating from the same sample.
Such analyses have been shown to find a large number of false positive
associations, as has been demonstrated by Squair et al., 2021, Zimmerman
et al., 2021, Junttila et al., 2022, and others. Pseudobulking can be
used to account for such within-sample correlation. To pseudobulk, we
will use AggregateExpression() to sum together gene counts of all the
cells from the same sample for each cell type. This results in one gene
expression profile per sample and cell type. We can then perform DE
analysis using DESeq2 on the sample level. This treats the samples,
rather than the individual cells, as independent observations. -\> for
average expression of gene counts, AverageExpression() -\> one embryo
for each time point but two embryos for 16hpf (embryo = dataset)

``` r
#Pseudobulking the counts from the same sample (fish) for each cell type 
#fil.seu.obj <- readRDS(file = "/home/fgsasse_lrs_1/Downloads/BA/BA_data/Rmd/seurat_obj.rds")
#pseudo_ct_seu.obj.fil <- AggregateExpression(fil.seu.obj,
                                             #group.by = c("zebrafish_anatomy_ontology_class"),
                                             #assays = "RNA",
                                             #slot = "data",
                                            # return.seurat = TRUE)

#Pseudobulking the counts from the same sample for each cell type across all developmental stages
#pseudo_all_seu.obj.fil <- AggregateExpression(fil.seu.obj, 
                              #group.by = c("zebrafish_anatomy_ontology_class", "timepoint"),
                              #assays = "RNA",
                              #slot = "data",
                              #return.seurat = FALSE)

#Idents(pseudo_all_seu.obj.fil) <- "zebrafish_anatomy_ontology_class"

#head(FindMarkers(pseudo_all_seu.obj.fil, ident.1 = "differentiating-neurons", ident.2 = NULL ),20)

#saving the pseudobulk matrices as MatrixMarket-format 
#bulk_ct_counts <- pseudo_ct_seu.obj.fil$RNA

#psudobulk_dir <- "/home/fgsasse_lrs_1/Downloads/BA/BA_data/Pseudobulks/RNA/celltypes_times"
                  #Matrix::writeMM(bulk_ct_counts, file.path(psudobulk_dir, "matrix.mtx"))
                  #write.table(colnames(bulk_ct_counts),
                             # file.path(psudobulk_dir, "barcodes.tsv"),
                              #quote = FALSE, row.names = FALSE, col.names = FALSE)
                 # write.table(rownames(bulk_ct_counts),
                              #file.path(psudobulk_dir, "features.tsv"),
                              #quote = FALSE, row.names = FALSE, col.names = FALSE)

#bulk_all_counts <- pseudo_all_seu.obj.fil$RNA
#psudobulk_dir2 <- "/home/fgsasse_lrs_1/Downloads/BA/BA_data/Pseudobulks/RNA/celltypes_times"
                  #Matrix::writeMM(bulk_ct_counts, file.path(psudobulk_dir2, "matrix.mtx"))
                  #write.table(colnames(bulk_ct_counts),
                             # file.path(psudobulk_dir2, "barcodes.tsv"),
                              #quote = FALSE, row.names = FALSE, col.names = FALSE)
                  #write.table(rownames(bulk_ct_counts),
                              #file.path(psudobulk_dir2, "features.tsv"),
                             # quote = FALSE, row.names = FALSE, col.names = FALSE)
```

New idea for pseudobulking: Geometric mean is more robust to outliers
than arithmetic mean. - Arithmetic mean is the sum of the values divided
by the total number of values = the average of the values - Geometric
mean is the product of the values raised to the multiplicative inverse
of the total number of values

``` r
mat <- matrix(1:12, nrow=4, ncol=3)
mat
```

    ##      [,1] [,2] [,3]
    ## [1,]    1    5    9
    ## [2,]    2    6   10
    ## [3,]    3    7   11
    ## [4,]    4    8   12

``` r
meta = data.frame(
  cells = paste("cell", 1:4),
  ct = rep(c("A", "B"), each = 2)
)
meta
```

    ##    cells ct
    ## 1 cell 1  A
    ## 2 cell 2  A
    ## 3 cell 3  B
    ## 4 cell 4  B

``` r
rownames(mat) <- meta$cells

split(meta$cells, meta$ct)
```

    ## $A
    ## [1] "cell 1" "cell 2"
    ## 
    ## $B
    ## [1] "cell 3" "cell 4"

``` r
agg_fun <- function(mat, meta, fun="mean") {
  if (fun == "mean") {
    agg_f <- function(x) {
      colMeans(mat[x, ])
    }
  } else if (fun == "geom_mean") {
    agg_f <- function(x) {
      exp(colMeans(log(mat[x,])))
    }
  }
  t(sapply(split(meta$cells, meta$ct), agg_f))
}
```

``` r
agg_geom <- agg_fun(mat, meta, fun = "geom_mean")
agg_geom
```

    ##       [,1]     [,2]      [,3]
    ## A 1.414214 5.477226  9.486833
    ## B 3.464102 7.483315 11.489125

``` r
agg_mean <- agg_fun(mat, meta, fun = "mean")
agg_mean
```

    ##   [,1] [,2] [,3]
    ## A  1.5  5.5  9.5
    ## B  3.5  7.5 11.5

``` r
abs(agg_mean - agg_geom)
```

    ##         [,1]       [,2]       [,3]
    ## A 0.08578644 0.02277442 0.01316702
    ## B 0.03589838 0.01668523 0.01087471

``` r
writeLines(cell_id_fil, "cell_id.txt")
```
