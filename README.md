# Linking scATAC-seq to scRNA-seq to identify CREs in zebrafish development
Multicellular organisms generate diverse cell types from a shared genome through differential gene regulation, but how cis-regulatory elements (CREs) coordinate this process at single-cell resolution remains poorly characterized. Using zebrafish, a vertebrate model with transparent, rapidly developing embryos, we integrate single-cell chromatin accessibility (scATAC-seq) and gene expression (scRNA-seq) data across hundreds of cell states and developmental stages. By aggregating similar cells into pseudobulks and metacells, we link accessible chromatin regions to nearby gene expression using linear models. This framework uncovers distinct classes of gene regulatory landscapes and offers a systematic approach to connect CREs with the expression dynamics underlying cellular diversification with implications for understanding development and gene regulation across vertebrates.

## Pearson r
### Correlation-based peak-gene pairs across cell types

<p align="center">
  <img width="400" height="250" src="https://github.com/user-attachments/assets/434139fb-0788-45db-8ec0-fbf59427678a" />
  <img width="400" height="250" src="https://github.com/user-attachments/assets/9f2c96bd-f0aa-486d-9a79-474aef8d2c34" />
</p>

**Figure A–B.** Correlation-based peak-gene links inferred across cell types.  
**A.** Distribution of identified peak-gene pairs across cell types.  
**B.** Summary of the genomic distance between correlated peaks and their target genes, highlighting the prevalence of distal regulatory interactions.

### Correlation-based peak-gene pairs across cell types + developmental stages

<p align="center">
  <img width="400" height="250" src="https://github.com/user-attachments/assets/84e93eab-e2cd-499d-8db2-9788a3f4e810" />
  <img width="400" height="250" src="https://github.com/user-attachments/assets/ed0d73a3-34e9-42c5-a2ff-e498cea7e9a6" />
</p>

**Figure C–D.** Correlation-based peak-gene links inferred across both cell types and developmental stages.  
**C.** Number and distribution of significant peak-gene pairs identified when developmental stage information is incorporated.  
**D.** Characteristics of the inferred regulatory interactions, including peak-to-gene distances and stage-specific patterns.
