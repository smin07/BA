# Linking scATAC-seq to scRNA-seq to identify CREs in zebrafish development
Multicellular organisms generate diverse cell types from a shared genome through differential gene regulation, but how cis-regulatory elements (CREs) coordinate this process at single-cell resolution remains poorly characterized. Using zebrafish, a vertebrate model with transparent, rapidly developing embryos, we integrate single-cell chromatin accessibility (scATAC-seq) and gene expression (scRNA-seq) data across hundreds of cell states and developmental stages. By aggregating similar cells into pseudobulks and metacells, we link accessible chromatin regions to nearby gene expression using linear models. This framework uncovers distinct classes of gene regulatory landscapes and offers a systematic approach to connect CREs with the expression dynamics underlying cellular diversification with implications for understanding development and gene regulation across vertebrates.

## Pearson r
<h5>Correlation-based peak-gene pairs across cell types</h5>

<table>
  <tr>
    <td align="center">
      <b>A</b><br>
      <img width="400" height="250" alt="first_pair_count_cor_ct" src="https://github.com/user-attachments/assets/7b60fc66-3dce-49fb-bcec-ad98f95ee661" />
    </td>
    <td align="center">
      <b>B</b><br>
      <img width="400" height="250" alt="first_pair_count_cor_ct2" src="https://github.com/user-attachments/assets/5a1aab16-7a98-4f94-8b2e-535a84b9751c" />
    </td>
  </tr>

**A.** Distribution of identified peak-gene pairs across cell types in different TSS-centered windows.  
**B.** Distribution of signigicant peak-gene pairs across cell types in different TSS-centered windows.

<h5>Correlation-based peak-gene pairs across cell types + developmental stages</h5>

<table>
  <tr>
    <td align="center">
      <b>C</b><br>
      <img width="400" height="250" alt="first_pair_count_cor_ct_time" src="https://github.com/user-attachments/assets/b2168e44-ef04-45a1-82db-70de74218add" />
    </td>
    <td align="center">
      <b>D</b><br>
      <img width="400" height="250" alt="first_pair_count_cor_ct_time2" src="https://github.com/user-attachments/assets/91c3958b-197c-47ea-9263-da079207cc2b" />
    </td>
  </tr>
</table>

**C.** Distribution of identified peak-gene pairs across cell types and developmental stages in different TSS-centered windows.  
**D.** Distribution of signigicant peak-gene pairs across cell types and developmental stages in different TSS-centered windows.

In both pseudobulk types, the number of non‑significant peaks per gene increases steadily with larger window sizes, reflecting the growing number of candidate regions captured at broader genomic ranges (Figure A,C). In cell‑type pseudobulks, most genes carry only a few significantly correlated peaks across all four TSS‑centered windows, with positively correlated peaks remaining consistently rare (Figure B). Negatively correlated peaks are present at all window sizes but become increasingly prominent at ±100 kb, suggesting that putative repressive regulatory elements are predominantly distal. In the cell type × developmental stage pseudobulks, the overall distribution shape is preserved, but genes exhibit higher peak counts across all window sizes compared to cell‑type‑only aggregation (Figure D). This increase is consistent with the capture of context‑specific accessible regions that are averaged out when collapsing across developmental stages. Notably, negatively correlated peaks become more prominent from ±50 kb onward, indicating that stage‑aware aggregation increases sensitivity to distal regulatory associations, including putative repressive elements that remain undetected in the coarser cell‑type pseudobulks.

## Ordinary Least Squares (OLS)
<h5>OLS-based peak-gene pairs across cell types</h5>
 <table>
  <tr>
    <td align="center">
      <b>A</b><br>
        <img width="1383" height="589" alt="image" src="https://github.com/user-attachments/assets/a278e810-6c8e-4994-8040-3fafe4c45085" />
    </td>
</table>

## Repository layout

The analysis code is being organized around a standard `src/` layout:

- `src/your_package/` for reusable data-processing, plotting, and modeling helpers
- `notebooks/` for exploratory and final analysis notebooks
- `scripts/` for runnable pipeline entry points
- `tests/` for unit tests
- `data/` for local data staging and generated artifacts


