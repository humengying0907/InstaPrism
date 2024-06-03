# InstaPrism
**InstaPrism** is an R package for cell type composition and gene expression deconvolution in bulk RNA-Seq data. Based on the same conceptual framework as [BayesPrism](https://github.com/Danko-Lab/BayesPrism), InstaPrism significantly increases calculation speed while providing nearly identical deconvolution results.

## Installation
```````
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("Biobase", quietly = TRUE))
    BiocManager::install("Biobase")

if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("humengying0907/InstaPrism")
```````
## Deconvolution with InstaPrism
Quick start with the following example
```````R
# bulk input with genes in rows and samples in columns; this example data is from ovarian cancer
bulk_expr = read.csv(system.file('extdata','example_bulk.csv',package = 'InstaPrism')) 

# load built-in reference; choose a reference that matches the tumor type of the bulk input
OV_ref = InstaPrism_reference('OV') 

# deconvolution with InstaPrism
deconv_res = InstaPrism(bulk_Expr = bulk_expr,refPhi_cs = OV_ref)

# The deconvoled cell type fraction is now accessible with
estimated_frac = t(deconv_res@Post.ini.ct@theta)
head(estimated_frac)

# InstaPrism also returns the deconvolved gene expression Z
Z = get_Z_array(deconv_res) # a sample by gene by cell-type array
head(Z[,1:10,'malignant'])
```````

## Deconvolution results comparison with BayesPrism
Using either scRNA-based reference (update = F) or updated reference (update = T), InstaPrism achieves nearly identical deconvolution results as BayesPrism.

<img src="https://github.com/humengying0907/InstaPrism/assets/54827603/36c6cfa1-308c-4a0b-adc2-bf4b6399139b" width=60% height=60%>

## Running time comparsion with BayesPrism
Below is a running time comparsion when running deconvolution on the tutorial [data](https://github.com/Danko-Lab/BayesPrism/tree/main/tutorial.dat) 
provided in BayesPrism. 

<img src="https://github.com/humengying0907/InstaPrism/assets/54827603/8e158249-9cc9-4f06-8e89-63867540bfc6" width=45% height=45%>

## Memory Comparison with BayesPrism
InstaPrism significantly reduced the memory usaged during deconvolution. Below is a memory usage comparsion when running deconvolution on the tutorial [data](https://github.com/Danko-Lab/BayesPrism/tree/main/tutorial.dat) 
provided in BayesPrism. 

<img src="https://github.com/humengying0907/InstaPrism/assets/54827603/3d9c8b8b-8aac-4c4b-b793-e64c33cac752" width=45% height=45%>

## Built-in Reference for InstaPrism
We have provided precompiled reference tailored for a wide range of cancer types. The reference can be downloaded directly using the download link below, or loaded directly using the following command
 ```````
BRCA_ref = InstaPrism_reference('BRCA')
```````



| reference name | tumor type                       | #cells used for reference construction | #cell types/cell states | umap | citation | download |
|----------------|----------------------------------|----------------------------------------|-------------------------|------|----------|----------|
| BRCA_refPhi    | breast cancer                    | 100,064                                    | 8/76                             |[UMAP](https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers)      | [Wu et al. 2021](https://www.nature.com/articles/s41588-021-00911-1) | [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/BRCA_refPhi.RDS) |
| CRC_refPhi     | colorectal cancer                | 371,223                                    | 15/98                           |[UMAP](https://singlecell.broadinstitute.org/single_cell/study/SCP1162/human-colon-cancer-atlas-c295?scpbr=human-cell-atlas-main-collection)      | [Pelka et al. 2021](https://www.cell.com/cell/fulltext/S0092-8674(21)00945-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421009454%3Fshowall%3Dtrue) |  [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/CRC_refPhi.RDS) |
| GBM_refPhi     | glioblastoma                     | 338,564                                    | 10/57                             |[cellxgeneLink](https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c), [UMAP](https://cellxgene.cziscience.com/e/56c4912d-2bae-4b64-98f2-af8a84389208.cxg/)      | [Ruiz et al. 2022](https://www.biorxiv.org/content/10.1101/2022.08.27.505439v1) |[↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/GBM_refPhi.RDS) |
| LUAD_refPhi    | lung adenocarcinomas             | 118,293                                    | 13/77                            |[UMAP](https://www.weizmann.ac.il/sites/3CA/study-data/umap/20115)      | [Xing et al. 2021](https://www.science.org/doi/10.1126/sciadv.abd9738) | [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/LUAD_refPhi.RDS) |
| OV_refPhi      | ovarian cancer                   | 929,690                                    | 9/40                             |[cellxgeneLink](https://cellxgene.cziscience.com/collections/4796c91c-9d8f-4692-be43-347b1727f9d8), [UMAP](https://cellxgene.cziscience.com/e/b252b015-b488-4d5c-b16e-968c13e48a2c.cxg/)      | [Vazquez et al. 2022](https://www.nature.com/articles/s41586-022-05496-1) | [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/OV_refPhi.RDS) |
| RCC_refPhi     | clear cell renal cell carcinoma  | 270,855                                    | 11/106                            |[cellxgeneLink](https://cellxgene.cziscience.com/collections/f7cecffa-00b4-4560-a29a-8ad626b8ee08), [UMAP](https://cellxgene.cziscience.com/e/5af90777-6760-4003-9dba-8f945fec6fdf.cxg/)      | [Li et al. 2022](https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00548-7) | [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/RCC_refPhi.RDS) |
| SKCM_refPhi    | skin cutaneous melanoma          | 4,645                                      | 8/23                             |[UMAP](https://www.weizmann.ac.il/sites/3CA/study-data/umap/20111)      | [Tirosh et al. 2016](https://www.science.org/doi/10.1126/science.aad0501) | [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/SKCM_refPhi.RDS) |

## Reference evaluation pipeline
To validate the performance of a reference, or to compare between different references, please refer to the [Reference_evalu_pipeline](https://github.com/humengying0907/InstaPrismSourceCode) we developed.

## Tutorial
Check [InstaPrism_tutorial](https://humengying0907.github.io/InstaPrism_tutorial.html) for detailed implementation of InstaPrism and compare its performance with BayesPrism.

## Reference
M. Hu and M. Chikina, “InstaPrism: an R package for fast implementation of BayesPrism.” bioRxiv, p. 2023.03.07.531579, Mar. 12, 2023.
doi: https://doi.org/10.1101/2023.03.07.531579
