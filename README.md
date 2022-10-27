Readme
=================
Nidia Barco Armengol

# Introduction 

The Deconv package contains seven functions that allow the user to generate graphs showing deconvolution results with various methods, and also has two functions that are useful in deconvolution. 
The Box_Deconv function generates a box plot comparing two conditions for each cell type and the Mouse_Human_converter function passes a matrix with mouse genes to a matrix with the equivalent human genes. 

# Deconvolution methods

| Method    | Description                                                                                                                                                           | Reference                                                                                                                                                                                                                                                                                                        |
|-----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| GEDIT     | - The method presents a great flexibility and precision. <br />- Has proven reference matrices.  <br />- Works well with different platforms, species and reference matrices.     | Nadel, B. B., Lopez, D., Montoya, D. J., Ma, F., Waddel, H., Khan, M. M., Mangul, S., & Pellegrini, M. (2021). The Gene Expression Deconvolution Interactive Tool (GEDIT): accurate cell type quantification from gene expression data. GigaScience, 10(2), giab002. https://doi.org/10.1093/gigascience/giab002 |
| QuanTIseq | - Method with very good results with tumor samples. <br />- Take into account the presence of type unspecified cells. <br />- It only allows you to use your own reference array. | Plattner, C., Finotello, F., & Rieder, D. (2020). Deconvoluting tumor-infiltrating immune cells from RNA-seq data using quanTIseq. Methods in enzymology, 636, 261–285. https://doi.org/10.1016/bs.mie.2019.05.056                                                                                               |
| EPIC      | - Method with good results with two own matrices but allows the use of other matrices. <br />- Take into account the presence of type unspecified cells.                    | Racle, J., & Gfeller, D. (2020). EPIC: A Tool to Estimate the Proportions of Different Cell Types from Bulk Gene Expression Data. Methods in molecular biology (Clifton, N.J.), 2120, 233–248. https://doi.org/10.1007/978-1-0716-0327-7_17                                                                      |
| FARDEEP   | - This method detects and removes outliers automatically.                                                                                                             | Hao, Y., Yan, M., Heath, B. R., Lei, Y. L., & Xie, Y. (2019). Fast and robust deconvolution of tumor infiltrating lymphocyte from expression profiles using least trimmed squares. PLoS computational biology, 15(5), e1006976. https://doi.org/10.1371/journal.pcbi.1006976                                     |
| CIBERSORT | - Method that allows a robust and reproducible deconvolution with different types of data, transcriptomics, epigenomics,...                                           | Chen, B., Khodadoust, M. S., Liu, C. L., Newman, A. M., & Alizadeh, A. A. (2018). Profiling Tumor Infiltrating Immune Cells with CIBERSORT. Methods in molecular biology (Clifton, N.J.), 1711, 243–259. https://doi.org/10.1007/978-1-4939-7493-1_12                                                            |

# Signature matrix

| Matrix                                      | Organism | Cell types  number | Reference                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
|---------------------------------------------|----------|--------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 10X immune (XImmune10)                      | Human    | 9                  | Nadel, B. B., Lopez, D., Montoya, D. J., Ma, F., Waddel, H., Khan, M. M., Mangul, S., & Pellegrini, M. (2021). The Gene Expression Deconvolution Interactive Tool (GEDIT): accurate cell type quantification from gene expression data. GigaScience, 10(2), giab002. https://doi.org/10.1093/gigascience/giab002                                                                                                                                                                                                 |
| BLUEPRINT                                   | Human    | 8                  | Martens, J. H., & Stunnenberg, H. G. (2013). BLUEPRINT: mapping human blood cell epigenomes. Haematologica, 98(10), 1487–1489. https://doi.org/10.3324/haematol.2013.094243                                                                                                                                                                                                                                                                                                                                      |
| BlueCode                                    | Human    | 35                 | Nadel, B. B., Lopez, D., Montoya, D. J., Ma, F., Waddel, H., Khan, M. M., Mangul, S., & Pellegrini, M. (2021). The Gene Expression Deconvolution Interactive Tool (GEDIT): accurate cell type quantification from gene expression data. GigaScience, 10(2), giab002.  https://doi.org/10.1093/gigascience/giab002                                                                                                                                                                                                |
| Human Primary Cell Atlas                    | Human    | 13                 | Mabbott, N. A., Baillie, J. K., Brown, H., Freeman, T. C., & Hume, D. A. (2013). An expression atlas of human primary cells: inference of gene function from coexpression networks. BMC genomics, 14, 632. https://doi.org/10.1186/1471-2164-14-632                                                                                                                                                                                                                                                              |
| ImmunoStates                                | Human    | 20                 | Vallania, F., Tam, A., Lofgren, S., Schaffert, S., Azad, T. D., Bongen, E., Haynes, W., Alsup, M., Alonso, M., Davis, M., Engleman, E., & Khatri, P. (2018). Leveraging heterogeneity across multiple datasets increases cell-mixture deconvolution accuracy and reduces biological and technical biases. Nature communications, 9(1), 4735. https://doi.org/10.1038/s41467-018-07242-6                                                                                                                          |
| LM22                                        | Human    | 22                 | Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., Hoang, C. D., Diehn, M., & Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature methods, 12(5), 453–457. https://doi.org/10.1038/nmeth.3337                                                                                                                                                                                                                                             |
| Skin Signatures                             | Human    | 21                 | Swindell, W. R., Johnston, A., Voorhees, J. J., Elder, J. T., & Gudjonsson, J. E. (2013). Dissecting the psoriasis transcriptome: inflammatory- and cytokine-driven gene expression in lesions from 163 patients. BMC genomics, 14, 527. https://doi.org/10.1186/1471-2164-14-527                                                                                                                                                                                                                                |
| Blood circulating immune  cells (EPIC_BCIC) | Human    | 6                  | Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. eLife, 6, e26476. https://doi.org/10.7554/eLife.26476                                                                                                                                                                                                                                                                            |
| Tumor infiltrating cells (EPIC_TIC)         | Human    | 6                  | Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. eLife, 6, e26476. https://doi.org/10.7554/eLife.26476                                                                                                                                                                                                                                                                            |
| TIL10                                       | Human    | 10                 | Finotello, F., Mayer, C., Plattner, C., Laschober, G., Rieder, D., Hackl, H., Krogsdam, A., Loncova, Z., Posch, W., Wilflingseder, D., Sopper, S., Ijsselsteijn, M., Brouwer, T. P., Johnson, D., Xu, Y., Wang, Y., Sanders, M. E., Estrada, M. V., Ericsson-Gonzalez, P., Charoentong, P., … Trajanoski, Z. (2019). Molecular and pharmacological modulators 21 of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome medicine, 11(1), 34. https://doi.org/10.1186/s13073-019-0638-6 |
| Wang2020                                    | Human    | 6                  | Wang, J., Devlin, B., & Roeder, K. (2020). Using multiple measurements of tissue to estimate subject-  and cell-type-specific gene expression. Bioinformatics (Oxford, England),36 (3), 782–788. https://doi.org/10.1093/bioinformatics/btz619                                                                                                                                                                                                                                                                   |
| Tabula Muris                                | Mouse    | 12                 | Tabula Muris Consortium, Overall coordination, Logistical coordination, Organ collection and processing, Library preparation and sequencing, Computational data analysis, Cell type annotation, Writing group, Supplemental text writing group, & Principal investigators (2018). Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. Nature, 562(7727), 367–372. https://doi.org/10.1038/s41586-018-0590-4                                                                                   |
| Mouse Body Atlas                            | Mouse    | 20                 | Lattin, J. E., Schroder, K., Su, A. I., Walker, J. R., Zhang, J., Wiltshire, T., Saijo, K., Glass, C. K., Hume, D. A., Kellie, S., & Sweet, M. J. (2008). Expression analysis of G Protein-Coupled Receptors in mouse macrophages. Immunome research, 4, 5. https://doi.org/10.1186/1745-7580-4-5                                                                                                                                                                                                                |
| ImmGen                                      | Mouse    | 137                | Heng, T. S., Painter, M. W., & Immunological Genome Project Consortium (2008). The Immunological Genome Project: networks of gene expression in immune cells. Nature immunology, 9(10), 1091–1094. https://doi.org/10.1038/ni1008-1091                                                                                                                                                                                                                                                                           |
| sig_matr_seqImmuCC (seqImmuCC)              | Mouse    | 10                 | Chen, Z., Quan, L., Huang, A., Zhao, Q., Yuan, Y., Yuan, X., Shen, Q., Shang, J., Ben, Y., Qin, F. X., & Wu, A. (2018).  seq-ImmuCC: Cell-Centric View of Tissue Transcriptome Measuring Cellular Compositions of Immune  Microenvironment From Mouse RNA-Seq Data.Frontiers in immunology, 9, 1286. https://doi.org/10.3389/fimmu.2018.01286                                                                                                                                                                    |
| Donovan_Brain.mouse                         | Mouse    | 7                  | Donovan, M., D'Antonio-Chronowska, A., D'Antonio, M., & Frazer, K. A. (2020).  Cellular deconvolution of GTEx tissues powers discovery of disease and cell-type associated regulatory variants. Nature communications,11 (1), 955. https://doi.org/10.1038/s41467-020-14561-0                                                                                                                                                                                                                                    |

# Functions

| Function              | Description                                                                                                                                                                                                                         |
|-----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| CIBERSORT_Deconv      | Function that generates the relative or absolute deconvolution plots with the CIBERSORT method,  with the non-fractional plot, returns a df resulting from the deconvolution and two graphs of deconvolution data.                  |
| FARDEEP_Deconv        | Function that generates the relative or absolute deconvolution plots with the FARDEEP method,  with the non-fractional plot, and returns a df resulting from the deconvolution.                                                     |
| EPIC_Deconv           | Function that generates the relative deconvolution plots with the EPIC method,  with the non-fractional plot, and returns a df resulting from the deconvolution.                                                                    |
| QuanTIseq_Deconv      | Function that generates the relative deconvolution graphs with the QuanTIseq method using the TIL10 reference matrix.  This function generates the two graphs without fractionation and returns a df with the deconvolution result. |
| GEDIT_Deconv          | Function that generates the relative deconvolution graphs with the GEDIT method.  This function generates the two graphs without fractionation and returns a df with the deconvolution result.                                      |
| Mouse_Human_converter | This function convert mouse expression data into human expression data, changing the mouse genes for her analog in human.                                                                                                           |
| Box_Deconv            | Function that generates a boxplot comparing two conditions for each cell type from a df and a condition vector.                                                                                                                     |

# Instalation 

``` r
library(devtools)
devtools::install_github("margenomics/Deconv")
library(Deconv)
```

# Usage

```r
# Charge signature matrix
sig.matrix<- Deconv::LM22
sig.matrix<- Deconv::EPIC_BCIC
sig.matrix<- Deconv::EPIC_TIC
sig.matrix<- Deconv::seqImmuCC
sig.matrix<- Deconv::TabulaMuris
sig.matrix<- Deconv::TIL10
sig.matrix<- Deconv::XImmune10
sig.matrix<- Deconv::Donovan_Brain.mouse
sig.matrix<- Deconv::Wang2020

# Function --> CIBERSORT_Devonv
c<- CIBERSORT_Deconv(matrix = matrix, sig.matrix = sig.matrix, results_dir = results_dir, 
cibersortpath = path_ciber, by_Cond = TRUE, cond = fractions, method = "abs")

# Function --> FARDEEP_Deconv
c<- FARDEEP_Deconv(matrix = matrix, sig.matrix = sig.matrix, results_dir = results_dir, byCond 
= TRUE, cond = fractions, method = "abs")

# Function --> EPIC_Deconv 
c<- EPIC_Deconv(matrix = matrix, sig.matrix = sig.matrix, results_dir = results_dir, byCond = 
TRUE, cond = fractions)

# Function --> QuanTIseq_Deconv 
c<- QuanTIseq_Deconv(matrix = matrix, results_dir = results_dir, byCond = TRUE, cond = 
fractions)

# Function --> GEDIT_Deconv
c<- GEDIT_Deconv(matrix = matrix, sig.matrix = sig.matrix, results_dir = results_dir, byCond = 
TRUE, cond = fractions

# Function --> Mouse_Human_converter 
M<- Mouse_Human_converter(matrix = matrix,resultant_file = "counts_mf_human.txt", results_dir
= results_dir)

# Function --> Box_Deconv
Box_Deconv(data= c, cond= fractions, results_dir = "/bicoh/nidia/Deconv/", f_name= 
"Box_CIBERSORT_abs")
```
# Examples
