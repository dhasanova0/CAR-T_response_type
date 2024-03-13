# Identifying biological factors influencing response type in CAR-T cell therapy

## Description
This study aims to discover biological factors associated with response type in CAR-T cell therapy from the transcriptome of PBMCs prior to treatment from patients diagnosed with Large B cell lymphoma. The study aims to understand the underlying biological factors in patients' PBMCs the could influence response type. This eventually would allow to improve patient stratification and future clinical trial designs. The study was conducted on publicly available scRNA-seq data originating from the study of [Haradhvala et. al](https://doi.org/10.1038/s41591-022-01959-0). Common approaches in scRNA-seq were applied to extract meaningful biological insights from the data. In addition, the novel statistical method, [Stator](https://doi.org/10.1101/2023.12.18.572232), was applied to further characterize the cellular heterogenity linked to a given response type (responder, non-responder). This repository contains the scripts that were used for the analysis. 

## Scripts
The folder **code** contains the scripts used for the analysis.

 ### 01_data_prep_stator: 
 - Contains scripts loading the count matrices and metadata to create Seurat objects. 
 - The data is quality filtered
 - The count matrices are split for Stator

### 02_data_exploration
- Contains scripts for metadata exploration
- Characterizes scRNA-seq data

### 03_create_md_stator
- Creates metadata for input in Stator of the cells from which interactions were inferred

### 04_projections
- Projection of sates to other cells
- Analysis of data with projected cells

### 05_downstream_analyses
- Comparison of various Stator runs
- Comparison of Stator results and initial analysis
- State enrichemnt analyses

get_state_function.R: Function which labels cells with states (cells already labeled with interactions)

qc_figures.R: Function which creates various figures for quality control

## Data
Raw count matrices obtained from: [GEO accession GSE197268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197268)

Metadata obtained from: [Haradhvala_et_al_2022](https://github.com/getzlab/Haradhvala_et_al_2022)

Input files for the scripts (except for the scripts in 01_data_prep_stator) are provided here: [Input files](https://polybox.ethz.ch/index.php/s/ehdTC9k7jRT7g5j)
- The input paths were preserved. Each script's working directory (wd) must be adjusted to the user's wd. 

Relevant outputs that were discussed in the thesis are provided here: [Results](https://polybox.ethz.ch/index.php/s/YVN7hLXujlEb24T)
