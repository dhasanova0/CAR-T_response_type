# CAR-T_response_type

## Scripts:
The folder code contains the scripts used for the analysis.

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

## Data:
Raw count matrices obtained from: [GEO accession GSE197268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197268)

Metadata obtained from: [Haradhvala_et_al_2022](https://github.com/getzlab/Haradhvala_et_al_2022)

Input files for the scripts (except for the scripts in 01_data_prep_stator) are provided here: [Input files](https://polybox.ethz.ch/index.php/s/ehdTC9k7jRT7g5j)

The input paths were preserved. Each script's working directory (wd) must be adjusted to the user's wd. 

Relevant outputs that were discussed in the thesis are provided here: [Results](https://polybox.ethz.ch/index.php/s/YVN7hLXujlEb24T)
