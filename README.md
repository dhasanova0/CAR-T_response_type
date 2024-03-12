# CAR-T_response_type

## Scripts:
The folder code contains the scripts used for the analysis.

 ### 01_data_prep_stator: 
 - Contains scripts loading the count matrices and metadata to create Seurat objects. 
 - The data is quality filtered
 - The count matrcies are split for Stator

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


