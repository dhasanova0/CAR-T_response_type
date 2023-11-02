# CAR-T_response_type

# Scripts:

- 00.1_samples_counts.R and 00.2_metadata_exploration: Initial description of data. Input:metadata, Ouput: plots
- 01_merge_add_metadata.R: Merge all patient's samples to one Seurat object. Add metadata to the Seurat object
- 02_cell_annotation_haradhvala.R: Annotate cells with existing annotations from the publication of Haradhvala et al.
- 03_qc_filtering.R: Perform QC and filter out low quality cells
- 04_doublet_removal.R: Identicfy doublets with DoubletFinder and remove them
- 05_normalize_batch_correction.R: Normalize data, perform batch correction on patient_id and remove TCR/BCR genes
- 06_subsample_20k.R: Create downsampled seurat objects containing 20k cells (10k R, 10k NR) preserving cell proportions per patient
- 07_export_rawcm_md.R: Export raw count matrices and metadata files as ".csv" from downsampled files 
