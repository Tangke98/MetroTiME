### Step 1. Batch effect remeval 
    ## For myeloid cells
    # Merge metacell  
    R CCA_corrected_batch_effect_for_metacell_Myeloid.R Seurat_TPM_myeloid Seurat_TPM_metadata_myeloid TF_Ligand_Metabolic integrated_object_path_myeloid


    ## For fibroblasts
    R CCA_corrected_batch_effect_for_metacell_Fibroblasts.R Seurat_TPM_fibroblasts Seurat_TPM_metadata_fibroblasts TF_Ligand_Metabolic integrated_object_path_fibroblasts

### Step 2. Identifying the optimal clustering parameters
    ## For myeloid cells
    # Silhouette score
    R Silhouette_ASW_Myeloid.R integrated_object_path_myeloid output_path

    # Clustree res
    R Clustree_Myeloid.R integrated_object_path_myeloid output_path


    ## For fibroblasts
    # Silhouette score
    R Silhouette_ASW_Fibroblasts.R integrated_object_path_fibroblasts output_path

    # Clustree res
    R Clustree_Fibroblasts.R integrated_object_path_fibroblasts output_path

### Step 3. Cell subtype annotation for Metacells
    ## For myeloid cells
    # Cell subtype annotation 
    R Cell_subtypes_annotation_Myeloid.R integrated_object_path_myeloid output_path


    ## For fibroblasts
    R Cell_subtypes_annotation_Fibroblasts.R integrated_object_path_fibroblasts output_path

### Step 4. Benchmark of cell subtype annotation for Metacells
    ## For myeloid cells
    R Rogue_score_Myeloid.R integrated_object_path_myeloid output_path
    

    ## For fibroblasts
    R Rogue_score_Fibroblasts.R integrated_object_path_fibroblasts output_path

