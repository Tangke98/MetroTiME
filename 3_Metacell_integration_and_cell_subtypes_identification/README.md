### Step 1. Batch effect remeval 
    # Merge Metacell  
    R CCA_corrected_batch_effect_for_Metacell.R cell_lineage Seurat_TPM Seurat_TPM_metadata TF_Ligand_Metabolic integrated_object_path

### Step 2. Identifying the optimal clustering parameters
    # Silhouette score
    R Silhouette_ASW.R cell_lineage integrated_object_path output_path

    # Clustree res
    R Clustree.R cell_lineage integrated_object_path output_path

### Step 3. Cell subtype annotation for Metacells
    # Cell subtype annotation 
    R Cell_subtypes_annotation.R cell_lineage integrated_object_path output_path

### Step 4. Benchmark of cell subtype annotation for Metacells
    R Rogue_score.R cell_lineage integrated_object_path output_path