### Step 1. Identifying the optimal clustering paramaters
    ## For myeloid cells
    # Silhouette score
    R Silhouette_ASW_metabolic_Myeloid.R integrated_object_path_myeloid Metabolic_gene output_path
    
    # Clustree res
    R Clustree_metabolic_Myeloid.R integrated_object_path_myeloid Metabolic_gene output_path


    ## For fibroblasts
    # Silhouette score
    R Silhouette_ASW_metabolic_Fibroblasts.R integrated_object_path_fibroblasts Metabolic_gene output_path

    # Clustree res
    R Clustree_metabolic_Fibroblasts.R integrated_object_path_fibroblasts Metabolic_gene output_path

### Step 2. Metabolic subtype identification for Metacells
    ## For myeloid cells
    # Metabolic subtypes annotation 
    R Metabolic_subtypes_annotation_Myeloid.R integrated_object_path_myeloid Metabolic_gene output_path MetaModule  MetaModule_info


    ## For fibroblasts
    # Metabolic subtypes annotation 
    R Metabolic_subtypes_annotation_Fibroblasts.R integrated_object_path_fibroblasts Metabolic_gene output_path MetaModule  MetaModule_info

### Step 3. Benchmark of metabolic subtype annotation for Metacells
    ## For myeloid cells
    R Rogue_score_metabolic_Myeloid.R integrated_object_path_myeloid metabolic_gene output_path


    ## For fibroblasts
    R Rogue_score_metabolic_Fibroblasts.R integrated_object_path_fibroblasts metabolic_gene output_path

