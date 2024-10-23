### Step 1. Metabolic subtypes correspondance to the cell subtypes
    ## For fibroblasts
    R Metabolic_subtypes_cell_subtypes_correspondance_Fibroblasts.R integrated_object_path_fibroblasts MetaModule_info metabolic_gene output_path
    

    ## For myeloid cells
    R Metabolic_subtypes_cell_subtypes_correspondance_Myeloid.R integrated_object_path_myeloid MetaModule_info metabolic_gene output_path

### Step 2. Metabolic subtypes correspondance to the cell function
    ## For fibroblasts
    ### Hallmarks score and functions score calculation
    R Metabolic_subtypes_hallmarks_functions_Fibroblasts.R integrated_object_path_fibroblasts metabolic_gene hallmark_path output_path

    ### MetaModule scores correlation with cell functions
    R Features_MetaModule_association_Fibroblasts.R  integrated_object_path_fibroblasts fibroblasts_metabolic_score fibroblasts_feature_score fibroblasts_DEM MetaModule_info output_path


    ## For myeloid cells
    ### Hallmarks score and functions score calculation
    R Metabolic_subtypes_hallmarks_functions_Myeloid.R integrated_object_path_myeloid metabolic_gene hallmark_path output_path

    ### MetaModule scores correlation with cell functions
    R Features_MetaModule_association_Myeloid.R  integrated_object_path_myeloid myeloid_metabolic_score myeloid_feature_score myeloid_DEM MetaModule_info output_path

### Step 3. Metabolic subtypes trajectory
    ## For fibroblasts
    ### Monocle analysis
    R Monocle_VAR_DEG_Fibroblasts.R integrated_object_path_fibroblasts fibroblasts_markers output_path

    ### Stream analysis
    R Stream_prepare_Fibroblasts.R integrated_object_path_fibroblasts fibroblasts_metabolic_score output_path
    Stream_fibroblasts.py firboblast_file work_dir
    R Stream_fibroblasts_analysis.R fibroblasts_DEM output_path


    ## For myeloid cells
    ### Monocle analysis
    R Monocle_VAR_DEG_Myeloid.R integrated_object_path_myeloid myeloid_markers output_path

    ### Stream analysis
    R Stream_prepare_Myeloid.R integrated_object_path_myeloid myeloid_metabolic_score output_path
    Stream_Myeloid.py Myeloid_file work_dir
    R Stream_myeloid_analysis.R myeloid_DEM output_path

### Step 4. Metabolic subtypes within a cell subtypes
    R Metabolic_subtypes_functions_within_cell_subtypes_Fibroblasts.R integrated_object_path_fibroblasts output_path
    