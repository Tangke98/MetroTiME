### Step 1. Metabolic subtypes correspondance to the cell subtypes
    R Metabolic_subtypes_cell_subtypes_correspondance.R cell_lineage integrated_object_path MetaModule_info metabolic_gene output_path
 
### Step 2. Metabolic subtypes correspondance to the cell function
    ### Hallmarks score and functions score calculation
    R Metabolic_subtypes_hallmarks_functions.R cell_lineage integrated_object_path metabolic_gene hallmark_path output_path

    ### MetaModule scores correlation with cell functions
    R Features_MetaModule_association.R  cell_lineage integrated_object_path metabolic_score feature_score DEM MetaModule_info output_path

### Step 3. Metabolic subtypes trajectory
    ### Monocle analysis
    R Monocle_VAR_DEG.R cell_lineage integrated_object_path markers output_path

    ### Stream analysis
    R Stream_prepare.R cell_lineage integrated_object_path metabolic_score output_path
    Stream.py cell_lineage file work_dir
    R Stream_analysis.R cell_lineage DEM output_path

### Step 4. Metabolic subtypes within a cell subtypes
    R Metabolic_subtypes_functions_within_cell_subtypes_Fibroblasts.R integrated_object_path_fibroblasts output_path
    