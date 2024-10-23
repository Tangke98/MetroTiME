### Step 1. Metabolic subtypes MetaRegulon inference
    # Fibroblasts
    R MetaRegulon_Fibroblasts.R integrated_object_path_fibroblasts fibroblasts_DEG fibroblasts_DEM MetaModule MetaModule_info output_path

    # Myeloid
    R MetaRegulon_Myeloid.R integrated_object_path_myeloid myeloid_DEG myeloid_DEM MetaModule MetaModule_info output_path

### Step 2. MetaRegulon analysis
    # Fibroblasts
    R MetaRegulon_analysis_Fibroblasts.R integrated_object_path_fibroblasts fibroblasts_DEG fibroblasts_DEM MetaModule MetaModule_info output_path

    # Myeloid
    R MetaRegulon_analysis_Myeloid.R integrated_object_path_myeloid myeloid_DEG myeloid_DEM MetaModule MetaModule_info output_path
    