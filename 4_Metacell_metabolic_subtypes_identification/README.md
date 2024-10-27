### Step 1. Identifying the optimal clustering paramaters
    # Silhouette score 
    R Silhouette_ASW_metabolic.R cell_lineage integrated_object_path Metabolic_gene output_path
    
    # Clustree res
    R Clustree_metabolic.R cell_lineage integrated_object_path Metabolic_gene output_path

### Step 2. Metabolic subtype identification for Metacells
    # Metabolic subtypes annotation 
    R Metabolic_subtypes_annotation.R cell_lineage integrated_object_path Metabolic_gene output_path MetaModule  MetaModule_info

### Step 3. Benchmark of metabolic subtype annotation for Metacells
    R Rogue_score_metabolic.R cell_lineage integrated_object_path metabolic_gene output_path

