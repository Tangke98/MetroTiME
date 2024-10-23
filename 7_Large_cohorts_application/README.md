### Step 1. Metabolic subtypes application in large cohort
    ## For clinical survival outcome
    R TCGA_MetaModule_survival.R MetaModule TPM_path clin_path output_path

    ## For immune function
    R TCGA_MetaModule_immune.R MetaModule MetaModule_info TPM_path immune_path output_path

    ## For ICB response
    R ICB_MetaModule.R ICB_data Tide_res MetaModule MetaModule_info output_path
    R ICB_MetaModule_analysis.R ICB_data MetaModule MetaModule_info output_path


### Step 2. Metabolic regulators application in large cohort
    ## For clinical survival outcome
    R TCGA_MetaModule_MetaRegulon_survival.R MM_survival factor_Fibroblasts factor_Myeloid gsva_path tpm_path clin_path output_path

    ## For drug response
    Drug_response.R myeloid_path fibroblast_path survival_path expression_path TF_list auc_value output_path



