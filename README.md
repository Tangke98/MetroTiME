# MetroTIME: Pan-cancer single-cell metabolism and regulatory analyses identifies microenvironment metabolic subtypes associated with tumor immunity

Metabolic reprogramming is a hallmark of cancer, yet the mechanisms underlying metabolic state heterogeneities and their regulation within the tumor microenvironment are still poorly understood. MetroSCREEN synergistically integrates intrinsic gene regulation with extrinsic cellular interactions to accurately predict potential regulators, demonstrating high consistency with results from single-cell perturbations. We developed a framework MetroTIME, which applied MetroSCREEN to a comprehensive pan-cancer single-cell atlas encompassing over 700,000 fibroblasts and myeloid cells across 36 cancer types. This framework consists of five main modules: tumor-related single-cell datasets collection, integrated data analysis using MetroSCREEN, pan-cancer MetaModule analysis, pan-cancer MetaRegulon analysis, and large-scale cohort validation.
<img src="https://github.com/yahan9416/TabulaTiME/blob/main/Image/TabulaTIME_workflow.png" alt="Image Description" width="100%" />
The MetroTIME can be explored and  visualization at [MetroTIME Website](http://wanglab-compbio.cn/MetroTIME/).

## 1. scRNA-seq data collection and preprocessing
To illuminate the intricate dynamics of the tumor microenvironment (TME) throughout the genesis and progression of tumors across a spectrum of cancer types, we aggregated single-cell RNA sequencing (scRNA-Seq) datasets from 746 donors representing 36 distinct cancer types. 
All collected datasets underwent rigorous pre-processing through the MAESTRO workflow, which encompassed quality control, elimination of doublets and batch effects, cell clustering, and precise cell type annotation. Ultimately, a total of 103 studies were collated, comprising an impressive 4,479,563 cells.

## 2. MetaCell identification
To mitigate the impact of technical disturbances and minimize computing resource expenses, TabulaTIME employed a strategy of grouping cells with similar expressions into MetaCells within each dataset, ensuring that each MetaCell encompassed approximately 30 cells. The average log TPM-transformed gene expression of all cells within each MetaCell was utilized as a representative measure of the MetaCell's expression.

## 3. Metacell integration and cell subtypes identification
For the entirety of metacells derived from all scRNA-seq datasets, we conducted an integrative analysis, subsequently assessing and correcting any prevailing batch effects.
To delve deeper into the metacell heterogeneity within specific cell types, we collected cells from fibroblasts and myeloid cells for downstream analysis. This pan-cancer single-cell atlas encompassing over 700,000 fibroblasts and myeloid cells across 36 cancer types. 


## 4. Metacell metabolic subtypes identification
To investigate the metabolic heterogeneity of fibroblasts and myeloid cells, metabolic genes were used to cluster the Metacells. To annotate the metabolic subtypes, MetroSCREEN was applied to the Metacells of these two lineages.

## 5. Characterization of cell metabolic subtypes
To facilitate the characterization of metabolic subtypes, we examined the correspondence between metabolic subtypes and cell subtypes, as well as between metabolic subtypes and cell functions. Additionally, we investigated the variation in metabolic states along lineage differentiation. 

## 6. Chracterization of MetaRegulon for MetaModule enriched in metabolic subtype
To explore the metabolic regulators of metabolic states, we applied MetroSCREEN. The regulators enriched in each metabolic subtype were identified and the metabolic regulatory networks for these states was constructed.

## 7. Large cohorts application
To explore the relationship between metabolic states, their regulators, and both clinical survival outcomes and treatment responses, MetroSCREEN was utilized to large-scale tumor profiles from TCGA, ICB-treated patients, and the Cancer Cell Line Encyclopedia (CCLE).

## Environment 
    Ubuntu 9.3.0
    R version 4.0.5	
    Python version 3.8.10	