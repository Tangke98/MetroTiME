#' Corrects the batch effects in datasets for specific cell lineages.
#' @param cell_lineage The cell lineage to analysis.
#' @param Seurat_TPM Path to the Metacell TPM file for cell lineage.
#' @param Seurat_TPM_metadata Path to the Metacell metadata file for cell lineage.
#' @param TF_Ligand_Metabolic Path to the genes file for integration.
#' @param integrated_object_path Path to the output file.
#' @author Ke Tang
#

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(stringr)
    library(stringi)
    library(MAESTRO)
})

args <- commandArgs(trailingOnly = TRUE)
cell_lineage<-args[1]
Seurat_TPM<- args[2]
Seurat_TPM_metadata<- args[3]
TF_Ligand_Metabolic<- args[4] ## integrated information
integrated_object_path<-args[5]

data<-readRDS(paste0(Seurat_TPM,cell_lineage,'/',cell_lineage,"_miniclusteR30_TPM.rds"))
data <- CreateSeuratObject(counts = data, project = cell_lineage, min.cells = 0, min.features = 0)
metadata<-readRDS(paste0(Seurat_TPM_metadata,cell_lineage,'/',cell_lineage,"_meta_info.rds"))

data@meta.data=metadata[,c('orig.ident','nCount_RNA','nFeature_RNA','percent.mito','percent.ercc','Dataset','Cancer_type','Tissue','batch')]
data.list <- SplitObject(data, split.by = "batch") ## split the file with the batch

nfeatures = 10000
dims.use = 1:30

TF_Ligand_Metabolic<-readRDS(TF_Ligand_Metabolic)
for(i in 1:length(data.list)){ ## use the tf ligand and metabolic gene and highly variable genes to integrate data
    data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = (length(nfeatures)-length(TF_Ligand_Metabolic)), verbose = FALSE)
    VariableFeatures(data.list[[i]])<-c(VariableFeatures(data.list[[i]]),TF_Ligand_Metabolic)
}

anchors <- FindIntegrationAnchors(object.list = data.list,
                                  dims = dims.use, 
                                  anchor.features = nfeatures)
saveRDS(anchors,paste0(integrated_object_path,cell_lineage,'/',cell_lineage,'_anchors.rds'))
RNA.integrated <- IntegrateData(anchorset = anchors, dims = dims.use)
saveRDS(RNA.integrated,paste0(integrated_object_path,cell_lineage,'/',cell_lineage,'_integration.rds'))