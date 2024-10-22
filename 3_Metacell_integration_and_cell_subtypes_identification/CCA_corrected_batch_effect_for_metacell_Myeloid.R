suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(stringr)
    library(stringi)
    library(MAESTRO)
})

args <- commandArgs(trailingOnly = TRUE)
Seurat_TPM_myeloid<- args[1]
Seurat_TPM_metadata_myeloid<- args[2]
TF_Ligand_Metabolic<- args[3] ## integrated information
integrated_object_path_myeloid<-args[4]

data<-readRDS(Seurat_TPM_myeloid)
data <- CreateSeuratObject(counts = data, project = "myeloid", min.cells = 0, min.features = 0)
metadata<-readRDS(Seurat_TPM_metadata_myeloid)

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
saveRDS(anchors,paste0(integrated_object_path_myeloid,'myeloid_anchors.rds'))
RNA.integrated <- IntegrateData(anchorset = anchors, dims = dims.use)
saveRDS(RNA.integrated,paste0(integrated_object_path_myeloid,'myeloid_integration.rds'))
saveRDS(RNA.integrated@meta.data,paste0(integrated_object_path_myeloid,'myeloid_integration_metadata.rds'))