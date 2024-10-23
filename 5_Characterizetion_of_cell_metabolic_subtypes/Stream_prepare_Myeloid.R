#' Prepares trajectory data for myeloid cells using STREAM.
#' @param integrated_object_path_myeloids Path to the myeloid cell object.
#' @param myeloids_metabolic_score Path to the MetaModule score of myeloid cells.
#' @param output_path Path to the output file.
#' @author Ke Tang
#
suppressMessages({
    library(Seurat)
    library(GSVA)
    data("human_gem")
    library(ggplot2)
    library(dplyr)
    library(stringr)
    library(stringi)
    library(ComplexHeatmap)
    library(BiocParallel)
    library(rPref)
    library(lazyeval)
    library(tidyr)
    library(SCopeLoomR)
    library(DropletUtils)
    library(Seurat)
    ## get the cell ..
    library(SeuratDisk)
    library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

integrated_object_path_myeloids<-args[1]
myeloids_metabolic_score<-args[2]
output_path<-args[3]

myeloid<-readRDS(paste0(integrated_object_path_myeloids,'myeloid_integration_annotation_celltype_subset_metatype.rds'))
sc=CreateSeuratObject(counts = myeloid@assays$integrated@data)
sc@meta.data=myeloid@meta.data
sdata.loom <- as.loom(x = sc, filename = paste0(output_path,"myeloid_integration_subset.loom", verbose = FALSE))
sdata.loom$close_all()

metabolic_gsva=readRDS(paste0(myeloids_metabolic_score,'myeloid_gsva.rds'))
metabolic_gsva.t=t(metabolic_gsva)
write.csv(metabolic_gsva.t,paste0(output_path,'myeloid_integration.gsva.csv'),
          row.names=T,quote=F,col.names=T)

