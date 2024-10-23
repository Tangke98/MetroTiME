#' Prepares trajectory data for fibroblasts using STREAM.
#' @param integrated_object_path_fibroblasts Path to the fibroblasts object.
#' @param fibroblasts_metabolic_score Path to the MetaModule score of fibroblasts.
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

integrated_object_path_fibroblasts<-args[1]
fibroblasts_metabolic_score<-args[2]
output_path<-args[3]

fibroblast<-readRDS(paste0(integrated_object_path_fibroblasts,'fibroblast_integration_annotation_celltype_metatype.rds'))
sc=CreateSeuratObject(counts = fibroblast@assays$integrated@data)
sc@meta.data=fibroblast@meta.data
sdata.loom <- as.loom(x = sc, filename = paste0(output_path,"fibroblast_integration_0618.loom", verbose = FALSE))
sdata.loom$close_all()

metabolic_gsva=readRDS(paste0(fibroblasts_metabolic_score,'fibroblast_gsva.rds'))
metabolic_gsva.t=t(metabolic_gsva)
write.csv(metabolic_gsva.t,paste0(output_path,'fibroblast_integration.gsva.csv'),
          row.names=T,quote=F,col.names=T)

