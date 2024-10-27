#' Prepares trajectory data for specific cell lineages using STREAM.
#' @param cell_lineage The cell lineage to analysis.
#' @param integrated_object_path Path to the specific cell lineages object.
#' @param metabolic_score Path to the MetaModule score of specific cell lineages.
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
cell_lineage<-args[1]
integrated_object_path<-args[2]
metabolic_score<-args[3]
output_path<-args[4]

integrated_object_path=paste0(integrated_object_path,cell_lineage,'/',cell_lineage)
output_path=paste0(output_path,cell_lineage,'/',cell_lineage)
metabolic_score=paste0(metabolic_score,cell_lineage,'/',cell_lineage)

seurat<-readRDS(paste0(integrated_object_path,'_integration_annotation_celltype_metatype.rds'))
sc=CreateSeuratObject(counts = seurat@assays$integrated@data)
sc@meta.data=seurat@meta.data
sdata.loom <- as.loom(x = sc, filename = paste0(output_path,"_integration.loom"), verbose = FALSE)
sdata.loom$close_all()

metabolic_gsva=readRDS(paste0(metabolic_score,'_gsva.rds'))
metabolic_gsva.t=t(metabolic_gsva)
write.csv(metabolic_gsva.t,paste0(output_path,'_integration.gsva.csv'),
          row.names=T,quote=F,col.names=T)

