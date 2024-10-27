#' Uses Clustree to evaluate the clustering of cell lineages based on highly variable genes.
#' @param cell_lineage The cell lineage to analysis.
#' @param integrated_object_path Path to the cell lineages Metacell TPM file.
#' @param output_path Path to the output file.
#' @author Ke Tang
# 
suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(stringr)
    library(stringi)
    library(MAESTRO)
    library(cluster)
    library(factoextra)
    library(clustree)
    library(tidyr)
    library(ComplexHeatmap)
    library(RColorBrewer)
})
args <- commandArgs(trailingOnly = TRUE)

cell_lineage<-args[1]
integrated_object_path<-args[2]
output_path<-args[3]


draw_cluster=function (data, cell_lineage,output_path, output_name, width, height) 
{
    clustree_plt <- clustree::clustree(data, prefix = paste0(DefaultAssay(data), 
        "_snn_res."))
    pdf(paste0(output_path, cell_lineage,"/",output_name, ".clutree.pdf"), width = width,height = height)
        print(clustree_plt)
    dev.off()
}
multi_res=function (resolution, inte_data) 
{
    for (res in resolution) { ## by the silhouette score, we have idenified the best dims and k.param to cluster the metacells
        inte_data <- RunUMAP(inte_data, dims = 1:40)
        inte_data <- FindNeighbors(inte_data, dims = 1:40, k.param = 80)
        inte_data <- FindClusters(inte_data, resolution = res)
    }
    return(inte_data)
}

inte_use<-readRDS(paste0(integrated_object_path,cell_lineage,"/",cell_lineage,'_integration.rds'))
inte_use <- FindVariableFeatures(inte_use, selection.method = "vst", nfeatures = 3000)
inte_use <- ScaleData(inte_use)
inte_use <- RunPCA(inte_use)
resolution=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
data=multi_res(resolution,inte_use) ## run multiple resolution

output_name=paste0(cell_lineage,'_celltype_pc40_cutree')
width=8
height=6
draw_cluster(data,output_path,cell_lineage,"/", output_name,width,height) ## double check the resolution to cluster the metacells
