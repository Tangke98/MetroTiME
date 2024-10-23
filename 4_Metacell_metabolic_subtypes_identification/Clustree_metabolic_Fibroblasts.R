# Uses Clustree to evaluate the clustering of fibroblasts based on metabolic genes.
# @param integrated_object_path_fibroblasts Path to the fibroblasts Metacell TPM file.
# @param metabolic_gene Path to the metabolic gene file.
# @param output_path Path to the output file.
# @author Ke Tang
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

integrated_object_path_fibroblasts<-args[1]
metabolic_gene<-args[2]
output_path<-args[3]

draw_cluster=function (data, output_path, output_name, width, height) 
{
    clustree_plt <- clustree::clustree(data, prefix = paste0(DefaultAssay(data), 
        "_snn_res."))
    pdf(paste0(output_path, output_name, ".clutree.pdf"), width = width, 
        height = height)
    print(clustree_plt)
    dev.off()
}
multi_res=function (resolution, inte_data) 
{
    for (res in resolution) { ## by the silhouette score, we have idenified the best dims and k.param to cluster the metacells
        inte_data <- RunUMAP(inte_data, dims = 1:45)
        inte_data <- FindNeighbors(inte_data, dims = 1:45, k.param = 70)
        inte_data <- FindClusters(inte_data, resolution = res)
    }
    return(inte_data)
}

inte_use<-readRDS(paste0(integrated_object_path_fibroblasts,'fibroblasts_integration.rds'))
metabolic=readRDS(metabolic_gene)
VariableFeatures(inte_use)=metabolic
inte_use <- ScaleData(inte_use)
inte_use <- RunPCA(inte_use)
resolution=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
data=multi_res(resolution,inte_use) ## run multiple resolution

output_name='fibroblasts_metabolic_pc45_cutree'
width=8
height=6
draw_cluster(data,output_path,output_name,width,height) ## double check the resolution to cluster the metacells

