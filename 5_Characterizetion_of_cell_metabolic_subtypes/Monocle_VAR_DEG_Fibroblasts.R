#' Analyzes the trajectory of fibroblasts.
#' @param integrated_object_path_fibroblasts Path to the fibroblasts Metacell TPM file.
#' @param fibroblasts_markers Path to the fibroblasts markers file.
#' @param output_path Path to the output file.
#' @author Ke Tang
#
suppressPackageStartupMessages({
    library(monocle)
    library(dplyr)
    library(tibble)
    library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

integrated_object_path_fibroblasts<-args[1]
fibroblasts_markers<-args[2]
output_path<-args[3]

fibroblast<-readRDS(integrated_object_path_fibroblasts,'fibroblast_integration_annotation_celltype_metatype.rds')

## highly variable genes
object=fibroblast
data <- as(as.matrix(object@assays$integrated@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = object@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily = negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

get_variable=function(object,gene_number,monocle_cds,output_path){
    test <- FindVariableFeatures(object, selection.method = "vst", nfeatures = gene_number) #
    ordering_genes <-VariableFeatures(test)
    monocle_cds <-
        setOrderingFilter(monocle_cds,
            ordering_genes = ordering_genes)
    monocle_cds <-
        reduceDimension(monocle_cds, 
                        max_components = 2)
    monocle_cds <- orderCells(monocle_cds) ## running step
    options(repr.plot.width =8, repr.plot.height = 6,repr.plot.res = 100)
    saveRDS(monocle_cds,paste0(output_path,'fibroblast_variable_',gene_number,'.rds'))
    
    p= plot_cell_trajectory(monocle_cds, color_by = "cell_type",cell_size = 0.5)+
    scale_color_manual(values=c('Fibro_SFRP1'="#97CADC",'Fibro_CCL5'='#B2B31D','Fibro_IL6'="#FBB065",'Fibro_CTHRC1'='#F87379',
                 'Fibro_GPX1'='#92C274','MyoFibro_RGS5'='#DB5DB8','MyoFibro_MYH11'='#984EA3','Fibro_CDC20'='#FED43B'))+
    theme(axis.title=element_text(size=25),axis.text=element_text(size=25),
                      legend.text=element_text(size=20),legend.title=element_text(size=25),
         legend.position = "right")+
    labs(title = 'Trajectory of fibroblast')
    pdf(paste0(output_path,'fibroblast_variable_',gene_number,'.pdf'),width=9,height=4)
        print(p) 
    dev.off()
}
gene_num=c(1500,2000,2500,3000,3500,4000,4500,5000)
for (num in gene_num){
    get_variable(object=object,gene_number=num,monocle_cds=monocle_cds,output_path=output_path)
}

## differentially expressed genes
object=fibroblast
fibroblast.markers<-readRDS(fibroblasts_markers)
data <- as(as.matrix(object@assays$integrated@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = object@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily = negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
markers.use=fibroblast.markers

get_variable=function(markers.use,gene_number,monocle_cds,output_path){
    markers.use %>%
    group_by(cluster) %>%
    slice_head(n = gene_number) %>%
    ungroup() -> top
    ordering_genes <-unique(top$'gene')
    
    monocle_cds <-
        setOrderingFilter(monocle_cds,
            ordering_genes = ordering_genes)
    monocle_cds <-
        reduceDimension(monocle_cds, 
                        max_components = 2)
    monocle_cds <- orderCells(monocle_cds) ## running step
    options(repr.plot.width =8, repr.plot.height = 6,repr.plot.res = 100)
    saveRDS(monocle_cds,
           paste0(output_path,'fibroblast_DEG_',gene_number,'.rds'))
    
    p= plot_cell_trajectory(monocle_cds, color_by = "cell_type",cell_size = 0.5)+
    scale_color_manual(values=c('Fibro_SFRP1'="#97CADC",'Fibro_CCL5'='#B2B31D','Fibro_IL6'="#FBB065",'Fibro_CTHRC1'='#F87379',
                 'Fibro_SAA1'='#92C274','MyoFibro_RGS5'='#984EA3','MyoFibro_MYH11'='#DB5DB8'))+
    theme(axis.title=element_text(size=25),axis.text=element_text(size=25),
                      legend.text=element_text(size=20),legend.title=element_text(size=25),
         legend.position = "right")+
    labs(title = 'Trajectory of fibroblast')
    pdf(paste0(output_path,'fibroblast_DEG_',gene_number,'.pdf'),width=9,height=4)
        print(p) 
    dev.off()
}
gene_num=c(50,100,150,200,250,300,350)
for (num in gene_num){
    get_variable(markers.use=fibroblast.markers,gene_number=num,monocle_cds=monocle_cds)
}





