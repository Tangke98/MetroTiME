#' Analyzes the trajectory of myeloid cells.
#' @param integrated_object_path_myeloid Path to the myeloid cells Metacell TPM file.
#' @param myeloids_markers Path to the myeloid cell markers file.
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

integrated_object_path_myeloids<-args[1]
myeloids_markers<-args[2]
output_path<-args[3]

myeloid<-readRDS(integrated_object_path_myeloids,'myeloid_integration_annotation_celltype_subset_metatype.rds')

## highly variable genes
object=myeloid
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
    saveRDS(monocle_cds,paste0(output_path,'myeloid_variable_',gene_number,'.rds'))
    
    p= plot_cell_trajectory(monocle_cds, color_by = "cell_type",cell_size = 0.5)+
    scale_color_manual(values=c('Mono_CD16'='#92C274','Mono_FCN1'='#1C7E76','Mono_THBS1'='#824DAE','Macro_SPP1'='#DB5DB8',
                 'Macro_SLPI'='#CA131F','Macro_C1QC'='#8ACDD4','Macro_IL1B'='#1965B0','Macro_IL32'='#FBB065',
                 'Macro_CDC20'='#FED43B','cDC2_CD1C'='#B2B31D','cDC1_CLEC9A'='#D2D77A','pDC_LILRA4'='#BBBBBB',
                 'Mast_CPA3'='#EDE2F6','Macro_CLEC10A'='#B8A1CD','Macro_CCL18'='#FD8586',
                 'cDC1_LAMP3'='#C9D79F'))+
    theme(axis.title=element_text(size=25),axis.text=element_text(size=25),
                      legend.text=element_text(size=20),legend.title=element_text(size=25),
         legend.position = "right")+
    labs(title = 'Trajectory of myeloid')
    pdf(paste0(output_path,'myeloid_variable_',gene_number,'.pdf'),width=9,height=4)
        print(p) 
    dev.off()
}
gene_num=c(1500,2000,2500,3000,3500,4000,4500,5000)
for (num in gene_num){
    get_variable(object=object,gene_number=num,monocle_cds=monocle_cds,output_path=output_path)
}

## differentially expressed genes
object=myeloid
myeloid.markers<-readRDS(myeloids_markers)
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
markers.use=myeloid.markers

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
           paste0(output_path,'myeloid_DEG_',gene_number,'.rds'))
    
    p= plot_cell_trajectory(monocle_cds, color_by = "cell_type",cell_size = 0.5)+
    scale_color_manual(values=c('Mono_CD16'='#92C274','Mono_FCN1'='#1C7E76','Mono_THBS1'='#824DAE','Macro_SPP1'='#DB5DB8',
                 'Macro_SLPI'='#CA131F','Macro_C1QC'='#8ACDD4','Macro_IL1B'='#1965B0','Macro_IL32'='#FBB065',
                 'Macro_CDC20'='#FED43B','cDC2_CD1C'='#B2B31D','cDC1_CLEC9A'='#D2D77A','pDC_LILRA4'='#BBBBBB',
                 'Mast_CPA3'='#EDE2F6','Macro_CLEC10A'='#B8A1CD','Macro_CCL18'='#FD8586',
                 'cDC1_LAMP3'='#C9D79F'))+
    theme(axis.title=element_text(size=25),axis.text=element_text(size=25),
                      legend.text=element_text(size=20),legend.title=element_text(size=25),
         legend.position = "right")+
    labs(title = 'Trajectory of myeloid')
    pdf(paste0(output_path,'myeloid_DEG_',gene_number,'.pdf'),width=9,height=4)
        print(p) 
    dev.off()
}
gene_num=c(50,100,150,200,250,300,350)
for (num in gene_num){
    get_variable(markers.use=myeloid.markers,gene_number=num,monocle_cds=monocle_cds)
}





