library(Seurat)
library(monocle)
library(ggplot2)
library(dplyr)

RNA.integrated<-readRDS('/fs/home/tangke/mouse_nsc/E13_E17/E13_E17_filter_annotation.rds')
markers <- FindAllMarkers(RNA.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

test=function(RNA.integrated,featureNum,markers_gene){
    data <- as(as.matrix(RNA.integrated@assays$integrated@data), 'sparseMatrix')
    pd <- new('AnnotatedDataFrame', data = RNA.integrated@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fData)
    monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily = negbinomial.size())
    monocle_cds <- estimateSizeFactors(monocle_cds)
    monocle_cds <- estimateDispersions(monocle_cds)

    gene_sle=markers_gene %>%
        dplyr::filter(p_val<0.01) %>%
        group_by(cluster) %>%
        slice_max(n = featureNum, order_by = avg_log2FC) %>%
        pull(gene) %>% unique()
        
    # gene_sle <-VariableFeatures(expression_matrix)
    ordering_genes <-gene_sle
    monocle_cds <-
    setOrderingFilter(monocle_cds,
        ordering_genes = ordering_genes)
    monocle_cds <-
    reduceDimension(monocle_cds,
                    max_components = 2)
    monocle_cds <- orderCells(monocle_cds) ## running step
    saveRDS(monocle_cds,paste0('/fs/home/tangke/mouse_nsc/E13_E17/',featureNum,'_topfeatures_annotation_monocle_trajectory.rds'))
}
 lapply(as.list(50,100,150,200),test,RNA.integrated=RNA.integrated,markers_gene=markers)
