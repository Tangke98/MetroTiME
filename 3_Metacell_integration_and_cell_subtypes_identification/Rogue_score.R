#' Uses ROGUE to evaluate the clustering of cell lineage based on highly variable genes.
#' @param cell_lineage The cell lineage to analysis.
#' @param integrated_object_path Path to the cell lineage Metacell TPM file.
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
    library(ROGUE)
    library(ComplexHeatmap)
    library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)

cell_lineage<-args[1]
integrated_object_path<-args[2]
output_path<-args[3]

integrated_object_path=paste0(integrated_object_path,cell_lineage,'/',cell_lineage)
output_path=paste0(output_path,cell_lineage,'/',cell_lineage)


cal_purity=function (seurat_rna, resolution, cluster, feature, metabolic_gene) 
{
    exp = seurat_rna@assays$integrated@data
    meta = seurat_rna@meta.data
    meta_use = meta[seurat_rna[[feature]] == cluster, ]
    exp_use = exp[, colnames(exp) %in% rownames(meta_use)]
    exp_use = exp_use[rownames(exp_use) %in% metabolic_gene, 
        ]
    exp.res <- SE_fun(exp_use)
    rogue.value <- CalculateRogue(exp.res, platform = "UMI")
    df = data.frame(cluster = cluster, resolution = resolution, 
        rogue = rogue.value)
    return(df)
}

inte_use<-readRDS(paste0(integrated_object_path,'_integration_annotation_celltype.rds'))
genes=VariableFeatures(inte_use)
if (cell_lineage=='Fibroblast'){
    levels_celltype=c('MyoFibro_MYH11','MyoFibro_RGS5','Fibro_CTHRC1',
                                                    'Fibro_SFRP1','Fibro_IL6',
                                                    'Fibro_CCL5','Fibro_SAA1','Fibro_CDC20')
    color_celltype=c('Fibro_SFRP1'="#97CADC",'Fibro_CCL5'='#B2B31D','Fibro_IL6'="#FBB065",'Fibro_CTHRC1'='#F87379',
                    'Fibro_SAA1'='#92C274','MyoFibro_RGS5'='#984EA3','MyoFibro_MYH11'='#DB5DB8')
}
if (cell_lineage=='Myeloid'){
    levels_celltype=c('Mono_FCN1','Mono_CD16','Macro_C1QC','Macro_SPP1','Mono_THBS1','Macro_SLPI',
                  'Macro_IL32','Macro_CDC20','cDC1_CLEC9A','cDC2_CD1C','pDC_LILRA4','Mast_CPA3')
    color_celltype=c('Mono_CD16'='#92C274','Mono_FCN1'='#1C7E76','Mono_THBS1'='#824DAE','Macro_SPP1'='#DB5DB8',
                    'Macro_SLPI'='#CA131F','Macro_C1QC'='#8ACDD4','Macro_IL1B'='#1965B0','Macro_IL32'='#FBB065',
                    'Macro_CDC20'='#FED43B','cDC2_CD1C'='#B2B31D','cDC1_CLEC9A'='#D2D77A','pDC_LILRA4'='#BBBBBB',
                    'Mast_CPA3'='#EDE2F6','Macro_CLEC10A'='#B8A1CD','Macro_CCL18'='#FD8586')
}

cal_purity_res=do.call(rbind,lapply(as.list(unique(as.character(inte_use$cell_type))),cal_purity,resolution=0.9,seurat_rna=inte_use,feature='cell_type',metabolic_gene=genes))
cal_purity_res$cluster=as.factor(cal_purity_res$cluster)
cal_purity_res$cluster=factor(cal_purity_res$cluster,levels=levels_celltype)
cal_purity_res=cal_purity_res[order(cal_purity_res$cluster),]

width=4
height=4

options(repr.plot.width = width, repr.plot.height = height+1,repr.plot.res = 100)
p=ggplot(cal_purity_res, aes(cluster, rogue)) +
    geom_segment(aes(x=cluster, 
               xend=cluster,
               y=0, 
               yend=rogue))+
    scale_y_continuous(limits = c(0.95,1),breaks = seq(0.95, 1, by = 0.01))+
    theme_bw()+
    theme(panel.grid =element_blank()) + 
    geom_point(shape=21,size=3,colour="black",aes(fill=cluster))+
    labs(x='Cluster',y='ROGUE Score')+
    scale_fill_manual(values=c(color_celltype))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(axis.title=element_text(size=15),axis.text=element_text(size=15),
          legend.title=element_text(size=15),legend.text=element_text(size=15),
          axis.text.x = element_text(angle = 90, hjust = 1))+NoLegend()
print(p)
pdf(paste0(output_path,"_celltype_rigue_score.pdf"),width=width,height=height+1)
    print(p)
dev.off()
