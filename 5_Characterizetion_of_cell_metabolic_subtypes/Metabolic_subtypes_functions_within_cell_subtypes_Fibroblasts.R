#' Corresponds metabolic subtypes and cell functions within CTHRC1+ eFibro.
#' @param integrated_object_path_fibroblasts Path to the fibroblasts Metacell TPM file.
#' @param output_path Path to the output file.
#' @author Ke Tang
#
suppressPackageStartupMessages({
    library(Seurat)
    library(GSVA)
    library(ggforce)
    library(ggplot2)
    library(ggalluvial)
    library(dplyr)
    library(stringr)
    library(stringi)
    library(ComplexHeatmap)
    library(BiocParallel)
    library(rPref)
    library(lazyeval)
    library(tidyr)
    library(RColorBrewer)
    library(tibble)
    library(reshape2)
    library(ggpubr)
    library(ggridges)
    library(ggrepel)
    suppressMessages(library(ROGUE))
    suppressMessages(library(tidyverse))
    library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

integrated_object_path_fibroblasts<-args[1]
output_path<-args[2]

singlecell_gene_test=function (SerautObj, genes.use, group.by = NULL, assay = "RNA", 
    comp = NULL, alpha_start = 0.05, Bonferroni = T, only_postive = F) 
{
    p_val.out <- c()
    stat.out <- c()
    condition.out <- c()
    gene.out <- c()
    if (only_postive == F) {
        for (gene in genes.use) {
            group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == 
                comp[1], ])
            group1_exp = SerautObj@assays[[assay]]@data[gene, 
                group1_cellname]
            group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == 
                comp[2], ])
            group2_exp = SerautObj@assays[[assay]]@data[gene, 
                group2_cellname]
            t_out = t.test(group1_exp, group2_exp)
            cond = paste(comp[1], comp[2], sep = "_")
            condition.out <- c(condition.out, cond)
            stat.out <- c(stat.out, t_out[["statistic"]])
            p_val.out <- c(p_val.out, t_out[["p.value"]])
            gene.out <- c(gene.out, gene)
        }
    }
    else {
        for (gene in genes.use) {
            group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == 
                comp[1], ])
            group1_exp = SerautObj@assays[[assay]]@data[gene, 
                group1_cellname]
            group1_exp <- group1_exp[which(group1_exp > 0)]
            group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == 
                comp[2], ])
            group2_exp = SerautObj@assays[[assay]]@data[gene, 
                group2_cellname]
            group2_exp <- group2_exp[which(group2_exp > 0)]
            t_out = t.test(group1_exp, group2_exp)
            cond = paste(comp[1], comp[2], sep = "_")
            condition.out <- c(condition.out, cond)
            stat.out <- c(stat.out, t_out[["statistic"]])
            p_val.out <- c(p_val.out, t_out[["p.value"]])
            gene.out <- c(gene.out, gene)
        }
    }
    if (Bonferroni == T) {
        new_alpha = alpha_start/(2 * length(genes.use))
        cat(paste("\n", "P-value for significance: p <", new_alpha, 
            "\n"))
        sig_out = p_val.out < new_alpha
        dfOUT <- data.frame(gene = gene.out, condition = condition.out, 
            p_val = p_val.out, statistic = stat.out, significant = sig_out)
        dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns", ifelse(dfOUT$p_val > 
            0.01, "*", ifelse(dfOUT$p_val > 0.001, "**", "****")))
    }
    else {
        dfOUT <- data.frame(gene = gene.out, condition = condition.out, 
            p_val = p_val.out, statistic = stat.out)
        dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns", ifelse(dfOUT$p_val > 
            0.01, "*", ifelse(dfOUT$p_val > 0.001, "**", "****")))
    }
    return(dfOUT)
}

metadata<-readRDS(paste0(integrated_object_path_fibroblasts,'fibroblast_integration_annotation_celltype_metatype_metadata.rds'))
fibroblast<-readRDS(paste0(integrated_object_path_fibroblasts,'fibroblast_integration_annotation_celltype_metatype.rds'))


## set the cell color
color_celltype=c('Fibro_SFRP1'="#97CADC",'Fibro_CCL5'='#B2B31D','Fibro_IL6'="#FBB065",'Fibro_CTHRC1'='#F87379',
                 'Fibro_SAA1'='#92C274','MyoFibro_RGS5'='#984EA3','MyoFibro_MYH11'='#DB5DB8')
color_metabolic=c('FAO'="#97CADC",'LYS'='#B2B31D','AA'="#92C274",
                                            'GLY'='#FBB2B4','GLYCAN'='#E47B7E',
                                            'PUFA'='#FBB065','IPM'='#984EA3',
                                            'OXP'='#DB5DB8')
color_cancer_type=c('BCC'="#4DAF4A",'BRCA'='#984EA3','CESC'="#FFFF33",'CHOL'='#A65628',
                    'CRC'='#F781BF','ESCA'='#999999','HNSC'='#66C2A5','KIRC'='#FC8D62',
                                                  'LIHC'='#8DA0CB','NSCLC'='#E78AC3','OV'='#A6D854','PAAD'='#FFD92F',
                                                  'SKCM'='#E5C494','BLCA'='#B3B3B3','LSCC'='#8DD3C7','OS'='#FFFFB3',
                                                  'SCC'='#BEBADA','SS'='#FB8072','STAD'='#80B1D3','THCA'='#FDB462',
                                                  'PRAD'='#B3DE69','LUAD'='#FCCDE5','MCC'='#D9D9D9','UCEC'='#BC80BD')
levels_celltype=c('MyoFibro_MYH11','MyoFibro_RGS5','Fibro_CTHRC1',
                                                 'Fibro_SFRP1','Fibro_IL6',
                                                 'Fibro_CCL5','Fibro_SAA1','Fibro_CDC20')
levels_metabolic=c('OXP','IPM','GLYCAN','GLY','FAO','PUFA','LYS','AA')

fibroblasts_score<-readRDS(fibroblasts_score)

width=3
height=3
output_name='fibroblast_metabolic_ACTA2.pdf'

object=subset(fibroblast,subset=metabolic_type%in% c('GLY','GLYCAN'))
A <- singlecell_gene_test(object, 
                    genes.use = c('ACTA2'),
                    group.by = 'metabolic_type', 
                    comp = c("GLY", "GLYCAN"))
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig
## draw the violin plot
options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=VlnPlot(object, features = 'ACTA2',pt.size = 0,cols =c('GLY'='#FBB2B4','GLYCAN'='#E47B7E'))+
    geom_signif()+theme_bw() + 
    theme(panel.grid =element_blank()) + 
    theme(axis.text = element_blank()) + 
    theme(axis.title=element_text(size=8),axis.text=element_text(size=8),
    legend.text=element_text(size=8),legend.title=element_text(size=8))+
    geom_boxplot(width=.2,col="black",fill="white")+
    geom_signif(annotations = anno_sig,
                    y_position = 1.3,
                    xmin = 1,
                    xmax = 2,
                    tip_length = 0)
print(p) 
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p) 
dev.off()





