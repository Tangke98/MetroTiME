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
MetaModule_info<-args[2]
metabolic_gene<-args[3]
output_path<-args[4]

metadata<-readRDS(paste0(integrated_object_path_fibroblasts,'fibroblasts_integration_annotation_celltype_metatype_metadata.rds'))
fibroblast<-readRDS(paste0(integrated_object_path_fibroblasts,'fibroblasts_integration_annotation_celltype_metatype.rds'))
MetaModule_info=readRDS(MetaModule_info)
metabolic_gene=readRDS(metabolic_gene)

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

## cell type proportion in metabolic states
Bar_prop=function (object, object_column1, object_column2, levels1, levels2, 
    color, title, width, height, output_path, output_name) 
{
    options(repr.plot.width = width + 2, repr.plot.height = height, 
        repr.plot.res = 100)
    prop <- prop.table(table(object[, object_column1], object[, 
        object_column2]), margin = 2)
    prop <- as.data.frame(prop)
    colnames(prop) <- c("cell_subtype", "metabolic_subtype", 
        "Freq")
    prop$cell_subtype = factor(prop$cell_subtype, levels = levels1)
    prop = prop[order(prop$cell_subtype), ]
    prop$metabolic_subtype = factor(prop$metabolic_subtype, levels = levels2)
    prop = prop[order(prop$metabolic_subtype), ]
    p = ggplot(prop, aes(x = metabolic_subtype, y = Freq, group = cell_subtype, 
        fill = cell_subtype)) + stat_summary(geom = "line", fun = "mean", 
        cex = 1, col = "white") + geom_bar(position = "fill", 
        stat = "identity") + scale_fill_manual(values = color_celltype) + 
        labs(x = NULL, y = NULL) + labs(title = title) + theme_bw() + 
        theme(plot.title = element_text(size = 8), axis.title.x = element_text(size = 8), 
            axis.title.y = element_text(size = 8), axis.text = element_text(size = 8), 
            axis.text.x = element_text(angle = 90, hjust = -0.5, 
                vjust = 0.5))
    print(p)
    pdf(paste0(output_path, output_name, ".prop.pdf"), width = width, 
        height = height)
    print(p)
    dev.off()
    pdf(paste0(output_path, output_name, ".prop.nolegend.pdf"), 
        width = width, height = height)
    print(p + NoLegend())
    dev.off()
}

alluvial_corres=function (object, levels1, levels2, color1, color2, width, height, 
    output_path, output_name) 
{
    options(repr.plot.width = width, repr.plot.height = height)
    object_use = to_lodes_form(as.data.frame(object), axes = 1:ncol(as.data.frame(object)), 
        id = "Cohort")
    object_use$stratum <- factor(object_use$stratum, levels = c(levels1, 
        levels2))
    p = ggplot(object_use, aes(x = x, stratum = stratum, alluvium = Cohort, 
        fill = stratum, label = stratum)) + scale_x_discrete(expand = c(0, 
        0)) + geom_flow(width = 1/15) + geom_stratum(alpha = 0.9, 
        width = 1/20) + geom_text(stat = "stratum", size = 0, 
        color = "black") + scale_fill_manual(values = c(color1, 
        color2)) + xlab("") + ylab("") + theme_bw() + theme(panel.grid = element_blank()) + 
        theme(panel.border = element_blank()) + theme(axis.line = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank()) + 
        theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6), 
            legend.text = element_text(size = 6), legend.title = element_text(size = 6))
    print(p)
    pdf(paste0(output_path, output_name, ".correspond.pdf"), 
        width = width, height = height)
    print(p)
    dev.off()
}


meta=fibroblast@meta.data
meta_use=meta[,c('cell_type','metabolic_type')]
Bar_prop(meta_use,'cell_type','metabolic_type',levels_celltype,levels_metabolic,
         color_celltype,'Cell Subtype Correspond to Metabolic Subtype',2,3,output_path,'celltype')

alluvial_corres(meta_use,levels_celltype,levels_metabolic,
                color_celltype,color_metabolic,3,5,output_path,'celltype_metabolic')