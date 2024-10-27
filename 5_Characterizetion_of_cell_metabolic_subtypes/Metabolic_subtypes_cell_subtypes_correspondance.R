#' Corresponds metabolic subtypes with cell subtypes in specific cell lineages.
#' @param cell_lineage The cell lineage to analysis.
#' @param integrated_object_path Path to the specific cell lineages Metacell TPM file.
#' @param metabolic_gene Path to the metabolic genes file.
#' @param MetaModule_info Path to the metabolic reaction information file.
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
cell_lineage<-args[1]
integrated_object_path<-args[2]
MetaModule_info<-args[3]
metabolic_gene<-args[4]
output_path<-args[5]

integrated_object_path=paste0(integrated_object_path,cell_lineage,'/',cell_lineage)
output_path=paste0(output_path,cell_lineage,'/',cell_lineage)

metadata<-readRDS(paste0(integrated_object_path,'_integration_annotation_celltype_metatype_metadata.rds'))
seurat<-readRDS(paste0(integrated_object_path,'_integration_annotation_celltype_metatype.rds'))
MetaModule_info=readRDS(MetaModule_info)
metabolic_gene=readRDS(metabolic_gene)

## set the cell color
if (cell_lineage=='Fibroblast'){
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
}
if (cell_lineage=='Myeloid'){
    ## set the cell color
    color_celltype=c('Mono_CD16'='#92C274','Mono_FCN1'='#1C7E76','Mono_THBS1'='#824DAE','Macro_SPP1'='#DB5DB8',
                    'Macro_SLPI'='#CA131F','Macro_C1QC'='#8ACDD4','Macro_IL1B'='#1965B0','Macro_IL32'='#FBB065',
                    'Macro_CDC20'='#FED43B','cDC2_CD1C'='#B2B31D','cDC1_CLEC9A'='#D2D77A','pDC_LILRA4'='#BBBBBB',
                    'Mast_CPA3'='#EDE2F6','Macro_CLEC10A'='#B8A1CD','Macro_CCL18'='#FD8586')
    # color_metabolic=c('PUFA+ Macro'="#DB5DB8",'GLY+ Mono'='#1C7E76','PURINE+ Macro'="#B8A1CD",'GAL+ Macro'='#8ACDD4',
    #                                               'TRP+ cDC2'='#B2B31D','AA+ Macro'='#689ECE','OXP+ Mono'='#92C274','SL+ Mast'='#EDE2F6',
    #                                               'GLY+ Macro'='#FD8586','GLYCAN+ Macro'='#CA131F','CHOL+ Macro'='#FBB065','OXP+ Macro'='#B3B3B3',
    #                                               'NMP+ Macro'='#FED43B','SL+ pDC'='#BBBBBB')

    color_metabolic=c('PUFA_LTM'='#1C7E76','PUFA_EM'='#92C274','OXP'='Thistle','FAO'='#824DAE',
                                                'SLM'='#8ACDD4','ARG'='#DB5DB8','FAS'='#608FBF','GST'='#CA131F',
                                                'PUFA'='#FBB065','AA'='#FED43B')
    # color_metabolic_type=c('Glycolysis'="#E31A1C",'Polysaccharide metabolism'='#F781BF','TCA cycle'="#D88867",
    #                        'Mitochondria'='#F1B05B','Lipid Metabolism'='#57998B','AA Metabolism'='#739FD6',
    #                        'Other AA Metabolism'='#20B2C5','Nucleotide Metabolism'='#BEAED4',
    #                        'Purine Metabolism'='#BEAED4','Pyrimidine metabolism'='#AD94C0',
    #                        'Transporters'='#B3B3B3')
    color_cancer_type=c('BCC'="#4DAF4A",'BRCA'='#984EA3','CESC'="#FFFF33",'CHOL'='#A65628',
                        'CRC'='#F781BF','ESCA'='#999999','HNSC'='#66C2A5','KIRC'='#FC8D62',
                                                    'LIHC'='#8DA0CB','NSCLC'='#E78AC3','OV'='#A6D854','PAAD'='#FFD92F',
                                                    'SKCM'='#E5C494','BLCA'='#B3B3B3','LSCC'='#8DD3C7','OS'='#FFFFB3',
                                                    'SCC'='#BEBADA','SS'='#FB8072','STAD'='#80B1D3','THCA'='#FDB462',
                                                    'PRAD'='#B3DE69','LUAD'='#FCCDE5','MCC'='#D9D9D9','UCEC'='#BC80BD')
    levels_celltype=c('Mono_FCN1','Mono_CD16','Macro_C1QC','Macro_SPP1','Mono_THBS1','Macro_SLPI',
                    'Macro_IL32','Macro_CDC20','cDC1_CLEC9A','cDC2_CD1C','pDC_LILRA4','Mast_CPA3')
    levels_metabolic=c("PUFA_LTM",'PUFA_EM','FAS','SLM','ARG','OXP','FAO','GST','PUFA','AA')
}


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


meta=seurat@meta.data
meta_use=meta[,c('cell_type','metabolic_type')]
Bar_prop(meta_use,'cell_type','metabolic_type',levels_celltype,levels_metabolic,
         color_celltype,'Cell Subtype Correspond to Metabolic Subtype',2,3,output_path,'_celltype')

alluvial_corres(meta_use,levels_celltype,levels_metabolic,
                color_celltype,color_metabolic,3,5,output_path,'_celltype_metabolic')