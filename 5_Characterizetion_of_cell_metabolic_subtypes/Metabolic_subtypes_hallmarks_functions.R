#' Corresponds metabolic subtypes and cell functions in specific cell lineages.
#' @param cell_lineage The cell lineage to analysis.
#' @param integrated_object_path Path to the specific cell lineages Metacell TPM file.
#' @param metabolic_gene Path to the metabolic genes file.
#' @param hallmark_path Path to the hallmarks file.
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
     library(MetroSCREEN)
})

args <- commandArgs(trailingOnly = TRUE)
cell_lineage<-args[1]
integrated_object_path<-args[2]
metabolic_gene<-args[3]
hallmark_path<-args[4]
output_path<-args[5]

integrated_object_path=paste0(integrated_object_path,cell_lineage,'/',cell_lineage)
output_path=paste0(output_path,cell_lineage,'/',cell_lineage)

metadata<-readRDS(paste0(integrated_object_path,'_integration_annotation_celltype_metatype_metadata.rds'))
seurat<-readRDS(paste0(integrated_object_path,'_integration_annotation_celltype_metatype.rds'))
metabolic_gene=readRDS(metabolic_gene)

if (cell_lineage=='Fibroblast'){
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

    ## hallmarks and cell features score calculation
    ### calculate feature scores for each metabolic states

    ECM_CAFs=c("ACTA2","TAGLN","BGN","COL8A1","COL15A1","IGFBP7","TPM1","TPM2","COL10A1","POSTN",
                        "MYL9,COL13A1","COL14A1","ACTA2","TAGLN","MYH11","MYLK","ACTG2","POSTN","FN1",
                        "LUM","DCN","VCAN","COL5A1","COL5A2","COL63A","ACTA2","TAGLN and PDGFA","MMP2",
                        "DCN","COL1A2","RGS5","MYL9","MYH11") ###collagens and contractile protein

    Immun_regulatory_CAFs=c("C5","C5AR1","TGFB1","TGFBR1","TGFB2","TGFBR1","CXCL12","CXCR4","CXCL13","CXCR5","C7","CFB",
                                        "CFH","CFI","BCAM","F11R","IRF5","CCL2","IL6","CXCL2","CXCL1","CXCL3",
                                    "VIM","FAP","COL3A1","DES","IL6","CXCL12","CXCL1","IGF1","FIGF","PDGFD","CXCL12","CXCL13",
                                        "PDPN","CXCL12","CXCL14","PDGFRA","CXCL12","CXCL14","CXCL1","CXCL2","PDGFRA")###marked by inflammatory and immune-regulatory genes
    Antigen_presenting_CAF=c("CD74","HLA-DQA1","CD83","HLA-DRA","HLA-DPA1","HLA-DQA","CD74","HLA-DRA","HLA-DRB1","H2-Q4") ###MHC class II-related genes
    features=list(ECM_CAFs,Immun_regulatory_CAFs,Antigen_presenting_CAF)
    names(features)=c('ECM','Immune_regulatory','Antigen_presenting')
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
    # levels_metabolic_type=c('Glycolysis', 'Polysaccharide metabolism',
    #                          'TCA cycle','Mitochondria',
    #                          'Lipid Metabolism','AA Metabolism','Other AA Metabolism',
    #                          'Purine Metabolism','Pyrimidine metabolism','Transporters')

    ## hallmarks and cell features score calculation
    ### calculate feature scores for each metabolic states

    Anti_inflammatory=c('IL1RN','IL10','IL4','IL11','IL13','TGFB1','TNFRSF1A','TNFRSF1B','IL1R2','IL18BP')
    Pro_inflammatory=c('IL1B','TNF','CCL2','CCL3','CCL5','CCL7','CCL8','CCL13','CCL17','CCL22','IL6','IL8')
    M1=c('IL23','TNF','CXCL9','CXCL10','CXCL11','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40','IDO1','KYNU','CCR7')
    M2=c('IL4R','CCL4','CCL13','CCL20','CCL17','CCL18','CCL22','CCL24','LYVE1','VEGFA',
                        'VEGFB','VEGFC','VEGFD','EGF','CTSA','CTSB','CTSC','CTSD','TGFB1','TGFB2','TGFB3',
                        'MMP14','MMP19','MMP9','CLEC7A','WNT7B','FASL','TNFSF12','TNFSF8','CD276','VTCN1',
                        'MSR1','FN1','IRF4')
    Antigen_Presenting=c("HLA-DQA1","CD83","HLA-DRA","HLA-DPA1","HLA-DQA","CD74","HLA-DRA","HLA-DRB1","H2-Q4")

    Angiogenesis=c('CCND2','CCNE1','CD44','CXCR4','E2F3','EDN1','EZH2','FGF18','FGFR1','FYN','HEY1',
                        'ITGAV','JAG1','JAG2','MMP9','NOTCH1','PDGFA','PTK2','SPP1','STC1','TNFAIP6',
                        'TYMP','VAV2','VCAN','VEGFA')
    Phagocytosis=c('MRC1','CD163','MERTK','C1QB','FCGR1A','FCGR2A','FCGR3A','CD68','TIMD4','SIRPA')
    ECM=c('ACTA2','TAGLN','BGN','COL8A1','COL15A1','IGFBP7','TPM1','TPM2','COL10A1','POSTN','MYL9',
            'COL13A1','COL14A1','MYH11','MYLK','ACTG2','FN1','LUM','DCN','VCAN','COL5A1','COL5A2',
            'COL63A','PDGFA','MMP2','COL1A2','RGS5')
    MDSC=c("CST3","LYZ","CD68","THBS1","S100A8","S100A9","TREM1")
    cell_cycle=c('CDK1', 'CDK5', 'CDC20', 'CCNA2', 'CCNB1','CCNB2')
    features=list(Anti_inflammatory,Pro_inflammatory,M1,M2,Antigen_Presenting,Angiogenesis,Phagocytosis,ECM,
            MDSC,cell_cycle,classical_mono,non_classical_monoc,intermediate_mono)
    names(features)=c('Anti_inflammatory','Pro_inflammatory','M1','M2','Antigen_Presenting','Angiogenesis',
                    'Phagocytosis','ECM', 'MDSC','cell_cycle')
}

DefaultAssay(seurat)='integrated'
metacell=as.matrix(GetAssayData(seurat))
cal_MetaModule(metacell,features,output_path,'_features_gsva')

### calculate hallmarks scores for each metabolic states
hallmarks=readRDS(hallmark_path)
cal_MetaModule(metacell,hallmarks,output_path,'_hallmarks_gsva')

## visualization of the feature scores and hallmarks scores
### for feature scores
output_name='_feature_ssgsea_res.pdf'
width=3
height=1.2

feature_score<-readRDS(paste0(output_path,'_features_gsva.rds'))
seurat@meta.data=seurat@meta.data[,1:13]
Idents(seurat)='metabolic_type'

object=seurat
if (cell_lineage=='Fibroblast'){
    features = c('ECM','Immune_regulatory','Antigen_presenting')
}
if (cell_lineage=='Myeloid'){
    features = c('Anti_inflammatory','Pro_inflammatory','M1','M2','Antigen_Presenting','Angiogenesis','Phagocytosis','ECM','cell_cycle')
}

levels_metabolic_use=levels_metabolic[levels_metabolic %in% unique(object$metabolic_type)]

object@meta.data=cbind(object@meta.data,t(feature_score[,rownames(object@meta.data)]))

p1 <- DotPlot(object, features = features,assay='RNA' ) + 
    coord_flip() + 
    theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
    labs(x=NULL,y=NULL) + 
    guides(size = guide_legend("Percent Expression") )+ 
    scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 
df<- p1$data
exp_mat<-df %>% 
    dplyr::select(-pct.exp, -avg.exp) %>%  
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
    as.data.frame()     
row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

# brewer_colors <- rev(c(brewer.pal(11, "PiYG")[2:4],'#F7F7F7',brewer.pal(11, "BrBG")[6:8]))
brewer_colors <- rev(c(brewer.pal(11, "PiYG")[2:4],'white','white','white'))
color_f <- colorRampPalette(brewer_colors)
colors_100 <- color_f(100)

options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=pheatmap(exp_mat[,levels_metabolic_use],fontsize = 8,treeheight_row =0,treeheight_col =0,color=colors_100,
        angle_col=90,cluster_rows=FALSE,cluster_cols=FALSE)

print(p)  
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p) 
dev.off()

### for hallmarks score
output_name='_hallmark_ssgsea_res.pdf'
width=4
height=4

hallmark<-readRDS(paste0(output_path,'_hallmark_gsva.rds'))
rownames(hallmark)=tolower(rownames(hallmark))
seurat@meta.data=seurat@meta.data[,1:13]
Idents(seurat)='metabolic_type'

object=seurat
object@meta.data=cbind(object@meta.data,t(hallmark[,rownames(object@meta.data)]))
feature=c('angiogenesis','apoptosis','hypoxia','il2_stat5_signaling','il6_jak_stat3_signaling','inflammatory_response','interferon_alpha_response','interferon_gamma_response','myc_targets_v1','myc_targets_v2','p53_pathway','tgf_beta_signaling','tnfa_signaling_via_nfkb','kras_signaling_up',
    'oxidative_phosphorylation','g2m_checkpoint','e2f_targets','epithelial_mesenchymal_transition','myogenesis','coagulation','glycolysis','bile_acid_metabolism','cholesterol_homeostasis')
levels_metabolic_use=levels_metabolic[levels_metabolic %in% unique(object$metabolic_type)]


p1 <- DotPlot(object, features = feature,assay='RNA' ) + 
    coord_flip() + 
    theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
    labs(x=NULL,y=NULL) + 
    guides(size = guide_legend("Percent Expression") )+ #legend
    scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 
df<- p1$data
exp_mat<-df %>% 
    dplyr::select(-pct.exp, -avg.exp) %>%  
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
    as.data.frame() 
row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

# brewer_colors <- rev(c(brewer.pal(11, "PiYG")[2:4],'#F7F7F7',brewer.pal(11, "BrBG")[6:8]))
brewer_colors <- rev(c(brewer.pal(11, "RdYlBu")[2:4],'white','white','white','white'))
color_f <- colorRampPalette(brewer_colors)
colors_100 <- color_f(100)

options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=pheatmap(exp_mat[,levels_metabolic_use],fontsize = 8,treeheight_row =0,treeheight_col =0,color=colors_100,
        angle_col=90,cluster_rows=FALSE,cluster_cols=FALSE)

print(p)  
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p) 
dev.off()

