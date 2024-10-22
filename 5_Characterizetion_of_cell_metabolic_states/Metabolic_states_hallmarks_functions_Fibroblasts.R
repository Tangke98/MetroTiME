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

integrated_object_path_fibroblasts<-args[1]
metabolic_gene<-args[2]
hallmark_path<-args[3]
output_path<-args[4]

metadata<-readRDS(paste0(integrated_object_path_fibroblasts,'fibroblasts_integration_annotation_celltype_metatype_metadata.rds'))
fibroblast<-readRDS(paste0(integrated_object_path_fibroblasts,'fibroblasts_integration_annotation_celltype_metatype.rds'))
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

DefaultAssay(fibroblast)='integrated'
metacell=as.matrix(GetAssayData(fibroblast))
cal_MetaModule(metacell,features,output_path,'fibroblast_features_gsva')

### calculate hallmarks scores for each metabolic states
hallmarks=readRDS(hallmark_path)
cal_MetaModule(metacell,hallmarks,output_path,'fibroblast_hallmarks_gsva')

## visualization of the feature scores and hallmarks scores
### for feature scores
output_name='fibroblast_feature_ssgsea_res.pdf'
width=3
height=1.2

feature_score<-readRDS(paste0(output_path,'fibroblast_features_gsva.rds'))
fibroblast@meta.data=fibroblast@meta.data[,1:13]
Idents(fibroblast)='metabolic_type'

object=fibroblast
features = c('ECM','Immune_regulatory','Antigen_presenting')
levels_metabolic_use=levels_metabolic[levels_metabolic %in% unique(object$metabolic_type)]

object@meta.data=cbind(object@meta.data,t(feature_score[,rownames(object@meta.data)]))

p1 <- DotPlot(object, features = c('ECM','Immune_regulatory','Antigen_presenting'),assay='RNA' ) + 
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
output_name='fibroblast_hallmark_ssgsea_res.pdf'
width=4
height=4

hallmark<-readRDS(paste0(output_path,'fibroblast_hallmark_gsva.rds'))
rownames(hallmark)=tolower(rownames(hallmark))
fibroblast@meta.data=fibroblast@meta.data[,1:13]
Idents(fibroblast)='metabolic_type'

object=fibroblast
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

