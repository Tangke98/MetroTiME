#' Associates function score with the MetaRegulon score in fibroblasts.
#' @param integrated_object_path_fibroblasts Path to the fibroblasts Metacell TPM file.
#' @param fibroblasts_metabolic_score Path to the metabolic score file.
#' @param fibroblasts_feature_score Path to the function score file.
#' @param fibroblasts_DEM Path to the differentially enriched MetaModule file.
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
    library(pheatmap)
    suppressMessages(library(ROGUE))
    suppressMessages(library(tidyverse))
})

args <- commandArgs(trailingOnly = TRUE)

integrated_object_path_fibroblasts<-args[1]
fibroblasts_metabolic_score<-args[2]
fibroblasts_feature_score<-args[3]
fibroblasts_DEM<-args[4]
MetaModule_info<-args[5]
output_path<-args[6]

get_cor=function(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,color){
    data_use=cor[,c(feature,'class')]
    data_use=data_use[order(data_use[,1],decreasing = TRUE),]
    data_use$order=1:nrow(data_use)
    module_use=unique(DEM[DEM$cluster %in% metabolic_type,'gene'])
    data_use$label=ifelse(rownames(data_use) %in% module_use,'1','0' )
    data_use$text = ifelse(rownames(data_use) %in% module_label, 
                           paste0(rownames(data_use), "|", MM.meta.nozero[rownames(data_use),'EQUATION']),'')
    colnames(data_use)[1]='feature'
    
    options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
    p=ggplot(data_use,aes(x=order,y=feature,fill=label))+geom_bar(stat = "identity")+
        theme_bw() + 
        theme(panel.grid =element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        scale_fill_manual(values=c('1'=color,'0'="grey"))+
        labs(x='MetaModule',y=paste0("Association with ",feature),title=)+
        theme(axis.title=element_text(size=12),
                axis.text=element_text(size=12),
                legend.text=element_text(size=12),
                legend.title=element_text(size=12)
                )+
        ggrepel::geom_text_repel(aes(label=text),max.overlaps = Inf)
    print(p)
    
    pdf(paste0(output_path,output_name),width=width,height=height)
        print(p) 
    dev.off()
}

metabolic_score=t(readRDS(fibroblasts_metabolic_score))
feature_score=t(readRDS(fibroblasts_feature_score))

cor=as.data.frame(cor(metabolic_score,feature_score))
MM.meta.nozero=readRDS(MetaModule_info) %>% as.data.frame()
cor$class=MM.meta.nozero[rownames(cor),'SUBSYSTEM']
DEM=readRDS(fibroblasts_DEM)

## For PUFA metabolic states
feature='Immune_regulatory'
metabolic_type=c('PUFA')
# metabolic_type=c('SLM')
# module_label=gsub("-","_",c('RE1514M','RE1514X','r1251','PROAKGOX1r','HMR-4776','HMR-4077','HMR-3838','HMR-6958','HMR-3816'))
# module_label=gsub("-","_",c(DEM[DEM$cluster %in% metabolic_type & DEM$class%in% c('Prostaglandin biosynthesis','Eicosanoid metabolism',
#                                                                                 'Leukotriene metabolism','Arachidonic acid metabolism',
#                                                                                 'Linoleate metabolism'),'gene']))
module_label=c('HMR-1330', ## prostaglandinD2
'HMR-0960', ## 5(S)-HETE
'HMR-0940'
# 'HMR_1129',
# 'HMR_1149'                          
                           )     ## 20-oxo-LTB4)
output_name='fibroblast_pufa_immune_regulatory.pdf'
width=5
height=5

color='#FBB065'
get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,'#FBB065')

## For GLYCAN metabolic states
feature='ECM'
metabolic_type=c('GLYCAN')
# metabolic_type=c('SLM')
# module_label=gsub("-","_",DEM[DEM$cluster %in% metabolic_type & DEM$length>1 & DEM$class%in% c('Chondroitin / heparan sulfate biosynthesis','Chondroitin sulfate degradation',
#                                                                                                'Heparan sulfate degradation','Keratan sulfate degradation',
#                                                                                               'Keratan sulfate biosynthesis'),'gene'])
module_label=c('HMR-7493', ## Chondroitin sulfate
'HMR-7494', ## Chondroitin sulfate
'HMR-3833') ## HMR_3833

output_name='fibroblast_glycan_cor_ecm.pdf'
width=5
height=5
get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,"#E47B7E")

## For LYS metabolic states
feature='Antigen_presenting'
metabolic_type=c('LYS')
# metabolic_type=c('SLM')
# module_label=gsub("-","_",DEM[DEM$cluster %in% metabolic_type & DEM$length>1 & DEM$class%in% c('Lysine metabolism'),'gene'])
module_label=c('HMR-6555','HMR-8026') ## HMR_3833

output_name='fibroblast_lys_cor_antigen_presenting.pdf'
width=5
height=5
get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,"#B2B31D")

