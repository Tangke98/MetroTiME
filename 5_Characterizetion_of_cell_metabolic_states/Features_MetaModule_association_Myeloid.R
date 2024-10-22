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

integrated_object_path_myeloid<-args[1]
myeloid_metabolic_score<-args[2]
myeloid_feature_score<-args[3]
myeloid_DEM<-args[4]
MetaModule_info<-args[5]
output_path<-args[6]

get_cor=function(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,color){
    data_use=cor[,c(feature,'class')]
    data_use=data_use[order(data_use[,1],decreasing = TRUE),]
    data_use$order=1:nrow(data_use)
    module_use=unique(DEM[DEM$cluster %in% metabolic_type & DEM$length>1,'gene'])
    
    data_use$label=ifelse(rownames(data_use) %in% module_use,'1','0' )
    data_use$text=ifelse(rownames(data_use) %in% module_label,paste0(rownames(data_use),"|",metabolic_metadata[rownames(data_use),'consume'],"->",metabolic_metadata[rownames(data_use),'produce']),'')
    colnames(data_use)[1]='feature'
    options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
    p=ggplot(data_use,aes(x=order,y=feature,fill=label))+geom_bar(stat = "identity")+
        theme_bw() + 
        theme(panel.grid =element_blank()) + #去除网格线
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

metabolic_score=t(readRDS(myeloid_metabolic_score))
feature_score=t(readRDS(myeloid_feature_score))

cor=as.data.frame(cor(metabolic_score,feature_score))

DEM=readRDS(myeloid_DEM)

## For ARG and FAO metabolic states
feature='Angiogenesis'
metabolic_type=c('ARG','FAO')
# metabolic_type=c('SLM')
# module_label=gsub("-","_",c('RE1514M','RE1514X','r1251','PROAKGOX1r','HMR-4776','HMR-4077','HMR-3838','HMR-6958','HMR-3816'))
# module_label=gsub("-","_",c(DEM[DEM$cluster %in% metabolic_type & DEM$class%in% c('Prostaglandin biosynthesis','Eicosanoid metabolism',
#                                                                                 'Leukotriene metabolism','Arachidonic acid metabolism',
#                                                                                 'Linoleate metabolism'),'gene']))
module_label=gsub("-","_",c('RE1514M','RE1514X','r1251','PROAKGOX1r','HMR-4776','HMR-4077'))
# 'HMR_1129',
# 'HMR_1149'                          
#                            )     ## 20-oxo-LTB4)
output_name='myeloid_arg_fao_cor_angiogenesis.pdf'
width=5
height=5

color='#FBB065'
get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,'#DB5DB8')

# For SLM metabolic states
feature='Antigen_Presenting'
metabolic_type=c('SLM')
# metabolic_type=c('SLM')
# module_label=gsub("-","_",DEM[DEM$cluster %in% metabolic_type & DEM$length>1 & DEM$class%in% c('Chondroitin / heparan sulfate biosynthesis','Chondroitin sulfate degradation',
#                                                                                                'Heparan sulfate degradation','Keratan sulfate degradation',
#                                                                                               'Keratan sulfate biosynthesis'),'gene'])
module_label=gsub("-","_",c('HMR-0642','HMR-0629','HEXA1l','HEXA2l'))

output_name='myeloid_slm_cor_antigen_presenting.pdf'
width=5
height=5
get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,"#1965B0")

## For GST metabolic states
feature='ECM'
metabolic_type=c('GST')
# metabolic_type=c('SLM')
# module_label=gsub("-","_",DEM[DEM$cluster %in% metabolic_type & DEM$length>1 & DEM$class%in% c('Lysine metabolism'),'gene'])
module_label=DEM[DEM$cluster %in% metabolic_type & DEM$length>1 & DEM$class%in% c('Glycine, serine and threonine metabolism'),'gene']

output_name='myeloid_gst_cor_ecm.pdf'
width=5
height=5
get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,"#CA131F")


