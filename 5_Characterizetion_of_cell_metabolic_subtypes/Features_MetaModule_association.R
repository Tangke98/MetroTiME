#' Associates function score with the MetaRegulon score in specific cell lineages.
#' @param integrated_object_path Path to the specific cell lineages Metacell TPM file.
#' @param metabolic_score Path to the metabolic score file.
#' @param feature_score Path to the function score file.
#' @param DEM Path to the differentially enriched MetaModule file.
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

cell_lineage<-args[1]
metabolic_score<-args[2]
feature_score<-args[3]
DEM<-args[4]
MetaModule_info<-args[5]
output_path<-args[6]

metabolic_score=paste0(metabolic_score,cell_lineage,'/',cell_lineage)
feature_score=paste0(feature_score,cell_lineage,'/',cell_lineage)
DEM=paste0(DEM,cell_lineage,'/',cell_lineage)
output_path=paste0(output_path,cell_lineage,'/',cell_lineage)

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

metabolic_score=t(readRDS(paste0(metabolic_score,"_gsva.rds")))
feature_score=t(readRDS(paste0(feature_score,'_features_gsva.rds')))

cor=as.data.frame(cor(metabolic_score,feature_score))
MM.meta.nozero=readRDS(MetaModule_info) %>% as.data.frame()
cor$class=MM.meta.nozero[rownames(cor),'SUBSYSTEM']
DEM=readRDS(paste0(DEM,'_metatype_type_DEM.rds'))

if (cell_lineage=='Fibroblast'){
    ## For PUFA metabolic states
    feature='Immune_regulatory'
    metabolic_type=c('PUFA')
    module_label=c('HMR-1330', ## prostaglandinD2
    'HMR-0960', ## 5(S)-HETE
    'HMR-0940')
    output_name='_pufa_immune_regulatory.pdf'
    width=5
    height=5

    color='#FBB065'
    get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,'#FBB065')

    ## For GLYCAN metabolic states
    feature='ECM'
    metabolic_type=c('GLYCAN')
    module_label=c('HMR-7493', ## Chondroitin sulfate
    'HMR-7494', ## Chondroitin sulfate
    'HMR-3833') ## HMR_3833

    output_name='_glycan_cor_ecm.pdf'
    width=5
    height=5
    get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,"#E47B7E")

    ## For LYS metabolic states
    feature='Antigen_presenting'
    metabolic_type=c('LYS')
    # metabolic_type=c('SLM')
    # module_label=gsub("-","_",DEM[DEM$cluster %in% metabolic_type & DEM$length>1 & DEM$class%in% c('Lysine metabolism'),'gene'])
    module_label=c('HMR-6555','HMR-8026') ## HMR_3833

    output_name='_lys_cor_antigen_presenting.pdf'
    width=5
    height=5
    get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,"#B2B31D")
}
if (cell_lineage=='Myeloid'){
    ## For ARG and FAO metabolic states
    feature='Angiogenesis'
    metabolic_type=c('ARG','FAO')
    
    module_label=gsub("-","_",c('RE1514M','RE1514X','r1251','PROAKGOX1r','HMR-4776','HMR-4077'))
    # 'HMR_1129',
    # 'HMR_1149'                          
    #                            )     ## 20-oxo-LTB4)
    output_name='_arg_fao_cor_angiogenesis.pdf'
    width=5
    height=5

    color='#FBB065'
    get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,'#DB5DB8')

    # For SLM metabolic states
    feature='Antigen_Presenting'
    metabolic_type=c('SLM')
    module_label=gsub("-","_",c('HMR-0642','HMR-0629','HEXA1l','HEXA2l'))

    output_name='_slm_cor_antigen_presenting.pdf'
    width=5
    height=5
    get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,"#1965B0")

    ## For GST metabolic states
    feature='ECM'
    metabolic_type=c('GST')
    module_label=DEM[DEM$cluster %in% metabolic_type & DEM$length>1 & DEM$class%in% c('Glycine, serine and threonine metabolism'),'gene']

    output_name='_gst_cor_ecm.pdf'
    width=5
    height=5
    get_cor(cor,feature,DEM,metabolic_type,module_label,output_path,output_name,width,height,"#CA131F")
}

