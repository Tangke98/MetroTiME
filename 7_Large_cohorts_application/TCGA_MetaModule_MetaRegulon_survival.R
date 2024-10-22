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
    library(GSVA)
    library(ComplexHeatmap)
    library(tidyr)
    library(rPref) 
    library(RColorBrewer)
    
    library("survival")
    library("survminer")
    library(stringr)
    library(BayesPrism)
    library(GSVA)
    library(survminer)
    library(ggforestplot)
    library(reshape2)
    library(dplyr)
    library(ggplot2)
    library(ComplexHeatmap)
    library(tidyr)
    library(circlize)
})

args <- commandArgs(trailingOnly = TRUE)

MM_survival= args[1]
factor_Fibroblasts= args[2]
factor_Myeloid= args[3]
gsva_path=args[4]
tpm_path=args[5]
clin_path=args[6]
output_path=args[7]

df_neg<-readRDS(paste0(MM_survival,'gsva_survival_bad.rds'))
df_pos<-readRDS(paste0(MM_survival,'gsva_survival_good.rds'))

ARG_factor<-readRDS(paste0(factor_Myeloid,'ARG:recom.rds'))
FAO_factor<-readRDS(paste0(factor_Ffactor_Myeloidibroblasts,'FAO:recom.rds'))
GST_factor<-readRDS(paste0(factor_Myeloid,'GST:recom.rds'))
GLYCAN_factor<-readRDS(paste0(factor_Fibroblasts,'GLYCAN:recom.rds'))
SLM_factor<-readRDS(paste0(factor_Myeloid,'SLM:recom.rds'))

get_module_factor_survival=function(ARG_factor,FAO_factor,GST_factor,GLYCAN_factor,module,
                                    tpm_path,gsva_path,clin_path,output_path,width,height){
    
    factor=get_factor(module,ARG_factor,FAO_factor,GST_factor,GLYCAN_factor)
    cancers=list.files(clin_path)
    lapply(cancers,draw_survival,tpm_path=tpm_path,gsva_path=gsva_path,clin_path=clin_path,module=module,factor=factor,output_path=output_path,width=width,height=height)
}

get_factor=function(module,ARG_factor,FAO_factor,GST_factor,GLYCAN_factor){
    if (module %in% names(ARG_factor)){
        factor1=ARG_factor[[module]]
        factor1=factor1[factor1$.level<51,]
        factor1=rownames(factor1)
    }else{
        factor1=NULL
    }
    if (module %in% names(FAO_factor)){
        factor2=FAO_factor[[module]]
        factor2=factor2[factor2$.level<51,]
        factor2=rownames(factor2)
    }else{
        factor2=NULL
    }
    if (module %in% names(GST_factor)){
        factor3=GST_factor[[module]]
        factor3=factor3[factor3$.level<51,]
        factor3=rownames(factor3)
    }else{
        factor3=NULL
    }
    if (module %in% names(GLYCAN_factor)){
        factor4=GLYCAN_factor[[module]]
        factor4=factor4[factor4$.level<51,]
        factor4=rownames(factor4)
    }else{
        factor4=NULL
    }
    factor=unique(c(factor1,factor2,factor3,factor4))
    return(factor)
}


draw_survival=function(cancer,tpm_path,gsva_path,clin_path,module,factor,output_path,width,height){
    tpm=as.matrix(readRDS(paste0(tpm_path,cancer,'.RNAseq.TPM.rds'))) ## read tpm
    gsva<-as.data.frame(t(readRDS(paste0(gsva_path,cancer,'.rds')))) ## read module res
    
    gsva<-as.data.frame(t(readRDS(paste0(gsva_path,cancer,'.rds'))))
    CLIN<-read.table(paste0(clin_path,cancer,'/TCGA_',cancer,'.clin.tsv'),sep='\t',header=T,row.names=1)
    CLIN$os<-ifelse(CLIN$vital_status=='Alive',CLIN$days_to_last_follow_up,CLIN$days_to_death)
    Character_sur=c('os', 'vital_status')
    SUR_use=CLIN[,Character_sur]
    SUR_use$vital_status=ifelse(SUR_use$vital_status=='Alive',0,1)
    cbind(SUR_use,t(as.data.frame(strsplit(rownames(SUR_use),'[.]'))))->SUR_use_new
    SUR_use_new[,'4']=as.numeric(SUR_use_new[,'4'])
    SUR_use_new$group=ifelse(SUR_use_new[,'4']< 10,'Tumor','Normal') ##定义normal和tumor
    SUR_use_new=SUR_use_new[order(SUR_use_new[,'4'],decreasing = FALSE),]
    SUR_use_new$patient=paste0(SUR_use_new[,'1'],'.',SUR_use_new[,'2'],'.',SUR_use_new[,'3'])
    SUR_use_new=SUR_use_new[!duplicated(SUR_use_new[,c('group','patient')]),] ##patient 和group重复的去除
    
    merge=merge(gsva,SUR_use_new,by='row.names')
    merge<-na.omit(merge)
    merge<-merge[merge$os>30,]  # save the patient whose os>30
    merge$os<-(merge$os)/365
    merge_use=merge[merge$os<10,] 
    merge_use=merge_use[merge_use$group=='Tumor',]
    rownames(merge_use)=merge_use[,1]
    merge_use=merge_use[,-1]
    
    tpm_use=t(tpm[rownames(tpm) %in% factor,rownames(merge_use)]) ##找到需要循环的factor
    merge_final=cbind(merge_use,tpm_use)
    
    merge_final$module_class=ifelse(merge_final[,module]> median(merge_final[,module]),'high_module','low_module')
    for (factor_use in intersect(factor,colnames(merge_final))){ ##对于循环的factor进行做图
        if (factor_use %in% colnames(merge_final)){
            output_name=paste0(cancer, ":", module, ":", factor_use,'.pdf')
            file_save=paste0(output_path,output_name)
            if (!file.exists(file_save)){
                merge_final$factor_class=ifelse(merge_final[,factor_use]> median(merge_final[,factor_use]),'high_factor','low_factor')
                merge_final$final_group=ifelse(merge_final$module_class=='high_module' & merge_final$factor_class=='high_factor','high_MM_high_MR',
                                ifelse(merge_final$module_class=='high_module' & merge_final$factor_class=='low_factor','high_MM_low_MR',
                                      ifelse(merge_final$module_class=='low_module' & merge_final$factor_class=='high_factor','low_MM_high_MR',
                                            'low_MM_low_MR')))
                fit <- survival::survfit(Surv(os, vital_status) ~ final_group, data = merge_final)
                merge_final <<- merge_final
                surv_diff <- NULL
                surv_diff <- survdiff(Surv(os, vital_status) ~ final_group,data = merge_final)
                p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
                if (p.value<0.05){
                    options(repr.plot.width = width, repr.plot.height = height, repr.plot.res = 100)
                    res_survival_plot = ggsurvplot(fit, pval = FALSE, pval.method = TRUE, 
                                                   conf.int = FALSE, palette = c("#FB8072", "#FDB462", "#608FBF", "#97CADC"),
                                                   legend.labs = c(paste0("high_MM_high_MR(",table(merge_final$final_group)['high_MM_high_MR'],")"),
                                                                   paste0("high_MM_low_MR(",table(merge_final$final_group)['high_MM_low_MR'],")"), 
                                                                   paste0("low_MM_high_MR(",table(merge_final$final_group)['low_MM_high_MR'],")"), 
                                                                   paste0("low_MM_low_MR(",table(merge_final$final_group)['low_MM_low_MR'],")")),
                                                   legend.title = "Group",
                                                   tables.theme = theme_cleantable(), 
                                                   ggtheme = theme_bw(), 
                                                   risk.table = FALSE, risk.table.y.text.col = FALSE,
                                                   font.main = c(20), font.x = c(20), font.y = c(20),
                                                   font.tickslab = c(20), legend = "right") 
                    p=res_survival_plot$plot+
                        theme(panel.grid =element_blank()) + 
                        ggtitle(paste0(cancer, ":", module, ":", factor_use, " P_value=", format(p.value, digits = 3)))
                    pdf(file_save,width=width,height=height)
                        print(p)
                    dev.off()
                    }
            }
    }
    }
}

width=7
height=5

lapply(as.list(rownames(df_neg)),get_module_factor_survival,ARG_factor=ARG_factor,FAO_factor=FAO_factor,GST_factor=GST_factor,GLYCAN_factor=GLYCAN_factor,
       tpm_path=tpm_path,gsva_path=gsva_path,clin_path=clin_path,output_path=output_path,width=width,height=height)

lapply(as.list(rownames(df_pos)),get_module_factor_survival,SLM_factor=SLM_factor,
       tpm_path=tpm_path,gsva_path=gsva_path,clin_path=clin_path,
       output_path=output_path,width=width,height=height)




















