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
    library(MetroSCREEN)
})

args <- commandArgs(trailingOnly = TRUE)

MetaModule=args[1]
TPM_path=args[2]
clin_path=args[3]
output_path=args[4]

files=list.files(TPM_path)
MM=readRDS(MetaModule)

for (file in files){
    cancer=strsplit(file,'[.]')[[1]][1]
    tpm<-as.matrix(readRDS(paste0(TPM_path,cancer,'.RNAseq.TPM.rds')))
    cal_MetaModule(tpm,MM,output_path,paste0(cancer,'_GSVA'))
}

each_module_survival=function(merge_use,cancer,module){
    merge_use=merge_use[order(merge_use[,module],decreasing = TRUE),]
    top_30_percent <- round(nrow(merge_use) * 0.3)
    top_30_df <- merge_use[1:top_30_percent,]
    # 计算tail 30%的观测数
    tail_30_percent <- round(nrow(merge_use) * 0.3)
    # 选择tail 30%的数据
    tail_30_df <- merge_use[(nrow(merge_use) - tail_30_percent + 1):nrow(merge_use),]
    top_30_df$group='high'
    tail_30_df$group='low'
    rbind=rbind(top_30_df,tail_30_df)
    fit<-coxph(Surv(os,vital_status)~group, data=rbind)
    res=summary(fit)
    outTab=data.frame()
    outTab=cbind(
        name=cancer,
        beta=log10(res$conf.int[,"exp(coef)"]),
        se=log10(res$coefficients[,"se(coef)"]),
        pvalue=res$coefficients[,"Pr(>|z|)"],
        N= res$n,
        module=module)
    return(outTab)
}
module_survival=function(gsva_path,clin_path,cancer,output_path){
    file_save=paste0(output_path,cancer,'.survival.rds')
    if (!file.exists(file_save)){
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
        each_module_survival_res=lapply(as.list(colnames(gsva)),each_module_survival,
                                        merge_use=merge_use,cancer=cancer)
        result <-Reduce(function(x,y){rbind(x,y)},each_module_survival_res)
        saveRDS(result,file_save)
    }
}

files=list.files(output_path)
cancers=lapply(strsplit(files,'[.]'),function(x){x[[1]]})

lapply(cancers[31:length(cancers)],module_survival,gsva_path=output_path,
       clin_path=clin_path,output_path=output_path)
       