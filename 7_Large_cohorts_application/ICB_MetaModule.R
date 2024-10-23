#' Analyzes ICB response of MetaModule and MetaRegulon.
#' @param ICB_data Path to ICB response patient data.
#' @param Tide_res Path to TIDE results for ICB response data.
#' @param MetaModule Path to metabolic genes in a metabolic reaction.
#' @param MetaModule_info Path to metabolic reaction information.
#' @param output_path Path to the output file.
#' @author Ke Tang
#
suppressMessages({
    library(harmony)
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library('nichenetr')
    library(ggpubr)
    library(GSVA)
    library(MAESTRO)
    library(ComplexHeatmap)
    library(limma)
    library(RColorBrewer)
    library(rPref) 
    library(kpcalg)
    library(igraph)
    library(RobustRankAggreg)
    library(parallel)
    library(pROC)
    library(MetroSCREEN)
})
args <- commandArgs(trailingOnly = TRUE)

ICB_data<-args[1]
Tide_res<-args[2]
MetaModule<-args[3]
MetaModule_info=args[4]
output_path=args[5]

data<-readRDS(ICB_data)
exp=data$expression
meta=as.data.frame(data$followup)
rownames(meta)=paste0(meta$dataset_name,"|",meta$Patient)
colnames(exp)=rownames(meta)
exp_norm <- log2(exp + 1)
exp_norm <- sweep(exp_norm, 1, rowMeans(exp_norm), "-")
write.table(exp,paste0(output_path,'icb.datasets.env.txt'),quote=F,sep='\t')
saveRDS(exp_norm,paste0(output_path,'icb.datasets.env.norm.rds'))

## calculate the metamodule score
MetaModule=readRDS(MetaModule)
MetaModule_info=readRDS(MetaModule_info)

cal_MetaModule(exp,MetaModule,output_path,'icb.datasets.GSVA')

get_auc=function(data,feature){
    roc_curve <- roc(data[,'Response'], data[,feature])
    auc_value <- auc(roc_curve)
    n=mean(data[data$Response==0,feature])
    r=mean(data[data$Response==1,feature])
    n_r=n-r
    df=data.frame(feature=feature,auc=auc_value,value=n_r)
    return(df)
}

gsva<-as.data.frame(t(readRDS(paste0(output_path,'icb.datasets.GSVA.rds'))))
## read the metadata
MM.meta.nozero=MetaModule_info %>% as.data.frame()
gsva_score=as.data.frame(readRDS(paste0(output_path,'icb.datasets.GSVA.rds')))

## enlarge our results
new_rows <- list()
for (i in rownames(gsva_score)){
    names=rownames(MM.meta.nozero[MM.meta.nozero$`GENE ASSOCIATION` == MM.meta.nozero[i,'GENE ASSOCIATION'],])
    if (length(names)>1){
        df=gsva_score[rownames(gsva_score) %in% names,]
        df_list <- replicate(length(names)-1, df, simplify = FALSE)
        all_dfs <- do.call(rbind, df_list)
        rownames(all_dfs)=names[-which(names==rownames(df))]
        new_rows[[i]]=all_dfs
    }  
}

new_rows_res=Reduce(function(x,y){rbind(x,y)},new_rows)
gsva_score_enlarge=rbind(gsva_score,new_rows_res)
saveRDS(gsva_score_enlarge,paste0(output_path,'icb.datasets.GSVA_enlarge.rds'))

gsva=as.data.frame(t(readRDS(paste0(output_path,'icb.datasets.GSVA_enlarge.rds'))))
gsva$Response=meta[rownames(gsva),'Response']
get_auc_gsva=lapply(as.list(colnames(gsva[1:(ncol(gsva)-1)])),get_auc,data=gsva)
get_auc_gsva_res=Reduce(function(x,y){rbind(x,y)},get_auc_gsva)

## get the auc values predic ted by tide
tide<-read.table(Tide_res,sep='\t',row.names = 1,header = T)
tide$Response=meta[rownames(tide),'Response']
get_auc_tide=lapply(as.list(colnames(tide[c(3:7,9:14)])),get_auc,data=tide)
get_auc_tide_res=Reduce(function(x,y){rbind(x,y)},get_auc_tide)

## get the auc from the immune signature
immune=as.data.frame(t(exp_norm)[,c('CCR5','CXCL13','CXCL9')])
immune$Response=meta[rownames(immune),'Response']
get_auc_immune=lapply(as.list(colnames(immune[,1:(ncol(immune)-1)])),get_auc,data=immune)
get_auc_immune_res=Reduce(function(x,y){rbind(x,y)},get_auc_immune)

## merge
rownames(get_auc_gsva_res)=get_auc_gsva_res[,1]
get_auc_gsva_res$class=MetaModule_info[rownames(get_auc_gsva_res),'SUBSYSTEM']
get_auc_gsva_res$length=MetaModule_info[rownames(get_auc_gsva_res),'length']
get_auc_gsva_res=na.omit(get_auc_gsva_res)
get_auc_gsva_res=get_auc_gsva_res[get_auc_gsva_res$length>1,]
get_auc_gsva_res=get_auc_gsva_res[!(get_auc_gsva_res$class)%in% c('Drug metabolism','Isolated','Xenobiotics metabolism','Pool reactions','Miscellaneous'),]
get_auc_gsva_res=get_auc_gsva_res[order(get_auc_gsva_res$auc,decreasing = TRUE),1:4]
get_auc_tide_res$class=get_auc_tide_res$feature
get_auc_immune_res$class=get_auc_immune_res$feature
rbind=rbind(rbind(get_auc_immune_res,get_auc_tide_res),get_auc_gsva_res)
rownames(rbind)=rbind[,1]
saveRDS(rbind,paste0(output_path,'AUC.rds'))

rbind=readRDS(paste0(output_path,'AUC.rds'))
rbind$name=paste0(rownames(rbind),"|",MetaModule_info[rownames(rbind),'EQUATION'])
rbind=rbind[!duplicated(rbind$auc),]
module_use=c('HMR-3398','FAOXC13C11m','HEXAHBl','HEXA1l','HMR-2281','HMR-2309','HMR-7220','HMR-7205','HMR-7207')
rbind$label=ifelse(rownames(rbind) %in% c('IFNG','CD274','CD8','CCR5','CXCL13','CXCL9','TIDE','MSI.Score','Dysfunction','Exclusion','MDSC','CAF','TAM.M2','CTL'),rbind$class,
                   ifelse(rownames(rbind) %in% module_use,rbind$name,''))
rbind$color=ifelse(rownames(rbind) %in% c('IFNG','CD274','CD8','CCR5','CXCL13','CXCL9','TIDE','MSI.Score','Dysfunction','Exclusion','MDSC','CAF','TAM.M2','CTL'),'1',
                           ifelse(rownames(rbind) %in% rbind[rbind$value<0,'feature'],'3','2')) 

output_name='ICB_all.pdf'
width=6
height=6
options(repr.plot.width = 10, repr.plot.height = 10,repr.plot.res = 100)
p=ggplot(rbind,aes(x=value,y=auc,color=color))+geom_point(aes(x=value,y=auc,color=color))+
    theme_bw() + 
    theme(panel.grid =element_blank()) + #去除网格线
    scale_color_manual(values=c('7'="grey",'1'="grey",'2'='#A10921','3'='#1A519C','4'='grey'))+
    labs(x='Different Score(NR-R)',y="AUC",title=)+
    theme(axis.title=element_text(size=40),
            axis.text=element_text(size=40),
            legend.text=element_text(size=40),
            legend.title=element_text(size=40),
          text=element_text(size=40),
            legend.position = "bottom")+
    ggrepel::geom_text_repel(aes(label=label,fontsize=20),max.overlaps = Inf)
p
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p)
dev.off()                           