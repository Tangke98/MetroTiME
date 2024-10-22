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
    library(missMDA)
    library(FactoMineR)
    library(remef)
    library(lme4)
    library(tibble)
    library(lmerTest)
    library(missMDA)
    library(MetroSCREEN)
})

args <- commandArgs(trailingOnly = TRUE)

MetaModule=args[1]
MetaModule_info=args[2]
TPM_path=args[3]
immune_path=args[4]
output_path=args[5]

get_tpm_integrated=function(input_path,file){
    setwd(input_path)
    cancer=strsplit(file,'[.]')[[1]][1]
    
    tpm<-as.data.frame(t(readRDS(paste0(input_path,file))))
    tpm_use=tpm[as.numeric(str_sub(rownames(tpm),14,15)) < 10,]
    tpm_use$cancer=cancer
    return(tpm_use)
}

files=list.files(TPM_path)
get_tpm_integrated_res=lapply(as.list(files),get_tpm_integrated,input_path=TPM_path)
get_tpm_integrated_res_rbind=Reduce(function(x,y){rbind(x,y)},get_tpm_integrated_res)
saveRDS(get_tpm_integrated_res_rbind,paste0(output_path,'tcga_integrated_data.rds'))

immune<-read.csv(immune_path,header = TRUE,row.names=1)

df=as.data.frame(t(as.data.frame(strsplit(rownames(get_tpm_integrated_res_rbind),'[.]'))))
rownames(df)=rownames(get_tpm_integrated_res_rbind)
df_use=df[order(df$V4),]
get_tpm_integrated_res_rbind$patient=paste0(df$V1,"-",df$V2,"-",df$V3)

get_tpm_integrated_res_rbind_use=get_tpm_integrated_res_rbind[rownames(df_use),]
get_tpm_integrated_res_rbind_use=get_tpm_integrated_res_rbind_use[get_tpm_integrated_res_rbind_use$patient %in% rownames(immune),]
get_tpm_integrated_res_rbind_use_nodup=get_tpm_integrated_res_rbind_use[!duplicated(get_tpm_integrated_res_rbind_use$patient),]
rownames(get_tpm_integrated_res_rbind_use_nodup)=get_tpm_integrated_res_rbind_use_nodup$patient

rownames(get_tpm_integrated_res_rbind_use_nodup)=get_tpm_integrated_res_rbind_use_nodup[,'patient']
get_tpm_integrated_res_rbind_use_nodup=get_tpm_integrated_res_rbind_use_nodup[,-18864:-18865]

MM=readRDS(MetaModule)

cal_MetaModule(t(get_tpm_integrated_res_rbind_use_nodup),MM,output_path,paste0('tcga_integrated_GSVA'))

immune<-read.csv(immune_path)
ck<- immune[,c(1:2,5:32)] %>%
    dplyr::rename(sample =`TCGA.Participant.Barcode`,
                  cancer = `TCGA.Study`) %>%
    mutate_at(colnames(.)[-c(1,2)],as.numeric) %>%
    mutate(sample = gsub(sample,pattern = "-",replacement = ".")) %>%
    column_to_rownames(var = "sample")

#impute features
rr.impute <- imputePCA(ck[,-1], ncp=2)
dataset.ck.impute<- as.data.frame(rr.impute$completeObs)
dataset.ck.impute$cancer =ck$cancer

#remove the cancertype effect from the actual value of each feature
feature.rev =list()
for (feature in colnames(dataset.ck.impute)[-ncol(dataset.ck.impute)]) {
    dataset.ck.impute[,c(feature,"cancer")] %>%
        lmer(as.formula(paste0("`",feature,"`","~","(1|cancer)")),data = .) ->tmp.mod
    y_partial <- remef(tmp.mod, ran = "all")
    feature.rev[[feature]]=y_partial
}
feature.rev.tcga = do.call(cbind,feature.rev)
rownames(feature.rev.tcga) = rownames(dataset.ck.impute)

saveRDS(feature.rev.tcga,paste0(output_path,"feature.rev.tcga.rds"))

feature.rev.tcga.use=feature.rev.tcga[,c(7,8,9,10)]
res.pca <- prcomp(feature.rev.tcga.use,scale. = T)
#Check the PC1 and PC2 of cancer types in tcga
res.ind <- get_pca_ind(res.pca)
res.ind.contri<- as.data.frame(res.ind$coord[,1:2])
res.ind.contri$object <- rownames(res.ind.contri)
res.ind.contri$cancer = dataset.ck.impute$cancer
saveRDS(res.ind.contri[,c(2,ncol(res.ind.contri))],paste0(output_path,"PC2_immune.rds"))

output_name='immune_pca.pdf'
width=6
height=6

res.pca <- prcomp(feature.rev.tcga[,c(7,8,9,10)],scale. = T)
res.ind <- get_pca_ind(res.pca)
res.ind.contri<- as.data.frame(res.ind$coord[,1:2])
res.ind.contri$object <- rownames(res.ind.contri)
res.ind.contri$cancer = dataset.ck.impute$cancer

fviz_pca_var(res.pca,axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)->p
print(p)

pdf(paste0(output_path,output_name),width=width,height=height)
    print(p) 
dev.off()

gsva<-as.data.frame(t(readRDS(paste0(output_path,'tcga_integrated_GSVA.rds'))))
rownames(gsva)=gsub('-',".",rownames(gsva))
pc<-readRDS(paste0(output_path,"PC2_immune.rds"))
merge=merge(gsva,pc,by='row.names')
rownames(merge)=merge[,1]
merge=merge[,-1]

#remove the cancertype effect from the actual value of each feature
feature.rev =list()
for (feature in colnames(merge)[1:(ncol(merge)-2)]) {
    merge[,c(feature,"cancer")] %>%
        lmer(as.formula(paste0("`",feature,"`","~","(1|cancer)")),data = .) ->tmp.mod
    y_partial <- remef(tmp.mod, ran = "all")
    feature.rev[[feature]]=y_partial
}
module.rev.tcga = do.call(cbind,feature.rev)
rownames(module.rev.tcga) = rownames(merge)

saveRDS(module.rev.tcga,paste0(output_path,'module.rev.tcga.rds'))

merge=merge(module.rev.tcga,pc,by='row.names')
rownames(merge)=merge[,1]
merge=merge[,-1]

get_cor=function(merge,feature){
    cor=cor.test(merge[,feature],merge[,'Dim.2'])
    df=data.frame(cor=cor$estimate,pvalue=cor$p.value,module=feature)
    return(df)
}

get_cor_res=lapply(as.list(colnames(merge)[1:(ncol(merge)-2)]),get_cor,merge=merge)
get_cor_res_rbind=Reduce(function(x,y){rbind(x,y)},get_cor_res)
saveRDS(get_cor_res_rbind,paste0(output_path,"module.feature.cor.tcga.rds"))

rownames(get_cor_res_rbind)=get_cor_res_rbind[,'module']
MetaModule_info<-readRDS(MetaModule_info)
get_cor_res_rbind$class=MetaModule_info[rownames(get_cor_res_rbind),'SUBSYSTEM']
get_cor_res_rbind$length=MetaModule_info[rownames(get_cor_res_rbind),'length']
get_cor_res_rbind_use=na.omit(get_cor_res_rbind[get_cor_res_rbind$length>1 & get_cor_res_rbind$pvalue<0.05,])
get_cor_res_rbind_use=get_cor_res_rbind_use[!get_cor_res_rbind_use$class%in% c('Protein degradatio','Isolated','Xenobiotics metabolism','Protein degradation','Protein modification'),]
get_cor_res_rbind_use=get_cor_res_rbind_use[order(get_cor_res_rbind_use$cor,decreasing = FALSE),]
get_cor_res_rbind_use$order=1:nrow(get_cor_res_rbind_use)

df_neg<-readRDS(paste0(output_path,'gsva_survival_bad.rds'))
df_pos<-readRDS(paste0(output_path,'gsva_survival_good.rds'))

module_bad=rownames(df_neg)
module_good=rownames(df_pos)

get_cor_res_rbind_use$name=paste0(rownames(get_cor_res_rbind_use),"|",MetaModule_info[rownames(get_cor_res_rbind_use),'EQUATION'])
get_cor_res_rbind_use$text=ifelse(get_cor_res_rbind_use$module %in% c('HMR-7493','HMR-7494','HMR-7333','HMR-4776','HMR-6710','HMR-6921',
                                                                     'HMR-2172'),get_cor_res_rbind_use$name,'')

get_cor_res_rbind_use$color=ifelse(get_cor_res_rbind_use$module %in% module_bad,'3',
                                   ifelse(get_cor_res_rbind_use$module %in% module_good,'4','5'))

output_name='immune_potential_module.pdf'
width=10
height=10

options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=ggplot(get_cor_res_rbind_use,aes(x=order,y=cor,fill=color,linewidth=2))+geom_bar(stat = "identity")+
    theme_bw() + 
    theme(panel.grid =element_blank()) + #去除网格线
    scale_fill_manual(values=c('1'="#1A519C",'2'='#A10921','3'='#A10921','4'='#1A519C','5'="grey"))+
    labs(x='MetaModule',y="Association with immunosuppressive signature",title=)+
    theme(axis.title=element_text(size=40),
            axis.text=element_text(size=40),
            legend.text=element_text(size=40),
            legend.title=element_text(size=40),
          legend.position = "bottom")+
    ggrepel::geom_text_repel(aes(label=text),max.overlaps = Inf)
p


