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
    library(circlize)
})

args <- commandArgs(trailingOnly = TRUE)

myeloid_path=args[1]
fibroblast_path=args[2]
survival_path=args[3]
expression_path=args[4]
TF_list=args[5]
auc_value=args[6]
output_path=args[7]

ARG_factor<-readRDS(paste0(myeloid_path,'/sc/regulon_v3/recom_res/ARG:recom.rds'))
FAO_factor<-readRDS(paste0(myeloid_path,'/sc/regulon_v3/recom_res/FAO:recom.rds'))
GST_factor<-readRDS(paste0(myeloid_path,'/sc/regulon_v3/recom_res/GST:recom.rds'))
GLYCAN_factor<-readRDS(paste0(fibroblast_path,'/sc/regulon_v3/recom_res/v1/GLYCAN:recom.rds'))
df_neg<-readRDS(paste0(survival_path,'gsva_survival_bad.rds'))

expression<-read.csv(paste0(expression_path,'Expression_Public_23Q4_subsetted.csv'))
rownames(expression)=expression[,1]

expression=expression[,9:19152]
expression=as.data.frame(t(expression))

output_name_use='CCLE_gsva_cancer'

li1 <- Filter(function(x) !is.null(x), ARG_factor[as.character(df_neg$Var1)])
li2 <- Filter(function(x) !is.null(x), FAO_factor[as.character(df_neg$Var1)])
li3 <- Filter(function(x) !is.null(x), GST_factor[as.character(df_neg$Var1)])
li4 <- Filter(function(x) !is.null(x), GLYCAN_factor[as.character(df_neg$Var1)])

for (i in 1:length(li1)){
    df=li1[[i]]
    df$module=names(li1[i])
    df$gene=rownames(df)
    li1[[i]]=df
}
for (i in 1:length(li2)){
    df=li2[[i]]
    df$module=names(li2[i])
    df$gene=rownames(df)
    li2[[i]]=df
}
for (i in 1:length(li3)){
    df=li3[[i]]
    df$module=names(li3[i])
    df$gene=rownames(df)
    li3[[i]]=df
}
for (i in 1:length(li4)){
    df=li4[[i]]
    df$module=names(li4[i])
    df$gene=rownames(df)
    li4[[i]]=df
}
              
df1=Reduce(function(x,y){rbind(x,y)},li1)
df2=Reduce(function(x,y){rbind(x,y)},li2)
df3=Reduce(function(x,y){rbind(x,y)},li3)
df4=Reduce(function(x,y){rbind(x,y)},li4)

rbind=rbind(df1,df2,df3,df4)
factor_use=rbind[rbind$rank<51,]
tf<-read.table(paste0(TF_list,'tf.txt'))
factor_use=factor_use[factor_use$gene %in% tf$V1,]
auc_tmp=read.csv(auc_value)

auc_tmp_use=auc_tmp[,c('depmap_id','broad_id','auc')]
auc_tmp_use_nodup <- auc_tmp_use %>%
    group_by(depmap_id, broad_id) %>%
    arrange(desc(auc), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    as.data.frame() 

auc_tmp_use_nodup_wide <- pivot_wider(auc_tmp_use_nodup, names_from = broad_id, values_from = auc) %>%
    as.data.frame()

auc_tmp_use_nodup_wide_nona <- auc_tmp_use_nodup_wide %>%
    filter(!is.na(depmap_id))

rownames(auc_tmp_use_nodup_wide_nona)=auc_tmp_use_nodup_wide_nona[,1]
auc_tmp_use_nodup_wide_nona=auc_tmp_use_nodup_wide_nona[,-1]
saveRDS(auc_tmp_use_nodup_wide_nona,
       paste0(auc_tmp,'drug_response_auc.rds'))

li=list()
index=1
for (factor in factor_use$gene){
    df=auc_tmp[grep(factor,auc_tmp$target),]
    if (nrow(df)>0){
        li[[index]]=df
        names(li)[index]=factor
        index=index+1
    }
}

MM<-readRDS(paste0(output_path,'CCLE_gsva_cancer.rds'))
AUC<-readRDS(paste0(auc_tmp,'drug_response_auc.rds'))

sample_use=intersect(rownames(AUC),colnames(MM))
AUC_use=AUC[sample_use,unique(unlist(drug))]
MM_use=t(MM[,sample_use]) 
li=list() 
index=1
for (i in colnames(AUC_use)){
    for (j in colnames(MM_use)){
        cor=cor(AUC_use[,i], MM_use[,j], use="pairwise.complete.obs")
        df=data.frame(drug=i,auc=j,value=cor)
        li[[index]]=df
        index=index+1
    }
}
saveRDS(li,paste0(output_path,'/drug_response/drug_factor_response.rds'))

drug_module_cor=Reduce(function(x,y){rbind(x,y)},li)
colnames(drug_module_cor)=c('drug','module','cor')
drug_factor_response=readRDS(paste0(output_path,'/drug_response/drug_factor_response.rds'))
li_drug_factor=list() 
for (i in 1:length(drug_factor_response)){
    df=data.frame(drug=unique(drug_factor_response[[i]]$broad_id),factor=names(drug_factor_response)[i])
    li_drug_factor[[i]]=df
}
li_drug_factor_res=Reduce(function(x,y){rbind(x,y)},li_drug_factor)
factor_use_res=factor_use[factor_use$gene %in% li_drug_factor_res$factor,c('gene','module')]
colnames(factor_use_res)[1]='factor'
drug_factor_module_info <- inner_join(li_drug_factor_res, factor_use_res, by = "factor",relationship = "many-to-many")
result <- inner_join(drug_module_cor, drug_factor_module_info, by = c("drug", "module"))

drug_screen=unique(result[result$cor<0,'drug'])
drug_module_cor_use=drug_module_cor[drug_module_cor$drug %in% drug_screen,]
drug_module_cor_use_wide <- pivot_wider(drug_module_cor_use, names_from = drug, values_from = cor) %>%
    as.data.frame() 
rownames(drug_module_cor_use_wide)=drug_module_cor_use_wide[,1]
drug_module_cor_use_wide=drug_module_cor_use_wide[,-1]
Heatmap(scale(t(drug_module_cor_use_wide)))

drug_factor=result[result$drug %in% drug_screen,c(1,4)]
drug_factor=drug_factor[!duplicated(drug_factor),]
li=list()
index=1
for (drug in unique(drug_factor$drug)){
    combined_string <- paste(drug_factor[drug_factor$drug==drug,'factor'], collapse = ",")
    df=data.frame(drug=drug,target=combined_string)
    li[[index]]=df
    index=index+1
}
li_res_drug_factor=Reduce(function(x,y){rbind(x,y)},li)
auc_tmp_use=auc_tmp[,c('broad_id','name')]
auc_tmp_use=auc_tmp_use[!duplicated(auc_tmp_use),]
rownames(auc_tmp_use)=auc_tmp_use[,1]
rownames(li_res_drug_factor)=li_res_drug_factor[,1]

li_res_drug_factor$name=auc_tmp_use[rownames(li_res_drug_factor),'name']
li_res_drug_factor$rowname=paste0(li_res_drug_factor$name,"|",li_res_drug_factor$target)  


drug_info<-read.csv(paste0(output_path,'Conpound_information.csv'))
drug_info=drug_info[!duplicated(drug_info[,c('broad_id','name')]),]
drug_info=drug_info[,c('broad_id','drug_category')]
drug_info=drug_info[!duplicated(drug_info),]
rownames(drug_info)=drug_info[,1]
li_res_drug_factor$drug_category=drug_info[rownames(li_res_drug_factor),'drug_category']  
module_factor=result[result$drug %in% drug_screen,c(2,4)]
module_factor=module_factor[!duplicated(module_factor),]
li=list()
index=1
for (module in unique(module_factor$module)){
    combined_string <- paste(module_factor[module_factor$module==module,'factor'], collapse = ",")
    df=data.frame(module=module,target=combined_string)
    li[[index]]=df
    index=index+1
}
li_res_module_factor=Reduce(function(x,y){rbind(x,y)},li)
li_res_module_factor$name=paste0(li_res_module_factor$module,"|",li_res_module_factor$target)

rownames(li_res_module_factor)=li_res_module_factor$module
draw=t(drug_module_cor_use_wide)[rownames(li_res_drug_factor),rownames(li_res_module_factor)]
rownames(draw)=li_res_drug_factor$rowname
colnames(draw)=li_res_module_factor$name

draw[draw>0]=0
row_anno=rowAnnotation(drug_category=li_res_drug_factor$drug_category,
                      col=list(drug_category=c( ##'noncancer'='#A02C76','targeted cancer'='#6B112A','chemo'='green',
                                              'chemo'="#FDB462", 'noncancer'="#608FBF", 'targeted cancer'="#FB8072")))
top_anno=HeatmapAnnotation(Macro_SPP1=li_res_module_factor_retain$Macro_SPP1,
                           Macro_THBS1=li_res_module_factor_retain$Macro_THBS1,
                           Macro_SLPI=li_res_module_factor_retain$Macro_SLPI,
                           Fibro_CTHRC1=li_res_module_factor_retain$Fibro_CTHRC1,
                      col=list(Macro_SPP1=c('0'='white','1'='#FB8072'),
                              Macro_THBS1=c('0'='white','1'='#FB8072'),
                              Macro_SLPI=c('0'='white','1'='#FB8072'),
                              Fibro_CTHRC1=c('0'='white','1'='#FB8072')))
drug_use=c('tyloxapol|NFKB2','TWS-119|MYC,JUN','nutlin-3|TP53','idasanutlin|TP53','BAY-87-2243|HIF1A','PAC-1|SP3','triptolide|RELA,REL','2-methoxyestradiol|HIF1A',
          '4-methylgenistein|ESRRA','desoxycortone|NR3C1','acetylcysteine|RELA,REL','carvedilol|HIF1A')
col_fun=colorRampPalette(rev(c('white',brewer.pal(9, "Reds"))))

width=8
height=9
options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=Heatmap(draw,
       col=col_fun(100),
          right_anno=label,
#           top_anno=top_anno,
       show_row_dend=FALSE,show_column_dend=FALSE,
          show_row_names=FALSE,
       border=TRUE,
       heatmap_legend_param=list(title="Pearson correlation"),
         left_anno=row_anno
         )
print(p)

pdf(paste0(output_path,'/drug_response/drug_response_heatmap2.pdf'),width=width,height=height)
    print(p) 
dev.off()

