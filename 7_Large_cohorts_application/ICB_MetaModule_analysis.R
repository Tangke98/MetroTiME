#' Analyzes ICB response based on MetaModule.
#' @param ICB_data Path to MetaModule and MetaRegulon in myeloid cells object.
#' @param MetaModule Path to metabolic genes in a metabolic reaction.
#' @param MetaModule_info Path to metabolic reaction information.
#' @param output_path Path to the output file.
#' @author Ke Tang
#
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
    
    library(pROC)
})

args <- commandArgs(trailingOnly = TRUE)

ICB_data<-args[1]
MetaModule=args[2]
MetaModule_info=args[3]
output_path=args[4]

data<-readRDS(ICB_data)
meta=as.data.frame(data$followup)
rownames(meta)=paste0(meta$dataset_name,"|",meta$Patient)

MetaModule=readRDS(MetaModule)
MetaModule_info=readRDS(MetaModule_info)

gsva<-as.data.frame(t(readRDS(paste0(output_path,'icb.datasets.GSVA_enlarge.rds'))))
gsva$Response=meta[rownames(gsva),'Response']
tide<-read.table(paste0(output_path,'tide_res.txt'),sep='\t',row.names = 1,header = T)
exp_norm<-readRDS(paste0(output_path,'icb.datasets.env.norm.rds'))
exp_norm=t(exp_norm)
cbind=cbind(tide[rownames(gsva),],gsva,exp_norm[rownames(gsva),])

get_auc_curve=function(data,feature){
    roc_curve <- roc(data[,'Response'], data[,feature])
    return(roc_curve)
}

feature=c('HMR-3398','HEXA1l','HEXAHBl','FAOXC13C11m','IFNG','CD274','CD8','CCR5','CXCL13','CXCL9','MSI.Score','CTL')
get_auc_curve_r=lapply(as.list(feature),get_auc_curve,data=cbind)
names(get_auc_curve_r)=c('HMR-3398','HEXA1l','HEXAHBl','FAOXC13C11m','IFNG','CD274','CD8','CCR5','CXCL13','CXCL9','MSI.Score','CTL')

output_name='ICB_auc_curve_r.pdf'
width=7
height=4
options(repr.plot.width =width, repr.plot.height = height,repr.plot.res = 100)
tf_auc_g <- ggroc(get_auc_curve_r)+
        annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
        theme_bw()+
        theme(panel.grid =element_blank()) + 
        theme(axis.text = element_blank()) + 
        theme(axis.title=element_text(size=20),axis.text=element_text(size=20),
              legend.text=element_text(size=20),legend.title=element_text(size=20))+
        labs(title = 'MMs ICB Response Prediction (R)')+
        scale_color_manual(values=c('HMR-3398'='#E31A1C','HEXA1l'='#FF7F00','HEXAHBl'='#FDBF6F','FAOXC13C11m'='#FB9A99','IFNG'='#33A02C','CD274'='#A6CEE3',
                                    'CD8'='#B2DF8A','CCR5'='#1F78B4','CXCL13'='#CAB2D6','CXCL9'='#6A3D9A','MSI.Score'='#FFFF99','CTL'='#B15928'),name='group')
print(tf_auc_g)
pdf(paste0(output_path,output_name),width=width,height=height)
    print(tf_auc_g)
dev.off()


feature=c('HMR-2309','HMR-7207','HMR-7220','HMR-2281','TIDE','Dysfunction','Exclusion','MDSC','CAF','TAM.M2')
get_auc_curve_nr=lapply(as.list(feature),get_auc_curve,data=cbind)
names(get_auc_curve_nr)=c('HMR-2309','HMR-7207','HMR-7220','HMR-2281','TIDE','Dysfunction','Exclusion','MDSC','CAF','TAM.M2')

output_name='ICB_auc_curve_nr.pdf'
width=6.5
height=4
options(repr.plot.width =width, repr.plot.height = height,repr.plot.res = 100)
tf_auc_g <- ggroc(get_auc_curve_nr)+
        annotate(geom = "segment", x = 1, y = 0, xend =0, yend = 1)+
        theme_bw()+
        theme(panel.grid =element_blank()) + 
        theme(axis.text = element_blank()) +
        theme(axis.title=element_text(size=20),axis.text=element_text(size=20),
              legend.text=element_text(size=20),legend.title=element_text(size=20))+
        labs(title = 'MMs ICB Response Prediction (NR)')+
        scale_color_manual(values=c('HMR-2309'='#E31A1C','HMR-7207'='#FF7F00','HMR-7220'='#FDBF6F','HMR-2281'='#FB9A99','TIDE'='#33A02C','Dysfunction'='#A6CEE3',
                                    'Exclusion'='#B2DF8A','MDSC'='#1F78B4','CAF'='#CAB2D6','TAM.M2'='#6A3D9A'),name='group')
print(tf_auc_g)
pdf(paste0(output_path,output_name),width=width,height=height)
    print(tf_auc_g)
dev.off()

## boxplot
rbind<-readRDS(paste0(output_path,'AUC.rds'))
rbind_NR=rbind[c('HMR-2309','HMR-7220','HMR-2281','HMR-7207','TIDE','Dysfunction','Exclusion','MDSC','CAF','TAM.M2'),]
rbind_NR$feature=factor(rbind_NR$feature,levels = c('HMR-2309','HMR-7220','HMR-2281','HMR-7207','TIDE','Dysfunction','Exclusion','MDSC','CAF','TAM.M2'))
rbind_NR=rbind_NR[order(rbind_NR$feature),]

output_name='ICB_auc_curve_nr_proportion.pdf'
width=6.5
height=4
options(repr.plot.width =width, repr.plot.height = height,repr.plot.res = 100)

p=ggplot(data=rbind_NR, mapping=aes(x=feature,y=auc,fill=feature))+
    geom_bar(stat="identity")+
    theme_bw()+
    theme(panel.grid =element_blank()) +
    theme(axis.text = element_blank()) + 
    theme(axis.title=element_text(size=20),axis.text=element_text(size=20),
          legend.text=element_text(size=20),legend.title=element_text(size=20),
          axis.text.x = element_text(angle = 90, hjust = 1))+
    coord_cartesian(ylim = c(0.4, 0.7))+
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")+
    scale_fill_manual(values=c('HMR-2309'='#E31A1C','HMR-7207'='#FF7F00','HMR-7220'='#FDBF6F','HMR-2281'='#FB9A99','TIDE'='#33A02C','Dysfunction'='#A6CEE3',
                                    'Exclusion'='#B2DF8A','MDSC'='#1F78B4','CAF'='#CAB2D6','TAM.M2'='#6A3D9A'),name='group')
print(p)
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p)
dev.off()

rbind_R=rbind[c('HMR-3398','FAOXC13C11m','HEXA1l','HEXAHBl','IFNG','CD274','CD8','CCR5','CXCL13','CXCL9','MSI.Score','CTL'),]
rbind_R$feature=factor(rbind_R$feature,levels = c('HMR-3398','FAOXC13C11m','HEXA1l','HEXAHBl','IFNG','CD274','CD8','CCR5','CXCL13','CXCL9','MSI.Score','CTL'))
rbind_R=rbind_R[order(rbind_R$feature),]

output_name='ICB_auc_curve_r_proportion.pdf'
width=6.5
height=4
options(repr.plot.width =width, repr.plot.height = height,repr.plot.res = 100)

p=ggplot(data=rbind_R, mapping=aes(x=feature,y=auc,fill=feature))+
    geom_bar(stat="identity")+
    theme_bw()+
    theme(panel.grid =element_blank()) + 
    theme(axis.text = element_blank()) + 
    theme(axis.title=element_text(size=20),axis.text=element_text(size=20),
          legend.text=element_text(size=20),legend.title=element_text(size=20),
          axis.text.x = element_text(angle = 90, hjust = 1))+
    coord_cartesian(ylim = c(0.4, 0.7))+
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")+
    scale_fill_manual(values=c('HMR-3398'='#E31A1C','HEXA1l'='#FF7F00','HEXAHBl'='#FDBF6F','FAOXC13C11m'='#FB9A99','IFNG'='#33A02C','CD274'='#A6CEE3',
                                    'CD8'='#B2DF8A','CCR5'='#1F78B4','CXCL13'='#CAB2D6','CXCL9'='#6A3D9A','MSI.Score'='#FFFF99','CTL'='#B15928'),name='group')
print(p)
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p)
dev.off()

