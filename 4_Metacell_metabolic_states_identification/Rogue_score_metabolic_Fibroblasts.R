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
    library(tidyr)
    library(ROGUE)
    library(ComplexHeatmap)
    library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)

integrated_object_path_fibroblasts<-args[1]
metabolic_gene<-args[2]
output_path<-args[3]

cal_purity=function (seurat_rna, resolution, cluster, feature, metabolic_gene) 
{
    exp = seurat_rna@assays$integrated@data
    meta = seurat_rna@meta.data
    meta_use = meta[seurat_rna[[feature]] == cluster, ]
    exp_use = exp[, colnames(exp) %in% rownames(meta_use)]
    exp_use = exp_use[rownames(exp_use) %in% metabolic_gene, 
        ]
    exp.res <- SE_fun(exp_use)
    rogue.value <- CalculateRogue(exp.res, platform = "UMI")
    df = data.frame(cluster = cluster, resolution = resolution, 
        rogue = rogue.value)
    return(df)
}

inte_use<-readRDS(paste0(integrated_object_path_fibroblasts,'fibroblasts_integration_annotation_celltype_metatype.rds'))
metabolic=readRDS(metabolic_gene)
genes=metabolic
levels_metabolic=c('OXP','IPM','GLYCAN','GLY','FAO','PUFA','LYS','AA')
color_metabolic=c('FAO'="#97CADC",'LYS'='#B2B31D','AA'="#92C274",
                                            'GLY'='#FBB2B4','GLYCAN'='#E47B7E',
                                            'PUFA'='#FBB065','IPM'='#984EA3',
                                            'OXP'='#DB5DB8')

cal_purity_res=do.call(rbind,lapply(as.list(unique(as.character(inte_use$metabolic_type))),cal_purity,resolution=1,seurat_rna=inte_use,feature='metabolic_type',metabolic_gene=genes))
cal_purity_res$cluster=as.factor(cal_purity_res$cluster)
cal_purity_res$cluster=factor(cal_purity_res$cluster,levels=levels_metabolic)
cal_purity_res=cal_purity_res[order(cal_purity_res$cluster),]

width=4
height=4

options(repr.plot.width = width, repr.plot.height = height+1,repr.plot.res = 100)
p=ggplot(cal_purity_res, aes(cluster, rogue)) +
    geom_segment(aes(x=cluster, 
               xend=cluster,
               y=0, 
               yend=rogue))+
    scale_y_continuous(limits = c(0.95,1),breaks = seq(0.95, 1, by = 0.01))+
    theme_bw()+
    theme(panel.grid =element_blank()) + 
    geom_point(shape=21,size=3,colour="black",aes(fill=cluster))+
    labs(x='Cluster',y='ROGUE Score')+
    scale_fill_manual(values=c(color_metabolic))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(axis.title=element_text(size=15),axis.text=element_text(size=15),
          legend.title=element_text(size=15),legend.text=element_text(size=15),
          axis.text.x = element_text(angle = 90, hjust = 1))+NoLegend()
print(p)
pdf(paste0(output_path,"rigue_score.pdf"),width=width,height=height+1)
    print(p)
dev.off()
