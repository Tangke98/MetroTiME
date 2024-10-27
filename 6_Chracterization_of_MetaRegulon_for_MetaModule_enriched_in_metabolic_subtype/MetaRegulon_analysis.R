#' Performs MetaRegulon analysis on specific cell lineages.
#' @param cell_lineage The cell lineage to analysis.
#' @param integrated_object_path Path to the specific cell lineages object.
#' @param MetaModule_info Path to the metabolic reaction information.
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
    library(circlize)
    library(RColorBrewer)
    suppressMessages(library(ROGUE))
    suppressMessages(library(tidyverse))
    library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

cell_lineage<-args[1]
integrated_object_path<-args[2]
MetaModule_info=args[3]
output_path=args[4]

integrated_object_path=paste0(integrated_object_path,cell_lineage,'/')
output_path=paste0(output_path,cell_lineage,'/')

seurat<-readRDS(paste0(integrated_object_path,cell_lineage,'_integration_annotation_celltype_metatype.rds'))

fontsize_label=8
fontsize_title=8
color_cancer_type=c('BCC'="#4DAF4A",'BRCA'='#984EA3','CESC'="#FFFF33",'GIST'="#FFFF33",'CHOL'='#A65628',
                    'CRC'='#F781BF','ESCA'='#999999','HNSC'='#66C2A5','HNSCC'='#66C2A5','NPC'='#66C2A5',
                    'KIRC'='#FC8D62',
                    'KIPAN'='#FC8D62','KICH'='#FC8D62',
                    'LIHC'='#8DA0CB','NSCLC'='#E78AC3','SCLC'='#E78AC3','OV'='#A6D854','PAAD'='#FFD92F',
                    'PRAD'='#B3DE69','SKCM'='#E5C494','STAD'='#80B1D3','BLCA'='#B3B3B3',
                    'LSCC'='#8DD3C7','OS'='#FFFFB3','OSCC'='#FFFFB3','SCC'='#BEBADA','SS'='#FB8072',
                    'THCA'='#FDB462','LUAD'='#FCCDE5','MCC'='#D9D9D9','UCEC'='#BC80BD','PBMC'='red',
                   'PPB'='red','UVM'='blue')

if (cell_lineage=='Fibroblast'){
    levels_metabolic=c('OXP','IPM','GLYCAN','GLY','FAO','PUFA','LYS','AA')
    color_metabolic=c('FAO'="#97CADC",'LYS'='#B2B31D','AA'="#92C274",
                                                'GLY'='#FBB2B4','GLYCAN'='#E47B7E',
                                                'PUFA'='#FBB065','IPM'='#984EA3',
                                                'OXP'='#DB5DB8')
    module=c(
        'HMR-6918', ## 4 OXP
        'HMR-6558',
    #     'HMR-6550',##'HMR-6549', ##3 IPM
        'HMR-7494',##'HMR-7493', ## 1 Chondroitin / heparan sulfate biosynthesis
        'HMR-4193', ## GLY
        'FAOXOHC16C16DCc' ,## 0 FAO
        'HMR-0940',##'HMR-1149' ,##2 PUFA  ##'RE3287R',
        'HMR-8025', ##'HMR-4241','HMR-6975', ##5 LYS
        'HMR-3912'
    )
    gene=c('CEBPD','YBX1',
    'NOTCH3','TCF4','EPAS1',
    'LOX','ZNF281',
    'ATF4','HIF1A',
    'KLF15','NFE2L2',
    'NFKB1','REL','SREBF1','MYC',
    'ZNF83','ZNF207',
    'TFDP1'

    #SLM2
    #STM
    #TRP
    )
    module_df=data.frame(module=module,metabolic_type=c(rep('OXP',1),rep('IPM',1),rep('GLYCAN',1),rep('GLY',1),rep('FAO',1),rep('PUFA',1),rep('LYS',1),rep('AA',1)))

}
if (cell_lineage=='Myeloid'){
    levels_metabolic=c("PUFA_LTM",'PUFA_EM','FAS','SLM','ARG','OXP','FAO','GST','PUFA','AA') ## ,'TRP+ cDC2','SL+ pDC','SL+ Mast'
    color_metabolic=c('PUFA_LTM'='#1C7E76','PUFA_EM'='#92C274','OXP'='Thistle','FAO'='#824DAE',
                                                'SLM'='#8ACDD4','ARG'='#DB5DB8','FAS'='#608FBF','GST'='#CA131F',
                                                'PUFA'='#FBB065','AA'='#FED43B')
    module=c('HMR-1080', ## PUFA_LTM A Leukotriene metabolism 
        'RE3434R',  ## PUFA_EM B Eicosanoid metabolism
                    'RE0344M', ## G FAS Fatty acid synthesis
                    'HMR-0642', ## SLM E Sphingolipid metabolism
                    'PROAKGOX1r', ## ARG F Arginine and proline metabolism
        'HMR-6914', ## OXP C Oxidative phosphorylation
        'RE1514M', ## FAO D Fatty acid oxidation
        'HMR-3750',## GST H Glycine, serine and threonine metabolism
        'FAEL183',#PUFA I Arachidonic acid metabolism
        'HMR-4786')##AA J Nucleotide metabolism
    gene=c('BACH1','KLF13',
        'CUX1',
        'NR1H3','PPARG','NFE2L1',
        'NFE2L2',
        'KLF7',
        'BNIP3',
        'CEBPD','YBX1',
        'RELB',
        'NR2F2','ELF3',
        'REL','NFKB1','MYC',
        'TFDP1','KLF5'
    )
    module_df=data.frame(module=module,metabolic_type=c(rep('PUFA_LTM',1),rep('PUFA_EM',1),rep('FAS',1),rep('SLM',1),rep('ARG',1),rep('OXP',1),rep('FAO',1),
                                                      rep('GST',1),rep('PUFA',1),rep('AA',1)))
}

MM.meta.nozero=readRDS(MetaModule_info)

for (i in 1:nrow(module_df)){
    file_res=paste0(integrated_object_path,module_df[i,2],'/MetaRegulon/',module_df[i,2],':',module_df[i,1],'.txt')
    if (!file.exists(file_res)){
        base=paste0(integrated_object_path,module_df[i,2],'/MetaRegulon/')
        names=rownames(MM.meta.nozero[MM.meta.nozero$`GENE ASSOCIATION` == MM.meta.nozero[module_df[i,1],'GENE ASSOCIATION'],])
        
        if (any(file.exists(paste0(base,module_df[i,2],':',names,'.txt')))){
            file=paste0(base,module_df[i,2],':',names,'.txt')[which(file.exists(paste0(base,module_df[i,2],':',names,'.txt')))]
            file.copy(file,file_res)
        }
    }
}

res=list()
for (i in 1:nrow(module_df)){
    file_res=paste0(integrated_object_path,module_df[i,2],'/MetaRegulon/',module_df[i,2],':',module_df[i,1],'.txt')
    file=read.csv(file_res,header=TRUE,row.names = 1)
    gene=rownames(file[file$rank<51 &file$Ligand_Receptor_interaction==0,])
    df=data.frame(gene=gene,metabolic_subtype=module_df[i,2])
    res[[i]]=df
}

res_rbind=Reduce(function(x,y){rbind(x,y)},res)
res_rbind=res_rbind[!duplicated(res_rbind$gene),]

data=GetAssayData(seurat)
data_use=as.data.frame(t(as.matrix(data[res_rbind$gene,])))
data_use$metabolic_type=seurat$metabolic_type
data_use$Cancer_type=seurat$Cancer_type
data_use_agg=aggregate(.~metabolic_type+Cancer_type,data_use,mean)
data_use_agg=data_use_agg[data_use_agg$metabolic_type %in% levels_metabolic,]
data_use_agg$metabolic_type=factor(data_use_agg$metabolic_type,levels=levels_metabolic)
data_use_agg=data_use_agg[order(data_use_agg$metabolic_type),]

draw_data=t(scale(data_use_agg[,3:ncol(data_use_agg)]))
draw_data_final=t(scale(data_use_agg[,3:ncol(data_use_agg)]))
draw_data_final[draw_data_final>1]=1
draw_data_final[draw_data_final<(-1)]=(-1)

li_factor=list()
li_sample=list()
index=1
for (mt in unique(res_rbind$metabolic_subtype)){
    use_factor=res_rbind[res_rbind$metabolic_subtype==mt,'gene']
    use_sample=draw_data_final[use_factor,rownames(data_use_agg[data_use_agg$metabolic_type==mt,])]
    row_sums <- rowSums(use_sample)
    col_sums <- colSums(use_sample)
    sorted_rows <- names(row_sums[order(row_sums, decreasing = TRUE)])
    sorted_cols <- names(col_sums[order(col_sums, decreasing = TRUE)])
    li_factor[[index]]=sorted_rows
    li_sample[[index]]=sorted_cols
    index=index+1
}

draw_data_final_order=draw_data_final[unlist(li_factor),unlist(li_sample)]
rownames(res_rbind)=res_rbind[,1]
res_rbind=res_rbind[rownames(draw_data_final_order),]
data_use_agg=data_use_agg[colnames(draw_data_final_order),]

row_anno=rowAnnotation(MetaRegulon=res_rbind$metabolic_subtype,
                       col=list(MetaRegulon=color_metabolic),
                       annotation_legend_param = list(
                           MetaRegulon=list(
                               title = "MetaRegulon",
                               labels_gp = gpar(fontsize = fontsize_label),
                               title_gp = gpar(fontsize = fontsize_title)
                                )))
top_anno=HeatmapAnnotation(Metabolic_type=data_use_agg$metabolic_type,
                           Cancer_type=data_use_agg$Cancer_type,
                           col=list(Metabolic_type=color_metabolic,
                                   Cancer_type=color_cancer_type),
                           annotation_legend_param = list(
                               Metabolic_type=list(
                                       title = "Metabolic Type",
                                       labels_gp = gpar(fontsize = fontsize_label),
                                       title_gp = gpar(fontsize = fontsize_title)),
                               Cancer_type=list(
                                       title = "Cancer type",
                                       labels_gp = gpar(fontsize = fontsize_label),
                                       title_gp = gpar(fontsize = fontsize_title))))

draw_data=t(scale(data_use_agg[,3:ncol(data_use_agg)]))



match=match(gene, rownames(draw_data_final_order), nomatch = 0)
match=match[!match==0]
label = rowAnnotation(foo = anno_mark(at = match,
                                      labels = rownames(draw_data_final_order[match,]),
                                      labels_gp = gpar(col = "black", fontsize = 6),
                                      padding =2,
                                      link_width = unit(3, "mm"), 
                                      link_height = unit(3, "mm"))
                     ) 

width=6
height=4
options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=Heatmap(draw_data_final_order,
        cluster_columns = FALSE,cluster_rows = FALSE,
        top_anno=top_anno,
        show_row_names=FALSE,
#         left_anno=row_anno,
        right_anno=label,
        show_column_names = FALSE,
        row_dend_width = unit(1, "mm"),
        col=c('#E0F3F8','white','#FDAE61'),
        column_dend_height = unit(1, "mm"),
        row_names_gp = gpar(fontsize = fontsize_label),
        column_names_gp = gpar(fontsize = fontsize_label),
        heatmap_legend_param=list(title="Scaled Expression",
                                 title_gp = gpar(fontsize = fontsize_title),
                                 labels_gp = gpar(fontsize = fontsize_label)))
print(p)

pdf(paste0(output_path,"_regulon_tf_expression_metabolictype.pdf"),width=width,height=height)
    print(p) 
dev.off()



