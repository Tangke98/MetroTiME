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
    library(GSVA)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(MetroSCREEN)
})
args <- commandArgs(trailingOnly = TRUE)
integrated_object_path_myeloid<-args[1]
metabolic_gene<-args[2]
output_path=args[3]
MetaModule=args[4]
MetaModule_info=args[5]

inte_use<-readRDS(paste0(integrated_object_path_myeloid,'myeloid_integration_annotation_celltype_subset.rds'))
metabolic=readRDS(metabolic_gene)
VariableFeatures(inte_use)=metabolic
inte_use <- ScaleData(inte_use)
inte_use <- RunPCA(inte_use)
## use the parameters which can best group the metacells
inte_use <- RunUMAP(inte_use, dims = 1:45)
inte_use <- FindNeighbors(inte_use, dims = 1:45,k.param =80)
inte_use <- FindClusters(inte_use, resolution = 1)

options(repr.plot.width = 8, repr.plot.height = 5,repr.plot.res = 100)
DimPlot(inte_use, reduction = "umap",cols=c('0'="#4DAF4A",'1'='#984EA3','2'="#FFFF33",'3'='#A65628',
                                              '4'='#F781BF','5'='#999999','6'='#66C2A5','7'='#FC8D62',
                                              '8'='#8DA0CB','9'='#E78AC3','10'='#A6D854','11'='#FFD92F',
                                              '12'='#E5C494','13'='#B3B3B3','14'='#BC80BD','15'='#FFFFB3',
                                              '16'='#BEBADA','17'='#D9D9D9','18'='#80B1D3','19'='#80B1D3',
                                            '20'='#FDB462','21'='#B3DE69'),label=TRUE)

Idents(inte_use)='cell_type'
options(repr.plot.width = 8, repr.plot.height = 5,repr.plot.res = 100)
DimPlot(inte_use, reduction = "umap",cols=c('Mono_CD16'='#92C274','Mono_FCN1'='#1C7E76','Mono_THBS1'='#824DAE','Macro_SPP1'='#DB5DB8',
                 'Macro_SLPI'='#CA131F','Macro_C1QC'='#8ACDD4','Macro_IL1B'='#1965B0','Macro_IL32'='#FBB065',
                 'Macro_CDC20'='#FED43B','cDC2_CD1C'='#B2B31D','cDC1_CLEC9A'='#D2D77A','pDC_LILRA4'='#BBBBBB',
                 'Mast_CPA3'='#EDE2F6','Macro_CLEC10A'='#B8A1CD','Macro_CCL18'='#FD8586',
                 'cDC1_LAMP3'='#C9D79F'))

saveRDS(inte_use,
        paste0(output_path,'myeloid_integration_metabolic_cluster_sub.rds'))

Idents(inte_use) <- 'seurat_clusters'
markers <- FindAllMarkers(inte_use, only.pos = TRUE)
markers_use=markers[markers$gene %in%metabolic, ] ##only focus metabolic genes
markers_use=markers_use[order(markers_use$avg_log2FC,decreasing = TRUE),]

cluster=c('19','20','6','0','16','18', ## A
          '4', # B
          '1','21', #C
          '3','14', #D
          '2','5','7','8','15','17', #E
          '9', #F
          '10', #G
          '11', #H
          '12', #I
          '13' #J
         )

## find the metabolic markers for each cluster
markers_use$cluster=factor(markers_use$cluster,levels=cluster)
markers_use=markers_use[order(markers_use$cluster),]
markers_use %>%
    group_by(cluster) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

features = unique(top10$gene)

p1 <- DotPlot(inte_use,assay='RNA',features = features) + 
    coord_flip() + 
    theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
    labs(x=NULL,y=NULL) + 
    guides(size = guide_legend("Percent Expression") )+ #legend
    scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 
df<- p1$data
exp_mat<-df %>% 
    dplyr::select(-pct.exp, -avg.exp) %>%  
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
    as.data.frame() 
row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

brewer_colors <- rev(c(brewer.pal(11, "PiYG")[2:4],'#F7F7F7',brewer.pal(11, "BrBG")[6:8]))
color_f <- colorRampPalette(brewer_colors)
colors_100 <- color_f(100)

output_name='myeloid_metatype_cluster_original_marker.pdf'
width=20
height=4
## use the top 10 metabolic markers to decide which metabolic clusters can group together
options(repr.plot.width = width, repr.plot.height =height,repr.plot.res = 100)
p=pheatmap(t(exp_mat)[cluster,],color=colors_100,cluster_rows = FALSE,cluster_cols = FALSE)
print(p)

pdf(paste0(output_path,output_name),width=width,height=height)
    print(p)
dev.off()

## metabolic states temporarily
Idents(inte_use)='seurat_clusters'
new.cluster.ids <- c("A",
                     'C','E','D','B','E',
                     'A','E','E','F','G',
                     'H','I','J','D','E',
                     'A','E','A','A','A',
                     'C')
names(new.cluster.ids) <- levels(inte_use)
inte_use <- RenameIdents(inte_use, new.cluster.ids)
inte_use@meta.data$metabolic_cluster=Idents(inte_use)
options(repr.plot.width =5, repr.plot.height = 5,repr.plot.res = 100)

DimPlot(inte_use, reduction = "umap",cols=c('A'='#92C274','B'='#1C7E76','C'='Thistle','D'='#F87379',
                                            'E'='#F781BF','F'='#97CADC','G'='#1965B0','H'='#6FB3A7',
                                            'I'='#A6D854','J'='#FDCDAC'))

Idents(inte_use) <- 'metabolic_cluster'
markers <- FindAllMarkers(inte_use, only.pos = TRUE)
markers_use=markers[markers$gene %in%metabolic, ] ## only focus on metabolic genes
markers_use=markers_use[order(markers_use$avg_log2FC,decreasing = TRUE),]

cluster=c('A','B','C','D','E','F','G','H','I','J')

markers_use$cluster=factor(markers_use$cluster,levels=cluster)
markers_use=markers_use[order(markers_use$cluster),]
markers_use %>%
    group_by(cluster) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

features = unique(top10$gene)

p1 <- DotPlot(inte_use,assay='RNA',features = features) + 
    coord_flip() + 
    theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
    labs(x=NULL,y=NULL) + 
    guides(size = guide_legend("Percent Expression") )+ #legend
    scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 
df<- p1$data
exp_mat<-df %>% 
    dplyr::select(-pct.exp, -avg.exp) %>%  
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
    as.data.frame() 
row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

brewer_colors <- rev(c(brewer.pal(11, "PiYG")[2:4],'#F7F7F7',brewer.pal(11, "BrBG")[6:8]))
color_f <- colorRampPalette(brewer_colors)
colors_100 <- color_f(100)    

output_name='myeloid_metatype_cluster_marker.pdf'
width=12
height=4

options(repr.plot.width = width, repr.plot.height =height,repr.plot.res = 100)
p=pheatmap(t(exp_mat)[cluster,],color=colors_100,cluster_rows = FALSE,cluster_cols = FALSE)
print(p)
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p)
dev.off()

## calculate MetaModule scores for each metabolic states
## MM information
MM=readRDS(MetaModule)

DefaultAssay(inte_use)='integrated'
metacell=as.matrix(GetAssayData(inte_use))
cal_MetaModule(metacell,MM,output_path,'myeloid_gsva_subset_nodup')

## read the metadata
MM.meta.nozero=readRDS(MetaModule_info) %>% as.data.frame()

## enlarge our results
gsva_score=readRDS(paste0(output_path,'myeloid_gsva_subset_nodup.rds'))
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
saveRDS(gsva_score_enlarge,paste0(output_path,'myeloid_subset_gsva.rds'))

## use the MetaModule scores to annotate cell metabolic states
# gsva_score=readRDS(paste0(output_path,'myeloid_gsva_subset_nodup'))
# gsva.seurat <- CreateSeuratObject(counts = gsva_score, project = "myeloid", min.cells = 0, min.features = 0)
# gsva.seurat$metabolic_cluster=as.factor(inte_use@meta.data[rownames(inte_use@meta.data),'metabolic_cluster'])

# Idents(gsva.seurat) <- 'metabolic_cluster'
# markers <- FindAllMarkers(gsva.seurat, only.pos = TRUE)

# MetaModule.markers$metabolic_type=MM.meta[MetaModule.markers$gene,'SUBSYSTEM']
# MetaModule.markers$reaction=MM.meta[MetaModule.markers$gene,'EQUATION']
# saveRDS(markers,paste0(output_path,'myeloid_metatype_cluster_DEM.rds'))

new.cluster.ids <- c("PUFA_LTM",'OXP','SLM','FAO','PUFA_EM','ARG','FAS','GST','PUFA','AA')
Idents(inte_use)='metabolic_cluster'
names(new.cluster.ids) <- levels(inte_use)
inte_use <- RenameIdents(inte_use, new.cluster.ids)
inte_use@meta.data$metabolic_type=Idents(inte_use)
saveRDS(inte_use,paste0(output_path,'myeloid_integration_annotation_celltype_subset_metatype.rds'))
saveRDS(inte_use@meta.data,paste0(output_path,'myeloid_integration_annotation_celltype_subset_metatype_metadata.rds'))

options(repr.plot.width =7, repr.plot.height = 6,repr.plot.res = 100)

DimPlot(inte_use, reduction = "umap",cols=c('PUFA_LTM'='#1C7E76','PUFA_EM'='#92C274','OXP'='Thistle','FAO'='#824DAE',
                                            'SLM'='#8ACDD4','ARG'='#DB5DB8','FAS'='#608FBF','GST'='#CA131F',
                                            'PUFA'='#FBB065','AA'='#FED43B'))

## dotplot to show the specific MetaModule in each metabolic states

features = c('HMR-1080', ## PUFA_LTM A Leukotriene metabolism 
'RE3434R',  ## PUFA_EM B Eicosanoid metabolism
'HMR-6911', ## OXP C Oxidative phosphorylation
'RE1514M', ## FAO D Fatty acid oxidation
'HMR-7580', ## SLM E Sphingolipid metabolism
'PROAKGOX1r','HMR-4776', ## ARG F Arginine and proline metabolism
'RE0344M', ## G FAS Fatty acid synthesis
'HMR-4699','HMR-7703',## GST H Glycine, serine and threonine metabolism
'HMR-0941',#PUFA I Arachidonic acid metabolism
'HMR-6780')##AA J Nucleotide metabolism

levels=rev(c("PUFA_LTM",'PUFA_EM','OXP','FAO','SLM','ARG','FAS','GST','PUFA','AA'))
output_name='myeloid_metatype_metatype_marker.pdf'
width=6
height=4

object=inte_use
# metabolic_gsva=t(gsva_score_enlarge)
object@meta.data=cbind(object@meta.data,metabolic_gsva[rownames(object@meta.data),])
options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=DotPlot_change(object,features,levels)+
    theme_bw() + 
    theme(panel.grid =element_blank()) + #去除网格线
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_colour_gradient2(low="#1A519C",mid='white',high="#A10921")+
    theme(axis.title=element_text(size=14),
            axis.text=element_text(size=14),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14)
            )
print(p) 
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p) 
dev.off()

## Find DEM for metabolic states
inte_use<-readRDS(paste0(output_path,'myeloid_integration_annotation_celltype_subset_metatype.rds'))
gsva_score_enlarge<-readRDS(paste0(output_path,'myeloid_subset_gsva.rds'))
MM.meta.nozero=readRDS(MetaModule_info) 
## use the MetaModule scores to annotate cell metabolic states
gsva.seurat <- CreateSeuratObject(counts = gsva_score_enlarge, project = "myeloid", min.cells = 0, min.features = 0)
gsva.seurat$metabolic_type=as.factor(inte_use@meta.data[rownames(inte_use@meta.data),'metabolic_type'])
Idents(gsva.seurat) <- 'metabolic_type' ## findall markers to decide the metabolic states
MetaModule.markers <- FindAllMarkers(gsva.seurat, only.pos = TRUE)
MetaModule.markers$metabolic_type=MM.meta.nozero[MetaModule.markers$gene,'SUBSYSTEM']
MetaModule.markers$reaction=MM.meta.nozero[MetaModule.markers$gene,'EQUATION']
MetaModule.markers$length=MM.meta.nozero[MetaModule.markers$gene,'length']
saveRDS(MetaModule.markers,paste0(output_path,'myeloid_metatype_type_subset_DEM.rds'))

## use the MetaModule scores to annotate cell metabolic states
inte_use<-readRDS(paste0(output_path,'myeloid_integration_annotation_celltype_subset_metatype.rds'))
gsva_score<-readRDS(paste0(output_path,'myeloid_gsva_subset_nodup.rds'))
MM.meta.nozero=readRDS(MetaModule_info) 
gsva.seurat <- CreateSeuratObject(counts = gsva_score, project = "myeloid", min.cells = 0, min.features = 0)
gsva.seurat$metabolic_type=as.factor(inte_use@meta.data[rownames(inte_use@meta.data),'metabolic_type'])
Idents(gsva.seurat) <- 'metabolic_type' ## findall markers to decide the metabolic states
MetaModule.markers <- FindAllMarkers(gsva.seurat, only.pos = TRUE)
MetaModule.markers$metabolic_type=MM.meta.nozero[MetaModule.markers$gene,'SUBSYSTEM']
MetaModule.markers$reaction=MM.meta.nozero[MetaModule.markers$gene,'EQUATION']
MetaModule.markers$length=MM.meta.nozero[MetaModule.markers$gene,'length']
saveRDS(MetaModule.markers,paste0(output_path,'myeloid_metatype_type_subset_DEM_nodup.rds')) ##修改名字

## Find DEG for metabolic states
inte_use<-readRDS(paste0(output_path,'myeloid_integration_annotation_celltype_subset_metatype.rds'))
myeloid.markers <- FindAllMarkers(inte_use, only.pos = TRUE)
myeloid.markers=myeloid.markers[myeloid.markers$p_val_adj<0.05,]
myeloid.markers=myeloid.markers[order(myeloid.markers$avg_log2FC,decreasing = TRUE),]
saveRDS(myeloid.markers,paste0(output_path,'myeloid_metatype_type_subset_DEG.rds'))