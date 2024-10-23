#' Annotates cell metabolism for the Metacells of fibroblasts.
#' @param integrated_object_path_fibroblasts Path to the fibroblasts object.
#' @param metabolic_gene Path to the metabolic gene file.
#' @param MetaModule Path to the file listing metabolic genes involved in a reaction.
#' @param MetaModule_info Path to the metabolic reaction information file.
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
    library(tidyr)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(MetroSCREEN)
})
args <- commandArgs(trailingOnly = TRUE)
integrated_object_path_fibroblasts<-args[1]
metabolic_gene<-args[2]
output_path=args[3]
MetaModule=args[4]
MetaModule_info=args[5]

DotPlot_change=function (object, features, levels) 
{
    cols = c("lightgrey", "blue")
    col.min = -2.5
    col.max = 2.5
    dot.min = 0
    dot.scale = 6
    idents = NULL
    group.by = NULL
    split.by = NULL
    cluster.idents = FALSE
    scale = TRUE
    scale.by = "radius"
    scale.min = NA
    scale.max = NA
    assay = NULL
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    split.colors <- !is.null(x = split.by) && !any(cols %in% 
        rownames(x = brewer.pal.info))
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    feature.groups <- NULL
    if (is.list(features) | any(!is.na(names(features)))) {
        feature.groups <- unlist(x = sapply(X = 1:length(features), 
            FUN = function(x) {
                return(rep(x = names(x = features)[x], each = length(features[[x]])))
            }))
        if (any(is.na(x = feature.groups))) {
            warning("Some feature groups are unnamed.", call. = FALSE, 
                immediate. = TRUE)
        }
        features <- unlist(x = features)
        names(x = feature.groups) <- features
    }
    cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
    data.features <- FetchData(object = object, vars = features, 
        cells = cells)
    data.features$id <- if (is.null(x = group.by)) {
        Idents(object = object)[cells, drop = TRUE]
    }
    else {
        object[[group.by, drop = TRUE]][cells, drop = TRUE]
    }
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
        if (split.colors) {
            if (length(x = unique(x = splits)) > length(x = cols)) {
                stop("Not enough colors for the number of groups")
            }
            cols <- cols[1:length(x = unique(x = splits))]
            names(x = cols) <- unique(x = splits)
        }
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident, 
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    if (cluster.idents) {
        mat <- do.call(what = rbind, args = lapply(X = data.plot, 
            FUN = unlist))
        mat <- scale(x = mat)
        id.levels <- id.levels[hclust(d = dist(x = mat))$order]
    }
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    ngroup <- length(x = levels(x = data.plot$id))
    if (ngroup == 1) {
        scale <- FALSE
        warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
    }
    else if (ngroup < 5 & scale) {
        warning("Scaling data with a low number of groups may produce misleading results", 
            call. = FALSE, immediate. = TRUE)
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == 
                x, "avg.exp"]
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min, 
                  max = col.max)
            }
            else {
                data.use <- log1p(x = data.use)
            }
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (split.colors) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, 
        levels = features)
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (split.colors) {
        splits.use <- vapply(X = as.character(x = data.plot$id), 
            FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                paste(sort(x = levels(x = object), decreasing = TRUE), 
                  collapse = "|"), ")_)"), replacement = "", 
            USE.NAMES = FALSE)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    if (!is.null(x = feature.groups)) {
        data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
            levels = unique(x = feature.groups))
    }
    data.plot$id = factor(data.plot$id, levels = levels)
    data.plot = data.plot[order(data.plot$id, decreasing = FALSE), 
        ]
    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
        color = color.by)) + scale.func(range = c(0, dot.scale), 
        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
        labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
            yes = "Identity", no = "Split Identity"))
}

inte_use<-readRDS(paste0(integrated_object_path_fibroblasts,'fibroblasts_integration_annotation_celltype.rds'))
metabolic=readRDS(metabolic_gene)
VariableFeatures(inte_use)=metabolic
inte_use <- ScaleData(inte_use)
inte_use <- RunPCA(inte_use)
## use the parameters which can best group the metacells
inte_use <- RunUMAP(inte_use, dims = 1:45)
inte_use <- FindNeighbors(inte_use, dims = 1:45,k.param =70)
inte_use <- FindClusters(inte_use, resolution = 0.5)

options(repr.plot.width = 8, repr.plot.height = 5,repr.plot.res = 100)
DimPlot(inte_use, reduction = "umap",cols=c('0'="#4DAF4A",'1'='#984EA3','2'="#FFFF33",'3'='#A65628',
                                              '4'='#F781BF','5'='#999999','6'='#66C2A5','7'='#FC8D62',
                                              '8'='#8DA0CB','9'='#E78AC3','10'='#A6D854','11'='#FFD92F',
                                              '12'='#E5C494','13'='#B3B3B3','14'='#BC80BD','15'='#FFFFB3',
                                              '16'='#BEBADA','17'='#D9D9D9','18'='#80B1D3','19'='#80B1D3',
                                            '20'='#FDB462','21'='#B3DE69'),label=TRUE)

Idents(inte_use)='cell_type'
options(repr.plot.width = 8, repr.plot.height = 5,repr.plot.res = 100)
DimPlot(inte_use, reduction = "umap",cols=c('Fibro_SFRP1'="#97CADC",'Fibro_CCL5'='#FFDBB0','Fibro_IL6'="#F7AA8A",'Fibro_CTHRC1'='#F87379',
                 'Fibro_SAA1'='#92C274','MyoFibro_RGS5'='Thistle','MyoFibro_MYH11'='#FBB2B4','Fibro_CDC20'='#57998B'))

saveRDS(inte_use,
        paste0(output_path,'fibroblasts_integration_metabolic_cluster.rds'))

Idents(inte_use) <- 'seurat_clusters'
markers <- FindAllMarkers(inte_use, only.pos = TRUE)
markers_use=markers[markers$gene %in%metabolic, ] ##only focus metabolic genes
markers_use=markers_use[order(markers_use$avg_log2FC,decreasing = TRUE),]

cluster=c('0','1','2','3','4','5','6','7')

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

output_name='fibroblasts_metatype_cluster_original_marker.pdf'
width=20
height=4
## use the top 10 metabolic markers to decide which metabolic clusters can group together
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
cal_MetaModule(metacell,MM,output_path,'fibroblasts_gsva_nodup')

## read the metadata
MM.meta.nozero=readRDS(MetaModule_info) %>% as.data.frame()

## enlarge our results
gsva_score=readRDS(paste0(output_path,'fibroblasts_gsva_nodup.rds'))
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
saveRDS(gsva_score_enlarge,paste0(output_path,'fibroblasts_gsva.rds'))

## use the MetaModule scores to annotate cell metabolic states
gsva.seurat <- CreateSeuratObject(counts = gsva_score, project = "fibroblasts", min.cells = 0, min.features = 0)
gsva.seurat$metabolic_cluster=as.factor(inte_use@meta.data[rownames(inte_use@meta.data),'seurat_clusters'])
Idents(gsva.seurat) <- 'metabolic_cluster'
markers <- FindAllMarkers(gsva.seurat, only.pos = TRUE)

MetaModule.markers$metabolic_type=MM.meta[MetaModule.markers$gene,'SUBSYSTEM']
MetaModule.markers$reaction=MM.meta[MetaModule.markers$gene,'EQUATION']
saveRDS(MetaModule.markers,paste0(output_path,'fibroblasts_metatype_cluster_DEM.rds'))

new.cluster.ids <- c("FAO",'GLYCAN','PUFA','IPM','OXP','LYS','AA','GLY')
names(new.cluster.ids) <- levels(inte_use)
inte_use <- RenameIdents(inte_use, new.cluster.ids)
inte_use@meta.data$metabolic_type=Idents(inte_use)
saveRDS(inte_use,paste0(output_path,'fibroblasts_integration_annotation_celltype_metatype.rds'))
saveRDS(inte_use@meta.data,paste0(output_path,'fibroblasts_integration_annotation_celltype_metatype_metadata.rds'))


options(repr.plot.width =5, repr.plot.height = 5,repr.plot.res = 100)

DimPlot(inte_use, reduction = "umap",cols=c('FAO'="#97CADC",'LYS'='#B2B31D','AA'="#92C274",
                                            'GLY'='#FBB2B4','GLYCAN'='#E47B7E',
                                            'PUFA'='#FBB065','IPM'='#984EA3',
                                            'OXP'='#DB5DB8'))

## dotplot to show the specific MetaModule in each metabolic states
metabolic_gsva=t(gsva_score_enlarge)
features=c(
    'HMR-6918', ## 4 OXP
    'HMR-6558',
#     'HMR-6550',##'HMR-6549', ##3 IPM
    'HMR-7494',##'HMR-7493', ## 1 Chondroitin / heparan sulfate biosynthesis
    'HMR-4193', ## GLY
    'FAOXOHC16C16DCc' ,## 0 FAO
     'HMR-0940',##'HMR-1149' ,##2 PUFA  ##'RE3287R',
    'HMR-0960',
    'HMR-1330',
    'HMR-8025',
    'HMR_6555', ##'HMR-4241','HMR-6975', ##5 LYS
    'HMR-3912'
 )

levels=rev(c('OXP','IPM','GLYCAN','GLY','FAO','PUFA','LYS','AA'))
output_name='fibroblasts_metatype_marker.pdf'
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
inte_use<-readRDS(paste0(output_path,'fibroblasts_integration_annotation_celltype_metatype.rds'))
gsva_score_enlarge<-readRDS(paste0(output_path,'fibroblasts_gsva.rds'))

MM.meta.nozero=readRDS(MetaModule_info) %>% as.data.frame()

## use the MetaModule scores to annotate cell metabolic states
gsva.seurat <- CreateSeuratObject(counts = gsva_score_enlarge, project = "fibroblasts", min.cells = 0, min.features = 0)
gsva.seurat$metabolic_type=as.factor(inte_use@meta.data[rownames(inte_use@meta.data),'metabolic_type'])
Idents(gsva.seurat) <- 'metabolic_type' ## findall markers to decide the metabolic states

MetaModule.markers <- FindAllMarkers(gsva.seurat, only.pos = TRUE)

MetaModule.markers$metabolic_type=MM.meta.nozero[MetaModule.markers$gene,'SUBSYSTEM']
MetaModule.markers$reaction=MM.meta.nozero[MetaModule.markers$gene,'EQUATION']
MetaModule.markers$length=MM.meta.nozero[MetaModule.markers$gene,'length']

saveRDS(MetaModule.markers,paste0(output_path,'fibroblasts_metatype_type_DEM.rds'))

## Find DEG for metabolic states
inte_use<-readRDS(paste0(output_path,'fibroblasts_integration_annotation_celltype_metatype.rds'))
fibroblasts.markers <- FindAllMarkers(inte_use, only.pos = TRUE)
fibroblasts.markers=fibroblasts.markers[fibroblasts.markers$p_val_adj<0.05,]
fibroblasts.markers=fibroblasts.markers[order(fibroblasts.markers$avg_log2FC,decreasing = TRUE),]
saveRDS(fibroblasts.markers,paste0(output_path,'fibroblasts_metatype_type_DEG.rds'))