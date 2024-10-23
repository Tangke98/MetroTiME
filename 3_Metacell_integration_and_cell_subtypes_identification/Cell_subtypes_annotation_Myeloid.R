#' Annotates cell subtypes for the Metacells of myeloid cells.
#' @param integrated_object_path_myeloid Path to the myeloid cells object.
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
})
args <- commandArgs(trailingOnly = TRUE)
integrated_object_path_myeloid<-args[1]
output_path=args[2]

Dimplot_umap=function (object, group, color, title, pt.size, width, height, 
    output_path, output_name) 
{
    options(repr.plot.width = width + 2, repr.plot.height = height, 
        repr.plot.res = 100)
    p = DimPlot(object, reduction = "umap", pt.size = pt.size, 
        group.by = group, cols = c(color)) + labs(title = title) + 
        theme(plot.title = element_text(size = 8), axis.title.x = element_text(size = 8), 
            axis.title.y = element_text(size = 8), axis.text = element_text(size = 8), 
            legend.text = element_text(size = 8), legend.title = element_text(size = 8))
    print(p)
    pdf(paste0(output_path, output_name, ".umap.pdf"), width = width, 
        height = height)
    print(p)
    dev.off()
    pdf(paste0(output_path, output_name, ".umap_nolegend.pdf"), 
        width = width, height = height)
    print(p + NoLegend())
    dev.off()
}
DotPlot_change=function (object, features, levels){
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

inte_use<-readRDS(paste0(integrated_object_path_myeloid,'myeloid_integration.rds'))
inte_use <- FindVariableFeatures(inte_use, selection.method = "vst", nfeatures = 3000)
inte_use <- ScaleData(inte_use)
inte_use <- RunPCA(inte_use)
## use the parameters which can best group the metacells
inte_use <- RunUMAP(inte_use, dims = 1:20)
inte_use <- FindNeighbors(inte_use, dims = 1:20,k.param =80)
inte_use <- FindClusters(inte_use, resolution = 0.9)

DimPlot(inte_use, reduction = "umap",cols=c("#4DAF4A",'#984EA3',"#FFFF33",'#A65628','#F781BF',
                    '#999999','#FFFFB3','#FC8D62','#8DA0CB','#E78AC3',
                    '#A6D854','#FFD92F','#E5C494','#B3B3B3','#8DD3C7',
                    '#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462',
                    '#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#E6F5C9',
                    '#F2F2F2','#B3E2CD','#FDCDAC',"#4DAF4A",'#984EA3',"#FFFF33",'#A65628','#F781BF',
                    '#999999','#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854'), label = TRUE)
inte_use.markers <- FindAllMarkers(inte_use, only.pos = TRUE)
saveRDS(inte_use.markers,paste0(,'myeloid_integration_DEG.rds'))

## pheatmap to show the markers of clusters
features = c('PTPRC','CSF1R','CD14','S100A9','FCN1','VCAN','FCGR3A','THBS1','APOE','C1QA','C1QB','C1QC','SLPI','SPP1',
            'IL1B','IL10','IL32','CDC20','FLT3','CCL22','CLEC9A','LAMP3','CD1C','GZMB','TPSAB1','CTSG')
cluster=c(
    '1','18','19','22','23', #Mono_FCN1
    '6', #Mono_CD16
    '4','0','16',#Mono_THBS1
    '21','2','7','10','3',#Macro_C1QC
    '11',#Macro_SLPI
    '9',#Macro_SPP1
    '12', #Macro_IL32
    '13',   #Macro_CDC20
    '17', #cDC1_CLEC9A
    '14', #cDC1_LAMP3
    '5','20', #cDC2_CD1C
    '15', #pDC_LILRA4
    '8') #Mast_CPA3
p1 <- DotPlot(inte_use,assay='RNA',features = features) + 
    coord_flip() + 
    theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
        labs(x=NULL,y=NULL) + 
        guides(size = guide_legend("Percent Expression") )+ #legend
        scale_color_gradientn(colours = rev(c('#C51B7D','#DE77AE','#F1B6DA','#F7F7F7','#F5F5F5','#C7EAE5','#80CDC1'))) #颜色
df<- p1$data
exp_mat<-df %>% 
    dplyr::select(-pct.exp, -avg.exp) %>%  
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
    as.data.frame() 
row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()
## set the color
brewer_colors <- rev(c(brewer.pal(11, "PiYG")[2:4],'#F7F7F7',brewer.pal(11, "BrBG")[6:8]))
color_f <- colorRampPalette(brewer_colors)
colors_100 <- color_f(100)

output_name='myeloid_celltype_cluster_marker.pdf'
width=6
height=7

options(repr.plot.width = width, repr.plot.height =height,repr.plot.res = 100)
p=pheatmap(t(exp_mat)[cluster,],color=colors_100,cluster_rows = FALSE,cluster_cols = FALSE)
print(p)
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p) 
dev.off()

## rename the cell subtypes
Idents(inte_use)='seurat_clusters'
new.cluster.ids <- c('Mono_THBS1',
                     'Mono_FCN1','Macro_C1QC','Macro_C1QC','Mono_THBS1','cDC2_CD1C', #5
                     'Mono_CD16','Macro_C1QC','Mast_CPA3','Macro_SPP1','Macro_C1QC', #10
                     'Macro_SLPI','Macro_IL32','Macro_CDC20','cDC1_LAMP3','pDC_LILRA4',#15
                     'Mono_THBS1','cDC1_CLEC9A','Mono_FCN1','Mono_FCN1','cDC2_CD1C',#20
                     'Macro_C1QC','Mono_FCN1','Mono_FCN1')

names(new.cluster.ids) <- levels(inte_use)
inte_use <- RenameIdents(inte_use, new.cluster.ids)
inte_use$cell_type=Idents(inte_use)
color_celltype=c('Mono_CD16'='#92C274','Mono_FCN1'='#1C7E76','Mono_THBS1'='#824DAE','Macro_SPP1'='#DB5DB8',
                 'Macro_SLPI'='#CA131F','Macro_C1QC'='#8ACDD4','Macro_IL1B'='#1965B0','Macro_IL32'='#FBB065',
                 'Macro_CDC20'='#FED43B','cDC2_CD1C'='#B2B31D','cDC1_CLEC9A'='#D2D77A','pDC_LILRA4'='#BBBBBB',
                 'Mast_CPA3'='#EDE2F6','Macro_CLEC10A'='#B8A1CD','Macro_CCL18'='#FD8586',
                 'cDC1_LAMP3'='#C9D79F')
Dimplot_umap(inte_use,'cell_type',color_celltype,'Cell Type',0.1,3,3,output_path,'cell_subtype')

# dotplot to show markers for the cell subtypes
object = inte_use
features = c('PTPRC','CSF1R','CD14','S100A9','FCN1','VCAN','FCGR3A','THBS1','APOE','C1QA','C1QB','C1QC','SLPI','SPP1',
            'IL1B','IL10','IL32','CDC20','FLT3','CCL22','CLEC9A','LAMP3','CD1C','GZMB','TPSAB1','CTSG')

levels=c('Mono_FCN1','Mono_CD16','Mono_THBS1','Macro_C1QC','Macro_SLPI','Macro_SPP1','Macro_IL32','Macro_CDC20','cDC1_CLEC9A',
            'cDC1_LAMP3','cDC2_CD1C','pDC_LILRA4','Mast_CPA3')

output_name='myeloid_celltype_marker.pdf'
width=9
height=4

options(repr.plot.width = 9, repr.plot.height = 4,repr.plot.res = 100)
p=DotPlot_change(object,features,rev(levels))+
    theme_bw() + 
    theme(panel.grid =element_blank()) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_colour_gradient2(low="#1A519C",mid='white',high="#A10921")+
    theme(axis.title=element_text(size=8),
            axis.text=element_text(size=8),
            legend.text=element_text(size=8),
            legend.title=element_text(size=8)
            )
p
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p) 
dev.off()

## find the cell type specific markers,remove the cell subtypes not interested
subset=subset(inte_use,subset=cell_type %in% c('cDC1_LAMP3','cDC2_CD1C','cDC1_CLEC9A','pDC_LILRA4','Mast_CPA3'),invert=TRUE)
DefaultAssay(subset)='integrated'
saveRDS(subset,paste0(output_path,'myeloid_integration_annotation_celltype_subset.rds'))
metadata=subset@meta.data
saveRDS(metadata,paste0(output_path,'myeloid_integration_annotation_celltype_subset_metadata.rds'))

Idents(subset)='cell_type'
inte_use.markers <- FindAllMarkers(subset, only.pos = TRUE)
saveRDS(inte_use.markers,paste0(output_path,'myeloid_integration_annotation_celltype_subset_DEG.rds'))