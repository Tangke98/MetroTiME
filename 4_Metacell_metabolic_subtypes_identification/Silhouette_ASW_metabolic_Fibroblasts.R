#' Uses Silhouette analysis to evaluate the clustering of fibroblasts based on metabolic genes.
#' @param integrated_object_path_fibroblasts Path to the fibroblasts Metacell TPM file.
#' @param metabolic_gene Path to the metabolic gene file.
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
})

args <- commandArgs(trailingOnly = TRUE)
integrated_object_path_fibroblasts<-args[1]
metabolic_gene<-args[2]
output_path=args[3]

inte_use<-readRDS(paste0(integrated_object_path_fibroblasts,'fibroblasts_integration_annotation_celltype.rds'))
metabolic=readRDS(metabolic_gene)

print_umap_metabolic=function(inte_use,metabolic,dims,k.param,res,output_path){
    file_res=paste0(output_path,'META_celltype_PC_',dims,"_Nei_",k.param,"_Res_",res,"_metadata",'.rds')
    if (!file.exists(file_res)){
        VariableFeatures(inte_use)=metabolic
        inte_use <- ScaleData(inte_use)
        inte_use <- RunPCA(inte_use)
        inte_use <- RunUMAP(inte_use, dims = 1:dims)
        inte_use <- FindNeighbors(inte_use, dims = 1:dims,k.param =k.param)
        inte_use <- FindClusters(inte_use, resolution = res)
        p=DimPlot(inte_use, reduction = "umap",label = TRUE)
        file_res=paste0(output_path,'META_cluster_PC_',dims,"_Nei_",k.param,"_Res_",res,'.pdf')
        pdf(file_res,width=7,height=5)
            print(p) 
        dev.off()

        p=DimPlot(inte_use, reduction = "umap",group.by='cell_type',label = TRUE)
        file_res=paste0(output_path,'META_celltype_PC_',dims,"_Nei_",k.param,"_Res_",res,'.pdf')
        pdf(file_res,width=7,height=5)
            print(p) 
        dev.off()

        file_res=paste0(output_path,'META_celltype_PC_',dims,"_Nei_",k.param,"_Res_",res,"_dist",'.rds')
        dist=inte_use@reductions$umap@cell.embeddings
        saveRDS(dist,file_res)

        file_res=paste0(output_path,'META_celltype_PC_',dims,"_Nei_",k.param,"_Res_",res,"_metadata",'.rds')
        metadata=inte_use@meta.data
        saveRDS(metadata,file_res)
    }
}
draw_si=function(data,PC,celltype,width,height,color,output_path){
    t1=data[data$PC==PC,]
    t1$Res=as.factor(t1$Res)
    
    output_name=paste0(celltype,"_",PC,'.pdf')
    options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
    
    p=ggplot(t1,aes(x=Nei,y=Si,color =Res)) + 
        geom_point(size = 2) + 
        scale_color_manual(values=color) + 
        theme(title = element_text(size = 15), text = element_text(size = 15)) + 
        theme_bw() + 
        geom_hline(yintercept = 0.15, linetype ="dashed") +
        labs(x='Neighbors',y='Silhouette',title=paste0('PC_',PC))+
        theme(panel.grid =element_blank()) + 
        theme(axis.text = element_blank()) + 
        theme(plot.title = element_text(size = 14),
        legend.title=element_text(size=14),
        axis.title=element_text(size=14),
        axis.text = element_text(size = 14),
        legend.text=element_text(size=14))+
        scale_x_continuous(breaks = seq(0, 80, by = 10),  
        limits = c(0, 80))
    print(p)
    pdf(paste0(output_path,output_name),width=width,height=height)
        print(p) 
    dev.off() 
}

dims_list <- c(10, 15, 20,25,30,35,40,45,50)
k_param_list <- c(10, 20, 30,40,50,60,70,80)
res_list <- c(0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

# test the parameters
for (dims in dims_list) {
  for (k.param in k_param_list) {
    for (res in res_list) {
      print_umap_metabolic(inte_use,metabolic,dims, k.param, res,output_path)
    }
  }
}

dims_list <- c(10, 15, 20,25,30,35,40,45,50)
k_param_list <- c(10, 20, 30,40,50,60,70,80)
res_list <- c(0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

li=list()
index <- 1 

for (dims in dims_list) {
  for (k.param in k_param_list) {
    for (res in res_list) {
        embedding=readRDS(paste0(input_path,'META_celltype_PC_',dims,"_Nei_",k.param,"_Res_",res,"_dist",'.rds'))
        metadata=readRDS(paste0(input_path,'META_celltype_PC_',dims,"_Nei_",k.param,"_Res_",res,"_metadata",'.rds'))
        cell_dists <- dist(embedding, method = "euclidean")
        cluster_info <- as.numeric(as.character(metadata[, grepl((paste0("integrated_snn_res.",res)), colnames(metadata))]))
        si <- silhouette(cluster_info, cell_dists)
        df=data.frame(PC=dims,Nei=k.param,Res=res,Si=mean(si[,3]))
        li[[index]] <- df
        index <- index + 1
    }
  }
}

rbind=do.call(rbind,li) ## rbind the results
color=c('0.1'='#8DD3C7','0.2'='#BEBADA','0.3'='#FB8072','0.4'='#80B1D3',
        '0.5'='#FDB462','0.6'='#B3DE69','0.7'='#FCCDE5','0.8'='#D9D9D9',
        '0.9'='#BC80BD','1'='#CCEBC5')
width=4
height=3
lapply(as.list(dims_list),draw_si,data=rbind,celltype='fibroblasts',
       width=width,height=height,color=color,output_path=output_path)

rbind_use=rbind[rbind$PC>15 & rbind$Res>0.4 & rbind$Nei==80,]
rbind_use$Res=as.character(rbind_use$Res)
color=c('0.5'='#E47B7E','0.6'='#B2B31D','0.7'='#92C274','0.8'='#FBB2B4',
        '0.9'='#97CADC','1'='#984EA3')

width=6
height=3

options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=ggplot(rbind_use, aes(x=PC, y=Si,color=Res) )+
    geom_point() +
    geom_line(aes(group = Res)) +
    xlab("PC")+
    ylab("Silhouette")+
    theme_bw() + 
    scale_color_manual(values=color) + 
    theme(title = element_text(size = 15), text = element_text(size = 15)) + 
            theme_bw() + 
            geom_hline(yintercept = 0.17, linetype ="dashed") +
            labs(x='PC',y='Silhouette')+
            theme(panel.grid =element_blank()) + 
            theme(axis.text = element_blank()) + 
            theme(plot.title = element_text(size = 15),
            legend.title=element_text(size=15),
            axis.title=element_text(size=15),
            axis.text = element_text(size = 15),
            legend.text=element_text(size=15))
print(p)
pdf(paste0(output_path,'meta_final_res.pdf'),width=width,height=height)
    print(p) 
dev.off() 
