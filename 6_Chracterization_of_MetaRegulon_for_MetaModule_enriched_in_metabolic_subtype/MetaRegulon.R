#' Identifies MetaRegulon in cell lineages.
#' @param integrated_object_path Path to the cell lineages object.
#' @param DEG Path to the differentially expressed genes of cell lineages.
#' @param DEM Path to the differentially enriched MetaModule score of cell lineages.
#' @param MetaModule Path to the metabolic reactions.
#' @param MetaModule_info Path to the metabolic reaction information.
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
    library(MetroSCREEN)
})

args <- commandArgs(trailingOnly = TRUE)

cell_lineage<-args[1]
integrated_object_path<-args[2]
DEG<-args[3]
DEM<-args[4]
MetaModule=args[5]
MetaModule_info=args[6]
output_path=args[7]

integrated_object_path=paste0(integrated_object_path,cell_lineage,'/')
output_path=paste0(output_path,cell_lineage,'/')
DEG=paste0(DEG,cell_lineage,'/')
DEM=paste0(DEM,cell_lineage,'/')

## Find the marker genes for each metabolic state of the metacell object
seurat<-readRDS(paste0(integrated_object_path,'_integration_annotation_celltype_metatype.rds'))
DEM <- readRDS(paste0(DEM,'_metatype_type_DEM_nodup.rds'))
DEG <- readRDS(paste0(DEG,'_metatype_type_DEG.rds'))

for (i in unique(DEG$cluster)){
  df=DEG[DEG$cluster==i,]
  if (nrow(df)>500){
     genes=df[,'gene'][1:500]
  } else{
     genes=df[,'gene']
  }

  write.table(genes,paste0(output_path,'lisa/',i,':marker.txt'),
        sep='\t',
        quote=F,
        row.names=FALSE,
        col.names=FALSE)
}

MetaModule=readRDS(MetaModule)
MetaModule_info=readRDS(MetaModule_info)

for (metabolic_state in unique(DEM$cluster)){
    ## set the parameters
    object=seurat
    feature='metabolic_type'
    state=metabolic_state
    ## Users can use the differentially enriched MetaModule
    # interested_MM=MetaModule.markers[MetaModule.markers$cluster=='COL11A1+ CAF','gene']
    interested_MM=DEM[DEM$cluster==metabolic_state & DEM$length>1,'gene']
    MM_list=MetaModule
    markers=DEG
    lisa_file=paste0(output_path,'lisa/',metabolic_state,':marker.txt.lisa.tsv')
    ligand_target_matrix='./ref/ligand_target_matrix.rds'
    lr_network='./ref/lr_network.rds'
    sample_tech='scRNA'
    output_path=output_path
    RP_path='/fs/data/cistrome/RP/'
    file_name=metabolic_state
    cal_MetaRegulon(object,feature,state,interested_MM,MM_list,markers,lisa_file,ligand_target_matrix,lr_network,sample_tech,output_path,RP_path,file_name)
}