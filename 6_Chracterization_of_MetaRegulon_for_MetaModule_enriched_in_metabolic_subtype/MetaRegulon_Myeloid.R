#' Identifies MetaRegulon in myeloid cells.
#' @param integrated_object_path_myeloid Path to the myeloid cells object.
#' @param myeloid_DEG Path to the differentially expressed genes of myeloid cells.
#' @param myeloid_DEM Path to the differentially enriched MetaModule score of myeloid cells.
#' @param MetaModule Path to the metabolic reactions
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

integrated_object_path_myeloid<-args[1]
myeloid_DEG<-args[2]
myeloid_DEM<-args[3]
MetaModule=args[4]
MetaModule_info=args[5]
output_path=args[6]

## Find the marker genes for each metabolic state of the metacell object
myeloid<-readRDS(paste0(integrated_object_path_myeloid,'myeloid_integration_annotation_celltype_subset_metatype.rds'))
myeloid_DEM <- readRDS(myeloid_DEM)
myeloid_DEG <- readRDS(myeloid_DEG)

for (i in unique(myeloid_DEG$cluster)){
  df=myeloid_DEG[myeloid_DEG$cluster==i,]
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

for (metabolic_state in unique(myeloid_DEM$cluster)){
    ## set the parameters
    object=myeloid
    feature='metabolic_type'
    state=metabolic_state
    ## Users can use the differentially enriched MetaModule
    # interested_MM=MetaModule.markers[MetaModule.markers$cluster=='COL11A1+ CAF','gene']
    interested_MM=myeloid_DEM[myeloid_DEM$cluster==metabolic_state & myeloid_DEM$length>1,'gene']
    MM_list=MetaModule
    markers=myeloid_DEG
    lisa_file=paste0(output_path,'lisa/',metabolic_state,':marker.txt.lisa.tsv')
    ligand_target_matrix='./ref/ligand_target_matrix.rds'
    lr_network='./ref/lr_network.rds'
    sample_tech='scRNA'
    output_path=output_path
    RP_path='/fs/data/cistrome/RP/'
    file_name=metabolic_state
    cal_MetaRegulon(object,feature,state,interested_MM,MM_list,markers,lisa_file,ligand_target_matrix,lr_network,sample_tech,output_path,RP_path,file_name)
}