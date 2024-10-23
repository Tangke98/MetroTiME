#' Performs trajectory analysis of myeloid cells using STREAM.
#' @param myeloid_DEM Path to the differentially enriched MetaModules in myeloid cells.
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
    suppressMessages(library(ROGUE))
    suppressMessages(library(tidyverse))
    library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

myeloid_DEM<-args[1]
output_path<-args[2]

color=c(    
'Arachidonic acid metabolism'='#6BAED6',
'Leukotriene metabolism'='#9ECAE1',
'Prostaglandin biosynthesis'='#DEEBF7',
'Eicosanoid metabolism'='#9ECAE1',
'Sphingolipid metabolism'='#1965B0',
'Glycosphingolipid metabolism'='#1965B0',
'Cysteine and methionine metabolism'='#9ECAE1',
    
'Chondroitin / heparan sulfate biosynthesis'='#41AB5D',
'Chondroitin sulfate degradation'='#74C476',
'Heparan sulfate degradation'='#A1D99B',
'Keratan sulfate degradation'='#C7E9C0',
'Keratan sulfate biosynthesis'='#E5F5E0',
'Pentose phosphate pathway'='#A1D99B',
'Fructose and mannose metabolism'='#C7E9C0',
    
'Glycine, serine and threonine metabolism'='#CA131F',
'Arginine and proline metabolism'='#DB5DB8',
    
'Glycolysis / Gluconeogenesis'='#FB8072',
'Oxidative phosphorylation'='#BC80BD',
'Fatty acid oxidation'='#FDB462',
    
'Glycerophospholipid metabolism'='#EFF3FF',
          
'Nicotinate and nicotinamide metabolism'='#EFF3FF',
'O-glycan metabolism'='#EFF3FF',

'Starch and sucrose metabolism'='#EFF3FF',
'Biopterin metabolism'='#EFF3FF',
'Transport reactions'='#EFF3FF',
'Xenobiotics metabolism'='#EFF3FF',
'Alanine, aspartate and glutamate metabolism'='#EFF3FF',
'Steroid metabolism'='#EFF3FF',
'N-glycan metabolism'='#EFF3FF',
'Inositol phosphate metabolism'='#EFF3FF',
'Porphyrin metabolism'='#EFF3FF',
'Androgen metabolism'='#EFF3FF',
'Bile acid biosynthesis'='#EFF3FF',
'Miscellaneous'='#EFF3FF',
'Linoleate metabolism'='#EFF3FF',
'Omega-3 fatty acid metabolism'='#EFF3FF',
'Tyrosine metabolism'='#EFF3FF',
'Valine, leucine, and isoleucine metabolism'='#EFF3FF',
'Serotonin and melatonin biosynthesis'='#EFF3FF',
'Vitamin B6 metabolism'='#EFF3FF',
'Isolated'='#EFF3FF',
'Amino sugar and nucleotide sugar metabolism'='#EFF3FF',
'Pyruvate metabolism'='#EFF3FF',
'Glutathione metabolism'='#EFF3FF',
'Protein assembly'='#EFF3FF',
'Galactose metabolism'='#EFF3FF',
'Glycosphingolipid biosynthesis-ganglio series'='#EFF3FF',
'Butanoate metabolism'='#EFF3FF',
'Drug metabolism'='#EFF3FF',
'Tryptophan metabolism'='#EFF3FF',
'Purine metabolism'='#EFF3FF',
'Pyrimidine metabolism'='#EFF3FF',
'Glycosphingolipid biosynthesis-globo series'='#EFF3FF',
'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism'='#EFF3FF',
'Glucocorticoid biosynthesis'='#EFF3FF',
'Nucleotide metabolism'='#EFF3FF',
'Retinol metabolism'='#EFF3FF',
'Estrogen metabolism'='#EFF3FF',
'Phenylalanine, tyrosine and tryptophan biosynthesis'='#EFF3FF',
'Ether lipid metabolism'='#EFF3FF',
'Cholesterol biosynthesis 1 (Bloch pathway)'='#EFF3FF',
'Metabolism of other amino acids'='#EFF3FF',
'Glycosphingolipid biosynthesis-lacto and neolacto series'='#EFF3FF',
'Folate metabolism'='#EFF3FF',
'Phosphatidylinositol phosphate metabolism'='#EFF3FF',
'Cholesterol metabolism'='#EFF3FF',
'Vitamin D metabolism'='#EFF3FF',
'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis'='#EFF3FF',
'Lysine metabolism'='#EFF3FF',
'Sulfur metabolism'='#EFF3FF',
'Pentose and glucuronate interconversions'='#EFF3FF',
'Riboflavin metabolism'='#EFF3FF',
'Histidine metabolism'='#EFF3FF',
'Fatty acid elongation (even-chain)'='#EFF3FF',
'Beta oxidation of poly-unsaturated fatty acids (mitochondrial)'='#EFF3FF',
'Fatty acid biosynthesis'='#EFF3FF',
'Terpenoid backbone biosynthesis'='#EFF3FF',
'Protein degradation'='#EFF3FF',
'Beta-alanine metabolism'='#EFF3FF',
'Fatty acid desaturation (even-chain)'='#EFF3FF',
'Fatty acid desaturation (odd-chain)'='#EFF3FF',
'Fatty acid biosynthesis (unsaturated)'='#EFF3FF',
'Fatty acid elongation (odd-chain)'='#EFF3FF',
'Omega-6 fatty acid metabolism'='#EFF3FF',
'Protein modification'='#EFF3FF',
'Ascorbate and aldarate metabolism'='#EFF3FF',
'Pantothenate and CoA biosynthesis',
'Carnitine shuttle (endoplasmic reticular)'='#EFF3FF',
'Carnitine shuttle (peroxisomal)'='#EFF3FF',
'Beta oxidation of di-unsaturated fatty acids (n-6) (peroxisomal)'='#EFF3FF',
'Fatty acid biosynthesis (even-chain)'='#EFF3FF',
'Fatty acid biosynthesis (odd-chain)'='#EFF3FF',
'Aminoacyl-tRNA biosynthesis'='#EFF3FF',
'Carnitine shuttle (mitochondrial)'='#EFF3FF',
'Heme degradation'='#EFF3FF',
'Blood group biosynthesis'='#EFF3FF',
'Vitamin B2 metabolism'='#EFF3FF',
'Acylglycerides metabolism'='#EFF3FF',
'Pool reactions'='#EFF3FF',
'Pantothenate and CoA biosynthesis'='#EFF3FF',
'Bile acid recycling'='#EFF3FF',
'Beta oxidation of phytanic acid (peroxisomal)'='#EFF3FF',
'Thiamine metabolism'='#EFF3FF',
'Urea cycle'='#EFF3FF',
'Propanoate metabolism'='#EFF3FF',
'Cholesterol biosynthesis 2'='#EFF3FF',
        'Cholesterol biosynthesis 3 (Kandustch-Russell pathway)'='#EFF3FF',
        'Lipoic acid metabolism'='#EFF3FF',
       'Formation and hydrolysis of cholesterol esters'='#EFF3FF',
        'Beta oxidation of di-unsaturated fatty acids (n-6) (mitochondrial)'='#EFF3FF',
        'Formation and hydrolysis'='#EFF3FF',
        'Glycerolipid metabolism'='#EFF3FF',
'Fatty acid activation (cytosolic)'='#EFF3FF',
'Beta oxidation of branched-chain fatty acids (mitochondrial)'='#EFF3FF',
'Beta oxidation of even-chain fatty acids (mitochondrial)'='#EFF3FF')

myeloid_DEM<-readRDS(myeloid_DEM)

arg.markers=myeloid_DEM[myeloid_DEM$cluster=='ARG',]
fao.markers=myeloid_DEM[myeloid_DEM$cluster=='FAO',]
gst.markers=myeloid_DEM[myeloid_DEM$cluster=='GST',]
slm.markers=myeloid_DEM[myeloid_DEM$cluster=='SLM',]

arg.markers$p_val_adj[arg.markers$p_val_adj==0]=1e-300
fao.markers$p_val_adj[fao.markers$p_val_adj==0]=1e-300
gst.markers$p_val_adj[gst.markers$p_val_adj==0]=1e-300
slm.markers$p_val_adj[slm.markers$p_val_adj==0]=1e-300

arg.markers$metabolic_type=factor(arg.markers$metabolic_type,levels = names(color))
arg.markers$shape=as.factor(ifelse(arg.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),'1',
                                    ifelse(arg.markers$metabolic_type %in% c('Arginine and proline metabolism','Glycine, serine and threonine metabolism',
                                                                   'Chondroitin / heparan sulfate biosynthesis','Chondroitin sulfate degradation',
                                                                    'Heparan sulfate degradation','Keratan sulfate degradation',
                                                                    'Keratan sulfate biosynthesis','Pentose phosphate pathway',
                                                                    'Fructose and mannose metabolism'),'2',
                                           ifelse(arg.markers$metabolic_type %in% c('Arachidonic acid metabolism','Leukotriene metabolism',
                                                                                   'Prostaglandin biosynthesis','Eicosanoid metabolism',
                                                                          'Cysteine and methionine metabolism','Sphingolipid metabolism','Glycosphingolipid metabolism'),'3','4'))))
arg.markers$size=ifelse(arg.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),5,5)
arg.markers=arg.markers[order(arg.markers$metabolic_type),]
arg.markers=arg.markers[order(arg.markers$metabolic_type,decreasing = TRUE),]
output_name='myeloid_arg_DEM_nolegend.pdf'
width=3
height=3
options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=ggplot(arg.markers,aes(x=avg_log2FC,y=-log10(p_val_adj),shape=shape,color = metabolic_type)) + 
  xlab("log2FC") + 
  geom_point(size = arg.markers$size) + 
  scale_shape_manual(values = c('1'=17,'2'=15,'3'=19,'4'=19) )+  
  scale_color_manual(values=color) + 
  geom_vline(xintercept = c(0.25), linetype ="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype ="dashed") + 
  theme(title = element_text(size = 15), text = element_text(size = 15)) + 
    theme_bw() + #去除背景色
            labs(x='Log2FC',y='-log10(Pvalue)')+
            theme(panel.grid =element_blank()) + #去除网格线
            theme(axis.text = element_blank()) + #去掉坐标轴
            #         scale_color_manual(values=color)+
            theme(plot.title = element_text(size = 8),
            legend.title=element_text(size=8),
            axis.title=element_text(size=8),
            axis.text = element_text(size = 8),
            legend.text=element_text(size=8))+NoLegend()

print(p)
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p) 
dev.off()


fao.markers$metabolic_type=factor(fao.markers$metabolic_type,levels = names(color))
fao.markers$shape=as.factor(ifelse(fao.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),'1',
                                    ifelse(fao.markers$metabolic_type %in% c('Arginine and proline metabolism','Glycine, serine and threonine metabolism',
                                                                   'Chondroitin / heparan sulfate biosynthesis','Chondroitin sulfate degradation',
                                                                    'Heparan sulfate degradation','Keratan sulfate degradation',
                                                                    'Keratan sulfate biosynthesis','Pentose phosphate pathway',
                                                                    'Fructose and mannose metabolism'),'2',
                                           ifelse(fao.markers$metabolic_type %in% c('Arachidonic acid metabolism','Leukotriene metabolism',
                                                                                   'Prostaglandin biosynthesis','Eicosanoid metabolism',
                                                                          'Cysteine and methionine metabolism','Sphingolipid metabolism'),'3','4'))))
fao.markers$size=ifelse(fao.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),5,5)
fao.markers=fao.markers[order(fao.markers$metabolic_type),]
fao.markers=fao.markers[order(fao.markers$metabolic_type,decreasing = TRUE),]
output_name='myeloid_fao_DEM_nolegend.pdf'
width=3
height=3
options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=ggplot(fao.markers,aes(x=avg_log2FC,y=-log10(p_val_adj),shape=shape,color = metabolic_type)) + 
  xlab("log2FC") + 
  geom_point(size = fao.markers$size) + 
  scale_color_manual(values=color) + 
 scale_shape_manual(values = c('1'=17,'2'=15,'3'=19,'4'=19) )+ 
  geom_vline(xintercept = c(0.25), linetype ="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype ="dashed") + 
  theme(title = element_text(size = 15), text = element_text(size = 15)) + 
    theme_bw() + #去除背景色
            labs(x='Log2FC',y='-log10(Pvalue)')+
            theme(panel.grid =element_blank()) + #去除网格线
            theme(axis.text = element_blank()) + #去掉坐标轴
            #         scale_color_manual(values=color)+
            theme(plot.title = element_text(size = 8),
            legend.title=element_text(size=8),
            axis.title=element_text(size=8),
            axis.text = element_text(size = 8),
            legend.text=element_text(size=8))+NoLegend()

print(p)
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p) 
dev.off()


gst.markers$metabolic_type=factor(gst.markers$metabolic_type,levels = names(color))
gst.markers$shape=as.factor(ifelse(gst.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),'1',
                                    ifelse(gst.markers$metabolic_type %in% c('Arginine and proline metabolism','Glycine, serine and threonine metabolism',
                                                                   'Chondroitin / heparan sulfate biosynthesis','Chondroitin sulfate degradation',
                                                                    'Heparan sulfate degradation','Keratan sulfate degradation',
                                                                    'Keratan sulfate biosynthesis','Pentose phosphate pathway',
                                                                    'Fructose and mannose metabolism'),'2',
                                           ifelse(gst.markers$metabolic_type %in% c('Arachidonic acid metabolism','Leukotriene metabolism',
                                                                                   'Prostaglandin biosynthesis','Eicosanoid metabolism',
                                                                          'Cysteine and methionine metabolism','Sphingolipid metabolism','Glycosphingolipid metabolism'),'3','4'))))
gst.markers$size=ifelse(gst.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),5,5)
gst.markers=gst.markers[order(gst.markers$metabolic_type),]
gst.markers=gst.markers[order(gst.markers$metabolic_type,decreasing = TRUE),]
gst.markers$p_val_adj[gst.markers$p_val_adj==0]=1e-300
output_name='myeloid_gst_DEM_nolegend.pdf'
width=3
height=3
options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=ggplot(gst.markers,aes(x=avg_log2FC,y=-log10(p_val_adj),shape=shape,color = metabolic_type)) + 
  xlab("log2FC") + 
  geom_point(size = gst.markers$size) + 
  scale_color_manual(values=color) + 
  scale_shape_manual(values = c('1'=17,'2'=15,'3'=19,'4'=19) )+ 
  geom_vline(xintercept = c(0.25), linetype ="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype ="dashed") + 
  theme(title = element_text(size = 15), text = element_text(size = 15)) + 
    theme_bw() + #去除背景色
            labs(x='Log2FC',y='-log10(Pvalue)')+
            theme(panel.grid =element_blank()) + #去除网格线
            theme(axis.text = element_blank()) + #去掉坐标轴
            #         scale_color_manual(values=color)+
            theme(plot.title = element_text(size = 8),
            legend.title=element_text(size=8),
            axis.title=element_text(size=8),
            axis.text = element_text(size = 8),
            legend.text=element_text(size=8))+NoLegend()

print(p)
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p) 
dev.off()


slm.markers$metabolic_type=factor(slm.markers$metabolic_type,levels = names(color))
slm.markers$shape=as.factor(ifelse(slm.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),'1',
                                    ifelse(slm.markers$metabolic_type %in% c('Arginine and proline metabolism','Glycine, serine and threonine metabolism',
                                                                   'Chondroitin / heparan sulfate biosynthesis','Chondroitin sulfate degradation',
                                                                    'Heparan sulfate degradation','Keratan sulfate degradation',
                                                                    'Keratan sulfate biosynthesis','Pentose phosphate pathway',
                                                                    'Fructose and mannose metabolism'),'2',
                                           ifelse(slm.markers$metabolic_type %in% c('Arachidonic acid metabolism','Leukotriene metabolism',
                                                                                   'Prostaglandin biosynthesis','Eicosanoid metabolism',
                                                                          'Cysteine and methionine metabolism','Sphingolipid metabolism'),'3','4'))))
slm.markers$size=ifelse(slm.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),5,5)
slm.markers=slm.markers[order(slm.markers$metabolic_type),]
slm.markers=slm.markers[order(slm.markers$metabolic_type,decreasing = TRUE),]
slm.markers$p_val_adj[slm.markers$p_val_adj==0]=1e-300
output_name='myeloid_slm_DEM_nolegend.pdf'
width=3
height=3
options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=ggplot(slm.markers,aes(x=avg_log2FC,y=-log10(p_val_adj),shape=shape,color = metabolic_type)) + 
  xlab("log2FC") + 
  geom_point(size = slm.markers$size) + 
  scale_color_manual(values=color) + 
  geom_vline(xintercept = c(0.25), linetype ="dashed") +
 scale_shape_manual(values = c('1'=17,'2'=15,'3'=19,'4'=19) )+  
  geom_hline(yintercept = -log10(0.05), linetype ="dashed") + 
  theme(title = element_text(size = 15), text = element_text(size = 15)) + 
    theme_bw() + #去除背景色
            labs(x='Log2FC',y='-log10(Pvalue)')+
            theme(panel.grid =element_blank()) + #去除网格线
            theme(axis.text = element_blank()) + #去掉坐标轴
            #         scale_color_manual(values=color)+
            theme(plot.title = element_text(size = 8),
            legend.title=element_text(size=8),
            axis.title=element_text(size=8),
            axis.text = element_text(size = 8),
            legend.text=element_text(size=8))+NoLegend()

print(p)
pdf(paste0(output_path,output_name),width=width,height=height)
    print(p) 
dev.off()