#' Performs trajectory analysis of fibroblast cells using STREAM.
#' @param fibroblasts_DEM Path to the differentially enriched MetaModules in fibroblasts.
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

fibroblasts_DEM<-args[1]
output_path<-args[2]

color=c('Arachidonic acid metabolism'='#6BAED6',
'Leukotriene metabolism'='#9ECAE1',
'Prostaglandin biosynthesis'='#DEEBF7',
'Eicosanoid metabolism'='#9ECAE1',
        
'Chondroitin / heparan sulfate biosynthesis'='#41AB5D',
'Chondroitin sulfate degradation'='#74C476',
'Heparan sulfate degradation'='#A1D99B',
'Keratan sulfate degradation'='#C7E9C0',
'Keratan sulfate biosynthesis'='#E5F5E0',
        
'Glycolysis / Gluconeogenesis'='#FB8072',
'Oxidative phosphorylation'='#BC80BD',
'Fatty acid oxidation'='#FDB462',
        
'Nicotinate and nicotinamide metabolism'='#EFF3FF',
'O-glycan metabolism'='#EFF3FF',
'Arginine and proline metabolism'='#EFF3FF',
'Starch and sucrose metabolism'='#EFF3FF',
'Biopterin metabolism'='#EFF3FF',
'Transport reactions'='#EFF3FF',
'Sphingolipid metabolism'='#EFF3FF',
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
'Glycine, serine and threonine metabolism'='#EFF3FF',
'Glycosphingolipid biosynthesis-ganglio series'='#EFF3FF',
'Butanoate metabolism'='#EFF3FF',
'Drug metabolism'='#EFF3FF',
'Tryptophan metabolism'='#EFF3FF',
'Purine metabolism'='#EFF3FF',
'Pyrimidine metabolism'='#EFF3FF',
'Glycosphingolipid biosynthesis-globo series'='#EFF3FF',
'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism'='#EFF3FF',
'Glycosphingolipid metabolism'='#EFF3FF',
'Pentose phosphate pathway'='#EFF3FF',
'Glucocorticoid biosynthesis'='#EFF3FF',
'Glycerophospholipid metabolism'='#EFF3FF',
'Nucleotide metabolism'='#EFF3FF',
'Retinol metabolism'='#EFF3FF',
'Estrogen metabolism'='#EFF3FF',
'Phenylalanine, tyrosine and tryptophan biosynthesis'='#EFF3FF',
'Ether lipid metabolism'='#EFF3FF',
'Cholesterol biosynthesis 1 (Bloch pathway)'='#EFF3FF',
'Cysteine and methionine metabolism'='#EFF3FF',
'Metabolism of other amino acids'='#EFF3FF',
'Glycosphingolipid biosynthesis-lacto and neolacto series'='#EFF3FF',
'Folate metabolism'='#EFF3FF',
'Phosphatidylinositol phosphate metabolism'='#EFF3FF',
'Cholesterol metabolism'='#EFF3FF',
'Vitamin D metabolism'='#EFF3FF',
'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis'='#EFF3FF',
'Lysine metabolism'='#EFF3FF',
'Sulfur metabolism'='#EFF3FF',
'Fructose and mannose metabolism'='#EFF3FF',
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
        'Lipoic acid metabolism'='#EFF3FF')

fibroblasts_DEM<-readRDS(fibroblasts_DEM)
pufa.markers=fibroblasts_DEM[fibroblasts_DEM$cluster=='PUFA',]
glycan.markers=fibroblasts_DEM[fibroblasts_DEM$cluster=='GLYCAN',]
oxp.markers=fibroblasts_DEM[fibroblasts_DEM$cluster=='IPM',]

pufa.markers$p_val_adj[pufa.markers$p_val_adj==0]=1e-300
glycan.markers$p_val_adj[glycan.markers$p_val_adj==0]=1e-300
oxp.markers$p_val_adj[oxp.markers$p_val_adj==0]=1e-300

pufa.markers$metabolic_type=factor(pufa.markers$metabolic_type,levels = names(color))
pufa.markers$shape=as.factor(ifelse(pufa.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),19,
                                    ifelse(pufa.markers$metabolic_type %in% c('Arachidonic acid metabolism','Leukotriene metabolism','Prostaglandin biosynthesis','Eicosanoid metabolism'),22,17)))
pufa.markers$size=ifelse(pufa.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),5,3)
pufa.markers=pufa.markers[order(pufa.markers$metabolic_type),]
pufa.markers=pufa.markers[order(pufa.markers$metabolic_type,decreasing = TRUE),]
output_name='fibroblast_pufa_DEM_nolegend.pdf'
width=3
height=3
options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=ggplot(pufa.markers,aes(x=avg_log2FC,y=-log10(p_val_adj),shape=shape,color = metabolic_type)) + 
  xlab("log2FC") + 
  geom_point(size = pufa.markers$size) + 
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

glycan.markers$metabolic_type=factor(glycan.markers$metabolic_type,levels = names(color))
glycan.markers$shape=as.factor(ifelse(glycan.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),19,
                                    ifelse(glycan.markers$metabolic_type %in% c('Arachidonic acid metabolism','Leukotriene metabolism','Prostaglandin biosynthesis','Eicosanoid metabolism'),22,17)))
glycan.markers$size=ifelse(glycan.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),5,3)
glycan.markers=glycan.markers[order(glycan.markers$metabolic_type),]
glycan.markers=glycan.markers[order(glycan.markers$metabolic_type,decreasing = TRUE),]
output_name='fibroblast_glycan_DEM.pdf'
width=3
height=3
options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=ggplot(glycan.markers,aes(x=avg_log2FC,y=-log10(p_val_adj),shape=shape,color = metabolic_type,size=)) + 
  xlab("log2FC") + 
  geom_point(size = glycan.markers$size) + 
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

# oxp.markers=na.omit(oxp.markers)
oxp.markers$metabolic_type=factor(oxp.markers$metabolic_type,levels = names(color))
oxp.markers$shape=as.factor(ifelse(oxp.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),19,
                                    ifelse(oxp.markers$metabolic_type %in% c('Arachidonic acid metabolism','Leukotriene metabolism','Prostaglandin biosynthesis','Eicosanoid metabolism'),22,17)))
oxp.markers$size=ifelse(oxp.markers$metabolic_type %in% c('Glycolysis / Gluconeogenesis','Oxidative phosphorylation','Fatty acid oxidation'),5,3)
oxp.markers=oxp.markers[order(oxp.markers$metabolic_type),]
oxp.markers=oxp.markers[order(oxp.markers$metabolic_type,decreasing = TRUE),]
output_name='fibroblast_oxp_DEM_nolegend.pdf'
width=3
height=3
options(repr.plot.width = width, repr.plot.height = height,repr.plot.res = 100)
p=ggplot(oxp.markers,aes(x=avg_log2FC,y=-log10(p_val_adj),shape=shape,color = metabolic_type,size=)) + 
  xlab("log2FC") + 
  geom_point(size = oxp.markers$size) + 
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
