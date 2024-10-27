#' Performs trajectory analysis of specific cell lineages using STREAM.
#' @param cell_lineage The cell lineage to analysis.
#' @param DEM Path to the differentially enriched MetaModules in specific cell lineages.
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
cell_lineage<-args[1]
DEM<-args[2]
output_path<-args[3]

DEM=paste0(DEM,cell_lineage,'/',cell_lineage)
output_path=paste0(output_path,cell_lineage,'/',cell_lineage)

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

plot_marker_data <- function(marker_data, output_name) {
    p <- ggplot(marker_data, aes(x = avg_log2FC, y = -log10(p_val_adj), shape = shape, color = metabolic_type)) +
        geom_point(size = marker_data$size) +
        scale_color_manual(values = color) +
        scale_shape_manual(values = c('1' = 17, '2' = 15, '3' = 19, '4' = 19)) +
        geom_vline(xintercept = 0.25, linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        labs(x = 'Log2FC', y = '-log10(Pvalue)') +
        theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), 
                           plot.title = element_text(size = 8), legend.title = element_text(size = 8),
                           axis.title = element_text(size = 8), legend.text = element_text(size = 8)) + NoLegend()
    
    pdf(paste0(output_path, output_name), width = 3, height = 3)
    print(p)
    dev.off()
}

# Helper function to set up markers
prepare_markers <- function(marker_data) {
    marker_data$metabolic_type <- factor(marker_data$metabolic_type, levels = names(color))
    marker_data$shape <- as.factor(ifelse(marker_data$metabolic_type %in% c('Glycolysis / Gluconeogenesis', 'Oxidative phosphorylation', 'Fatty acid oxidation'), '1', '2'))
    marker_data$size <- ifelse(marker_data$metabolic_type %in% c('Glycolysis / Gluconeogenesis', 'Oxidative phosphorylation', 'Fatty acid oxidation'), 5, 3)
    marker_data$p_val_adj[marker_data$p_val_adj == 0] <- 1e-300
    marker_data[order(marker_data$metabolic_type, decreasing = TRUE), ]
}

# Load data and apply functions
DEM <- readRDS(paste0(DEM, '_metatype_type_DEM.rds'))
if (cell_lineage == 'Fibroblast') {
    pufa.markers <- prepare_markers(DEM[DEM$cluster == 'PUFA', ])
    glycan.markers <- prepare_markers(DEM[DEM$cluster == 'GLYCAN', ])
    oxp.markers <- prepare_markers(DEM[DEM$cluster == 'IPM', ])

    plot_marker_data(pufa.markers, '_pufa_DEM_nolegend.pdf')
    plot_marker_data(glycan.markers, '_glycan_DEM.pdf')
    plot_marker_data(oxp.markers, '_oxp_DEM_nolegend.pdf')
} else if (cell_lineage == 'Myeloid') {
    arg.markers <- prepare_markers(DEM[DEM$cluster == 'ARG', ])
    fao.markers <- prepare_markers(DEM[DEM$cluster == 'FAO', ])
    gst.markers <- prepare_markers(DEM[DEM$cluster == 'GST', ])
    slm.markers <- prepare_markers(DEM[DEM$cluster == 'SLM', ])

    plot_marker_data(arg.markers, '_arg_DEM_nolegend.pdf')
    plot_marker_data(fao.markers, '_fao_DEM_nolegend.pdf')
    plot_marker_data(gst.markers, '_gst_DEM_nolegend.pdf')
    plot_marker_data(slm.markers, '_slm_DEM_nolegend.pdf')
}

