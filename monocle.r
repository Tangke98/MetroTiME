library(Seurat)
library(monocle)
library(ggplot2)
library(dplyr)
library(stringr)

RNA.integrated<-readRDS('/fs/home/tangke/mouse_nsc/E13_E17/E13_E17_filter_GE_annotation.rds')
DefaultAssay(RNA.integrated)='integrated'
RNA.integrated$cell=rownames(RNA.integrated@meta.data)

test=function(RNA.integrated,gene){
    # expression_matrix <- FindVariableFeatures(RNA.integrated, selection.method = "vst", nfeatures = featureNum)
    expression_matrix=RNA.integrated
    data <- as(as.matrix(expression_matrix@assays$integrated@data), 'sparseMatrix')
    pd <- new('AnnotatedDataFrame', data = expression_matrix@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fData)
    monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily = negbinomial.size())
    monocle_cds <- estimateSizeFactors(monocle_cds)
    monocle_cds <- estimateDispersions(monocle_cds)
    
    # gene_sle <-VariableFeatures(expression_matrix)
    gene_sle<-gene
    
    ordering_genes <-gene_sle
    monocle_cds <-
    setOrderingFilter(monocle_cds,
        ordering_genes = ordering_genes)
    monocle_cds <-
    reduceDimension(monocle_cds, 
                    max_components = 2)
    monocle_cds <- orderCells(monocle_cds) ## running step
    # saveRDS(monocle_cds,paste0('/fs/home/tangke/mouse_nsc/E13_E17/',featureNum,'_GE_monocle_trajectory.rds'))
    saveRDS(monocle_cds,paste0('/fs/home/tangke/mouse_nsc/E13_E17/','choosed_gene_GE_monocle_trajectory.rds'))
}



rgc=list(c('APOE','BMPR1B','CD38','CDK6','CDON','CRYAB','CTGF','DACH1','DBI','DDIT3','DMRT3',
           'DMRTA2','DNAJC1','EMX1','EMX2','ETNK1','FABP7','FBXO32','FGF1','FGFBP3','FGFR1','FGFR2',
           'FGFR3','FGFR4','GLI3','HES1','HES5','HOPX','KIAA1217','LDHA','LECT1','LRIG3','LYN','MAFF',
           'MFGE8','NAAA','NAMPT','NAPEPLD','NCKAP5','NES','NKAIN4','NR4A1','PALLD','PARD3B','PDGFD',
           'PDGFRB','PLCH1','PPARGC1A','ROBO1','ROBO2','ROBO3','ROBO4','SALL1','SAMD4A','SHISA2','SHROOM3',
           'SLC1A2','SLC1A3','SLIT2','SMARCA4','SNTG1','SOX1','SOX2','SOX6','SOX8','SOX9','SPARCL1 ','STOX1',
           'TBC1D1','TLE4','TMEM33','TMEM47','TNFRSF19','TSPAN12','VCAM1','VIM','BLBP','GLAST'))
rgc <- lapply(lapply(rgc, tolower),str_to_title)

quiescent_RGC=list(c('Mt3','Fabp7','Vim','Aldoc','HOPX','VCAM1'))
quiescent_RGC <- lapply(lapply(quiescent_RGC, tolower),str_to_title)


dividing_apical_RGC=list(c('Top2a','Cdk1','Cenpf','Cenpa','Cenpf','Ccnb1','Lockd','Ccnb2','Mdk','Pcna','Aldoc',
              'GP1BA','COL6A2','COL6A3','GP1BB','COL5A2','COL6A1','LAMA1','VWF','HSPG2','TNN','FN1','ITGA9','GP9',
           'COMP','IBSP','CD36','CHAD','GP5','VTN','THBS4','ITGA4','ITGA3','ITGA2B','ITGA7','ITGA5','COL5A1',
           'COL4A6','ITGA11','SV2C','COL2A1','COL3A1','COL4A1','AGRN','COL4A2','COL4A4','ITGB3','ITGB4','RELN',
           'ITGB5','ITGB6','ITGB7','LAMC2','ITGAV','ITGB1','LAMB2','SPP1','LAMB3','LAMC1','COL1A1','LAMA4','LAMA5',
           'LAMB1','COL1A2','ITGA10','GP6','ITGA8','LAMB4','TNR','CD47','SV2A','CD44','DAG1','TNXB','LAMA3','LAMA2',
           'SDC3','ITGB8','ITGA6','ITGA2','ITGA1','SV2B','TNC','COL11A1','LAMC3','COL11A2','HMMR','SDC2','SDC4','COL5A3',
           'THBS3','COL6A6','THBS2','SDC1','THBS1','APOE','BMPR1B','CD38','CDK6','CDON','CRYAB','CTGF','DACH1','DBI','DDIT3','DMRT3',
           'DMRTA2','DNAJC1','EMX1','EMX2','ETNK1','FABP7','FBXO32','FGF1','FGFBP3','FGFR1','FGFR2',
           'FGFR3','FGFR4','GLI3','HES1','HES5','HOPX','KIAA1217','LDHA','LECT1','LRIG3','LYN','MAFF',
           'MFGE8','NAAA','NAMPT','NAPEPLD','NCKAP5','NES','NKAIN4','NR4A1','PALLD','PARD3B','PDGFD',
           'PDGFRB','PLCH1','PPARGC1A','ROBO1','ROBO2','ROBO3','ROBO4','SALL1','SAMD4A','SHISA2','SHROOM3',
           'SLC1A2','SLC1A3','SLIT2','SMARCA4','SNTG1','SOX1','SOX2','SOX6','SOX8','SOX9','SPARCL1 ','STOX1',
           'TBC1D1','TLE4','TMEM33','TMEM47','TNFRSF19','TSPAN12','VCAM1','VIM'))
dividing_apical_RGC <- lapply(lapply(dividing_apical_RGC, tolower),str_to_title)


ipc=list(c('ASCL1','BCAN','EGFR','HES6','MYCN','NOTUM','OLIG1','OLIG2','S100A6','Eomes','BTG2'))
ipc <- lapply(lapply(ipc, tolower),str_to_title)
         
IPC_immature_neuron=list(c('Neurod6','Dll1','Neurod2','Dcx','Neurog2','Gadd45g','Sox4','Eomes'))
         
neuron=list(c('Stmn2','Ptprd','Neurod6','Tmsb10','Neurod2'))
         
obin=list(c('ARX','DCX','DLX1','DLX2','DLX5','DLX6','ETV1','ETV4','ETV5','GAD1','GAD2','GSX1',
'GSX2','HTR3A','MEIS1','MEIS2','PBX1','PBX2','PBX3','PROKR2','SLC32A1','SP8','SP9',
'TSHZ1','TUBB3','VAX1'))
obin <- lapply(lapply(obin, tolower),str_to_title)
         
pyn=list(c('BCL11A','BCL11B','BHLHE22','CDK4','CDK5','CUX1','CUX2','EMX1','EMX2','EOMES','FEZF2',
           'FOXG1','INSM1','LHX2','NDN','NEUROD1','NEUROD2','NEUROD6','NEUROG1','NEUROG2','NR2F1',
           'NR2F2','OTX1','PAX6','POU3F1','POU3F2','POU3F3','PPP1R17','RORB','SATB2','SOX5','SSTR2',
           'TBR1','TTF1','ZEB2'))
pyn <- lapply(lapply(pyn, tolower),str_to_title)
         
opc=list(c('APC','APOD','C1QL1','CHD7','CHD8','CLDN11','CNP','CSPG4','DLL1','DLL3','ERBB4','GALC',
           'MAG','MBP','MOG','NF1','NKX2-2','OLIG1','OLIG2','PCDH15','PDGFRA','PLP1','PPP1R14B',
           'QK','S100B','SALL3','SOX10','SOX6','SOX8'))
opc <- lapply(lapply(opc, tolower),str_to_title)
         
RGC_primed=list(c('Ntm','Prdm16','Fgfr1','Nckap5','Fgfr2','Adgrv1','Qk','Tnc','Gli3','Creb5'))
         
RGC_late_activated=list(c('Diaph3','Aspm','Sgo2a','Ect2','Smc4','Ncapd3','Cenpp','Cenpe','Creb5','Ncam1'))

bmp=list(c('ACVR1','ACVR1B','ACVR2A','ACVR2B','BMP2','BMP4','BMP5','BMP6','BMP7',
           'BMPR1A','BMPR1B','BMPR2','GDF5','GDF6','GDF7','ID1','ID2','ID3','ID4',
           'MSX1','NOG','SMAD1','SMAD2','SMAD3','SMAD4','SMAD5','SMAD6','SMAD7',
           'SMAD8','SMAD9','TTR'))
bmp <- lapply(lapply(bmp, tolower),str_to_title)

# mge=list(c('DLX1','DLX2','LHX6','LHX8','MAF','MAFB','NPY','SST'))
notch=list(c('HES1','HES5','HES6','HEY1','JAG1','NOTCH1','NOTCH2','NOTCH3','RBPJ'))
notch <- lapply(lapply(notch, tolower),str_to_title)

wnt<-list(c('JUN','LRP5','LRP6','PPP3R2','SFRP2','SFRP1','PPP3CC','VANGL1','PPP3R1','FZD1','FZD4','APC2','FZD6','FZD7','SENP2','FZD8',
            'LEF1','CREBBP','FZD9','PRICKLE1','CTBP2','ROCK1','CTBP1','WNT9B','WNT9A','CTNNBIP1','DAAM2','TBL1XR1','MMP7','CER1','MAP3K7',
            'VANGL2','WNT2B','WNT11','WNT10B','DKK2','SKP1P2','CHP2','AXIN1','AXIN2','DKK4','NFAT5','MYC','SOX17','CSNK2A1','CSNK2A2',
            'NFATC4','CSNK1A1','NFATC3','CSNK1E','BTRC','PRKX','SKP1','FBXW11','RBX1','CSNK2B','SIAH1','TBL1Y','WNT5B','CCND1','CAMK2A',
            'NLK','CAMK2B','CAMK2D','CAMK2G','PRKACA','APC','PRKACB','PRKACG','WNT16','DAAM1','CHD8','FRAT1','CACYBP','CCND2','NFATC2',
            'NFATC1','CCND3','PLCB2','PLCB1','CSNK1A1L','PRKCB','PLCB3','PRKCA','PLCB4','WIF1','PRICKLE2','PORCN','RHOA','FRAT2','PRKCG',
            'MAPK9','MAPK10','WNT3A','DVL3','RAC2','DVL2','RAC3','FZD3','DKK1','CXXC4','DVL1','FOSL1','CUL1','WNT10A','WNT4','SMAD3','TCF7',
            'SMAD4','RAC1','TCF7L2','SMAD2','WNT1','MAPK8','EP300','WNT7A','GSK3B','WNT7B','PSEN1','WNT8A','WNT8B','WNT2','WNT3','WNT5A','WNT6',
            'CTNNB1','PPP2CB','PPP2CA','PPP2R1A','TBL1X','PPP2R1B','ROCK2','NKD1','FZD10','FZD5','NKD2','TCF7L1','RUVBL1','PPARD','PPP3CB','TP53',
            'PPP3CA','PPP2R5A','PPP2R5E','PPP2R5D','PPP2R5C','PPP2R5B','FZD2','SFRP5','SFRP4','CHP1'))
wnt <- lapply(lapply(wnt, tolower),str_to_title)

shh=list(c('GAS1','GLI1','GLI2','GLI3','HHIP','PTCH1','SMO','SUFU'))
shh <- lapply(lapply(shh, tolower),str_to_title)

tgfb=list(c('TFDP1','NOG','TNF','GDF7','INHBB','INHBC','COMP','INHBA','THBS4','RHOA','CREBBP','ROCK1','ID1','ID2','RPS6KB1','RPS6KB2',
            'CUL1','SKP1P2','ID4','SMAD3','MAPK3','RBL2','SMAD4','RBL1','NODAL','SMAD1','MYC','SMAD2','MAPK1','SMURF2','SMURF1','EP300',
            'BMP8A','GDF5','SKP1','CHRD','TGFB2','TGFB1','IFNG','CDKN2B','PPP2CB','PPP2CA','PPP2R1A','ID3','SMAD5','RBX1','FST','PITX2',
            'PPP2R1B','TGFBR2','AMHR2','LTBP1','LEFTY1','AMH','TGFBR1','SMAD9','LEFTY2','SMAD7','ROCK2','TGFB3','SMAD6','BMPR2','GDF6',
            'BMPR1A','BMPR1B','ACVRL1','ACVR2B','ACVR2A','ACVR1','BMP4','E2F5','BMP2','ACVR1C','E2F4','SP1','BMP7','BMP8B','ZFYVE9',
            'BMP5','BMP6','ZFYVE16','THBS3','INHBE','THBS2','DCN','THBS1'))
tgfb <- lapply(lapply(tgfb, tolower),str_to_title)

ecm=list(c('GP1BA','COL6A2','COL6A3','GP1BB','COL5A2','COL6A1','LAMA1','VWF','HSPG2','TNN','FN1','ITGA9','GP9',
           'COMP','IBSP','CD36','CHAD','GP5','VTN','THBS4','ITGA4','ITGA3','ITGA2B','ITGA7','ITGA5','COL5A1',
           'COL4A6','ITGA11','SV2C','COL2A1','COL3A1','COL4A1','AGRN','COL4A2','COL4A4','ITGB3','ITGB4','RELN',
           'ITGB5','ITGB6','ITGB7','LAMC2','ITGAV','ITGB1','LAMB2','SPP1','LAMB3','LAMC1','COL1A1','LAMA4','LAMA5',
           'LAMB1','COL1A2','ITGA10','GP6','ITGA8','LAMB4','TNR','CD47','SV2A','CD44','DAG1','TNXB','LAMA3','LAMA2',
           'SDC3','ITGB8','ITGA6','ITGA2','ITGA1','SV2B','TNC','COL11A1','LAMC3','COL11A2','HMMR','SDC2','SDC4','COL5A3',
           'THBS3','COL6A6','THBS2','SDC1','THBS1'))
ecm <- lapply(lapply(ecm, tolower),str_to_title)

cell_cycle=list(c('CDC16','CDC7','CDC45','GADD45B','DBF4','ANAPC1','CREBBP','MDM2','ABL1',
                  'SMC1B','SKP1P2','GADD45G','ATM','ATR','ANAPC7','RBL2','ANAPC5','RBL1',
                  'MYC','CDC14B','SMC1A','CDC14A','SKP1','TGFB2','TGFB1','GADD45A','STAG1',
                  'PLK1','RBX1','STAG2','TGFB3','MCM4','ORC6','CCND1','MAD1L1','MCM3','MCM6',
                  'MCM5','YWHAB','CCNA2','MCM7','BUB1','CHEK1','WEE2','MCM2','PTTG1','CDC27',
                  'CDC25B','CDC25C','CDC25A','CDC6','CDC20','BUB3','YWHAZ','CCND2','YWHAH','CCNB1',
                  'YWHAG','YWHAE','CCNE1','ORC3','CCND3','PTTG2','SFN','E2F1','TFDP1','ZBTB17','CDK1',
                  'ESPL1','ANAPC10','RAD21','BUB1B','ANAPC11','RB1','SKP2','CUL1','SMAD3','SMAD4','ANAPC2',
                  'TFDP2','PRKDC','MAD2L1','ANAPC4','SMAD2','YWHAQ','CHEK2','CDC23','EP300','GSK3B','CDKN2A',
                  'CDKN1C','CDKN1B','CDKN1A','CDKN2D','CCNA1','CDKN2B','CDKN2C','FZR1','SMC3','ANAPC13','PCNA',
                  'TTK','PKMYT1','CDK2','CDC26','E2F5','CDK4','WEE1','E2F4','E2F3','TP53','E2F2','ORC1','ORC2',
                  'CCNE2','CDK6','ORC4','CCNB2','CDK7','MAD2L2','ORC5','HDAC1','HDAC2','CCNH','CCNB3'))
cell_cycle <- lapply(lapply(cell_cycle, tolower),str_to_title)

cell_cycle_s<-list(c('MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7','DTL','PRIM1','UHRF1',
                     'MLF1IP','HELLS','RFC2','RPA2','NASP','RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2',
                     'ATAD2','RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1','CLSPN','POLA1',
                     'CHAF1B','BRIP1','E2F8'))
cell_cycle_s <- lapply(lapply(cell_cycle_s, tolower),str_to_title)

cell_cycle_g2m<-list(c('HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2','NUF2','CKS1B','MKI67','TMPO',
                     'CENPF','TACC3','FAM64A','SMC4','CCNB2','CKAP2L','CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1',
                     'KIF20B','HJURP','CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1','NCAPD2','DLGAP5','CDCA2','CDCA8',
                     'ECT2','KIF23','HMMR','AURKA','PSRC1','ANLN','LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA'))
cell_cycle_g2m <- lapply(lapply(cell_cycle_g2m, tolower),str_to_title)

centromere=list(c('CENPA','CENPH','CENPL','CENPM','CENPN','CENPQ','CENPW','INCENP'))
centromere <- lapply(lapply(centromere, tolower),str_to_title)

# ep=list(c('AMOTL2','ANXA2','ATF3','BTG2','COL4A4','CRYAB','DUSP1','ENKUR','FAM13A',
#           'FBXO32','FOXB','FOXJ1','GPC4','IFIM3','ITIH5','JUNB','LGALS3','MAFF','NR4A1',
#           'OGN','PBXIP1','PDLIM3','PGX3','PROM1','RHOB','S100A11','SPARC','TMEM47','TPM1',
#           'TSPAN15','TTC12','ZFP36','ZFP36L1'))
fgf=list(c('FGF1','FGF10','FGF12','FGF13','FGF17','FGF18','FGF19','FGF2','FGF8','FGFBP3',
           'FGFR1','FGFR2','FGFR3','FGFR4','SPRY1','SPRY2'))
fgf <- lapply(lapply(fgf, tolower),str_to_title)

kinase=list(c('ARAF','BRAF','EGR1','EGR3','ETV1','ETV4','ETV5','HRAS','KRAS','MAP2K1',
             'MAP2K2','MAPK1','MAPK3','NRAS','RAF1'))
kinase <- lapply(lapply(kinase, tolower),str_to_title)
gene=c(rgc,quiescent_RGC,dividing_apical_RGC,ipc,
 IPC_immature_neuron,neuron,obin,pyn,opc,
 RGC_primed,RGC_late_activated)
