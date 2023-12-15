library(stringr)
library(RColorBrewer)
library(foreach)
library(ggplot2)
library(gridExtra)
library(data.table)
library(dplyr)
library(beepr)
library(parallel)
metadatafile = args[1]
pthresholds =  "p0001"
PANCANCER = fread("sup_files/mc3.v0.2.8.PUBLIC.maf", header = T,  stringsAsFactors = F, sep = "\t")
PANCANCER = data.frame(PANCANCER, stringsAsFactors = F)
PANCANCER = PANCANCER[PANCANCER$Variant_Classification != "3'UTR",]
PANCANCER = PANCANCER[PANCANCER$Variant_Classification != "5'UTR",]
PANCANCER = PANCANCER[PANCANCER$Variant_Classification != "Intron",]
PANCANCER = PANCANCER[PANCANCER$Variant_Classification != "Silent",]
sampleID = str_split_fixed(PANCANCER$Tumor_Sample_Barcode, "-", 4)
sampleID = apply(sampleID[,1:3], 1, paste0, collapse = "-")
PANCANCER$TCGA_ID_shorten = sampleID

npvalue = 1
PANCANCERNFE2L2 =unique( data.frame(NRF2 = PANCANCER$Hugo_Symbol == "NFE2L2" , PANCANCER[c("TCGA_ID_shorten", "Start_Position", "End_Position", "Protein_position", "Variant_Classification")]))
PANCANCERNFE2L2$Variant_Classification[PANCANCERNFE2L2$Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del")] <- "Frame_Shift_InDel"
PANCANCERNFE2L2$Variant_Classification[PANCANCERNFE2L2$Variant_Classification %in% c("In_Frame_Ins", "In_Frame_Del")] <- "In_Frame_InDel"
PANCANCERNFE2L2$Variant_Classification  = gsub("_", " ", PANCANCERNFE2L2$Variant_Classification)
saveRDS(PANCANCERNFE2L2,"sup_files/mc3.v0.2.8.PUBLIC_NRF2.maf.RDS")