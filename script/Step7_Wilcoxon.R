#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(stringr)
library(data.table)
library(Matrix)
print(getwd())
print(args[1])
print(args[2])
print(args[3])

print("namewritten")
print("TCGA_mut_metadata_mutburden2023_1000_exskipupdated.tsv")
filesfolder  =args[1]
print(args[1])
cancers = list.files(filesfolder)
mutmetadata  =fread(args[2])
print(cancers)
cont=as.numeric(args[3])
print(cancers[cont])
dir.create("tmp/SJselection")
mutmetadata_mutated = mutmetadata$TCGA_ID[mutmetadata$relevance != "Not_relevant"]
mutmetadata_notmutated = mutmetadata$TCGA_ID[mutmetadata$relevance == "Not_relevant"]

#for(cont in 1:length(cancers)){
  print(cont)
  caso  = readRDS(paste0(filesfolder,cancers[cont]))
  samples = rownames(caso)
  samplesshort = substring(samples,1,16)
  mut= caso[samplesshort %in% mutmetadata_mutated ,]
  notmut= caso[samplesshort %in% mutmetadata_notmutated ,]
  pvalues = c()
  if(!is.null(dim(mut))){
    if(dim(mut)[1]!= 0 ){  
      for(nSJ in 1: dim(caso)[2]){
        
        pvalues = c(pvalues,  wilcox.test(notmut[,nSJ],mut[,nSJ],"less")$p.value)
      }
      
      
      df =list(data.frame(samples, samplesshort) ,data.frame(SJ = colnames(caso),  pvalues, Bonferroni = p.adjust(pvalues, method = "bonferroni"), FDR =p.adjust(pvalues, method = "fdr") ))
      
      saveRDS(df, paste0("tmp/SJselection/sampleandpvaluesperSJ", cancers[cont]))
    }
  }
  
#}

