#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(stringr)
library(data.table)
library(Matrix)
Bonferronip0001 = c()

cancers = list.files("tmp/SJselection")
for(cont in 1:length(cancers)){
  
  df =readRDS( paste0("tmp/SJselection/", cancers[cont]))
  Bonferronip0001 = c(Bonferronip0001, as.character(df[[2]]$SJ[df[[2]]$Bonferroni <= 0.001]))
}

Bonferronip0001 = Bonferronip0001[str_split_fixed(Bonferronip0001, ":",2)[,1] %in% apply(cbind("chr", c(1:22,"X")),1, paste, collapse = "")]

saveRDS(unique(Bonferronip0001), "tmp/SJselection.RDS")