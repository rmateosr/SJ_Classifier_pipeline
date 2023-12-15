#!/usr/bin/Rscript
list.of.packages <- c( "dplyr", "stringr","data.table","ggplot2", "Matrix" )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
library(dplyr)
library(stringr)
library(data.table)
library(Matrix)
dir.create("tmp/TAN")
dirs = list.files("tmp/Normal")
ndir = 1
allcombined = data.frame(fread(paste0("tmp/Normal/", dirs[ndir])),check.names = F)
for(ndir in 2:length(dirs)){
    print(ndir)
    allcombined = rbind(allcombined, data.frame(fread(paste0("tmp/Normal/", dirs[ndir])),check.names = F))
}
colnames(allcombined)[1] = "Sample_name"
write.table(allcombined,paste0("tmp/TAN/Normal.tsv") , sep = "\t", quote = F, row.names = F, col.names = T)

dirs = list.files("tmp/Altered")
ndir = 1
allcombined = data.frame(fread(paste0("tmp/Altered/", dirs[ndir])),check.names = F)
for(ndir in 2:length(dirs)){
    print(ndir)
    allcombined = rbind(allcombined, data.frame(fread(paste0("tmp/Altered/", dirs[ndir])),check.names = F))
}
colnames(allcombined)[1] = "Sample_name"
write.table(allcombined,paste0("tmp/TAN/Altered.tsv") , sep = "\t", quote = F, row.names = F, col.names = T)


Altered = data.frame(fread(paste0("tmp/TAN/Altered.tsv")),check.names = F)
Normal = data.frame(fread(paste0("tmp/TAN/Normal.tsv")),check.names = F)
Total = cbind(Normal[,1] , Altered[,-1] + Normal[,-1])
colnames(Total)[1] = "Sample_name"
write.table(Total,paste0("tmp/TAN/Total.tsv") , sep = "\t", quote = F, row.names = F, col.names = T)

