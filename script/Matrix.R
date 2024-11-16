#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
# List of required packages
required_packages <- c("stringr", "tools", "data.table")

# Function to install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# Install and load required packages
install_if_missing(required_packages)
library(stringr)
library(tools)
library(Matrix)
library(data.table)

alteredSJnames <- readRDS("resources/SJnames.RDS")
sample <- as.data.frame(fread(args[1]))
# Extract the file name without extension
file_name <- file_path_sans_ext(basename(args[1]))


rownames(sample) = sample[,1]
SJinfo = sample[,2]
sample = sample[,-(1:2)]
sampleinfo = colnames(sample)


alteredSJnames = alteredSJnames[alteredSJnames!= "Sample_name"]
alteredSJcoords =  str_split_fixed(alteredSJnames, "[:-]", 4)
alteredstarts = apply(cbind(alteredSJcoords[,1], ":", alteredSJcoords[,2]),1,paste0, collapse = "")
alteredends = apply(cbind(alteredSJcoords[,1], ":", alteredSJcoords[,3]),1,paste0, collapse = "")
notannotatedarrayfromjxn =  sample[SJinfo==0, ,drop = FALSE]
alteredSJappearing = as.matrix(sample[rownames(sample) %in% alteredSJnames & SJinfo==0,] )


AlteredthisProject = matrix(0, ncol = length(alteredSJnames), nrow = dim(notannotatedarrayfromjxn)[2])
colnames(AlteredthisProject) = alteredSJnames
rownames(AlteredthisProject) = colnames(notannotatedarrayfromjxn)
for(cont in 1:dim(AlteredthisProject)[2]){
  alteredSJcountsinthisproject=as.matrix(alteredSJappearing[rownames(alteredSJappearing) %in% alteredSJnames[cont],])
  if(dim(alteredSJcountsinthisproject)[1]> 0 & dim(alteredSJcountsinthisproject)[2]> 0){
    AlteredthisProject[,cont] =  t(alteredSJcountsinthisproject)
  }
}
annotatedarrayfromjxn =  sample[SJinfo==1, ,drop = FALSE]
if(dim(annotatedarrayfromjxn)[1] > 0){
  SJcoords =  str_split_fixed(rownames(annotatedarrayfromjxn), "[:-]", 4)
  SJstarts = apply(cbind(SJcoords[,1], ":", SJcoords[,2]),1,paste0, collapse = "")
  SJends = apply(cbind(SJcoords[,1], ":", SJcoords[,3]),1,paste0, collapse = "")
  
  whichstarts= SJstarts %in% alteredstarts
  arraystarts = annotatedarrayfromjxn[whichstarts,,drop  =FALSE]
  startscoordsofarraystarts = SJstarts[whichstarts]
  
  whichends= SJends %in% alteredends
  arrayends = annotatedarrayfromjxn[whichends,,drop  =FALSE]
  endscoordsofarrayends = SJends[whichends]
  
  
  startsthisProject = matrix(0, ncol = length(alteredstarts), nrow = dim(annotatedarrayfromjxn)[2])
  colnames(startsthisProject) = alteredSJnames
  rownames(startsthisProject) = colnames(annotatedarrayfromjxn)
  
  
  for(cont in 1:dim(startsthisProject)[2]){
    alteredSJstartcountsinthisproject=as.matrix(arraystarts[(startscoordsofarraystarts) %in% alteredstarts[cont],,drop  =FALSE])
    if(dim(alteredSJstartcountsinthisproject)[1]> 0 & dim(alteredSJstartcountsinthisproject)[2]> 0 ){
      startsthisProject[,cont] =  colSums(alteredSJstartcountsinthisproject)
    }
  }
  endsthisProject = matrix(0, ncol = length(alteredends), nrow = dim(annotatedarrayfromjxn)[2])
  colnames(endsthisProject) = alteredSJnames
  rownames(endsthisProject) = colnames(annotatedarrayfromjxn)
  
  
  for(cont in 1:dim(endsthisProject)[2]){
    alteredSJendcountsinthisproject=as.matrix(arrayends[(endscoordsofarrayends) %in% alteredends[cont],,drop=FALSE])
    if(dim(alteredSJendcountsinthisproject)[1]> 0 & dim(alteredSJendcountsinthisproject)[2]> 0){
      endsthisProject[,cont] =  colSums(alteredSJendcountsinthisproject)
    }
  }
  
}  else { 
  startsthisProject = 0
  endsthisProject = 0
}
Normal = startsthisProject + endsthisProject
Altered = AlteredthisProject
Total = Normal + Altered
dataset = list(Altered, Total,sampleinfo)
saveRDS(dataset, paste0("tmp/", file_name,  ".RDS"))
