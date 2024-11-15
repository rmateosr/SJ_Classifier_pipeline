#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
list.of.packages <- c( "dplyr", "stringr","data.table","ggplot2", "Matrix" )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
library(dplyr)
library(stringr)
library(data.table)
library(Matrix)
print("hola")
print(args)
print(getwd())
samplefolder= args[1]
ndir= as.numeric(args[2])
print(ndir)
dir.create("tmp/Altered")
SJandpvalues = readRDS("tmp/SJselection.RDS")
adjustment = c("bonferroni")
pthres = c(0.001)
SJs = (SJandpvalues)
chromosomes  =apply(cbind("chr", c(1:22,"X")), 1, paste0, collapse="")
for(nadjustment in 1:length(adjustment)){
    for(npthres in 1:length(pthres)){
        SJlabelselection = SJandpvalues
        dirs = list.files(samplefolder)
        cancersample =   readRDS(paste0(samplefolder, dirs[as.numeric(ndir)]))
        SJs = rownames(cancersample[[1]])
        SJsofcancer =  matrix(0,ncol = length(SJlabelselection), nrow = dim(cancersample[[1]])[2])
        cancersample = (as.matrix(cancersample[[1]][SJs %in% SJlabelselection,]))
        SJs = rownames(cancersample)
        for(selectedSJ in 1:length(SJlabelselection)){
            SJfound = cancersample[SJs %in% SJlabelselection[selectedSJ],,drop=FALSE]
            if(dim(SJfound)[1] > 0){
            SJsofcancer[,selectedSJ] = SJfound
        }
    }
    colnames(SJsofcancer) <- SJlabelselection
    rownames(SJsofcancer) <- colnames(cancersample)
    write.table(SJsofcancer,paste0("tmp/Altered/Altered",str_split_fixed(dirs[as.numeric(ndir)], "\\.", 2)[1], ".tsv")  , sep = "\t", quote = F, row.names = T, col.names = T)
    gc()
    }
}