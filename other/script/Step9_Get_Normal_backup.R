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
dir.create("tmp/Normal")
SJandpvalues = readRDS("tmp/SJselection.RDS")
adjustment = c("bonferroni")
pthres = c(0.001)
SJs = str_split_fixed(SJandpvalues, "[-:]", 4)
SJstarts = apply(SJs[,1:2], 1, paste0, collapse = ":")
SJends =  apply(SJs[,c(1,3)], 1, paste0, collapse = ":")
chromosomes  =apply(cbind("chr", c(1:22,"X")), 1, paste0, collapse="")
for(nadjustment in 1:length(adjustment)){
    for(npthres in 1:length(pthres)){
        SJlabelselection = SJandpvalues
        SJlabelselectionsplit =    str_split_fixed(SJlabelselection, "[-:]", 4)
        SJselectionstarts = apply(SJlabelselectionsplit[,1:2], 1, paste0, collapse = ":")
        SJselectionends =  apply(SJlabelselectionsplit[,c(1,3)], 1, paste0, collapse = ":")
        dirs = list.files(samplefolder)
        cancersample =   readRDS(paste0(samplefolder, dirs[as.numeric(ndir)]))
        preselection = cancersample[[1]] [SJstarts %in% SJselectionstarts| SJends %in% SJselectionends,,drop=FALSE]
        SJsofcancer =  matrix(0,ncol = length(SJlabelselection), nrow = dim(preselection)[2])
        SJs = rownames(preselection)
        SJs = str_split_fixed(SJs, "[-:]", 4)
        SJstarts = apply(SJs[,1:2], 1, paste0, collapse = ":")
        SJends =  apply(SJs[,c(1,3)], 1, paste0, collapse = ":")
        for(selectedSJ in 1:length(SJlabelselection)){
            normalcounts  = preselection[SJstarts %in% SJselectionstarts[selectedSJ] | SJends %in% SJselectionends[selectedSJ],,drop=FALSE]
            SJsofcancer[,selectedSJ]  = colSums(normalcounts)
        }
        colnames(SJsofcancer) <- SJlabelselection
        rownames(SJsofcancer) <- colnames(preselection)
        write.table(SJsofcancer,paste0("tmp/Normal/Normal_",str_split_fixed(dirs[as.numeric(ndir)], "\\.", 2)[1], ".tsv") , sep = "\t", quote = F, row.names = T, col.names = T)
        #}
        gc()
    }
}