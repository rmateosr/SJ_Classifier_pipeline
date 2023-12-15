#!/usr/bin/Rscript
set.seed(12)   
args = commandArgs(trailingOnly=TRUE)
library(stringr)
library(RColorBrewer)
library(foreach)
#library(ggsci)
library(ggplot2)
library(gridExtra)
#library(viridis)
#library(gplots)
library(data.table)
library(dplyr)
library(beepr)
library(parallel)
#set.seed(10) 
betabinomial = function(alphabeta, total, altered) { - sum(
    lgamma(total + 1) +
        lgamma(altered + alphabeta[1]) +
        lgamma(total - altered + alphabeta[2]) +
        lgamma(alphabeta[1] + alphabeta[2])-
        lgamma(altered+1) -
        lgamma(total - altered + 1 ) -
        lgamma(total + alphabeta[1] + alphabeta[2]) -
        lgamma(alphabeta[1]) - 
        lgamma(alphabeta[2])
    )
    
}


# Sjscorethreshold = 10
# scorediffthreshold=10
Sjscorethreshold = 10
scorediffthreshold= 10
c = 0.001


#compilation of all cancer types


metadatafile = args[1]
dataforbetabin  =data.frame(fread("tmp/TAN/Altered.tsv"), check.names = F)
totalreadsperSJPANCANCER=data.frame(fread("tmp/TAN/Total.tsv"), check.names = F)
mutmetadata  =data.frame(fread(metadatafile))
dataforbetabin$TCGA_ID = substring(dataforbetabin$Sample_name,1,16)
totalreadsperSJPANCANCER$TCGA_ID = dataforbetabin$TCGA_ID 
dataforbetabin = inner_join(dataforbetabin, mutmetadata)
dataforbetabin = dataforbetabin[dataforbetabin$relevance != "Ambiguous",]
dataforbetabin = dataforbetabin[,colnames(dataforbetabin)!= "Sample_name"]
totalreadsperSJPANCANCER = inner_join(totalreadsperSJPANCANCER, mutmetadata)

totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[totalreadsperSJPANCANCER$relevance != "Ambiguous",]
totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[,grep("chr", colnames(totalreadsperSJPANCANCER))]
dataforbetabin_SJ= dataforbetabin[,grep("chr", colnames(dataforbetabin))]


cancers =  unique(dataforbetabin$Cancer_Type) 
resultspercancer = data.frame()

totalreadsperSJPANCANCERonecancer = totalreadsperSJPANCANCER
dataforbetabinonecancer = dataforbetabin



totalreadsperSJPANCANCERonecancer = totalreadsperSJPANCANCER
dataforbetabinonecancer = dataforbetabin[,colnames(dataforbetabin)!= "Sample_name"]
if( sum(dataforbetabinonecancer$mutlabels != "None") > 1){
    
    
    dataforbetabinonecancer$mutornot =  dataforbetabinonecancer$mutlabels != "None"
    dataforbetabinonecancer$mutornot[dataforbetabinonecancer$mutlabels != "None"] <- "Mut"
    dataforbetabinonecancer$mutornot[dataforbetabinonecancer$mutlabels == "None"] <- "Not_Mut"
    dataforbetabinonecancer$mutornot= factor(dataforbetabinonecancer$mutornot)
    
    nonmutset =  dataforbetabinonecancer[dataforbetabinonecancer$mutornot == "Not_Mut",]
    totalnonmutset = totalreadsperSJPANCANCERonecancer[dataforbetabinonecancer$mutornot == "Not_Mut",]
    mutset =  dataforbetabinonecancer[dataforbetabinonecancer$mutornot != "Not_Mut",]
    totalmutset = totalreadsperSJPANCANCERonecancer[dataforbetabinonecancer$mutornot != "Not_Mut",]
    
    
    whichsampletestNotmut = sample(1:dim(nonmutset)[1], (dim(nonmutset)[1]) / 2)
    
    trainnonmutset =  nonmutset
    totaltrainnonmutset =  totalnonmutset
    testnonmutset =  nonmutset
    totaltestnonmutset =  totalnonmutset
    
    
    whichsampletestmut = sample(1:dim(mutset)[1], (dim(mutset)[1]) / 2)
    
    trainmutset =  mutset
    totaltrainmutset = totalmutset
    testmutset =  mutset
    totaltestmutset =  totalmutset
    
    
    
    trainset = rbind(trainnonmutset, trainmutset)
    totaltrainset = rbind(totaltrainnonmutset, totaltrainmutset)
    testset = rbind(testnonmutset, testmutset)
    totaltestset = rbind(totaltestnonmutset, totaltestmutset)
    
    
    SJtrainset = trainset[,substring(colnames(trainset), 1,3) == "chr"]
    SJtestset = testset[,substring(colnames(testset), 1,3) == "chr"]
    alpha = matrix(0, nrow = 2, ncol = dim(SJtrainset)[2])
    beta = alpha
    
    
    for(feature in 1:dim(SJtrainset)[2] ){
        print(feature)
        sj =  SJtrainset[,feature]
        totalsj =  totaltrainset[,feature]
        
        
        K = sj[trainset$mutlabels == "None"]
        N = totalsj[trainset$mutlabels == "None"]
        
        
        
        
        o = optim(par = c(0.5,0.5), fn = betabinomial, method = "L-BFGS-B", lower=c(0+c,0+c), upper=c(10000,10000) , control = list(maxit = 1000000), total = N , altered = K)
        alpha[1, feature] = o$par[1]
        beta[1,feature] = o$par[2]
        
        
        
        K = sj[trainset$mutlabels != "None"]
        N = totalsj[trainset$mutlabels != "None"]
        
        
        o = optim(par = c(0.5,0.5), fn = betabinomial, method = "L-BFGS-B", lower=c(0+c,0+c), upper=c(10000,10000) , control = list(maxit = 1000000), total = N , altered = K)
        alpha[2, feature] = o$par[1]
        beta[2,feature] = o$par[2]
        
        
    }
    
    

    #popt = 1-popt
    lambda = alpha/(alpha + beta)
    colnames(lambda) <- NULL
    rownames(lambda) <- c("Not_Mut", "Mut")

    
    ismutornot = c("Not_Mut", "Mut")
    #SJtestset[SJtestset == 0] <- exp(-10)
    #ropt = mpfr(ropt, 100)
    
    py = table(trainset$mutornot) / sum(table(trainset$mutornot))
    py =  py[order(names(py), decreasing = T)]
    lpy = log(py)
    
    trainedparameters = list(alpha, beta, lpy)
}

saveRDS(trainedparameters,  "output/Trained_Parameters_Full.RDS")


