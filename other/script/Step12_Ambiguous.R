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
  
  alphabetamodel = function(SJandtotalcolumn,SJortotal) {
    c = 0.001
    SJcolumn =  SJandtotalcolumn[SJortotal== "SJ"]
    totalcolumn =  SJandtotalcolumn[SJortotal== "Total"]
    
    K = SJcolumn
    N = totalcolumn
    
    o = optim(par = c(0.5,0.5), fn = betabinomial, method = "L-BFGS-B", lower=c(0+c,0+c), upper=c(Inf,Inf) , control = list(maxit = 1000000), total = N , altered = K)
    alpha = o$par[1]
    beta = o$par[2]
    return(data.frame(alpha, beta))
  }


  Sjscorethreshold = 10
  scorediffthreshold= 20
  c = 0.001

  metadatafile = args[1]
  dataforbetabin  =data.frame(fread("tmp/TAN/Altered.tsv"), check.names = F)
  totalreadsperSJPANCANCER=data.frame(fread("tmp/TAN/Total.tsv"), check.names = F)
  mutmetadata  =data.frame(fread(metadatafile))
  dataforbetabin$TCGA_ID = substring(dataforbetabin$Sample_name,1,16)
  totalreadsperSJPANCANCER$TCGA_ID = dataforbetabin$TCGA_ID 
  splitcolnames = str_split_fixed(colnames(dataforbetabin), ":", 2)[,1]
  thechromosomes = apply(cbind("chr", c(1:22,"X")), 1, paste, collapse = "")
  dataforbetabin = dataforbetabin[,substring(splitcolnames,1,3) != "chr" |  splitcolnames %in% thechromosomes]
  
  splitcolnames = str_split_fixed(colnames(totalreadsperSJPANCANCER), ":", 2)[,1]
  totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[,substring(splitcolnames,1,3) != "chr" |  splitcolnames %in% thechromosomes]
  
  dataforbetabin = inner_join(dataforbetabin, mutmetadata)
  if(sum(dataforbetabin$relevance == "Ambiguous")> 0){
  dataforbetabin = dataforbetabin[dataforbetabin$relevance == "Ambiguous",]
  dataforbetabin = dataforbetabin[,colnames(dataforbetabin)!= "Sample_name"]
  totalreadsperSJPANCANCER = inner_join(totalreadsperSJPANCANCER, mutmetadata)
  
  totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[totalreadsperSJPANCANCER$relevance == "Ambiguous",]
  totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[,grep("chr", colnames(totalreadsperSJPANCANCER))]
  dataforbetabin_SJ= dataforbetabin[,grep("chr", colnames(dataforbetabin))]
  testset  =dataforbetabin
  SJtestset = dataforbetabin_SJ
  totaltestset = totalreadsperSJPANCANCER
  # if(nonormalthres != 1){
  #  SJlist = readRDS("C:/Users/Raul/Documents/NRF2/classifier_pipeline_SRA_based_individual_cancers/individualGATHEREDSJbypvalueasin02032000.RDS")
  pthresholds =  "p0001"
  
  npvalue = 1
  set.seed(12) 
  
  
  #dataforbetabin = dataforbetabin_initial
  #totalreadsperSJPANCANCER = totalreadsperSJPANCANCER_initial
   # dataforbetabin$TCGA_ID = substring(dataforbetabin$Sample_name,1,16)
  
  ###########################################
  ######## Removing duplicates ##############
  ###########################################
  
  
  theduplicates =unique( dataforbetabin$TCGA_ID[duplicated(dataforbetabin$TCGA_ID)])
  duppositionsOUT  = c()
  for(cont in 1:length(theduplicates) ){
    duppositions =   which(dataforbetabin$TCGA_ID == theduplicates[cont])
    duppositionsOUT = c(duppositionsOUT, duppositions[-which.max(rowSums(totalreadsperSJPANCANCER[duppositions,-1] ))])
    
  }
  
  dataforbetabin = dataforbetabin[-duppositionsOUT,]
  totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[-duppositionsOUT,]
  ###########################################
  ###########################################
  ###########################################
  
  
  
  
  
  
  
  ##########################################
  
  set.seed(12) 
  
  trainedparameters = readRDS("output/Trained_Parameters_Full.RDS")
  
  alpha = trainedparameters[[1]]
  beta = trainedparameters[[2]]
  lpy = trainedparameters[[3]]
    ismutornot = c("Not_Mut", "Mut")
    
    ###positive values
    #####
    lgammatotalplusone = lgamma(totaltestset + 1)
    
    #ln(r(k+alpha))
    #alpha[1,] -> notmut alpha
    #alpha[2, ] -> mut alpha
    lgammakplusalphanotmut = lgamma(SJtestset + rep(alpha[1,], each = nrow(SJtestset)))
    lgammakplusalphamut = lgamma(SJtestset + rep(alpha[2,], each = nrow(SJtestset)))
    
    #ln(r(n-k+b))
    lgammatotalminuskplusbetanotmut = lgamma(totaltestset - SJtestset + rep(beta[1,], each = nrow(totaltestset)))
    lgammatotalminuskplusbetamut = lgamma(totaltestset - SJtestset + rep(beta[2,], each = nrow(totaltestset)))
    
    #this will later be separated into not mut and mut 
    lgammaalphaplusbeta = lgamma(alpha + beta)
    #####
    
    
    #negative values
    #####
    lgammakplusone = lgamma(SJtestset + 1)
    lgammatotalminuskplusone = lgamma(totaltestset - SJtestset + 1)
    
    #ln(r(n + alpha + beta))
    lgammatotalplusalphaplusbetanotmut = lgamma(totaltestset + rep(alpha[1,], each = nrow(totaltestset)) + rep(beta[1,], each = nrow(totaltestset)) )
    lgammatotalplusalphaplusbetamut = lgamma(totaltestset + rep(alpha[2,], each = nrow(totaltestset)) + rep(beta[2,], each = nrow(totaltestset)) )
    
    #this will later be separated into not mut and mut 
    lgammaalpha = lgamma(alpha)
    
    #this will later be separated into not mut and mut 
    lgammabeta = lgamma(beta)
    #####
    
    
    
    
    notmutscorematrix = 
      lgammatotalplusone + 
      lgammakplusalphanotmut + 
      lgammatotalminuskplusbetanotmut + 
      rep(lgammaalphaplusbeta[1,], each = nrow(totaltestset)) -
      
      lgammakplusone - 
      lgammatotalminuskplusone -
      lgammatotalplusalphaplusbetanotmut - 
      rep(lgammaalpha[1,], each = nrow(totaltestset)) -
      rep(lgammabeta[1,], each = nrow(totaltestset)) 
    
    
    
    
    
    
    mutscorematrix = 
      lgammatotalplusone + 
      lgammakplusalphamut + 
      lgammatotalminuskplusbetamut + 
      rep(lgammaalphaplusbeta[2,], each = nrow(totaltestset)) -
      
      lgammakplusone - 
      lgammatotalminuskplusone -
      lgammatotalplusalphaplusbetamut - 
      rep(lgammaalpha[2,], each = nrow(totaltestset)) -
      rep(lgammabeta[2,], each = nrow(totaltestset)) 
    
    
    
    diffscoremutnotmut = mutscorematrix - notmutscorematrix
    diffscoremutnotmut[diffscoremutnotmut > Sjscorethreshold] <- Sjscorethreshold
    diffscoremutnotmut[diffscoremutnotmut < -Sjscorethreshold] <- -Sjscorethreshold
    ddiffscoremutnotmut = diffscoremutnotmut
    #pdf(paste0("C:/Users/Raul/Documents/NRF2/Features added version 201912/retake_January_2020/retakeonpipeline-pancancer/inputsoutputs/sdconditionalscore", cancers[cont],"bonferronitop1000exonskipupdatedindividualGATHERED.pdf"))
    varperfeature = (apply(ddiffscoremutnotmut, 2,sd))
    #boxplot(varperfeature, main = paste0("Standard deviation in conditional score per sample: ", cancers[cont]))
    #dev.off()
    #my_palette = viridis(256)
    my_palette = colorRampPalette(c("blue","white", "red"))(256)
    
    
    
    finalscore = lpy[2] - lpy[1] +  rowSums(diffscoremutnotmut)
    resultadosnumericPANCANCER= data.frame(mutlabels = testset$mutlabels, score = finalscore)
    resultadosnumericPANCANCER= data.frame(mutlabels = testset$mutlabels, score = finalscore)
    ###
    
    
    resultadosnumericPANCANCER= data.frame(TCGA_ID = testset$TCGA_ID, TCGA_ID_shorten = testset$TCGA_ID_shorten, mutlabels = testset$mutlabels, score = finalscore, cancer = testset$Cancer_Type)
    ###
    resultadosnumericPANCANCERNFE2L2 = inner_join(resultadosnumericPANCANCER, mutmetadata)
    
    write.table(resultadosnumericPANCANCERNFE2L2, paste0("output/Scores_Ambiguous.tsv"), sep = "\t", quote = F, row.names = F)
    resultadosnumericPANCANCER_CSample_name= resultadosnumericPANCANCER
    
  }