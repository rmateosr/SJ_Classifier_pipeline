#rm(list = ls())
set.seed(12)   
#set.seed(27)   
#set.seed(32)   
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


nonormalthres=0.9
nonormalthreslist = c(1,0.9,0.8, 0.5, 0.2, 0.1, 0)
nonormalthreslist = c(0, 0.1, 0.2, 0.5, 0.8, 0.9, 1)
nonormalthreslist = c(0.3, 0.4,  0.6, 0.7)
nonormalthreslist = c("notused")
#nonormalthreslist = c(0.5)
# nonormalthreslist = c(1,0.9,0.8, 0.5, 0.2, 0.1, 0)
#for(t in 1:length(nonormalthreslist)){
t  =1 
print(t)
nonormalthres=nonormalthreslist[t]
# Sjscorethreshold = 10
# scorediffthreshold=10
Sjscorethreshold = 10
scorediffthreshold= 20
c = 0.001
#N = rep(1500, 10000)
#P = rbeta(10000, 10000, 2.8)
# K = c()
# 
#for(cont in 1:length(P)){
#      K = c(K, rbinom(1, 1500, P[cont]))
#}
# 
# 
# o = optim(par = c(0.5,0.5), fn = betabinomial, method = "L-BFGS-B", lower=c(0+c,0+c), upper=c(Inf,Inf) , control = list(maxit = 1000000), total = N , altered = K)
# 
# o 
#set.seed(40)
library(stringr)
library(RColorBrewer)
#library(Rmpfr)
#library(doParallel)
library(foreach)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(viridis)
library(gplots)
library(data.table)
library(dplyr)
#set.seed(40)
library(beepr)
library(parallel)
thechromosomes = apply(cbind("chr", c(1:22, "X")), 1, paste0, collapse="")
#compilation of all cancer types
setwd("C:/Users/Raul/Documents/NRF2/SJclassifier_newdesign_072023_includingambiguouswithburdenandexonskipping//")

dataforbetabin_initial  =data.frame(fread("C:/Users/Raul/Documents/NRF2/SJclassifier_newdesign_072023_includingambiguouswithburdenandexonskipping//TotalAlteredNormalCounts_pipelinemutburden_includingambiguous2023/Alteredcounts_bonferroni0001pipelinemutburden_includingambiguous2023.tsv"), check.names = F)
totalreadsperSJPANCANCER_initial=data.frame(fread("C:/Users/Raul/Documents/NRF2/SJclassifier_newdesign_072023_includingambiguouswithburdenandexonskipping/TotalAlteredNormalCounts_pipelinemutburden_includingambiguous2023/Totalcounts_bonferroni0001pipelinemutburden_includingambiguous2023.tsv"), check.names = F)

splitcolnames = str_split_fixed(colnames(dataforbetabin_initial), ":", 2)[,1]
dataforbetabin_initial = dataforbetabin_initial[,substring(splitcolnames,1,3) != "chr" |  splitcolnames %in% thechromosomes]

splitcolnames = str_split_fixed(colnames(totalreadsperSJPANCANCER_initial), ":", 2)[,1]
totalreadsperSJPANCANCER_initial = totalreadsperSJPANCANCER_initial[,substring(splitcolnames,1,3) != "chr" |  splitcolnames %in% thechromosomes]


#  SJlist = readRDS("C:/Users/Raul/Documents/NRF2/classifier_pipeline_SRA_based_individual_cancers/individualGATHEREDSJbypvalueasin02032023.RDS")
pthresholds =  "p0001"

npvalue = 1
set.seed(12) 

dataforbetabin = dataforbetabin_initial
totalreadsperSJPANCANCER = totalreadsperSJPANCANCER_initial
dataforbetabin =  dataforbetabin_initial
totalreadsperSJPANCANCER =  totalreadsperSJPANCANCER_initial












mutmetadata  =fread("C:/Users/Raul/Documents/NRF2/SJclassifier_newdesign_072023_includingambiguouswithburdenandexonskipping/TCGA_mut_metadata_mutburden2023_1000_exskipupdated.tsv")

dataforbetabin$TCGA_ID = substring(dataforbetabin$Sample_name,1,16)

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



totalreadsperSJPANCANCER$TCGA_ID = dataforbetabin$TCGA_ID 
dataforbetabin = inner_join(dataforbetabin, mutmetadata)
dataforbetabin = dataforbetabin[,colnames(dataforbetabin)!= "Sample_name"]
#dataforbetabin = dataforbetabin[dataforbetabin$relevance != "Ambiguous",]
totalreadsperSJPANCANCER = inner_join(totalreadsperSJPANCANCER, mutmetadata)
#totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[totalreadsperSJPANCANCER$relevance != "Ambiguous",]
totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[,grep("chr", colnames(totalreadsperSJPANCANCER))]

dataforbetabin_SJ= dataforbetabin[,grep("chr", colnames(dataforbetabin))]

totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[,str_split_fixed(colnames(totalreadsperSJPANCANCER),":",2) %in% thechromosomes]

dataforbetabin_SJ = dataforbetabin_SJ[,str_split_fixed(colnames(dataforbetabin_SJ),":",2) %in% thechromosomes]


  cancers =  unique(dataforbetabin$Cancer_Type) 
#cancers = "LUAD"
resultspercancer = data.frame()

totalreadsperSJPANCANCERonecancer = totalreadsperSJPANCANCER
dataforbetabinonecancer = dataforbetabin

resultspercancer = data.frame()
# PANCANCER = fread("C:/Users/Raul/Documents/NRF2/mc3.v0.2.8.PUBLIC.maf", header = T,  stringsAsFactors = F, sep = "\t")
# PANCANCER = data.frame(PANCANCER, stringsAsFactors = F)
# #PANCANCER = PANCANCER[PANCANCER$Variant_Classification != "Silent",]
# PANCANCER = PANCANCER[PANCANCER$Variant_Classification != "3'UTR",]
# PANCANCER = PANCANCER[PANCANCER$Variant_Classification != "5'UTR",]
# PANCANCER = PANCANCER[PANCANCER$Variant_Classification != "Intron",]
# PANCANCER = PANCANCER[PANCANCER$Variant_Classification != "Silent",]
# sampleID = str_split_fixed(PANCANCER$Tumor_Sample_Barcode, "-", 4)
# sampleID = apply(sampleID[,1:3], 1, paste0, collapse = "-")
# PANCANCER$TCGA_ID_shorten = sampleID
# PANCANCERNFE2L2 =unique( data.frame(PANCANCER$Hugo_Symbol == "NFE2L2" , PANCANCER[c("TCGA_ID_shorten", "Protein_position", "Variant_Classification")]))
## PANCANCERNFE2L2 =unique( data.frame(PANCANCER$Hugo_Symbol == "NFE2L2" , PANCANCER[c("TCGA_ID_shorten","Hugo_Symbol","Chromosome" ,"Start_Position" , "End_Position", "Strand" ,  "Reference_Allele" ,"Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Codons",  "Protein_position", "Variant_Classification")]))
# PANCANCERNFE2L2$Variant_Classification[PANCANCERNFE2L2$Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del")] <- "Frame_Shift_InDel"
# PANCANCERNFE2L2$Variant_Classification[PANCANCERNFE2L2$Variant_Classification %in% c("In_Frame_Ins", "In_Frame_Del")] <- "In_Frame_InDel"
#for(cont in 1:length(cancers)){
PANCANCERNFE2L2 =  readRDS("C:/Users/Raul/Documents/NRF2/mc3.v0.2.8.PUBLIC_NFE2L2.RDS")
totalreadsperSJPANCANCERonecancer = totalreadsperSJPANCANCER#[dataforbetabin$Cancer_Type==cancers[cont],]
dataforbetabinonecancer = dataforbetabin[,colnames(dataforbetabin)!= "Sample_name"]#[dataforbetabin$Cancer_Type==cancers[cont],]
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
  
  trainnonmutset =  nonmutset[whichsampletestNotmut,]
  totaltrainnonmutset =  totalnonmutset[whichsampletestNotmut,]
  testnonmutset =  nonmutset[-whichsampletestNotmut,]
  totaltestnonmutset =  totalnonmutset[-whichsampletestNotmut,]
  
  
  
  whichsampletestmut = sample(1:dim(mutset)[1], (dim(mutset)[1]) / 2)
  
  trainmutset =  mutset[whichsampletestmut,]
  totaltrainmutset = totalmutset[whichsampletestmut,]
  testmutset =  mutset[-whichsampletestmut,]
  totaltestmutset =  totalmutset[-whichsampletestmut,]
  
  
  
  trainset = rbind(trainnonmutset, trainmutset)
  totaltrainset = rbind(totaltrainnonmutset, totaltrainmutset)
  testset = rbind(testnonmutset, testmutset)
  totaltestset = rbind(totaltestnonmutset, totaltestmutset)
  
  # 
  # a = testset
  # testset =trainset
  # trainset = a
  # 
  # a = totaltestset
  # totaltestset = totaltrainset
  # totaltrainset = a
  # 
  
  SJtrainset = trainset[,substring(colnames(trainset), 1,3) == "chr"]
  SJtestset = testset[,substring(colnames(testset), 1,3) == "chr"]
  
  SJandtotal = rbind(SJtrainset,totaltrainset)
  SJortotal= c(rep("SJ", dim(SJtrainset)[1]),rep("Total", dim(totaltrainset)[1]))
  mutlabels = trainset$mutlabels
  
  notmutresults = mclapply(as.list(data.frame(SJandtotal[mutlabels == "None",])), alphabetamodel,SJortotal = SJortotal[mutlabels == "None"])
  mutresults = mclapply(as.list(data.frame(SJandtotal[mutlabels != "None",])), alphabetamodel,SJortotal = SJortotal[mutlabels != "None"] )
  
  abnotmut = matrix(unlist(notmutresults), nrow=2)
  abmut = matrix(unlist(mutresults), nrow=2)
  
  alpha = rbind(abnotmut[1,], abmut[1,])
  beta = rbind(abnotmut[2,], abmut[2,])
  
  
  
  
  
  my_palette = viridis(256)
  #popt = 1-popt
  lambda = alpha/(alpha + beta)
  colnames(lambda) <- NULL
  rownames(lambda) <- c("Not_Mut", "Mut")
  #pdf("C:/Users/Raul/Documents/NRF2/Features added version 201912/retake_January_2020/retakeonpipeline-pancancer/inputsoutputs/myscript_betabin_MLE_viridis_trainset_PANCANCER_TEST_WEIGHTEDNOTINVERSE_totalSJwithNAbonferronitop1000exonskipupdatedindividualGATHERED.pdf")
  #heatmap.2(
  #   lambda,
  #   na.rm = T,
  #   Rowv = NULL,
  #   Colv = NULL,
  #   cexCol = 0.5,
  #   cexRow = 1.5,
  #   margins = c(8, 8),
  #   trace = "none",
  #   col = my_palette,
  #   main = "PANCANCER lambda",
  #   labCol = FALSE
  # )
  #dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ismutornot = c("Not_Mut", "Mut")
  #SJtestset[SJtestset == 0] <- exp(-10)
  #ropt = mpfr(ropt, 100)
  
  py = table(trainset$mutornot) / sum(table(trainset$mutornot))
  py =  py[order(names(py), decreasing = T)]
  lpy = log(py)
  
  
  #Score items
  
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
  resultadosnumericPANCANCER= data.frame(testset$TCGA_ID, mutlabels = testset$mutlabels, score = finalscore)
  differenceSJPANCANCER= data.frame(mutlabels = testset$mutlabels, diff = diffscoremutnotmut)
  write.table(resultadosnumericPANCANCER, paste0("diffscoremutnotmutp0001parallel_CV1_name_chrremoved.tsv"), sep = "\t", quote = F, row.names = F)
  
  ###
  
  
  resultadosnumericPANCANCER= data.frame(TCGA_ID = testset$TCGA_ID, TCGA_ID_shorten = testset$TCGA_ID_shorten, mutlabels = testset$mutlabels, score = finalscore, cancer = testset$Cancer_Type)
  ###
  resultadosnumericPANCANCERNFE2L2 = inner_join(resultadosnumericPANCANCER, mutmetadata)
  
  write.table(resultadosnumericPANCANCER, paste0("resultsindividualGATHERED_directresultclassifier_",pthresholds[npvalue],"parallel_CV1_name_otherchrremoved.tsv"), sep = "\t", quote = F, row.names = F)
  resultadosnumericPANCANCER_CSample_name= resultadosnumericPANCANCER
  
  
  
  
  #########  2  #################
  pivot = trainset
  trainset = testset
  testset = pivot 
  
  pivot = totaltrainset
  totaltrainset = totaltestset
  totaltestset = pivot 
  
  SJtrainset = trainset[,substring(colnames(trainset), 1,3) == "chr"]
  SJtestset = testset[,substring(colnames(testset), 1,3) == "chr"]
  
  SJandtotal = rbind(SJtrainset,totaltrainset)
  SJortotal= c(rep("SJ", dim(SJtrainset)[1]),rep("Total", dim(totaltrainset)[1]))
  mutlabels = trainset$mutlabels
  
  notmutresults = mclapply(as.list(data.frame(SJandtotal[mutlabels == "None",])), alphabetamodel,SJortotal = SJortotal[mutlabels == "None"])
  mutresults = mclapply(as.list(data.frame(SJandtotal[mutlabels != "None",])), alphabetamodel,SJortotal = SJortotal[mutlabels != "None"] )
  
  abnotmut = matrix(unlist(notmutresults), nrow=2)
  abmut = matrix(unlist(mutresults), nrow=2)
  
  alpha = rbind(abnotmut[1,], abmut[1,])
  beta = rbind(abnotmut[2,], abmut[2,])
  
  
  
  
  
  my_palette = viridis(256)
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
  
  
  #Score items
  
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
  resultadosnumericPANCANCER= data.frame(testset$TCGA_ID, mutlabels = testset$mutlabels, score = finalscore)
  differenceSJPANCANCER= data.frame(mutlabels = testset$mutlabels, diff = diffscoremutnotmut)
  write.table(resultadosnumericPANCANCER, paste0("diffscoremutnotmutp0001parallel_CV2_name_chrremoved.tsv"), sep = "\t", quote = F, row.names = F)
  
  ###
  
  
  resultadosnumericPANCANCER= data.frame(TCGA_ID = testset$TCGA_ID, TCGA_ID_shorten = testset$TCGA_ID_shorten, mutlabels = testset$mutlabels, score = finalscore, cancer = testset$Cancer_Type)
  ###
  resultadosnumericPANCANCERNFE2L2 = inner_join(resultadosnumericPANCANCER, PANCANCERNFE2L2)
  resultadosnumericPANCANCERNFE2L2 = inner_join(resultadosnumericPANCANCER, mutmetadata)
  
  write.table(resultadosnumericPANCANCER, paste0("resultsindividualGATHERED_directresultclassifier_",pthresholds[npvalue],"parallel_CV2_otherchrremoved.tsv"), sep = "\t", quote = F, row.names = F)
  
  
  resultadosnumericPANCANCER_CV2  = resultadosnumericPANCANCER
  
  resultadosnumericPANCANCER_CVfull =  rbind(resultadosnumericPANCANCER_CSample_name, resultadosnumericPANCANCER_CV2)
  
  dim(resultadosnumericPANCANCER_CVfull)
  
}

resultadosnumericPANCANCER_CVfull = unique(resultadosnumericPANCANCER_CVfull)
write.table(resultadosnumericPANCANCER_CVfull, paste0("resultsindividualGATHERED_directresultclassifier_",pthresholds[npvalue],"parallel_CVfull_otherchrremoved.tsv"), sep = "\t", quote = F, row.names = F)




