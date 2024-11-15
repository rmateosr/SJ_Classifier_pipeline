#!/usr/bin/Rscript
#args is the SRA file
args = commandArgs(trailingOnly = TRUE)
set.seed(12)

library(stringr)
library(foreach)
library(ggsci)
library(data.table)
library(dplyr)

samplefilefullname = args[1]
SRAID = tools::file_path_sans_ext(basename(samplefilefullname))

trainedparameters <- readRDS("resources/Trained_Parameters.RDS")
alpha = as.matrix(data.frame(trainedparameters[1]))
beta = as.matrix(data.frame(trainedparameters[2]))
lpy = unlist(trainedparameters[3])

Sjscorethreshold = 10
scorediffthreshold = 10
c = 0.001


SRAproject <- readRDS(paste0("input/", SRAID , ".RDS"))
dataforbetabin  = data.frame(SRAproject[[1]], check.names = F)
totalreadsperSJPANCANCER = data.frame(SRAproject[[2]], check.names = F)


thechromosomes = apply(cbind("chr", c(1:22, "X")), 1, paste0, collapse =
                         "")
splitcolnames = str_split_fixed(colnames(dataforbetabin), ":", 2)[, 1]
dataforbetabin = dataforbetabin[, substring(splitcolnames, 1, 3) != "chr" |
                                  splitcolnames %in% thechromosomes]

splitcolnames = str_split_fixed(colnames(totalreadsperSJPANCANCER), ":", 2)[, 1]
totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[, substring(splitcolnames, 1, 3) != "chr" |
                                                      splitcolnames %in% thechromosomes]

whichones = colnames(dataforbetabin) %in% colnames(dataforbetabin)
dataforbetabin = dataforbetabin[, whichones]
totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[, whichones]
alpha = alpha[, whichones]
beta = beta[, whichones]
if (length(SRAproject)  >  2) {
  dataforbetabin$TCGA_ID = SRAproject[[3]]$external_id
  totalreadsperSJPANCANCER$TCGA_ID = SRAproject[[3]]$external_id
} else {
  dataforbetabin$TCGA_ID = SRAID[nSRA]
  totalreadsperSJPANCANCER$TCGA_ID = SRAID[nSRA]
  
}

dataforbetabin = dataforbetabin[, colnames(dataforbetabin) != "Sample_name"]
totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[, grep("chr", colnames(totalreadsperSJPANCANCER))]
dataforbetabin_SJ = dataforbetabin[, grep("chr", colnames(dataforbetabin))]


cancers =  unique(dataforbetabin$Cancer_Type)

resultspercancer = data.frame()

totalreadsperSJPANCANCERonecancer = totalreadsperSJPANCANCER
dataforbetabinonecancer = dataforbetabin


resultspercancer = data.frame()

totalreadsperSJPANCANCERonecancer = totalreadsperSJPANCANCER
dataforbetabinonecancer = dataforbetabin[, colnames(dataforbetabin) != "Sample_name"]#

testset = dataforbetabinonecancer
totaltestset = totalreadsperSJPANCANCERonecancer

SJtestset = testset[, substring(colnames(testset), 1, 3) == "chr"]

lgammatotalplusone = lgamma(totaltestset + 1)

lgammakplusalphanotmut = lgamma(SJtestset + rep(alpha[1, ], each = nrow(SJtestset)))
lgammakplusalphamut = lgamma(SJtestset + rep(alpha[2, ], each = nrow(SJtestset)))


lgammatotalminuskplusbetanotmut = lgamma(totaltestset - SJtestset + rep(beta[1, ], each = nrow(totaltestset)))
lgammatotalminuskplusbetamut = lgamma(totaltestset - SJtestset + rep(beta[2, ], each = nrow(totaltestset)))


lgammaalphaplusbeta = lgamma(alpha + beta)


lgammakplusone = lgamma(SJtestset + 1)
lgammatotalminuskplusone = lgamma(totaltestset - SJtestset + 1)


lgammatotalplusalphaplusbetanotmut = lgamma(totaltestset + rep(alpha[1, ], each = nrow(totaltestset)) + rep(beta[1, ], each = nrow(totaltestset)))
lgammatotalplusalphaplusbetamut = lgamma(totaltestset + rep(alpha[2, ], each = nrow(totaltestset)) + rep(beta[2, ], each = nrow(totaltestset)))

lgammaalpha = lgamma(alpha)

lgammabeta = lgamma(beta)


notmutscorematrix =
  lgammatotalplusone +
  lgammakplusalphanotmut +
  lgammatotalminuskplusbetanotmut +
  rep(lgammaalphaplusbeta[1, ], each = nrow(totaltestset)) -
  
  lgammakplusone -
  lgammatotalminuskplusone -
  lgammatotalplusalphaplusbetanotmut -
  rep(lgammaalpha[1, ], each = nrow(totaltestset)) -
  rep(lgammabeta[1, ], each = nrow(totaltestset))


mutscorematrix =
  lgammatotalplusone +
  lgammakplusalphamut +
  lgammatotalminuskplusbetamut +
  rep(lgammaalphaplusbeta[2, ], each = nrow(totaltestset)) -
  
  lgammakplusone -
  lgammatotalminuskplusone -
  lgammatotalplusalphaplusbetamut -
  rep(lgammaalpha[2, ], each = nrow(totaltestset)) -
  rep(lgammabeta[2, ], each = nrow(totaltestset))



diffscoremutnotmut = mutscorematrix - notmutscorematrix
diffscoremutnotmut[diffscoremutnotmut > Sjscorethreshold] <-
  Sjscorethreshold
diffscoremutnotmut[diffscoremutnotmut < -Sjscorethreshold] <-
  -Sjscorethreshold
ddiffscoremutnotmut = diffscoremutnotmut
varperfeature = (apply(ddiffscoremutnotmut, 2, sd))

finalscore = lpy[2] - lpy[1] +  rowSums(diffscoremutnotmut)

if (length(SRAproject)  >  2) {
  resultadosnumericPANCANCER = data.frame(
    Sample_ID = SRAproject[[3]]$external_id,
    Score = finalscore
  )
  
} else {
  resultadosnumericPANCANCER = data.frame(
    Sample_ID =  SRAID[nSRA],
    Score = finalscore
  )
  
  
}

output_file <- paste0("output/", SRAID, "_Score.tsv")

write.table(resultadosnumericPANCANCER, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

print("Analysis complete")


