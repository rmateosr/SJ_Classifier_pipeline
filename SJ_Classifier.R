#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)


if (length(args) != 2) {
  stop("Error: Please provide both the input file and output folder.\nUsage: Rscript SJ_Classifier.R <input_file> <output_folder>", call. = FALSE)
}


input_file <- args[1]
output_file <- args[2]

if (!file.exists(input_file)) {
  stop(paste("Error: The input file does not exist:", input_file), call. = FALSE)
}

output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  message("Output folder does not exist. Creating the folder...")
  dir.create(output_dir, recursive = TRUE)
}

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
sample <- as.data.frame(fread(input_file))
file_name <- file_path_sans_ext(basename(input_file))


rownames(sample) = sample[,1]
SJinfo = sample[,2]
sample = sample[,-(1:2),drop=FALSE]
sampleinfo = colnames(sample)


alteredSJnames = alteredSJnames[alteredSJnames!= "Sample_name"]
alteredSJcoords =  str_split_fixed(alteredSJnames, "[:-]", 4)
alteredstarts = apply(cbind(alteredSJcoords[,1], ":", alteredSJcoords[,2]),1,paste0, collapse = "")
alteredends = apply(cbind(alteredSJcoords[,1], ":", alteredSJcoords[,3]),1,paste0, collapse = "")
notannotatedarrayfromjxn =  sample[SJinfo==0, ,drop = FALSE]
alteredSJappearing = as.matrix(sample[rownames(sample) %in% alteredSJnames & SJinfo==0,,drop = FALSE] )


AlteredthisProject = matrix(0, ncol = length(alteredSJnames), nrow = dim(notannotatedarrayfromjxn)[2])
colnames(AlteredthisProject) = alteredSJnames
rownames(AlteredthisProject) = colnames(notannotatedarrayfromjxn)
for(cont in 1:dim(AlteredthisProject)[2]){
  alteredSJcountsinthisproject=as.matrix(alteredSJappearing[rownames(alteredSJappearing) %in% alteredSJnames[cont],,drop = FALSE])
  if(dim(alteredSJcountsinthisproject)[1]> 0 & dim(alteredSJcountsinthisproject)[2]> 0){
    AlteredthisProject[,cont] =  t(alteredSJcountsinthisproject)
  }
}
annotatedarrayfromjxn =  sample[SJinfo==1, ,drop = FALSE]
if(dim(annotatedarrayfromjxn)[1] > 0){
  SJcoords =  str_split_fixed(rownames(annotatedarrayfromjxn), "[:-]", 4)
  SJstarts = apply(cbind(SJcoords[,1,drop  =FALSE], ":", SJcoords[,2]),1,paste0, collapse = "")
  SJends = apply(cbind(SJcoords[,1,drop  =FALSE], ":", SJcoords[,3]),1,paste0, collapse = "")
  
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

samplefilefullname = input_file
SRAID =file_path_sans_ext(basename(samplefilefullname))

trainedparameters <- readRDS("resources/Trained_Parameters.RDS")
alpha = as.matrix(data.frame(trainedparameters[1]))
beta = as.matrix(data.frame(trainedparameters[2]))
lpy = unlist(trainedparameters[3])

Sjscorethreshold = 10
scorediffthreshold = 10
c = 0.001


dataforbetabin  = Altered
totalreadsperSJPANCANCER = Total


thechromosomes = apply(cbind("chr", c(1:22, "X")), 1, paste0, collapse =
                         "")
splitcolnames = str_split_fixed(colnames(dataforbetabin), ":", 2)[, 1]
dataforbetabin = dataforbetabin[, substring(splitcolnames, 1, 3) != "chr" |
                                  splitcolnames %in% thechromosomes, drop = FALSE]

splitcolnames = str_split_fixed(colnames(totalreadsperSJPANCANCER), ":", 2)[, 1]
totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[, substring(splitcolnames, 1, 3) != "chr" |
                                                      splitcolnames %in% thechromosomes, drop = FALSE]

whichones = colnames(dataforbetabin) %in% colnames(dataforbetabin)
dataforbetabin = dataforbetabin[, whichones, drop = FALSE]
totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[, whichones, drop = FALSE]
alpha = alpha[, whichones]
beta = beta[, whichones]


dataforbetabin = dataforbetabin[, colnames(dataforbetabin) != "Sample_name", drop = FALSE]
totalreadsperSJPANCANCER = totalreadsperSJPANCANCER[, grep("chr", colnames(totalreadsperSJPANCANCER)), drop = FALSE]
dataforbetabin_SJ = dataforbetabin[, grep("chr", colnames(dataforbetabin)), drop = FALSE]


resultspercancer = data.frame()

totalreadsperSJPANCANCERonecancer = totalreadsperSJPANCANCER
dataforbetabinonecancer = dataforbetabin


resultspercancer = data.frame()

totalreadsperSJPANCANCERonecancer = totalreadsperSJPANCANCER
dataforbetabinonecancer = dataforbetabin[, colnames(dataforbetabin) != "Sample_name", drop = FALSE]#

testset = dataforbetabinonecancer
totaltestset = totalreadsperSJPANCANCERonecancer

SJtestset = testset[, substring(colnames(testset), 1, 3) == "chr", drop = FALSE]

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

resultadosnumericPANCANCER = data.frame(
  Sample_ID = sampleinfo,
  Score = finalscore
)


output_dir <- dirname(output_file)

write.table(resultadosnumericPANCANCER, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

print("Analysis complete")


