#!/bin/sh
#PBS -e $HOME/SJproject/error.txt
module use /usr/local/package/modulefiles
module load R
Rscript  script/Step12_Ambiguous.R $1


