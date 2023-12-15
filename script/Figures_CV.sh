#!/bin/sh
#PBS -e $HOME/SJproject/error.txt
module use /usr/local/package/modulefiles
module load R
Rscript  script/Figures_CV.R $1


