#!/bin/sh
#PBS -e $HOME/SJproject/error.txt
module use /usr/local/package/modulefiles
module load R
echo $2
echo holi
Rscript script/Step7_Wilcoxon.R $1 $2 $3


