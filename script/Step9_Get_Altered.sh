#!/bin/sh
#PBS -e $HOME/SJproject/error.txt
module use /usr/local/package/modulefiles
module load R
Rscript  script/Step9_Get_Altered.R $1 $2


