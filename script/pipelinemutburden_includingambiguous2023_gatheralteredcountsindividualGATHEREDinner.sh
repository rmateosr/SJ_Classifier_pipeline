#!/bin/sh
#PBS -e $HOME/SJproject/error.txt
module use /usr/local/package/modulefiles
module load R
Rscript  pipelinemutburden_includingambiguous2023_Step8_getalteredcountsfromSJselection_individualGATHEREDparallel.R $1


