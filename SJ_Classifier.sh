#!/bin/bash
mkdir -p input
mkdir -p output
INPUT=$1
Rscript script/Matrix.R $INPUT
Rscript script/Test.R $INPUT
