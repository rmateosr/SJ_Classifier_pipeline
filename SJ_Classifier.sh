#!/bin/bash
INPUT=$1
if [ -z "$INPUT" ]; then
    echo "Error: No input file provided."
    echo "Usage: $0 <input_file> [output_folder]"
    exit 1
fi
OUTPUT_FOLDER=${2:-output}
mkdir -p tmp
mkdir -p $OUTPUT_FOLDER
Rscript script/Matrix.R $INPUT
Rscript script/Test.R $INPUT $OUTPUT_FOLDER
rm -r tmp
