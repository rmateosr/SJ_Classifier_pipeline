# SJ-Classifier
**A splicing junction-based classifier for detecting abnormal constitutive activation of the KEAP1-NRF2 pathway.**
## Introduction 

SJ-Classifier is a splicing Junction-based classifier trained on The Cancer Genome Atlas datasets for detection of abnormal constitutive activation of the KEAP1-NRF2 pathway.

----
## How to use 
To run the classifier, use the following command from the directory containing the script:
```bash
./SJ_Classifier.sh File_to_analyze.RDS
```

----
## Input Format

The input to this script is an `.RDS` file containing three key tables, structured as described below. These tables must adhere to the same design as those produced by the `create_rse()` function from the `recount3` R package with the parameter `type = "jxn"`, though the data itself does not need to originate from that specific function.

1. **Sample Overview**  
    - An overview detailing the samples for each cancer type.

2. **Splice Junction (SJ) List**  
    -	A comprehensive list of SJ, complete with coordinates and indications of whether the junction is annotated in reference genomes.

3. **SJ Count Matrix**  
    -	A matrix capturing the counts of each SJ, organized with junctions as rows and samples as columns..


----
## Output Format

The script generates a tab-separated values (TSV) file as output with the sample IDs as well as the Score obtained from the classifier. Below is an example of the output format:

Sample_ID | Score
--------- | -----
Sample1   | -64.3326178927976
Sample2   | 14.9246427665612
Sample3   | 37.5658065396102


The output will be saved in the output/ folder with a filename matching the input file's base name followed by \_Score.tsv. 

----
## Example

To demonstrate how to use the script, you can run the following command:

```bash
./SJ_Classifier.sh example_data/SRP066737.RDS
```

