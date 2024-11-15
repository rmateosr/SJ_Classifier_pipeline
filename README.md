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

For context, we applied the `recount3` package to obtain SJ data from the TCGA dataset. Using the `available_projects()` function, we identified projects labeled with `project_type = "tcga"`, which yielded 10,507 tumoral samples from 33 distinct cancer types. The SJ information was downloaded by setting the `type = "jxn"` parameter in the `create_rse()` function. Through this process, we extracted the following tables:

1. **Sample Overview**  
   A table providing metadata about the samples, such as:
   - Cancer type or condition.
   - Unique identifiers for each sample.

2. **Splice Junction (SJ) List**  
   A comprehensive list of splice junctions (SJs), including:
   - **Coordinates**: Genomic positions of the junctions.
   - **Annotation Status**: Indication of whether each junction is annotated in reference genomes.

3. **SJ Count Matrix**  
   A matrix capturing the counts of each SJ across samples, structured as follows:
   - **Rows**: Individual splice junctions.
   - **Columns**: Individual samples.
   - **Entries**: Observed counts for each SJ-sample pair.

Although these tables are commonly generated using the `recount3` package, the script only requires that the input `.RDS` file contains data in this format.

----
## Output Format

The script generates a tab-separated values (TSV) file as output with the sample IDs as well as the Score obtained from the classifier. Below is an example of the output format:

```tsv
Sample_ID	Score
Sample1	-64.3326178927976
Sample2	14.9246427665612
Sample3	37.5658065396102
```

The output will be saved in the output/ folder with a filename matching the input file's base name followed by \_Score.tsv. 

----
## Example

To demonstrate how to use the script, you can run the following command:

```bash
./SJ_Classifier.sh example_data/SRP066737.RDS
```

