# SJ-Classifier
**A splicing junction-based classifier for detecting abnormal constitutive activation of the KEAP1-NRF2 pathway.**
## Introduction 

SJ-Classifier is a splicing junction-based classifier trained on The Cancer Genome Atlas datasets, designed to detect abnormal constitutive activation of the KEAP1-NRF2 pathway.

----
## How to use 
To run the classifier, exacute the following command from the directory containing the script:
```bash
Rscript SJ_Classifier.R /path/to/input_file.tsv /path/to/output_file.tsv
```

----
## Input Format

The input to this script is a `.tsv` file containing multiple columns, structured as follows:

1. **1st Column: Splice Junction (SJ)**
    - Contains unique identifiers for each splice junction in the dataset.
    - Format: chromosome:start-end:strand
        - chromosome: The chromosome where the splice junction is located (e.g., chr1, chrX).
        - start: The genomic position of the start of the splice junction (1-based).
        - end: The genomic position of the end of the splice junction (1-based).
        - strand: Indicates the transcription strand where the splice junction occurs (+ for forward strand, - for reverse strand).
    - Example format: `chr1:12228-12612:+`

2. **2nd Column: Annotation**
    - A binary indicator (0 or 1) specifying whether the junction is annotated in reference genomes, similar to the annotations reported in datasets accessed via the recount3 resource:
        - `1`: Annotated in reference genomes.
        - `0`: Not annotated in reference genomes.

3. **Columns 3 and Beyond: SJ Count Matrix**
    - Represents a matrix capturing the counts of each splice junction (SJ).
    - Each row corresponds to a splice junction.
    - Each additional column represents the counts for a specific sample.
    - Column headers must reflect the sample names.

### Additional Notes
- Input files must have a consistent structure with tab-delimited columns.
- Ensure the input file includes all required columns in the correct order to avoid errors during processing.

----
## Output Format

The script generates a tab-separated values (TSV) file as output with the sample IDs as well as the score obtained from the classifier. Below is an example of the output format:

Sample_ID | Score
--------- | -----
Sample1   | -64.3326178927976
Sample2   | 14.9246427665612
Sample3   | 37.5658065396102

Note: The score represents the classifier's confidence level, where positive scores indicate a higher likelihood of abnormal activation of the KEAP1-NRF2 pathway, and negative scores indicate a lower likelihood.

----
## Example

To demonstrate how to use the script, the following command can be executed:

```bash
Rscript SJ_Classifier.R example_data/SRP066737_test.tsv output/SJ_Classifier_results.tsv
```
