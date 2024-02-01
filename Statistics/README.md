## Step 7
# Quantification and statistical analysis

### Purpose
The purpose of this analysis was to assess the efficacy of the mapping process by employing `SAMtools stats` for a comprehensive statistical evaluation of read characteristics within a specific BAM file. 
Additionally, the mean coverage across the entire exome region was determined using `mosdepth`. 
Visual representations in the Results section, including boxplots, probability density function (PDF) curves, and column plots, were crafted using the R programming language. 

# Prerequisites

Ensure the following prerequisites are met before conducting this analysis:

- SAMtools is installed for statistical analysis.
- Mosdepth is available for calculating mean coverage.
- The R programming language is installed along with required packages such as tidyverse, ggpubr, conflicted, vcfR, and writexl.

# Usage
```bash
samtools stats --threads $THREAD $SAMPLE_NAME.bam > $SAMPLE_NAME.stats

mosdepth $SAMPLE_NAME $SAMPLE_NAME.bam
```
The statistical R script is available at [Statistics/allele_count.Rmd](https://github.com/LuciaNhuNguyen/Masterarbeit/blob/main/Statistics/allele_count.Rmd).
