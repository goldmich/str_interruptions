# str_interruptions
The selective dynamics of interruptions at short tandem repeats
[bioRxiv](https://www.biorxiv.org/content/10.1101/2025.06.09.658724v1.abstract)

# Data processing
The manuscripts largely relies on STR genotypes from 1KGP and H3A as described in Jam et al., 2023. Please refer to the manuscript for the full list of databases analyzed and for the respective sources.

The majority of the data processing for all main and supplementary figures occurs in the `interruptions.snake` file. The environment to run the snakemake can be found in the `environment.yaml` file. Some steps in the data processing are not necessary for producing the manuscript figures or models (i.e. the vep analysis) though I am leaving them in for posterity.

The figures and models for the two main results sections are in separate R scripts. As their names would suggest, the `analyze_jam_coding_strs.R` includes the coding STR analysis and the `analyze_jam_strs.R` includes the noncoding analysis. The output from the snakemake files is input for the R code. I ran R version 4.2.2, largely in interactive mode.
