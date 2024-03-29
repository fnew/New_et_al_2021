# snp_metagenomes
Repo for code for associating SNPs and Metagenomic genes.

/qc 
----
This directory contains scripts and documentation for performing quality control on the metagenomic data.

* Removal of human reads
* Removal of duplicate reads
* Quality trimming

/metagenomic_analysis
---------------------
This directory contains scripts, submission scripts, and documentation for processing and analysing metagenomic data.

* Alignments
* Gene abundance calculations
* Creating the abundance table
* and more

/taxonomic_data_processing
--------------------------
This directory contains scripts, submission scripts, and R markdowns for processing and visualizing the taxonomic assignments of the metagenomic data.

* Kraken2 and Bracken 
* Species visualizations in R

/snp_processing
---------------
This directory contains scripts and documentation for processing of the SNP data.

* Preprocessing and reformatting
* Filtering: HWE, LD, missingness, MAF

/statistical_analysis
---------------------
This directory contains scripts, submission scripts, and documentation for the statistical analysis of associating metagenomic data and SNPs.

* PCA, linear mixed model to extract residuals
* SCCA functions
* SCCA
* SCCA post-processing
