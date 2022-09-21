# Multi-trait rare variant gene burden test pipeline
This repository is for a pipeline that runs a multi-trait rare variant burden test  on a cohort for which you have rare variant data for genes across the human  genome. It designed to include  genome-wide common variant PCs, genome-wide rare variant PCs, and local variants around each gene PCs, as well as any other covariates you may wish to include.

GaMuT is the software used to run a gene-based association testing of rare variants with multiple phenotypes in a sample of unrelated subjects. Original paper is here: Broadaway etl al. 2016 http://www.cell.com/ajhg/abstract/S0002-9297(16)00052-5. This code uses a burden based genotype matrix for each gene.

- GaMuT_pipeline.R is an R script that performs the GaMuT functions for each gene

- mergeCovarMats.R is an R script that merges the covariant inputs together. Aall the PC files and the user defined covariates(e.g. sex and age)

- GaMuT_bsub_template_v4.bsub is the bash script template that performs the preprocessing of the data (primarily the genetic data) for each gene and coordinates the running of the GaMuT code for each gene. The GaMuT_pipeline_driver_v4.sh file will fill in the gene specific data for each gene and fire off each gene as a separate job.

- GaMuT_pipeline_driver_v4.sh is the bash script that loops through the list of genes you would like to analyze and submits a separate job for each gene based on editing GaMuT_bsub_template_v4.bsub

- GaMuT_pipeline_driver.bsub is just a bsub file to execute and track the success of GaMuT_pipeline_driver_v4.sh


- GaMuT_pipeline_config.sh is a config file for providing the paths to the input files needed for performing your analysis. For example you need to provide paths to the in VCF file, the input phenotype data, the input covariates, and the genotype PCs. There eare also a numbeer of parameters the user can set including the minimum number of carriers for a gene to be analyzed and the size of the window around the gene to use for defining local PCs.

## How to run the pipeline
After configuring your GaMuT_pipeline_config.sh to include all your paths to input files, paths to required scripts (GaMuT scripts, bsub template file, etc.) and other run specific parameters (memory allocation, submission queue, minimum number of carriers, etc. ) The pipeline can launched by navigating to the analysis directory and running:

```
bsub < GaMuT_pipeline_driver.bsub
```
or

```
./GaMuT_pipeline_driver_v4.sh
```
