#!/usr/bin/env Rscript
#mergeCovarMats.R
#The purpose of this script is to combine all the features to be used as covariates into one matrix for input into GaMuT. It is intended to be run as part of the GaMuT_pipeline_driver.sh script.
#NOTE: The code currently assumes the unique IDs are in a column named "PMBB_ID"

library(tidyverse)

#read in the command line arguments
args = commandArgs(TRUE)

#analysis specific covartiates matrix (should only have sample IDs for this analysis)
analysis_covar_mat_path = args[1]
#path to common PCs file
common_PCs_path = args[2]
#path to rare PCs file
rare_PCs_path = args[3]
#path to local variants PCs file
local_PCs_path = args[4]
#path to the output file
output_path = args[5]

#read in the input files
covar_mat <- read_tsv(analysis_covar_mat_path,col_names = TRUE)
common_vars_PCs <- read_tsv(common_PCs_path,col_names = TRUE)
rare_vars_PCs <- read_tsv(rare_PCs_path,col_names = TRUE)
local_vars_PCs <- read_tsv(local_PCs_path,col_names = TRUE)

#perform merges on the PMBB_ID
covar_mat <- merge(covar_mat,common_vars_PCs, by = "PMBB_ID")
covar_mat <- merge(covar_mat,rare_vars_PCs, by = "PMBB_ID")
covar_mat <- merge(covar_mat,local_vars_PCs, by = "PMBB_ID")

#write the merged covar_mat to file
write.table(covar_mat, file = output_path , sep="\t", row.names=FALSE, quote=FALSE)
