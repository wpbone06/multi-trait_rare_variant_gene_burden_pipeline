#!/usr/bin/env Rscript
#GaMuT_pipeline.R
#The purpose of this script is to run the burden GaMuT test given the gene specific genotype matrix, phenotype matrix, and covariates matrix (provided in a config file) it is designed to be run as a part of the GaMuT_pipeline_driver.sh script.

#read in required libraries
library(CompQuadForm)
library(tidyverse)

#parse command line arguments
args = commandArgs(TRUE)
GaMuT_config = args[1]

#source R code and config file
#confing file for this run
source(GaMuT_config)
#source published GaMuT functions
source(GaMuT_functions)

#function to prep the genotype data for a burden test. Convert to burden vector and mean center (input is a genotype matrix)
prep_burden_geno <- function(X){
  #make a single vector of either variant in gene or none
  X <-rowSums(X) > 0 
  X <- as.numeric(X)
    
  #mean center
  X  = as.matrix(scale(X,center=T,scale=F))
  return(X)
}

#Read in the input files from the config file and convert the IDs to rownames

#making a note here as a reminder for applying code to other datasets
print("The sample IDs here are PMBB specific")
geno <- read_tsv(gt_mat ,col_names = TRUE)
geno <- as.data.frame(geno)
rownames(geno) <- geno$IID
drop_geno <- c("IID")
geno = geno[,!(names(geno) %in% drop_geno)]

pheno <- read_tsv(pheno_mat ,col_names = TRUE)
pheno <- as.data.frame(pheno)
rownames(pheno) <- pheno$PMBB_ID
drop_pheno <- c("PMBB_ID")
pheno = pheno[,!(names(pheno) %in% drop_pheno)]

covar <- read_tsv(covar_mat ,col_names = TRUE)
covar <- as.data.frame(covar)
rownames(covar) <- covar$PMBB_ID
drop_covar <- c("PMBB_ID")
covar = covar[,!(names(covar) %in% drop_covar)]


#NOTE: need to make the sample IDs row names using rownames(data_frame) <- data_frame$col1

#NOTE: could filter out samples with no array data here if need be

##################################################################
#Adjust for covariates
##################################################################

#initialize the matrix to be filled with the residualized phenotypes from pheno
res_pheno = matrix(NA, nrow=nrow(pheno), ncol=ncol(pheno))

#for each phenotype (t) in pheno residualize it for all the covariate in covar code based on example from https://epstein-software.github.io/GAMuT/GAMuT-example.html
for(t in 1:ncol(pheno)){

    model.residual = lm(pheno[,t] ~ as.matrix(covar))
    res_pheno[,t] = resid(model.residual)

}

##################################################################
#Prep phenotype data for GaMuT test
##################################################################

#prep the residualized phenotype matrix for GaMuT
res_pheno <- apply(res_pheno, 2, scale, scale=T, center=T)
res_pheno  <- as.matrix(res_pheno)
proj_pheno = proj_GAMuT_pheno(res_pheno)

# projection matrix
pheno_Yc = proj_pheno$Kc
# eigenvalues of Yc
pheno_lambda_Y = proj_pheno$ev_Kc
  
##################################################################
#Prep genotype burden for default GaMuT test
##################################################################
geno = as.matrix(geno)
G = prep_burden_geno(geno)
#note this should be the normal GaMuT function, the prep step is the important one
burden_geno <- linear_GAMuT_geno(G)

# linear kernel similarity matrix
burden_linear_Xc <- burden_geno$Lc
# eigenvalues of Xc
burden_lambda_X <- burden_geno$ev_Lc

##################################################################
#Perform the burden GaMuT Test
##################################################################
GaMuT_burden_pvalue <- TestGAMuT(pheno_Yc,pheno_lambda_Y,burden_linear_Xc,burden_lambda_X)

##################################################################
#Write the results to file
##################################################################

#NOTE: need to add gene to the GaMuT_config file
result = as.data.frame(t(c(gene, GaMuT_burden_pvalue)))

colnames(result) = c("ENSG_ID","GaMuT_Burden_Test_P-value" )

#write the results to file
write.table(result, file = output , sep="\t", row.names=FALSE, quote=FALSE)


