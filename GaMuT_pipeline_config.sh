# This is the config file for the GaMuT_pipeline_driver.sh for the PMBB data. It is designed to provide the paths to files such as the paths to the VCFs for the exome and array data, the path to the phenotype files for the analysis, the path the gene list file, etc.

###################################################################################################
#Pipeline parameters (paths to scripts and BSUB parameters)
###################################################################################################

#path to the GaMuT R code provided at https://github.com/epstein-software/GAMuT from Broadaway et al. (2016)
GaMuT_functions="/home/wbone/group/personal/wbone/Rare_Var_Multi-trait_GWAS/epstein-software-GAMuT-e0da206/GAMuT-functions.R"

#path to the R script that runs GaMuT
GaMuT_Rscript="/home/wbone/group/personal/wbone/Rare_Var_Multi-trait_GWAS/GaMuT_pipeline_code/GaMuT_pipeline.R"

#path to the R script that combines the covariates into one file
merge_Covars="/home/wbone/group/personal/wbone/Rare_Var_Multi-trait_GWAS/GaMuT_pipeline_code/mergeCovarMats.R"

#path to the template bsub for each GaMuT run
bsub_template="/home/wbone/group/personal/wbone/Rare_Var_Multi-trait_GWAS/GaMuT_pipeline_code/GaMuT_bsub_template_v4.bsub"

#requested memory for each GaMuT job
mem_req="20"

#queue jobs should be submitted to
bsub_queue="epistasis_long"


###################################################################################################
#GENOTYPE DATA
###################################################################################################

#path to a file with all unrelated individuals in the cohort
unrelated_sample_IDs="/project/PMBB/PMBB-Release-2020-2.0/Exome/IBD/PMBB-Release-2020-2.0_genetic_exome_3rd_degree_unrelated"

#path to file that lists each ancestry specific samples file and provides an ancestry specific string tag
#DONT FORGET TO SWITCH TO THE BETTER FILE!
#ancestry_files="/home/wbone/group/personal/wbone/Rare_Var_Multi-trait_GWAS/GaMuT_pipeline_code/GaMuT_dependency_files/test_ancestry_file.txt"

#path to a hg38 fasta file (with fai file)
fasta_file="/project/ritchie02/projects/PMBB_VEP_Annotations/Homo_sapiens_assembly38.fasta"

#Regenie.anno file, this is the Regenie output file that contains all the variants you want to perform the analyses on
regenie_anno="/project/ritchie02/projects/PMBB_VEP_Annotations/FinalAnnotations/Regenie_Annotations/PMBB.Regenie.RareOnly.anno"

#\n delimited list of all the genes you wish to perform your analyses on. Example all the genes that have at least one variant that meet your criteria 
gene_list="/home/wbone/group/personal/wbone/Rare_Var_Multi-trait_GWAS/GaMuT_pipeline_code/PMBB.Regenie.RareOnly.ensembl_ID.45k.test.anno"

#Ensembl gene BED file. BED file of the ensembl ID of the chr, start, stop, gene IDi
gene_bed="/home/wbone/group/datasets/Ensembl/Homo_sapiens.GRCh38.Ensemble.GeneCoordinates.txt"

#Path to VCF files that have the gene variants that you wish to use for your analysis (path should include all the string for the chr specific . For PMBB this should be the exome data.
gene_vcf_input_path="/project/PMBB/PMBB-Release-2020-2.0/Exome/pVCF/GL_by_chrom/PMBB-Release-2020-2.0_genetic_exome_chr"
#str that follows the chromosome in the input vcf file names
gene_vcf_input_path_str_end="_GL.vcf.gz"

#Path to the VCF files that have the array data (these files are used for correcting for local ancestry around each gene
array_vcf_input_path="/project/PMBB/PMBB-Release-2020-2.0/Genotype/PMBB-Release-2020-2.0_genetic_genotype.vcf.gz"
#NOTE: Might switch back in the imputed data TBI files become available
#array_vcf_input_path="/project/PMBB/PMBB-Release-2020-2.0/Imputed/vcf/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_chr"
#string that follows the chromosome number in the array VCF file names
#array_vcf_input_path_str_end=".vcf.gz"

#gene level variant carrier frequency cut off (supplied as a float representing the prop of the sample size for each analysis) NOTE this is not the same as the case_carrier_prop_min which needs to be set separately if you have a case-control trait
carrier_prop_min=0.00075

#local window is the number of base pairs around the gene for which you wish to include variants for the local ancestry PCA for correcting for local ancestry. NOTE: it will go this many bps both upstream and downstream of the gene.
local_window=5000000

#path to SNP array data PCs for common variants
common_vars_PCs="/home/wbone/group/personal/wbone/Rare_Var_Multi-trait_GWAS/GaMuT_pipeline_code/CAD_phenos/PMBB_CAD_lipid_array_PCs_20220919_complete.txt"

#path to rare variant PCs
rare_vars_PCs="/home/wbone/group/personal/wbone/Rare_Var_Multi-trait_GWAS/GaMuT_pipeline_code/lipids_test_phenos/PMBB_test_lipid_rare_PCs_20220513.txt"
#NOTE: MAKE SURE TO CHANGE THIS. Switching it to Common PCs to test the code

###################################################################################################
#PHENOTYPE DATA
###################################################################################################

#case control trait included boolean
case_control="TRUE"

#Path to the phenotype matrix (should ONLY include samples with complete phenotype data MAKE SURE this includes pertinent covariates)
#CHANGE THIS TO BE THE NORMALIZED ONLY DATA
pheno_mat="/home/wbone/group/personal/wbone/Rare_Var_Multi-trait_GWAS/GaMuT_pipeline_code/CAD_phenos/PMBB_CAD_HDL_LDL_TG_medians_normalized_pheno_mat_20220908_complete.txt"

#gene level variant case-carrier frequency cut off (supplied as a float representing the prop of the sample size for each analysis) NOTE this is not the same as the carrier_prop_min which needs to be set separately
case_carrier_prop_min=0.000625

#Path to the covariates matrix (should ONLY include samples with complete data MAKE SURE this includes the phenotypes of interest)
covar_mat="/home/wbone/group/personal/wbone/Rare_Var_Multi-trait_GWAS/GaMuT_pipeline_code/CAD_phenos/PMBB_CAD_lipid_covar_mat_20220908_complete.txt"

###################################################################################################
#RUN LOG FILES
###################################################################################################
#log of the genes that were not run and the reason why they were not run
filtered_genes="/home/wbone/group/scratch/wbone/CAD_and_lipids_GaMuT_20220919/filtered_genes.log"
