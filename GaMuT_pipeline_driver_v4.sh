#!/bin/bash
#GaMuT_pipeline_driver_v4.sh
#This script is designed to be the bash script cooardinating the GaMuT Rare Variant Multi-trait GWAS pipeline.
#This version is reconfigured to include the local variants PCA in the gene level bsub (same bsub that handles the GaMuT burden test)

##################################################################
#Prep to run GaMuT pipeline 
##################################################################

#load modules
module load plink/1.90Beta6.18
module load htslib
module load bcftools/1.9

#save the starting directory as the analysis directory
analysis_dir=`pwd`

#source the analysis config file
source GaMuT_pipeline_config.sh

#generate a list of the samples that have complete phenotype data from the phenotype matrix (the first column contains all of the sample IDs. Need to remove the PMBB_ID at the top
cat $pheno_mat | cut -f 1 | tail -n +2 > complete_phenos_sample_IDs.txt

##################################################################
#Loop through each gene
##################################################################

#loop through the list of genes "gene_list" (ensembl IDs) that are represented in the data
while read gene
do

    echo $gene

    #move into a gene specific directory 
    mkdir $gene

    cd $gene
    gene_dir=`pwd`

    #find the genomic positions of the gene from a BED file of all the Ensembl IDs
    gene_chr=`grep $gene $gene_bed | cut -f 1`

    #continue to next gene if the gene is on the X or the Y chormosome
    #if [[ "$gene_chr" == "X" ]] || [[ "$gene_chr" == "x" ]] || [[ "$gene_chr" == "Y" ]] || [[ "$gene_chr" == "y" ]]
    #check that the chromosome contains letters
    if [[ $gene_chr =~ [a-zA-Z] ]]
    then

        echo "Skipping " $gene " because it is on the " $gene_chr " chromosome."
        echo
        echo "Skipping " $gene " because it is on the " $gene_chr " chromosome." >> $filtered_genes
        cd $analysis_dir
        continue

    fi

    #NOTE: for v4 plan to make everything after this part of the bsub
    gene_start=`grep $gene $gene_bed | cut -f 2`

    gene_stop=`grep $gene $gene_bed | cut -f 3`

    ##################################################################
    #Generate GaMuT Rscript config file
    ##################################################################
    echo GaMuT_functions=\"$GaMuT_functions\" > GaMuT_config.R

    echo gene=\"$gene\" >> GaMuT_config.R

    echo gt_mat=\"$gene_dir/$gene"_vars_GT_mat.txt"\" >> GaMuT_config.R
    echo pheno_mat=\"$gene_dir/$gene"_pheno_mat.txt"\" >> GaMuT_config.R
    echo covar_mat=\"$gene_dir/$gene"_all_covar_mat.txt"\" >> GaMuT_config.R
    echo output=\"$gene_dir/$gene"_GaMuT_burden_test_result.txt"\" >> GaMuT_config.R

    ##################################################################
    #Prep the bsub file
    ################################################################## 

    #copy bsub template to this directory
    cp $bsub_template ./$gene"_GaMuT.bsub"

    #add path to the analysis directory (overall starting directory)
    sed -i "s|ANALYSIS_DIR|$analysis_dir|g" $gene"_GaMuT.bsub"

    #add path to GaMuT Rscript to the bsub
    sed -i "s|GAMUT_RSCRIPT|$GaMuT_Rscript|" $gene"_GaMuT.bsub"

    #add path to this directory to the bsub (pipe delimiters due to "/" in path
    sed -i "s|GENE_DIR|$gene_dir|g" $gene"_GaMuT.bsub"

    #populate variables specific to this run in the bsub
    sed -i "s/GENE_ID/$gene/" $gene"_GaMuT.bsub"
    sed -i "s/GENE_CHR/$gene_chr/" $gene"_GaMuT.bsub"
    sed -i "s/LOCAL_REGION_START/$local_region_start/" $gene"_GaMuT.bsub"
    sed -i "s/LOCAL_REGION_STOP/$local_region_stop/" $gene"_GaMuT.bsub"
    #sed -i "s|ARRAY_CHR_VCF|$array_chr_vcf|" $gene"_GaMuT.bsub"

    echo "Submitting job for:" $gene $gene"_GaMuT.bsub"

    #submit the GaMuT bsub
    bsub < $gene"_GaMuT.bsub" -R \"rusage[mem="$mem_req"GB]\" -M "$mem_req"G -q $bsub_queue

    #cd back to the top level analysis directory
    cd $analysis_dir

done < $gene_list


echo "All GaMuT association jobs submitted!" 
