#!/bin/bash

# This script is intended to create a consensus SNP list for the datasets (Chen) we'll be using for the basis plus the largest databanks we're using: UKBB (Neale + PanUKBB), and FinnGen.
# To achieve this, we'll download some example files, process them if necessary, and create the intersection.

# Download the files

# Neale UKBB Asthma
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz -O Asthma_Neale.tsv.gz

# PanUKBB Asthma
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/categorical-20002-both_sexes-1111.tsv.bgz -O Asthma_PanUKBB.tsv.gz

# FinnGen Asthma
wget https://storage.googleapis.com/finngen-public-data-r6/summary_stats/finngen_R6_J10_ASTHMA.gz -O Asthma_FinnGen.tsv.gz


# Prepare UKBB files prior to liftover
Rscript 00-Preprocess_cons_files.R
# This will output *_processed.tsv.gz files

# PanUKBB and Neale are in hg19, let's use the pipeline to liftover them
~/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/pipeline_v5.3.2_beta.sh -f Asthma_Neale_processed.tsv.gz
~/rds/rds-cew54-basis/GWAS_tools/01-Pipeline/pipeline_v5.3.2_beta.sh -f Asthma_PanUKBB_processed.tsv.gz

# Create the consensus list of SNPs and assign LD blocks to it using MacDonald LDdetect hg38 LD blocks. This results in 9800900 SNPs.
Rscript 00-create_BCB4_consensus_SNP_list.R

# Once we have everything, we can delete large sumstats files used to create the manifest.
rm *.tsv.gz

