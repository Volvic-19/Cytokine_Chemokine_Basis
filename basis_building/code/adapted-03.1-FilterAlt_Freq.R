# This script is an intermediate step. It filters out ALT_FREQ < 1% or >99% to help size down gwas.DT.tsv.gz 

# 2022-06-20

# Load libraries
library(data.table)
setDTthreads(15)






#########################################################
gwas.DT <- fread("../manifest/gwas.DT.tsv.gz", showProgress = TRUE, tmpdir = "tmp")# 25 GB

message("gwas.DT loaded!")

summary(gwas.DT)

# Filter out SNPs with allele frequency < 1% or >99%
gwas.DT <- gwas.DT[ALT_FREQ > 0.01]
gwas.DT <- gwas.DT[ALT_FREQ < 0.99]

message("After filtering ALT_FREQ ~ (0.01, 0.99),we get:")

# Check gwas.DT
dim(gwas.DT)

summary(gwas.DT)

# Write overlaped SNPs 

message("writing filtered SNPs")

fwrite(gwas.DT, "../manifest/filtered.1.99.gwas.DT.tsv.gz", sep="\t")

message("finished writing!")

