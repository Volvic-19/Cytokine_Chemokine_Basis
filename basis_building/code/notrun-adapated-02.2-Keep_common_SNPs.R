# Creating a common-SNP Ferkingstad dataset

# We're now selecting 104 immune-related Ferkingstad traits to create the basis.
# Once we have filtered and aligned Ferkingstad datasets, we'll create a single Ferkingstad file and remove all SNPs not common to all 104 datasets.
# MHC has been removed on the previous step 


# Load libraries
library(data.table)
setDTthreads(15)

#########################################################
gwas.DT <- fread("../manifest/gwas.DT.tsv.gz")# 25 GB


# Check for dupilcated SNPs in any of the datasets

gwas.DT[ , .(count = .N), by = c("pid","Trait")][ count > 1]

# Count SNP presence and keep only SNPs that are present in all datasets

SNP.counter <- gwas.DT[.(count := uniqueN(Trait)), by = pid]
SNP.counter <- SNP.counter[count == max(count)]
gwas.DT <- gwas.DT[pid %in% SNP.counter$pid]


# Check gwas.DT
summary(gwas.DT)


# Check if there is any column missing any value (e.g. LD block)
gwas.DT[is.na(ld.block)]
gwas.DT <- na.omit(gwas.DT)

# Save gwas.DT
fwrite(gwas.DT, "../data/gwas.DT.tsv.gz", sep="\t")

# Update manifest
SNP.manifest <- fread("../data/IL5_consensus_manifest_8M.tsv") #update file name
SNP.manifest <- SNP.manifest[pid %in% unique(gwas.DT$pid)]

fwrite(SNP.manifest, "../data/IL5_consensus_manifest_8M.tsv", sep="\t") #update file name


