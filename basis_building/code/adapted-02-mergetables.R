# This script is the intermediate step to merge all "manifested" Ferkingstad datasets into one table. Otherwise fread() fails to load 104 files in slurm job.

# Once we have filtered and aligned Ferkingstad datasets, we'll create a single Ferkingstad file and remove all SNPs not common to all 104 datasets.
# MHC has been removed on the previous step 


# Load libraries
library(data.table)
setDTthreads(18)

#########################################################
# read all 104 manifested_data (intersect with first manifest) into one file

files <- dir(path = "../manifested_data")
ppath <- "../manifested_data/"

gwas.DT <- data.table()

for (i in files){
	message("loading" , i)
	ds <- fread(paste0(ppath , i), tmpdir = "tmp")
	ds[, Trait := gsub("-ft.tsv.gz", "", i)]
	gwas.DT <- rbindlist(list(gwas.DT, ds))
	rm(ds)
}


# Save gwas.DT
fwrite(gwas.DT, "../manifest/gwas.DT.tsv.gz", sep="\t")

