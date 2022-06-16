# Creating a common-SNP Chen dataset

# We're now selecting 8 immune-related Chen traits to create the basis.
# Once we have filtered and aligned Chen datasets, we'll create a single Chen file and remove all SNPs not common to all 8 datasets.
# We'll also remove the MHC 

# 2022-06-08
# Guillermo Reales

# Load libraries
library(data.table)
setDTthreads(18)

#########################################################

# Chen trait selection

seltraits <- c("BASC", "EOSC", "LYMC", "MONC", "PLAV", "NEUC", "PLAC", "LEUC")
chends <- paste0("../data/",seltraits, "_Chen_6.tsv.gz")


# Create empty file and append to it
gwas.DT <- data.table()

for (i in chends){
	message("Loading ", i)
	ds <- fread(i, tmpdir="tmp")
	gwas.DT <- rbindlist(list(gwas.DT, ds))
	rm(ds)
}

# Check for duplicated SNPs in any of the datasets
gwas.DT[ , .(count = .N), by=c("pid","Trait") ][ count > 1 ]
# There seems to be one SNP (19:39892560) that's duplicated in some files
# This is because two SNPs in hg19 are mapped to the same pid in hg38. To prevent errors, we'll simply remove this SNP 
gwas.DT  <- gwas.DT[ pid != "19:39892560" ]

# Count SNP presence and keep only those present in all datasets
snp.counter <- gwas.DT[, .(count = uniqueN(Trait)), by =pid]
snp.counter <- snp.counter[count == max(count)] ## 8 copies
gwas.DT <- gwas.DT[pid %in% snp.counter$pid]
# Check gwas.DT
summary(gwas.DT)

# There are some SNPs missing LD.block, we need to remove them
gwas.DT[is.na(ld.block)]
gwas.DT <- na.omit(gwas.DT)

# We'll remove the MHC, too. We'll use the coordinates from Dilthey, 2021 (https://www.sciencedirect.com/science/article/pii/S1357272520301990), plus a window, so we'll remove blocks containing 20-40Mb
mhc.blocks <- gwas.DT[  CHR38 == 6 & BP38 > 20000000 & BP38 < 40000000, unique(ld.block) ]
gwas.DT  <- gwas.DT[ !ld.block %in% mhc.blocks ]

# Save gwas.DT
fwrite(gwas.DT, "../data/gwas.DT.tsv.gz", sep="\t")

# Update manifest

SNP.manifest <- fread("../data/BCB4_consensus_manifest_8M.tsv")
SNP.manifest <- SNP.manifest[pid %in% unique(gwas.DT$pid)]

fwrite(SNP.manifest, "../data/BCB4_SNP_manifest_after_Chen_8M.tsv", sep="\t")
