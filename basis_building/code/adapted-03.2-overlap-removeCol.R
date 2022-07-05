# This is an intermediate step to keep only SNPs present in all 104 Ferkingstad traits.
# Unnecessary columns were removed to save computing memory. Otherwise Slurm job could easily be shut down.

# 2022-06-20

# Load libraries
library(data.table)
setDTthreads(15)



#########################################################
gwas.DT <- fread("../manifest/filtered.1.99.gwas.DT.tsv.gz", showProgress = TRUE, tmpdir = "tmp")# 

message("gwas.DT loaded!")

# Check for dupilcated SNPs in any of the datasets

message("Is there duplicated pid-trait?")

gwas.DT[ , .(count = .N), by = c("pid","Trait")][ count > 1]
# gwas.DT <- gwas.DT[,("Trait","pid","ld.block","Name")]
# rm(gwas.DT)

message("remove unwantted columns~")

gwas.DT[, c("CHR38","BP38","REF","ALT","Name","rsids","ImpMAF"):=NULL]


# Count SNP presence and keep only SNPs that are present in all datasets
message("start counter!")

SNP.counter <- gwas.DT[, .(count = uniqueN(Trait)), by = pid]
SNP.counter <- SNP.counter[count == max(count)]
gwas.DT <- gwas.DT[pid %in% SNP.counter$pid]


# Check if there is any column missing any value (e.g. LD block)
message("check ld.block NA")
gwas.DT[is.na(ld.block)]
gwas.DT <- na.omit(gwas.DT)

# Check gwas.DT
message("After keeping overlap & remove unwanted columns, we get:")
dim(gwas.DT)
summary(gwas.DT)

# Write overlaped SNPs 

message("writing overlapped SNPs")

fwrite(gwas.DT, "../manifest/overlap.gwas.DT.tsv.gz", sep="\t")

message("finished writing!")

