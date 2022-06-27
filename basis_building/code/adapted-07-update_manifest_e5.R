# This script is to update the IL5 <- 6M manifest with overlap results. After updating, only SNPs with shrinkage value > 1e-5 will be kept in IL5 manifest.
# Load libraries
library(data.table)
setDTthreads(18)

#########################################################


# Update manifest

SNP.manifest <- fread("../manifest/IL5_consensus_manifest_6M.tsv", tmpdir = "tmp")
e5.DT <- fread("../manifest/shrinkage.e5.DT.tsv.gz", tmpdir = "tmp")

dim(SNP.manifest)
# [1] 6161058       9

SNP.manifest <- SNP.manifest[pid %in% unique(e5.DT$pid)]

dim(SNP.manifest)
# [1] 673318      9

fwrite(SNP.manifest, "../manifest/IL5_consensus_manifest_6M_e5.tsv", sep="\t")









