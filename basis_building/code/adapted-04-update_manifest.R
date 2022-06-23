# This script is to update the IL5 manifest with overlap results. After updating, only SNPs kept in all 104 traits will be kept in IL5 manifest.
# Load libraries
library(data.table)
setDTthreads(18)

#########################################################


# Update manifest

SNP.manifest <- fread("../manifest/IL5_consensus_manifest_8M.tsv", tmpdir = "tmp")
gwas.DT <- fread("../manifest/overlap.gwas.DT.tsv.gz", tmpdir = "tmp")

message("before filtering, we have:")
dim(SNP.manifest)

SNP.manifest <- SNP.manifest[pid %in% unique(gwas.DT$pid)]

message("start writing!")
fwrite(SNP.manifest, "../manifest/IL5_consensus_manifest_6M.tsv", sep="\t")
message("finished writing!")


message("after filtering, we have:")
dim(SNP.manifest)

warnings()







