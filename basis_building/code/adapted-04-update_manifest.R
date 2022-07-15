# This script is to update the CCL8 initial manifest with overlap results. After update, only SNPs shared by all 40 traits will be kept in CCL8 manifest.
# Slurm: sbatch slurm_04_update_manifest_FK

# Load libraries
library(data.table)
setDTthreads(18)

#########################################################


# Update manifest

SNP.manifest <- fread("../manifest/CCL8_consensus_manifest_8M.tsv", tmpdir = "tmp")
gwas.DT <- fread("../manifest/overlap.gwas.DT.tsv.gz", tmpdir = "tmp")

message("before filtering, we have:") #[1] 7833926
dim(SNP.manifest)

SNP.manifest <- SNP.manifest[pid %in% unique(gwas.DT$pid)]

message("start writing!")
fwrite(SNP.manifest, "../manifest/CCL8_consensus_manifest_6M.tsv", sep="\t")
message("finished writing!")


message("after filtering, we have:") #[1] 6147256       8
dim(SNP.manifest)

warnings()







