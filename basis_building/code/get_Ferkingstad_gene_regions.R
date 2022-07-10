##### Get 100kb regions around genes in Ferkingstad

# Guillermo Reales
# 2022-07-08

# Background: We want to explore the presence/absence of eQTLs near the genes that encode
# the selected proteins in Ferkingstad, as that may be an indicator for aptamer accuracy.
# We'll search for GWAS hits within 100kb from the limits of the genes, so to do that we'll
# extract gene coordinates from biomaRt

######################
### Load packages  ###
######################


library(data.table)
setDTthreads(30)
library(biomaRt)
library(magrittr)

######################
### Process data  ####
######################

# Load files
f <- fread("../../Pre-work/Selected_traits/Ferkingstad_basistraits_withduplicates.tsv")

gid <- unique(f$HGNC.ID)

# Retrieve data from biomaRt
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)


res <- getBM(attributes = c("hgnc_id", "ensembl_gene_id", "hgnc_symbol","chromosome_name", "start_position", "end_position"),
            filters = c("hgnc_id"), values = gid, mart = mart)
res <- as.data.table(res)

# Some checks
setdiff(gid, res$hgnc_id) # No missing values
missing <- setdiff(f$Ensembl.Gene.ID, res$ensembl_gene_id)
f[ Ensembl.Gene.ID %in% missing, .(FileName, Gene, HGNC.ID)]
chm <- f[ Ensembl.Gene.ID %in% missing, HGNC.ID] # Cross-check
res[ hgnc_id %in% chm]

chrs <- c(paste0(1:22), "X") # Keep valid chromosomes only
res <- res[chromosome_name %in% chrs]

m <- merge(res[,.(hgnc_id, hgnc_symbol, chromosome_name, start_position, end_position)], f, by.x=c("hgnc_id", "hgnc_symbol"), by.y=c("HGNC.ID", "Gene"))

# Calculate regions
m[, region_start:= start_position - 10^5][ region_start < 0, region_start:= 0][, region_end:=end_position + 10^5]
all(m$region_start < m$region_end) # TRUE

# Save
fwrite(m, "../../Pre-work/Selected_traits/Ferkingstad_basistraits_withregions.tsv", sep="\t")


## Now we'll loop over all 104 ferkingstad files to search for peaks in these regions
files <- m$FileName
path <- "../raw_fk_data/"


peaks <- lapply(files, function(file) {
    message("Working on ", file, " (", which(file == files), "/", length(files),").")
    ff <- fread(paste0(path, file), tmpdir="tmp")
    mfile <- m[FileName == file, .(chromosome_name, region_start, region_end)]
    s <- ff[Chrom == paste0("chr", mfile$chromosome_name) & Pos >= mfile$region_start & Pos <= mfile$region_end & minus_log10_pval >= 8]
    s[, FileName:= file]
    s
})
    
ps <- rbindlist(peaks)
unique(ps$FileName) # 40 unique files

m[, `Protein (short name)`] %>% unique # 100 unique proteins
m[FileName %in% unique(ps$FileName), `Protein (short name)`] %>% unique # 38 unique proteins
    
fwrite(ps, "../../Pre-work/Selected_traits/Ferkingstad_peaks.tsv", sep="\t")

