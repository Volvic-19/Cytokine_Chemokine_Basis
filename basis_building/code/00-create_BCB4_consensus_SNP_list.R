# This script will use the processed Chen, UKBB (Neale + PanUKBB), and FinnGen and create an intersection SNP list 


library(data.table)
setDTthreads(15)
library(GenomicRanges)

# Load files
ch  <- fread("~/rds/rds-cew54-basis/02-Processed/EOSC_Chen_32888493_6-hg38.tsv.gz") # Example Chen file
pu  <- fread("Asthma_PanUKBB_processed-hg38.tsv.gz") # PanUKBB 
ne  <- fread("Asthma_Neale_processed-hg38.tsv.gz") # Neale
fg  <- fread("Asthma_FinnGen.tsv.gz") # FinnGen R6 example file

# Let's remove strange positions
stchr <- c(1:22, "X")

unique(ch$CHR38)
ch <- ch[ CHR38 %in% stchr ]

unique(pu$CHR38)
pu <- pu[CHR38 %in% stchr ]

unique(ne$CHR38)
ne <- ne[CHR38 %in% stchr ]

unique(fg$`#chrom`)
fg[,`#chrom`:=as.character(`#chrom`)][ `#chrom` == "23", `#chrom`:="X"]


ch[, pid:=paste(CHR38, BP38, sep=":")]
pu[, pid:=paste(CHR38, BP38, sep=":")]
ne[, pid:=paste(CHR38, BP38, sep=":")]
fg[, pid:=paste(`#chrom`, pos, sep=":")]

summary(ch)
summary(pu)
summary(ne)
summary(fg)

cor <- ch$pid
cor <- intersect(cor, pu$pid)
cor <- intersect(cor, ne$pid)
cor <- intersect(cor, fg$pid)

cor <- unique(cor)
length(cor)
# [1] 9800900
head(cor)

manifest <- ch[pid %in% cor][, .(pid, CHR38, BP38, REF, ALT)]

# Remove duplicates
manifest <- unique(manifest)
dups <- manifest[duplicated(manifest$pid), pid]
manifest <- manifest[!pid %in% dups]

# Remove SNPs with unallignable alleles
manifest[, alleles:=paste(REF,ALT, sep="/")]
manifest <- manifest[!alleles %in% c("A/T","T/A","G/C","C/G")]
manifest[,alleles:=NULL]
table(manifest$CHR38)

# Order by chr and bp
manifest <- manifest[, CHR38:=as.numeric(CHR38)][order(CHR38, BP38)]


# Bring MacDonald et al. liftover blocks, available at
# https://raw.githubusercontent.com/jmacdon/LDblocks_GRCh38/master/data/EUR_LD_blocks.bed
bl  <- fread("~/rds/rds-cew54-wallace-share/Data/reference/lddetect/EUR_hg38_MacDonald/EUR_LD_blocks.bed")

blranges <- GRanges(seqnames=bl$chr, ranges=IRanges(start=bl$start, end=bl$stop), strand="*")

snpranges <- GRanges(seqnames=paste("chr",manifest$CHR38, sep=""), ranges=IRanges(start=manifest$BP38, end=manifest$BP38), strand="*")

manifest[, ld.block:=findOverlaps(snpranges, blranges, select='last')] # Assing blocks. This will respect the [start, stop) block nomenclature from lddetect.

# Save the new consensus manifest
fwrite(manifest, "../data/BCB4_consensus_manifest_8M.tsv", sep="\t")






