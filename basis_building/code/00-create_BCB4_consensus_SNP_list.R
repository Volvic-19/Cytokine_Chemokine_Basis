# This script will use the processed Astle, UKBB (Neale + PanUKBB), and FinnGen and create an intersection SNP list.
# It will also align the manifest to a 1000GP reference manifest


library(data.table)
setDTthreads(20)
library(GenomicRanges)
library(IMDtools)


# Load files
ch  <- fread("~/rds/rds-cew54-basis/02-Processed/EOSC_Astle_27863252_1-hg38.tsv.gz") # Example Astle file
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
# [1] 10227870
head(cor)

manifest <- ch[pid %in% cor][, .(pid, CHR38, BP38, REF, ALT)]
manifest[, CHR38:=as.character(CHR38)]
hg19 <- pu[pid %in% cor, .(CHR19, BP19, CHR38, BP38)] # Get hg19 coords
manifest <- merge(manifest, hg19, by=c("CHR38", "BP38"))



# Remove duplicates
manifest <- unique(manifest)
dups <- manifest[duplicated(manifest$pid), pid]
# Many of these duplicates correspond to multiple alleles, so we'll keep nchar = 1 only
manifest <- manifest[ nchar(REF) == 1 & nchar(ALT) == 1]
dups <- manifest[duplicated(manifest$pid), pid] 
# The remaining 28 are either mapped to different hg19 or have different alleles. Remove them
#manifest[pid %in% dups]
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

# There are some SNPs missing LD.block, we need to remove them
manifest[is.na(ld.block)]
manifest <- na.omit(manifest)

# We'll remove the MHC, too. We'll use the coordinates from Dilthey, 2021 (https://www.sciencedirect.com/science/article/pii/S1357272520301990), plus a window, so we'll remove blocks containing 20-40Mb
mhc.blocks <- manifest[  CHR38 == 6 & BP38 > 20000000 & BP38 < 40000000, unique(ld.block) ]
manifest  <- manifest[ !ld.block %in% mhc.blocks ]
nrow(manifest)
# [1] 8390337

# Import refmanifest and use those alleles
refman <- fread("../data/Reference_1000GP_manifest.tsv.gz")

mm <- merge(manifest, refman[,.(pid, REF, ALT)], by="pid", suffixes=c(".ast", ".ref"))
all(mm$REF.ast == mm$REF.ref)
mm[ REF.ast != REF.ref] # not aligned
rm(mm) # Very big, let's remove it for now

# Since our aligner is meant for summary statistics datasets, we need to fool it by adding some of the "minimum columns required"
manifest[, c("BETA", "SE", "P"):= 0]
refman[,alleles:=NULL]
d <- copy(manifest)
rm(manifest)

m.al <- g.align(ds = d, manifest = refman) # here refman will be our manifest, and the manifest will be the dataset to align

m.al[,c("BETA", "SE", "P"):=NULL]
m.al <- m.al[order(CHR38,BP38)]

# Save the new consensus manifest
fwrite(m.al, "../data/BCB4_consensus_manifest_8M_aligned.tsv", sep="\t")






