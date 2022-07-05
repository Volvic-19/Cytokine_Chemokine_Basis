
# This script uses the processed Chen, UKBB (Neale + PanUKBB),FinnGen and IL5 from Ferkingstad to create an intersection SNP list 
# Compared to the old version script, note that CHR19 and BP19 were kept for further mapping

library(data.table)
setDTthreads(15)
library(GenomicRanges)



setwd("~/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_building/")

# Load files

fk  <- fread("data/11071_1_IL5_IL_5.txt.gz") # Example Ferkingstad file (IL5)
pu  <- fread("Asthma_PanUKBB_processed-hg38.tsv.gz") # PanUKBB 
ne  <- fread("Asthma_Neale_processed-hg38.tsv.gz") # Neale
fg  <- fread("Asthma_FinnGen.tsv.gz") # FinnGen R6 example file

# change fk chrosome to character numbers
setnames(fk, old = c("Chrom","Pos","effectAllele","otherAllele"),new = c("CHR38","BP38","ALT","REF"))

fk[ , CHR38 := gsub("chr", "", CHR38)]

# Let's remove strange positions
stchr <- c(1:22, "X")


unique(fk$CHR38)
fk <- fk[ CHR38 %in% stchr ]

unique(pu$CHR38)
pu <- pu[CHR38 %in% stchr ]

unique(ne$CHR38)
ne <- ne[CHR38 %in% stchr ]

unique(fg$`#chrom`)
fg[,`#chrom`:=as.character(`#chrom`)][ `#chrom` == "23", `#chrom`:="X"]


fk[, pid:=paste(CHR38, BP38, sep=":")]
pu[, pid:=paste(CHR38, BP38, sep=":")]
ne[, pid:=paste(CHR38, BP38, sep=":")]
fg[, pid:=paste(`#chrom`, pos, sep=":")]

summary(fk)
summary(pu)
summary(ne)
summary(fg)

cor <- fk$pid
cor <- intersect(cor, pu$pid)
cor <- intersect(cor, ne$pid)
cor <- intersect(cor, fg$pid)

cor <- unique(cor)
length(cor)
# [1] 10429723
head(cor)

manifest <- fk[pid %in% cor][, .(pid, CHR38, BP38, REF, ALT)]

# Remove duplicates
manifest <- unique(manifest)
dups <- manifest[duplicated(manifest$pid), pid]
manifest <- manifest[!pid %in% dups] # from 10957466 to 10062388

# Remove SNPs with unallignable alleles
manifest[, alleles:=paste(REF,ALT, sep="/")]
manifest <- manifest[!alleles %in% c("A/T","T/A","G/C","C/G")]
manifest[,alleles:=NULL]
table(manifest$CHR38)

# Order by chr and bp
manifest <- manifest[, CHR38:=as.numeric(CHR38)][order(CHR38, BP38)] #8638423

# Remove NA caused by chromosome X 
manifest <- manifest[which(!is.na(manifest$CHR38)), ]

# Bring MacDonald et al. liftover blocks, available at
# https://raw.githubusercontent.com/jmacdon/LDblocks_GRCh38/master/data/EUR_LD_blocks.bed
bl  <- fread("~/rds/rds-cew54-wallace-share/Data/reference/lddetect/EUR_hg38_MacDonald/EUR_LD_blocks.bed")

blranges <- GRanges(seqnames=bl$chr, ranges=IRanges(start=bl$start, end=bl$stop), strand="*")

snpranges <- GRanges(seqnames=paste("chr",manifest$CHR38, sep=""), ranges=IRanges(start=manifest$BP38, end=manifest$BP38), strand="*")

manifest[, ld.block:=findOverlaps(snpranges, blranges, select='last')] # Assign blocks. This will respect the [start, stop) block nomenclature from lddetect.

# Add 'CHR19' 'BP19' by mapping to pu
map <- pu[,c('pid','CHR19','BP19')] 
dups <- map[duplicated(map$pid), pid]
map <- map[!pid %in% dups]

output <- merge(manifest, map, by = "pid", all.x = TRUE)
output <- output[order(CHR38, BP38)]
output <- output[nchar(ALT)==1 & nchar(REF)==1] 

# Remove "A/T" "T/A" "C/G" "G/C"
output[, REF_ALT:= paste(REF,ALT,sep='/')] 
rmREF_ALT <- c("A/T","T/A","C/G","G/C","A/A","T/T","C/C","G/G")
output <- output[!REF_ALT %in% rmREF_ALT] # 8023580 SNPs

# Remove MHC and associated LD blocks
mhc.blocks <- output[  CHR38 == 6 & BP38 > 20000000 & BP38 < 40000000, unique(ld.block) ]
output  <- output[ !ld.block %in% mhc.blocks ]

# 7938504 SNPs



# Save the new consensus manifest
fwrite(output, "./IL5_consensus_manifest_8M.tsv", sep="\t") 






