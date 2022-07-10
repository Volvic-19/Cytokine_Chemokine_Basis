
# This script uses the processed  UKBB (Neale + PanUKBB),FinnGen, IL5 from Ferkingstad and reference_1000GP_manifest to create an intersection SNP list 
# Compared to the old version script, note that CHR19 and BP19 were kept for further mapping
# 2022/07/07 : added the reference_1000 genome alignment step

library(data.table)
setDTthreads(15)
library(GenomicRanges)
library(IMDtools)


setwd("~/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_building/")

# Load files

fk  <- fread("raw_fk_data/11071_1_IL5_IL_5.txt.gz") # Example Ferkingstad file (IL5)
pu  <- fread("Asthma_PanUKBB_processed-hg38.tsv.gz") # PanUKBB 
ne  <- fread("Asthma_Neale_processed-hg38.tsv.gz") # Neale
fg  <- fread("Asthma_FinnGen.tsv.gz") # FinnGen R6 example file

# change fk chrosome to character numbers
setnames(fk, old = c("Chrom","Pos","effectAllele","otherAllele"),new = c("CHR38","BP38","ALT","REF"))

fk[ , CHR38 := gsub("chr", "", CHR38)]

# Let's remove strange positions and X chromosome
stchr <- as.character(c(1:22))


unique(fk$CHR38)
fk <- fk[ CHR38 %in% stchr ]

unique(pu$CHR38)
pu <- pu[CHR38 %in% stchr ]

unique(ne$CHR38)
ne <- ne[CHR38 %in% stchr ]

unique(fg$`#chrom`)
fg[,`#chrom`:=as.character(`#chrom`)][ `#chrom` == "23", `#chrom`:="X"]
fg <- fg[ `#chrom` %in% stchr ]


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
# [1] 10138616


# Add hg19 from pu
manifest <- fk[pid %in% cor][, .(pid, CHR38, BP38, REF, ALT)]
hg19 <- pu[pid %in% cor, .(CHR19, BP19, CHR38, BP38)] # Get hg19 coords
manifest <- merge(manifest, hg19, by=c("CHR38", "BP38"))

# Remove duplicated pid from manifest
manifest <- unique(manifest) #10656457
dups <- manifest[duplicated(manifest$pid), pid]
manifest <- manifest[!pid %in% dups] # 9778559


# Remove SNPs with unallignable alleles
manifest[, alleles:=paste(REF,ALT, sep="/")]
manifest <- manifest[nchar(alleles) == 3]
manifest <- manifest[!alleles %in% c("A/T","T/A","G/C","C/G","A/A","T/T","C/C","G/G")]
manifest[,alleles:=NULL]
table(manifest$CHR38) # 8023573

# Order by chr and bp
manifest <- manifest[, CHR38:=as.numeric(CHR38)][order(CHR38, BP38)] # 8023573

# Bring MacDonald et al. liftover blocks, available at
# https://raw.githubusercontent.com/jmacdon/LDblocks_GRCh38/master/data/EUR_LD_blocks.bed
bl  <- fread("~/rds/rds-cew54-wallace-share/Data/reference/lddetect/EUR_hg38_MacDonald/EUR_LD_blocks.bed")

blranges <- GRanges(seqnames=bl$chr, ranges=IRanges(start=bl$start, end=bl$stop), strand="*")

snpranges <- GRanges(seqnames=paste("chr",manifest$CHR38, sep=""), ranges=IRanges(start=manifest$BP38, end=manifest$BP38), strand="*")

manifest[, ld.block:=findOverlaps(snpranges, blranges, select='last')] # Assign blocks. This will respect the [start, stop) block nomenclature from lddetect.

# There are some SNPs missing LD.block, we need to remove them
manifest[is.na(ld.block)]
manifest <- na.omit(manifest) # 8023572

# Remove MHC and associated LD blocks
mhc.blocks <- manifest[  CHR38 == 6 & BP38 > 20000000 & BP38 < 40000000, unique(ld.block) ]
manifest  <- manifest[ !ld.block %in% mhc.blocks ] # 7938496 SNPs

# Import refmanifest and use those alleles. This step is to make sure
# our manifest is in the same REF and ALT as reference_100GP, as later we need to 
# use LD metrices in computation of significance in projection
# Note: refman can have different alleles at the same pid! so when loading LD metrices later, 
#       pid, REF, ALT have to be all taken into consideration to map!

refman <- fread("manifest/Reference_1000GP_manifest.tsv.gz")

# Since our aligner is meant for summary statistics datasets, we need to fool it by adding some of the "minimum columns required"
manifest[, REF_ALT := alleles][, alleles:=NULL]
manifest[, c("BETA", "SE", "P"):= 0]
refman[,alleles:=NULL]
d <- copy(manifest) #7938496
rm(manifest)

m.al <- g.align(ds = d, manifest = refman) # here refman will be our manifest, and the manifest will be the dataset to align
#  1160494  to flip,  0  to find their complement, and  0  to find their reverse complement.
#  Unfortunately,  0  SNPs are ambiguous, and  38664  were impossible to align. These will be removed now.

m.al[,c("BETA", "SE", "P"):=NULL]
m.al <- m.al[order(CHR38,BP38)] # FOR example, 1:818802 flipped from "A/G" to "G/A"
m.al[, REF_ALT:=NULL]

# Sanity check: if all REF, ALT in m.al is the same as what is in Reference_1000GP_manifest
#mm <- merge(m.al, refman[,.(pid, REF, ALT)], by="pid", suffixes=c(".ast", ".ref"))
#all(mm$REF.ast == mm$REF.ref)
#mm[ REF.ast != REF.ref] # not aligned
#rm(mm) # Very big, let's remove it for now
any(duplicated(m.al$pid))
# [1] FALSE

# Save the new consensus manifest
fwrite(m.al, "manifest/IL5_consensus_manifest_8M.tsv", sep="\t") 

# 7833926 SNP





