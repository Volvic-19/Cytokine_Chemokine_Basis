# This script is to make side-by-side manhattan plot of IL1B from Ferkingstad and Ahola datasets. And M-CSF from FK, AH, Folkersen datasets
# We noticed that IL1B is an outlier in many PCs and intend to see if there is any difference in the raw summary
# statistics data.

library(data.table)
setDTthreads(10)
library(karyoploteR)
library(GenomicRanges)

################################  IL1-B  ##################################################

# load IL1-B data from Ferkingstad and Ahola
IL1B_fk <- fread("../../basis_building/raw_fk_data/3037_62_IL1B_IL_1b.txt.gz")
IL1B_ah <- fread("~/rds/rds-cew54-basis/02-Processed/IL1B_AholaOlli_27989323_1-hg38.tsv.gz")

# get ride of "chrx" from IL1B_fk
# get ride of odd CHR38 (such as "11_KI270927v1_alt" instead of "11") in IL1B_ah
IL1B_fk <- IL1B_fk[Chrom != "chrX"]
stchr <- as.character(1:22)
IL1B_ah <- IL1B_ah[CHR38 %in% stchr]

# plot data has to be in Grange form
s1 <- GRanges(seqnames=IL1B_fk$Chrom, ranges=IRanges(start=IL1B_fk$Pos, end=IL1B_fk$Pos), strand="*")
s1$pval <- IL1B_fk$Pval
s2 <- GRanges(seqnames=paste0("chr",IL1B_ah$CHR38), ranges=IRanges(start=IL1B_ah$BP38, end=IL1B_ah$BP38), strand="*")
s2$pval <- IL1B_ah$P

# make manhattan plot
ids <- c("Ferkinstad_IL1B","Ahola_IL1B")

png(filename="IL1B.FKvsAH.manhattanplot.png", width=12, height=8, units="in", res=300)
	kp <- plotKaryotype(plot.type=4)
	kpAddLabels(kp, labels = ids[1], srt=90, pos=3, r0=0.5, r1=1, cex=1.8, label.margin = 0.025)
	kpAxis(kp, ymin=0, ymax=15, r0=0.5)
	kp <- kpPlotManhattan(kp, data=s1,  r0=0.5, r1=1, ymax=10, points.col = "brewer.set1")
	kpAddLabels(kp, labels = ids[2], srt=90, pos=3, r0=0, r1=0.5, cex=1.8, label.margin = 0.025)
	kpAxis(kp, ymin=0, ymax=15, r0=0.5, r1=0, tick.pos = c(5,10))
	kp <- kpPlotManhattan(kp, data=s2, r0=0.5, r1=0, ymax=10, points.col = "2blues")
dev.off()	


################################  MCS-F  ##################################################

# load MCSF data from Ferkingstad , Ahola , Folkersen, Hoglund
MCSF_fk <- fread("../../basis_building/raw_fk_data/3738_54_CSF1_CSF_1.txt.gz")
MCSF_ah <- fread("~/rds/rds-cew54-basis/02-Processed/MCSF_AholaOlli_27989323_1-hg38.tsv.gz")
MCSF_hoglund <- fread("~/rds/rds-cew54-basis/02-Processed/MCSF_Hoglund_31727947_1-hg38.tsv.gz")
MCSF_folkersen <- fread("~/rds/rds-cew54-basis/02-Processed/MCSF_Folkersen_33067605_1-hg38.tsv.gz")

# get ride of "chrx" from fk
# get ride of odd CHR38 (such as "11_KI270927v1_alt" instead of "11") in ah, folkersen, hoglund
MCSF_fk <- MCSF_fk[Chrom != "chrX"]
stchr <- as.character(1:22)
MCSF_ah <- MCSF_ah[CHR38 %in% stchr]
MCSF_hoglund <- MCSF_hoglund[CHR38 %in% stchr]
MCSF_folkersen <- MCSF_folkersen[CHR38 %in% stchr]

# plot data has to be in Grange form
s1 <- GRanges(seqnames=MCSF_fk$Chrom, ranges=IRanges(start=MCSF_fk$Pos, end=MCSF_fk$Pos), strand="*")
s1$pval <- MCSF_fk$Pval
s2 <- GRanges(seqnames=paste0("chr",MCSF_ah$CHR38), ranges=IRanges(start=MCSF_ah$BP38, end=MCSF_ah$BP38), strand="*")
s2$pval <- MCSF_ah$P
s3 <- GRanges(seqnames=paste0("chr",MCSF_hoglund$CHR38), ranges=IRanges(start=MCSF_hoglund$BP38, end=MCSF_hoglund$BP38), strand="*")
s3$pval <- MCSF_hoglund$P
s4 <- GRanges(seqnames=paste0("chr",MCSF_folkersen$CHR38), ranges=IRanges(start=MCSF_folkersen$BP38, end=MCSF_folkersen$BP38), strand="*")
s4$pval <- MCSF_folkersen$P


# make manhattan plot
ids <- c("Ferkingstad_MCSF","Folkersen_MCSF")

png(filename="MCSF.FKvsFolkersen.manhattanplot.png", width=12, height=8, units="in", res=300)
	kp <- plotKaryotype(plot.type=4)
	kpAddLabels(kp, labels = ids[1], srt=90, pos=3, r0=0.5, r1=1, cex=1.8, label.margin = 0.025)
	kpAxis(kp, ymin=0, ymax=15, r0=0.5)
	kp <- kpPlotManhattan(kp, data=s1,  r0=0.5, r1=1, ymax=10, points.col = "brewer.set1")
	kpAddLabels(kp, labels = ids[2], srt=90, pos=3, r0=0, r1=0.5, cex=1.8, label.margin = 0.025)
	kpAxis(kp, ymin=0, ymax=15, r0=0.5, r1=0, tick.pos = c(5,10))
	kp <- kpPlotManhattan(kp, data=s4, r0=0.5, r1=0, ymax=10, points.col = "2blues")
dev.off()	




################################  CXCL12  ##################################################

# load CXCL12 data from Ferkingstad , Ahola
CXCL12_fk <- fread("../../basis_building/raw_fk_data/3516_60_CXCL12_SDF_1.txt.gz")
CXCL12_ah <- fread("~/rds/rds-cew54-basis/02-Processed/SDF1a_AholaOlli_27989323_1-hg38.tsv.gz")

# get ride of "chrx" from IL1B_fk
# get ride of odd CHR38 (such as "11_KI270927v1_alt" instead of "11") in IL1B_ah
CXCL12_fk <- CXCL12_fk[Chrom != "chrX"]
stchr <- as.character(1:22)
CXCL12_ah <- CXCL12_ah[CHR38 %in% stchr]

# plot data has to be in Grange form
s1 <- GRanges(seqnames=CXCL12_fk$Chrom, ranges=IRanges(start=CXCL12_fk$Pos, end=CXCL12_fk$Pos), strand="*")
s1$pval <- CXCL12_fk$Pval
s2 <- GRanges(seqnames=paste0("chr",CXCL12_ah$CHR38), ranges=IRanges(start=CXCL12_ah$BP38, end=CXCL12_ah$BP38), strand="*")
s2$pval <- CXCL12_ah$P

# make manhattan plot
ids <- c("Ferkinstad_CXCL12","Ahola_CXCL12")

png(filename="../Plots/manhattan/CXCL12.FKvsAH.manhattanplot.png", width=12, height=8, units="in", res=300)
	kp <- plotKaryotype(plot.type=4)
	kpAddLabels(kp, labels = ids[1], srt=90, pos=3, r0=0.5, r1=1, cex=1.8, label.margin = 0.025)
	kpAxis(kp, ymin=0, ymax=15, r0=0.5)
	kp <- kpPlotManhattan(kp, data=s1,  r0=0.5, r1=1, ymax=10, points.col = "brewer.set1")
	kpAddLabels(kp, labels = ids[2], srt=90, pos=3, r0=0, r1=0.5, cex=1.8, label.margin = 0.025)
	kpAxis(kp, ymin=0, ymax=15, r0=0.5, r1=0, tick.pos = c(5,10))
	kp <- kpPlotManhattan(kp, data=s2, r0=0.5, r1=0, ymax=10, points.col = "2blues")
dev.off()	

# CXCL9


