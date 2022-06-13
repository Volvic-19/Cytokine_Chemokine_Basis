# Preprocess consensus files 
# This will prepare UKBB files for liftover. FinnGen is already in hg38 so no need to liftover them.

library(data.table)
setDTthreads(10)

pu <- fread("Asthma_PanUKBB.tsv.gz")
ne <- fread("Asthma_Neale.tsv.gz")

pu  <- pu[,.(chr, pos, ref,alt, beta_meta, se_meta,pval_meta)] 
setnames(pu, c("beta_meta", "se_meta", "pval_meta"), c("BETA","SE", "P"))
pu

ne[, c("CHR", "BP", "REF", "ALT") := tstrsplit(variant, ":", fixed=TRUE)][, variant :=NULL]
ne <- ne[, c("CHR","BP", "REF","ALT","beta",  "se","pval")]

fwrite(pu, "Asthma_PanUKBB_processed.tsv.gz", sep="\t", na=NA, quote=FALSE, row.names=FALSE)
fwrite(ne, "Asthma_Neale_processed.tsv.gz", sep="\t", na=NA, quote=FALSE, row.names=FALSE)
