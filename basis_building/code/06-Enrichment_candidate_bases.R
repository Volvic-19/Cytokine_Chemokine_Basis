## Perform enrichment on driver SNPs from bases with several parameters.


## Background: After taking a look at several statistics of bases created using different parameters, we chose to perform enrichment for the nearest genes associated with this SNPs to see if this can help us make an informed decision on the parameters.
# We created basis using 4 parameters:1e-9, 3e-8, 1e-7, 1e-5. Then we extracted the driver SNPs for each PC and saved them. We then applied a custom script to extract the OpenTargetsGenetics (OTG) annotations for these SNPs.

#############################################
######## Load libraries     #################
#############################################

library(data.table)
setDTthreads(8)
library(magrittr)
library(gprofiler2)

#############################################
########     Load files     #################
#############################################

snp <- fread("../data/drivers_8t_1e5s_9to5.tsv") # SNP driver SNPs, by PC and param
length(unique(snp$pid))
# [1] 5112
otg.o <- readRDS("~/rds/rds-cew54-basis/People/CHRIS/basis-blood/bloodbasisv4+1-driver-snps-otg-annotations.RDS")
length(unique(otg.o[,pid]))
# [1] 4986
# Most SNPs have annotations, but not all

# Load manifest
man <- fread("../data/BCB4_SNP_manifest_after_Chen_8M.tsv")
mm  <- man[pid %in% unique(snp$pid)]

otg <- otg.o[ otg.o[, .I[tss_score == max(tss_score)], by=pid]$V1]  # We're interested in the nearest genes, so we'll choose those annotations with highest tss_score.
length(unique(otg[,pid]))
# [1] 4986


#############################################
########     Process data       #############
#############################################

# Filter out the SNPs in the OTG dataset that do not correspond to those in the manifest
mm[, "alleles":=paste(REF,ALT, sep="/")] 
otg.test <- merge(otg, mm[, .(pid,alleles)], by=c("pid")) 
length(unique(otg.test[,pid]))
# [1] 4986
otg.test[, c("allele1", "allele2"):=list(paste(ref_allele, alt_allele, sep="/"), paste(alt_allele, ref_allele, sep="/"))]
otg.test2  <- otg.test[allele1 %in% alleles | allele2 %in% alleles]
length(unique(otg.test2[,pid]))
# [1] 4971
dffsnp <- setdiff(otg.test$pid, otg.test2$pid)
mm[pid %in% dffsnp]
otg[pid %in% dffsnp, .(pid, ref_allele, alt_allele)] # Some SNPs are non biallelic, and hence don't match. Not too many, so we can exclude them.

otg2 <- otg.test2
otg2[, c("alleles", "allele1","allele2"):=NULL]


# Function for gprofiler2
run_gprofiler2 <- function(x, annot){

	params <- unique(x$Param)
	flist <- lapply(params, function(i){
		message("Working on parameter: ",i)
		p.ds <- x[Param == i]
		PCs <- unique(p.ds$PC)
		pclist <- lapply(PCs, function(y){
			message("Working on ", y)
			pids <- p.ds[PC == y, pid] %>% unique
			geneids <- annot[pid %in% pids, gene_id] %>% unique
			gpr <- gost(query = geneids, organism = "hsapiens", ordered_query = FALSE, significant = TRUE, exclude_iea = TRUE)
			gpr <- gpr$result %>% as.data.table
			gpr[, c("Param", "PC"):=list(i, y)]
		})
		pclist <- rbindlist(pclist, fill=TRUE)
	})
	flist <- rbindlist(flist, fill=TRUE)
}

enr <- run_gprofiler2(snp, annot =otg2)

enr

saveRDS(enr, "../data/Enrichment_analysis_results_8t.RDS")

