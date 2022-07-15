### Create LD matrices from 1000GP Phase III

## Background: To project datasets with confidence intervals we need to use information on SNP correlation, which we'll
## obtain from LD matrices. To generate those, we'll prepare LD matrices using 1000GP Phase III data on our 4483 non-zero SNPs.

# Guillermo Reales
# 2022-07-05

#############################################
######## Load libraries     #################
#############################################

library(data.table)
setDTthreads(20)
library(snpStats)
library(annotSnpStats)
#devtools::install_github("chr1swallace/annotSnpStats") 
library(magrittr)
library(snpStats)

#############################################
######## Helper functions   #################
#############################################

bdiag_with_dimnames <- function(lst) {
    ret <- Matrix::bdiag(lst)
    rownames(ret) <- colnames(ret) <- unlist(lapply(lst,rownames))
    ret
}



#############################################
######## Load files         #################
#############################################

# Prepare SNPs and manifest

load("../data/Astle_basis_1e4_p1e8.RData")
rot.pca <- basis.spca$loadings 
dim(rot.pca)# Note: 110646 SNPs in original, but only 4472 with at least one non-zero row entry
lite.rot.pca <- rot.pca[rowSums(rot.pca) != 0,]
bpids <- rownames(lite.rot.pca)

SNP.manifest <- fread("../data/BCB4_SNP_manifest_after_Astle_6M.tsv")
lSNP <- SNP.manifest[pid %in% bpids]

### Get reference files

rpath <- "/home/gr440/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3_GRCh38_30X/"
pops <- fread(paste0(rpath,"20130606_g1k_3202_samples_ped_population.txt"))
pops <- pops[Superpopulation == "EUR" & FatherID == 0 & MotherID == 0] # Select Europeans, remove children from trios

# Loop zone
# Apply only if not already available
if(!file.exists("../data/LD_bychr/LDmatrix_1e4s_1e8p_chr1.RDS")){

    lapply(1:22, function(i){
        SNPi <- lSNP[CHR38 == i, .(pid, CHR38, BP38, REF, ALT, ld.block)]
        bim <- fread(paste0(rpath,"chr",i,".bim"))
        bim[, idx:=1:nrow(bim)]
        bim  <- bim[V4 %in% SNPi$BP38,]
        merged  <- merge(bim, SNPi, 
        		      by.x=c("V4","V5","V6"),
        		      by.y=c("BP38","ALT", "REF"))
        if(any(duplicated(merged$V2))){
        		message("Removing ", nrow(merged[duplicated(merged$V2),]), " SNPs with duplicated names")
        		merged  <- merged[!duplicated(merged$V2),]
        	}
        idx.snps <- merged$idx
        
        fam <- fread(paste0(rpath,"chr",i,".fam"))
        # We'll use the following to filter SNPs at load time.
        # There's a bug in read.plink, so we'll filter by individuals later
        message("Loading panel file for chr", i, ".")
        panelfile <- read.plink(paste0(rpath, "chr", i), select.snps = idx.snps)
        sm <- annot.plink(panelfile)
        sm <- sm[ pops$SampleID, ] # is this the step to filter by individuals?
        colnames(sm) <- with(snps(sm), paste(chromosome,position,sep=":"))
        pids <- colnames(sm)
        
        snps.i <- SNPi[pid %in% pids]
        ss <- split(snps.i, snps.i$ld.block)
        LDm <- 	lapply(ss, function(snpsub){
        			sm.map <- match(snpsub$pid, pids)
        			r <- ld(sm[,sm.map],sm[,sm.map], stats="R")
        			r[is.na(r)] <- 0
        			r
        		      }) %>% bdiag_with_dimnames(.)
        message("Saving LD matrix for ", i, "...")
        saveRDS(LDm,file=paste0("../data/LD_bychr/LDmatrix_1e4s_1e8p_chr",i, ".RDS"))
        message("Done. LD matrix for ", i," ready!")
    
    })
}


# Join in a single matrix
LD <- lapply(1:22, function(x){
    t <- readRDS(paste0("../data/LD_bychr/LDmatrix_1e4s_1e8p_chr",x,".RDS"))
}) %>% bdiag_with_dimnames 

saveRDS(LD, file="../data/LD_bychr/LDmatrix_1e4s_1e8p_ALL.RDS")
