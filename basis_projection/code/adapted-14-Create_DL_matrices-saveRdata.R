# This script is to compute LD matrices of sparse SNPs in our Ferkingstad-40-basis. Then save updated projection function
#   together with key stats into Rdata
# it is adapted from 10-Create_LD_matrices.R. 

## Background: To project datasets with confidence intervals we need to use information on SNP correlation, which we'll
## obtain from LD matrices. To generate those, we'll prepare LD matrices using 1000GP Phase III data on our 5510 non-zero SNPs.

# 2022-07-13

# https://genomicsbootcamp.github.io/book/genotype-files-in-practice.html#fam-file---info-on-individuals
# Read above link for introductions to .bim, .fam files

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
    ret <- Matrix::bdiag(lst) # bdiag() function Construct a Block Diagonal Matrix
    rownames(ret) <- colnames(ret) <- unlist(lapply(lst,rownames))
    ret
}


#############################################
########     Load files     #################
#############################################

# Load sparse manifest and sparse rotation matrix
load("../Rdata/cytokine_sparseSNP_40basis-noLD-noProjFun.RData")
README

#                                                                                 rot.pca 
#         "sparse rotation matrix (5519 * 40)(SNPs with at least one non-zero row entry)" 
#                                                                                 use.pca 
#                         "logical value matrix. hard thresholded rot.pca!=0 (5519 * 40)" 
#                                                                            beta.centers 
#              "vector of column means of 40 Ferkingstad +1 control at sparse basis snps" 
#                                                                               shrinkage 
#                                            "vector of w/sigma_maf at sparse basis snps" 
#                                                                            SNP.manifest 
#                                                       "data.table of sparse basis snps" 
#                                                                                   stats 
# "info on the threshold quantile and number of SNPs per PC in the sparse basis (40 * 3)" 

bpids <- rownames(rot.pca) # 5519 pids of sparse SNPs

# Load reference files
rpath <- "~/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3_GRCh38_30X/"
pops <- fread(paste0(rpath,"20130606_g1k_3202_samples_ped_population.txt")) # 3202 * 7
pops <- pops[Superpopulation == "EUR" & FatherID == 0 & MotherID == 0] # Select Europeans, remove children from trios
#Africans (AFR), Admixed Americans (AMR), East Asians (EAS), Europeans (EUR) and South Asians (SAS)
# 525 * 7

# Loop to 
setwd(rpath) # not sure why read.plink() can not work if I pass paste0(rpath, "chr", i, ".bed") to it.

# Apply only if not already available
lSNP <- SNP.manifest

    lapply(1:22, function(i){
        SNPi <- lSNP[CHR38 == i, .(pid, CHR38, BP38, REF, ALT, ld.block)]
        bim <- fread(paste0(rpath,"chr",i,".bim"))
        bim[, idx:=1:nrow(bim)]
        bim  <- bim[V4 %in% SNPi$BP38,]
        merged  <- merge(bim, SNPi,  #431
        		      by.x=c("V4","V5","V6"),
        		      by.y=c("BP38","ALT", "REF"))
        if(any(duplicated(merged$V2))){
        		message("Removing ", nrow(merged[duplicated(merged$V2),]), " SNPs with duplicated names")
        		merged  <- merged[!duplicated(merged$V2),]
        	}
        idx.snps <- merged$idx #original snp row position in the bim file
        
        fam <- fread(paste0(rpath,"chr",i,".fam")) # 6 cols: family id, sample id(or individual id)
        # We'll use the following to filter SNPs at load time.
        # There's a bug in read.plink, so we'll filter by individuals later
        message("Loading panel file for chr", i, ".")
        panelfile <- read.plink(paste0("chr", i, ".bed"), select.snps = idx.snps) # add .bed extension
        sm <- annot.plink(panelfile)
        sm <- sm[ pops$SampleID, ] #filter by individuals?
        colnames(sm) <- with(snps(sm), paste(chromosome,position,sep=":")) #snps is in sm
        pids <- colnames(sm)
        
        snps.i <- SNPi[pid %in% pids] #431 rows, when i=1
        ss <- split(snps.i, snps.i$ld.block)
        LDm <- 	lapply(ss, function(snpsub){
        			sm.map <- match(snpsub$pid, pids)
        			r <- ld(sm[,sm.map],sm[,sm.map], stats="R")
        			r[is.na(r)] <- 0
        			r
        		      }) %>% bdiag_with_dimnames(.)
        message("Saving LD matrix for ", i, "...")
        saveRDS(LDm,file=paste0("~/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_projection/LDmetrices/LDmetrices_Ferkingstad40_PC1-40sparse_chr",i, ".RDS"))
        message("Done. LD matrix for ", i," ready!")
    
    })

setwd("/home/qz284/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_projection/code/")

# Join in a single matrix
LD <- lapply(1:22, function(x){
    t <- readRDS(paste0("../LDmetrices/LDmetrices_Ferkingstad40_PC1-40sparse_chr",x,".RDS"))
}) %>% bdiag_with_dimnames 

saveRDS(LD, file="../LDmetrices/LDmetrices_FK40_PC1-40sparse_ALL.RDS") # 5519 * 5519



############### Update project_sparse function. Save it together with other essential stats for projection into a Rdata #####


#############################################
######## Helper functions   #################
#############################################

# Update projection function. This time we include the parts to compute confidence interval using LD metrices
# This project_sparse function is identical to that in cupcake.

project_sparse <- function (beta, seb, pids) {
    if (length(beta) != length(seb) || length(beta) != length(pids) || 
        !length(beta)) 
        stop("arguments must be equal length vectors > 0")
    if (!all(pids %in% SNP.manifest$pid)) 
        stop("all pids must be members of sparse basis (SNP.manifest$pid)")
    if (length(pids) < 0.95 * nrow(rot.pca)) 
        warning("more than 5% sparse basis snps missing")
    b <- beta * shrinkage[pids] - beta.centers[pids]
    proj <- b %*% rot.pca[pids, ]
    # ATTENTION: seb needs to be NA free! otherwise will introduce NA into the matrix computation
    v <- seb * shrinkage[pids] * rot.pca[pids, ] # v is not variance, it is an intermediate term, like a projection of seb, parallel to b and proj
    var.proj <- t(v) %*% LD[pids, pids] %*% v # why NA?
    ctl <- (-beta.centers[pids]) %*% rot.pca[pids, ]
    delta <- (proj - ctl)[1, ]
    chi2 <- (t(delta) %*% solve(var.proj) %*% delta)[1, 1]
    ret <- data.table::data.table(PC = colnames(proj), proj = proj[1, 
        ], var.proj = Matrix::diag(var.proj), delta = delta, 
        p.overall = stats::pchisq(chi2, df = 14, lower.tail = FALSE)) # NOTE: remember to update df as number of PCs!
    ret$z = ret$delta/sqrt(ret$var.proj)
    ret$p = stats::pnorm(abs(ret$z), lower.tail = FALSE) * 2
    copy(ret)
}

# Load LD
LD <- readRDS("../LDmetrices/LDmetrices_FK40_PC1-40sparse_ALL.RDS")

# Reorder rows and columns of LD to match rot.pca and SNP.manifest
all(rownames(rot.pca)== SNP.manifest$pid)
#[1] TRUE
LD <- LD[SNP.manifest$pid, SNP.manifest$pid] # order by SNP.manifest pid
all(rownames(LD)== SNP.manifest$pid)
#[1] TRUE
all(colnames(LD)== SNP.manifest$pid)
#[1] TRUE


# Save all the stats useful for projection into one Rdata 

README <- c(rot.pca = "sparse rotation matrix (5519 * 40)",
            use.pca = "logical value matrix. hard thresholded rot.pca!=0 (5519 * 40)",
            beta.centers = "vector of column means of 40 Ferkingstad +1 control at sparse basis snps (5519)",
            shrinkage = "vector of w/sigma_maf at sparse basis snps (5519)",
            LD = "matrix of genotype correlations in 1000 Genomes EUR at sparse basis snps (5519 * 5519)",
            SNP.manifest = "data.table of sparse basis snps (5519 * 8)",
            stats = "info on the threshold quantile and number of SNPs per PC in the sparse basis (40 * 3)"
            project_sparse="function to project new beta + seb into sparse basis")


save(rot.pca,use.pca,beta.centers,shrinkage,LD,SNP.manifest,
     project_sparse,README,stats,
     file="../Rdata/cytokine_sparseSNP_40basis-LD-ProjFun.RData")



# We noticed a lot more LD==1 (15587)[5519*5519] in FK basis than in Cell basis (4826)[4472*4472]. what does it mean?
# LD gives us cor(X,Y)
# Sanity check!
# mannually check the snp_x, snp_y (where snp_x and snp_y have an LD of 1) in the basis rotation matrix
snp_x <- c("10:104628428", "10:104640883","10:122341927")
snp_y <- c("10:104628931", "10:104640933","10:122346943")


(rot.pca[snp_x[2],] - rot.pca[snp_y[2],])/rot.pca[snp_x[2],]

# Broadly SNPs with LD=1 share similar coefficients in our pca loading matrix. some disagreements occurred at the later 20 PCs