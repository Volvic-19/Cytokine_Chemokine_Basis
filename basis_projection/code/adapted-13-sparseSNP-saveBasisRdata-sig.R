# This script is adapted from 'adapted-09-saveBasisRdata.R', instead of using the complete e5 SNPs, this time we 
#    filter SNPs to find a minimum number of SNPs to keep while keeping a high correlation (>0.999) with the rotation matrix.
# The new Ferkingstad-40-sparse basis will be saved along with the new projection function (project_sparse) as
#    Rdata for further usage.

# 1. filter sparse SNP, build sparse 40 basis.
# 2. update projection function that computes confidence interval.

# 2022-07-12

library(data.table)
library(Matrix)
library(magrittr)

setDTthreads(10)

## Key parameters

nc = 40 # we start with keep the first 40 (out of 41 PCs)

##############################################
### Helper functions to search sparse SNP
##############################################

# Functions from https://grealesm.github.io/Bases/cell_basis_v2/basis_building/building_cell_basis_v2.html

## Formalise search
##' calculate use vector for given quantile and PC
##' This function will compare if each value in the rotation vector of a specific PC is bigger (more negative, taking negative absolute values) than the value at the proposed quantile.
##' @param q quantile
##' @param j PC
##' @param X PCA rotation matrix vector at PC j
##' @return use vector. A `logical` vector that states whether a SNP has a rotation value bigger in absolute terms than the value at the quantile. 
f0 <- function(q,j, X) { #j is PC, but not used in this function
  quants <- quantile(-abs(X),q)
  (-abs(X)) < quants
}

##' correlation with full data for given quantile and PC
##' This function will compute the use vector for a specified (j) PC and quantile (q), then calculate beta0, which is the product of multiplying the centred matrix (zm) for each `TRUE` rotation value in the use vector. Then outputs the correlation between beta0 and the projected basis matrix at the specified PC (j)
##' @param q quantile
##' @param j PC
##' @param zm Centred basis.matrix
##' @param pc.emp Basis object containing rotation matrix
##' @return correlation
f <- function(q, j, zm, pc.emp) {
  use <- f0(q, j, pc.emp$rotation[,j])
  beta0 = zm %*% ifelse(use, pc.emp$rotation[,j], 0)
  cor(beta0, pc.emp$x[,j])
}

##' Do search for given PC
##' This funtion takes a PC (j), the desired minimum correlation (mincor), the centred matrix (zm), and the basis object. Then it will perform a search to find the minimum quantile of SNPs that we can go (starting from 0.5 downwards) while keeping a correlation above the desired threshold (mincor). In other words, for a given PC, we'll try to find the minimum amount of SNPs that we can, with the biggest rotation values, while keeping the desired minimum correlation with the values in the transformed matrix (pc.emp$x).
##' @param j PC
##' @param mincor Minimum correlation threshold. 
##' @return quantile
g <- function(j, mincor, zm, pc.emp) {
  n <- nrow(pc.emp$rotation) # max number
  qi <- c(1/n,1) # test quantile
  # cat("Initial qi (Test quantile) is ", paste(qi, collapse = " "), ".\n", sep="")
  ni <- floor(qi * n) # test number
  # cat("Initial ni (Test number) is ", paste(ni, collapse = " "), ".\n", sep="")
  iteration = 1
  while(ni[2]-ni[1]>1) {
    ## check midpoint
    newq <- mean(qi)
    newc <- f(newq, j, zm, pc.emp)
    if(newc > mincor) {
      qi[2] <- newq
    } else {
      qi[1] <- newq
    }
    ni <- floor(qi * n) # test number
    # cat("After iteration ", iteration, " newq was ", paste(newq, collapse = " "),", correlation was ", paste(newc, collapse = " "), ", qi was updated to ", paste(qi, collapse = " "), " and ni to ", paste(ni, collapse = " "), ".\n",sep="")
    # iteration <- iteration +1
    ## cat(qi,"\t",ni,"\n")
  }
  return(qi[2])
}

##' Find quantiles for sparse basis
##' This function will take a covariance matrix and a PCA object created from it, as well as the number of retained components and the minimum desired correlation for quantile search, and perform a by-PC search for the minumim quantile of SNPs that we can keep while maintaining correlation above a desired threshold the values in the transformed matrix (pc.emp$x)
##' @param basis.mat The covariance matrix used for PCA.
##' @param basis The `princomp` object result of PCA.
##' @param nc Number of maximum PCs to keep.
##' @param mincor Minimum correlation threshold for quantile search. 
##' @return A list containing a summary matrix `summ` containing the minimum quantile and the number of SNPs in it for each PC, and a `user` logical matrix of SNPs in each PC in the sparse matrix.
find_sparse_q <- function(basis.mat.emp, pc.emp, nc=NULL, mincor=0.999) {
  if(is.null(nc))
    nc <- ncol(pc.emp$rotation)-1
  zm <- scale(basis.mat.emp,center=TRUE,scale=FALSE) # centered input
  q <- sapply(1:nc, g, mincor=mincor, zm=zm, pc.emp=pc.emp)
  summ <- cbind(PC=1:nc, q=q, n=floor(q * nrow(pc.emp$rotation)))
  user <- mapply(f0, q=q, j=1:nc, X=lapply(1:nc, function(j) pc.emp$rotation[,j]))
  print(summ)
  invisible(list(summ=summ, user=user))
}



################################################
### Load basis and Filter objects by sparse SNPs
################################################

basis <- readRDS("../../basis_building/PCA/cytokine_basis_40.RDS") # # 641,079 unique SNPs
basis.mat <- readRDS("../../basis_building/PCA/cytokine_basis_matrix_40.RDS")
shrinkage.DT <- fread("../../basis_building/manifest/shrinkage.e5.DT.tsv.gz", tmpdir = "tmp")
SNP.manifest <- fread("../../basis_building/manifest/CCL8_consensus_manifest_6M_e5.tsv")


ret <- find_sparse_q(basis.mat, basis, nc)

# We keep these stats to save them together with the rest of the basis files
stats <- ret$summ 

## We'll keep all SNPs with at least one TRUE in a PC
sparseSNPs <- apply(ret$user,1,any) # named vector of logical value
table(sparseSNPs)
# sparseSNPs
# FALSE   TRUE 
# 635560   5519

# Check the frequency of sparse SNPs
hist(rowSums(ret$use)[sparseSNPs],main=paste("SNP use in sparse PCA"))

# filter the “use” matrix, the rotation matrix, and the centred matrix to contain only our sparse SNPs. 
# In the rotation matrix, we’ll keep only values for SNPs that have a significant effect, and thus are 
# listed in each PC, turning the rest into zeros.

# We reduce the use matrix to contain only our sparse SNPs
use.pca <- ret$use[sparseSNPs, 1:nc]
# Same with the rotation matrix. Plus, we assign the value of each sparse SNP not in the specific PC zero.
rot.pca <- basis$rotation[sparseSNPs, 1:nc]  # new sparse rotation matrix
rot.pca[!use.pca]=0

# Filter centred covariance matrix : We remove every not sparse SNPs from the centred covariance matrix
M.centred <- scale(basis.mat,center=TRUE,scale=FALSE) # centered input
M.centre <- attr(M.centred,"scaled:center") # centre

Z <- M.centred[1:nc,which(sparseSNPs)] #num [1:40, 1:5519] , 40 PCs, 5519 pids, 40 * 5519 metrices


# We compute the vector of column means of 40 traits + one control. We modify beta.centers as we only need colmeans of sparse SNP columns
beta.centers <- structure(attr(M.centred,"scaled:center"),names=colnames(M.centred))[rownames(rot.pca)]# named vector length 5519. dnsized from M.centred-"scaled:center" 641079/

# Filter the shrinkage metrics :
ishrink <- shrinkage.DT[match(rownames(rot.pca),pid),c("pid", "shrinkage")]
shrinkage <- structure(ishrink$shrinkage,names=ishrink$pid) # named vector of length 5519. name: pid

# Filter manifest : remove non-sparse SNPs
SNP.manifest.cell.sparse <- SNP.manifest[match(rownames(rot.pca), pid), ] # 641079 dn to 5519
SNP.manifest <- SNP.manifest.cell.sparse

# Filter shrinkage e5 table: remove non-sparse SNPs
shrinkage.DT <- shrinkage.DT[pid %in% rownames(rot.pca)] # 25643160 rows dn to 220760 rows

# Save filtered shrinkage table( basis matrix to project as reference later), updated manifest
fwrite(shrinkage.DT, file = "../../basis_building/manifest/shrinkage.sparse.e5.tsv",sep = "\t") # 5519 * 40 = 220760
fwrite(SNP.manifest.cell.sparse, file = "../../basis_building/manifest/CCL8_consensus_manifest_6M_e5_sparse.tsv", sep = "\t") #5519

# Save filtered rotation matrix, filtered beta.centers, filtered sparse manifest, 
README <- c(rot.pca = "sparse rotation matrix (5519 * 40)",
            use.pca = "logical value matrix. hard thresholded rot.pca!=0 (5519 * 40)",
            beta.centers = "vector of column means of 40 Ferkingstad +1 control at sparse basis snps",
            shrinkage = "vector of w/sigma_maf at sparse basis snps",
            #LD = "matrix of genotype correlations in 1000 Genomes EUR at sparse basis snps",
            SNP.manifest = "data.table of sparse basis snps",
            stats = "info on the threshold quantile and number of SNPs per PC in the sparse basis (40 * 3)")
            #project_sparse="function to project new beta + seb into sparse basis"

save(rot.pca,use.pca,beta.centers,shrinkage,SNP.manifest,
     README,stats,#project_sparse
     file="../Rdata/cytokine_sparseSNP_40basis-noLD-noProjFun.RData")







