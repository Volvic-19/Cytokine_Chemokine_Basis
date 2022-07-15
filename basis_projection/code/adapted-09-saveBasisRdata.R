# This script and the following adapted-10 is to apply projection of Ohola.olli datasets (reduced by adapted-08-reduction.R) onto
# our Ferkingstad-40 basis (complete, not sparse yet) for a quick visualisation test. Signficance is not computed. 

# Before running projection, we need to prepare Rdata to avoid repeated computing of basis, 

# 2022-07-11

##############################################
### LOAD LIBRARIES AND WRITE Rdata
##############################################

library(data.table)
library(Matrix)
library(magrittr)

setDTthreads(10)

# 641,079 unique SNPs

basis <- readRDS("../../basis_building/PCA/cytokine_basis_40.RDS")
basis.mat <- readRDS("../../basis_building/PCA/cytokine_basis_matrix_40.RDS")
shrinkage.DT <- fread("../../basis_building/manifest/shrinkage.e5.DT.tsv.gz", tmpdir = "tmp")
SNP.manifest <- fread("../../basis_building/manifest/CCL8_consensus_manifest_6M_e5.tsv")

M.centred <- scale(basis.mat,center=TRUE,scale=FALSE) # centered input
#M.centre <- attr(M.centred,"scaled:center") 
#M.projected <- M.centred %*% basis$rotation # projection # 105*105

# rot.pca: As we did not remove any sparse SNPs at this stage, it will be just complete rotation matrix
rot.pca <- basis$rotation

# beta.centers: named vector, length of complete unique SNPs
beta.centers <- structure(attr(M.centred,"scaled:center"),names=colnames(M.centred))[rownames(rot.pca)]

# shrinkage: named vector, length of complete unique SNPs
ishrink <- shrinkage.DT[match(rownames(rot.pca),pid),c("pid", "shrinkage")] #641079 * 2
shrinkage <- structure(ishrink$shrinkage,names=ishrink$pid)

# Adapted from projection_sparse function. It does not compute significance, only coordinates.
project_nosig <- function (beta, seb, pids) { #seb is SE
    if (length(beta) != length(seb) || length(beta) != length(pids) || 
        !length(beta)) 
        stop("arguments must be equal length vectors > 0")
    if (!all(pids %in% SNP.manifest$pid)) 
        stop("all pids must be members of basis (SNP.manifest$pid)")
    if (length(pids) < 0.95 * nrow(rot.pca)) #original 0.95 threshold. But as we are not using sparse SNP, harder to meet 0.95 request
        warning("more than 5% sparse basis snps missing")
    b <- beta * shrinkage[pids] - beta.centers[pids]
    proj <- b %*% rot.pca[pids, ]

   # v <- seb * shrinkage[pids] * rot.pca[pids, ]
   # var.proj <- t(v) %*% LD[pids, pids] %*% v
    ctl <- (-beta.centers[pids]) %*% rot.pca[pids, ]
    delta <- (proj - ctl)[1, ]
   # chi2 <- (t(delta) %*% solve(var.proj) %*% delta)[1, 1]
    ret <- data.table::data.table(PC = colnames(proj), proj = proj[1, 
        ],  delta = delta) 
        # NOTE: remember to update df as number of PCs!
   # ret$z = ret$delta/sqrt(ret$var.proj)
   # ret$p = stats::pnorm(abs(ret$z), lower.tail = FALSE) * 2
    copy(ret)
}


README <- c(rot.pca = "rotation matrix (not sparsed yet, complete SNPs)",
			beta.centers = "named vector of column means of 40 + 1 at complete SNPs",
			shrinkage = "named vector of w/sigma_maf at complete SNPs",
			SNP.manifest = "SNP.manifest, 6M e5",
			project_nosig = "adapted projection function, no sig was computed"
	) # check if in the original script, shrinkage was recomputed ?


save(rot.pca, beta.centers, shrinkage, SNP.manifest, project_nosig, README, file = "../Rdata/cytokine_completeSNP_40.Rdata")

## preview the number of non-zero SNPs in PC 1-20
#summary(rot.pca)
#colSums(!rot.pca == 0)

#min(abs(rot.pca[,"PC1"]))
#[1] 1.692257e-11

