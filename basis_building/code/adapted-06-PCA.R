## This script is to be applied after shrinkage. We first input shrinkage metrics, create matrix of shrinked values,
## and use PCA to create the basis


library(data.table)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(cowplot)

setDTthreads(10)

# 21-06-2022




shrinkage.DT <- fread("../manifest/shrinkage.e5.DT.tsv.gz")

# create basis using PCA

B <- data.table::dcast(shrinkage.DT,pid ~ Trait, value.var="metric")
snames <- B[,1]$pid
tmp.mat <- t(as.matrix(B[,-1]))
colnames(tmp.mat) <- snames
basis.mat <- rbind(tmp.mat, control = rep(0, ncol(tmp.mat)))

basis <- stats::prcomp(basis.mat,center = TRUE, scale = FALSE)


# Visualising the basis(draw PCA biplot)

plot.DT <- data.table(trait=rownames(basis$x),basis$x)
biplot <- ggplot(plot.DT, aes(PC1,PC2, label = trait))+
  geom_point()+
  geom_label_repel(size=3,force = 20, seed = 1) + 
  ggtitle("Cytokine_chemokine basis biplot: e5")+
  theme_bw()
# xlim <- with(plot.DT,c(min(PC1),max(PC1)) * 1.2)
# ylim <- with(plot.DT,c(min(PC2),max(PC2)) * 1.2)
# with(plot.DT,plot(PC1,PC2,type='n',xlim=xlim,ylim=ylim,main="Blood cell basis biplot"))
# with(plot.DT,text(PC1,PC2,labels=trait, cex= 0.8, adj=c(0,0)))
# with(plot.DT,points(PC1,PC2,cex=0.5,pch=19))
biplot
ggsave("../plots/e5_PC1_PC2.png", biplot)


# Visualising the basis(hierarchical clustering)

#cluster <- basis$x %>% dist %>% hclust
#plot(cluster)



# Save files

saveRDS(basis, "../PCA/cytokine_basis.RDS")
saveRDS(basis.mat, "../PCA/cytokine_basis_matrix.RDS")


# Compute a sparse basis (we have  673,318_ SNPs at the moment, but can further filter to get _ sparse SNPs that recreate the projected matrix )

M.centred <- scale(basis.mat,center=TRUE,scale=FALSE) # centered input
M.centre <- attr(M.centred,"scaled:center") # centre
    M.projected <- M.centred %*% basis$rotation # projection
cor(as.vector(M.projected),as.vector(basis$x)) # check that our covariance matrix's correlation with basis$x is 1.
#[1] 1
sum(abs(as.vector(M.projected)-as.vector(basis$x))) # and the substraction of both vectors is 0
#[1] 0


# Plotting the reconstruction error
# The following function will help us identify the number of PCs that we’ll need to take into account in our basis, by plotting the Mean Square Error for each PC, comparing our centred matrix with a reconstructed one. This example is displayed both linear and log scale.

reconstruction.error.plot <- function(basis, basis.mat, scale="linear"){
  M.centred <- scale(basis.mat,center=TRUE,scale=FALSE) # centered input
  M.centre <- attr(M.centred,"scaled:center") # centre
  M.projected <- M.centred %*% basis$rotation # projection # each column of basis$rotation is an eigenvector
  p <- ncol(M.centred) # 265701
  n <- nrow(M.centred) # 36
  u <- ncol(M.projected) # 36
  if (n != u) stop("Something's wrong, nrow(M.centred) and ncol(M.projected) are not the same!")
  err <- sapply(1:u, function(k) {
    reconstructed <- matrix(M.centre,n,p) +
      M.projected[,1:k,drop=FALSE] %*% t(basis$rotation[,1:k,drop=FALSE])
    residuals = M.centred - reconstructed
    error = mean(residuals^2)
  })
  if(scale=="log"){
    plot(1:u,log(err),ylim=c(min(log(err)),max(log(err))),ylab="Log Mean square error",xlab="Number of components"); abline(h=min(log(err)),lty=2)
  } else {
    plot(1:u,err,ylim=c(0,max(err)),ylab="Mean square error",xlab="Number of components"); abline(h=min(err),lty=2)
  }
}



rec.plot <- reconstruction.error.plot(basis = basis, basis.mat = basis.mat)

rec.plot.log <- reconstruction.error.plot(basis = basis, basis.mat = basis.mat, scale = "log")

# Is our basis naturally sparse?
# Plot rotation vectors by PC to see what their distribution looks like. We expect most elements in the rotation vectors to be zero, or close to zero, meaning that only a few SNPs will have an effect in the rotation of basis.mat during PCA, and thus, we can make our basis sparse by shedding off most SNPs.

sparse.plot <- function(basis, pcs){
  sX <- basis$rotation[,1:pcs]
  m14 <- reshape2::melt(sX)  %>% as.data.table()
  theme_set(theme_cowplot())
  ggplot(m14,aes(x=value)) + 
    geom_histogram(binwidth=0.01) +
    facet_wrap(~Var2) +
    scale_y_sqrt("Count (sqrt scale)",
                 breaks=c(1,1e+3,1e+4,1e+5,#2e+5,
                          3e+5),labels=c("0","1,000","10,000","100,000",#"200,000",
                                         "300,000"),
                 limits=c(0,3e+5)) +
    background_grid()
}

# sparse plot
# nc =2
spplot <- sparse.plot(basis, nc)
spplot

# Find the SNPs for our sparse basis
# The aim now is to find a minimum number of SNPs to keep while keeping a high correlation (>0.999) with the rotation matrix. To do that we’ll use the functions below to perform a quantile search. In a nutshell, what these functions will do is to take each PC in the rotation matrix, quantilize its values in terms of how far from zero they are, and then search for the lowest quantile (ie. the minimum amount of SNPs starting from the farthest from zero) that keeps a correlation with the full data above 0.999. These SNPs will make up the SNPs in the sparse matrix

## Formalise search
##' calculate use vector for given quantile and PC
##' This function will compare if each value in the rotation vector of a specific PC is bigger (more negative, taking negative absolute values) than the value at the proposed quantile.
##' @param q quantile
##' @param j PC
##' @param X PCA rotation matrix vector at PC j
##' @return use vector. A `logical` vector that states whether a SNP has a rotation value bigger in absolute terms than the value at the quantile. 

f0 <- function(q,j, X) {
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





