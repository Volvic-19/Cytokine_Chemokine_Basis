## This script is to be applied after shrinkage. We first input shrinkage metrics, create matrix of shrinked values,
## and then use PCA to create the basis

## 2022-07-11

library(data.table)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(cowplot)

setDTthreads(10)




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
biplot <- ggplot(plot.DT, aes(PC5,PC6, label = trait))+
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
ggsave("../plots/basis40/basis40_e5_PC5_PC6.png", biplot)



# Save files

saveRDS(basis, "../PCA/cytokine_basis_40.RDS")
saveRDS(basis.mat, "../PCA/cytokine_basis_matrix_40.RDS")

# Read files
basis <- readRDS("../PCA/cytokine_basis_40.RDS")
basis.mat <- readRDS("../PCA/cytokine_basis_matrix_40.RDS")


# Inspect PCs

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
  p <- ncol(M.centred)  # 641079
  n <- nrow(M.centred)  # 41
  u <- ncol(M.projected)   # 41
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


# Plot MSE reconstruction error
png(file = "../plots/basis40/basis40_MSE.png")
reconstruction.error.plot(basis = basis, basis.mat = basis.mat)
dev.off()

png(file = "../plots/basis40/basis40_MSE_log.png")
reconstruction.error.plot(basis = basis, basis.mat = basis.mat, scale = "log")
dev.off()


# Plot scree plot (eigenvalue)
# compute total variance

plot_scree <- function(basis){

	variance = basis$sdev^2 / sum(basis$sdev^2)

	plot <- ggplot(data=data.frame(x=1:length(basis$sdev),y=variance), aes(x=x,y=y))+
	geom_line()+
	xlab("Principle Component")+
	ylab("Variance Explained")+
	ggtitle("Scree Plot")

	return(plot)
}

scree.plot <- plot_scree(basis)
ggsave("../plots/basis40/basis40_scree_plot_e5.png", scree.plot)

