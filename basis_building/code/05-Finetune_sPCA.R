### Fine-tuning sPCA parameters.

## Background: We want to use sPCA for the new basis. However, sPCA requires setting sparsity parameters. Since we don't know them, we'll project the same traits as those in the basis but from European Chen (aka. Chen_2) and Astle using different parameters. In summary:

# 1. select some sparsity params
# 2. run spca with those params
# 3. while at it, save the driver SNPs per PC
# 3. project Chen-EUR, Astle and calculate the squared difference
# 4. store params + mean difference per PC (mean taken over the cell types)
# 5. adjust the params (slightly up / slightly down)
# 6. goto 1
# eventually you have a vector of params and mean sq diffs for each PC, and can select the param that minimises the difference

# Guillermo Reales
# 2022-06-08

#############################################
######## Load libraries     #################
#############################################

library(data.table)
setDTthreads(0)
library(magrittr)
library(elasticnet)

#library(ggplot2)
#library(cowplot)
#library(RColorBrewer)

#############################################
######## Helper functions   #################
#############################################

prepare.datasets <- function(paths, snames){
	message("Preparing datasets for projection.")
	ds <- lapply(paths, function(x){
		message("Working on ", x)
		tname <- paste(sapply(strsplit(x, split="/|-|_"), `[`,3:4), collapse="_")
		d <- fread(x, tmpdir="tmp")
		d[, pid:=paste(CHR38,BP38, sep=":")]
		d <- d[pid %in% snames]
		d <- d[, .(pid, BETA, SE)]
		d[, Trait:=tname]
		d
	})
	ds <- rbindlist(ds)
 }

project.ds <- function(dspr, rot.spca, shrinkage, beta.centers){
	
   	ts <- unique(dspr$Trait)
	proj.table <- lapply(ts, function(i){
	
	    tt <- dspr[Trait == i]
	    beta <- tt$BETA
	    pids <- tt$pid

	    if ( length(beta) != length(pids) ) 
	         stop("arguments must be equal length vectors > 0")
	    if (length(pids) < 0.95 * nrow(rot.spca)) 
	         warning("more than 5% sparse basis snps missing")
	    
	    b <- beta * shrinkage[pids] - beta.centers[pids]
	    proj <- b %*% rot.spca[pids, ]
	    ctl <- (-beta.centers[pids]) %*% rot.spca[pids, ]
            delta <- (proj - ctl)[1, ]
	    pdt <- data.table(PC = colnames(proj), Projection=as.numeric(proj), Delta=as.numeric(delta), Trait=i)
       })
       proj.table <- rbindlist(proj.table)
}


create.basis.project <- function(prds, param, basis.mat, shrinkage.DT, destfile, snps.destfile){
	
	message("Creating bases using sPCA with provided parameters and projecting datasets.")
	message("Computing centred input and beta.centers.")
	# We compute the vector of column means of 15 blood cell trait GWAS + one control
	M.centred <- scale(basis.mat,center=TRUE,scale=FALSE) # centered input
	beta.centers <- structure(attr(M.centred,"scaled:center"),names=colnames(M.centred))	
	# We do the same with the shrinkage metrics
	ishrink <- shrinkage.DT[match(colnames(basis.mat),pid),c("pid", "shrinkage")]
	shrinkage <- structure(ishrink$shrinkage,names=ishrink$pid)

	p.tables <- lapply(param, function(i){
		message("Creating basis using param: ", i,".")
		basis.spca <- spca(basis.mat, K = 8, para = rep(i, 8), trace=TRUE)
		message("Done!")
		rot.spca <- basis.spca$loadings
		var_snps <- data.table(PC = paste0("PC",1:8), Var.exp=basis.spca$pev, SNP.PC= colSums(rot.spca !=0), SNP.total= sum(rowSums(rot.spca !=0) != 0))
		
		message("Extracting driver SNPs for basis param ", i, ".")
		snp <- rbindlist(lapply(1:8, function(x){
			data.table(pid = rownames(rot.spca[rot.spca[, x] != 0,]), PC = paste0("PC", x), Param= i)
		}))
		if(file.exists(snps.destfile)){
			fwrite(snp, snps.destfile , sep="\t", append=TRUE)
		} else {
			fwrite(snp, snps.destfile , sep="\t")
		}

		message("Projecting datasets onto the basis...")
		proj.table <- project.ds(prds, rot.spca = rot.spca, shrinkage = shrinkage, beta.centers = beta.centers) # Project datasets onto the basis
		proj.table[, Trait2:=sapply(strsplit(Trait, "_"), `[`, 1)][, Param:=i]
		proj.table <- merge(proj.table, var_snps, by="PC")
		if(file.exists(destfile)){
			fwrite(proj.table, destfile , sep="\t", append=TRUE)
		} else {
			fwrite(proj.table, destfile , sep="\t")
		}

       })
	p.tables <- rbindlist(p.tables)
	p.tables
}

projection.stats <- function(p.table){

	message("Creating summary statistics from projection tables...")
	var_snps <- unique(p.table[,.(PC, Param, Var.exp, SNP.PC, SNP.total)])
	d1 <- p.table[, .(sqdiff = diff(Delta)^2), by=c("Trait2", "PC", "Param")][, .(meansqdiff.Delta = mean(sqdiff)), by=c("PC", "Param")]
	d2 <- p.table[, .(mean.Delta = mean(abs(Delta)), sd.Delta = sd(abs(Delta))), by=c("PC","Param")]
	diff.table <- merge(d1, d2, by=c("PC", "Param"))
	diff.table <- merge(diff.table, var_snps, by=c("PC", "Param"))

}

# create.basis.extract.drivers <- function(prds, param, basis.mat, shrinkage.DT, destfile){
# 
# 	message("Creating bases using sPCA with provided parameters and projecting datasets.")
# 	message("Computing centred input and beta.centers.")
# 	# We compute the vector of column means of 15 blood cell trait GWAS + one control
# 	M.centred <- scale(basis.mat,center=TRUE,scale=FALSE) # centered input
# 	beta.centers <- structure(attr(M.centred,"scaled:center"),names=colnames(M.centred))	
# 	# We do the same with the shrinkage metrics
# 	ishrink <- shrinkage.DT[match(colnames(basis.mat),pid),c("pid", "shrinkage")]
# 	shrinkage <- structure(ishrink$shrinkage,names=ishrink$pid)
# 
# 	p.tables <- lapply(param, function(i){
# 		message("Creating basis using param: ", i,".")
# 		basis.spca <- spca(basis.mat, K = 15, para = rep(i, 15), trace=TRUE)
# 		message("Done!")
# 		rot.spca <- basis.spca$loadings
# 		message("Extracting driver SNPs...")
# 		snp <- rbindlist(lapply(1:15, function(x){
# 			data.table(pid = rownames(rot.spca[rot.spca[, x] != 0,]), PC = paste0("PC", x), Param= i)
# 		}))
# 		if(file.exists(destfile)){
# 			fwrite(snp, destfile , sep="\t", append=TRUE)
# 		} else {
# 			fwrite(snp, destfile , sep="\t")
# 		}
# 
#        })
# 	p.tables <- rbindlist(p.tables)
# 	p.tables
# }


#############################################
######## Load files         #################
#############################################

 # We'll run it in an array, so take the argument from command line
 # args <- commandArgs(trailingOnly = TRUE)
 # args <- as.numeric(args)

# Prepare necessary files, or simply load them if already there.
 # Prepare matrix 
 shrinkage.DT <- fread("../data/shrinkage.DT.1e5.tsv.gz", tmpdir="tmp") # Prepare basis matrix containing all the SNPs. Don't comment this line, since we'll need it anyway.
 #  B <- dcast(shrinkage.DT,pid ~ Trait, value.var="metric")
 #  snames <- B[,1]$pid
 #  tmp.mat <- t(as.matrix(B[,-1]))
 #  colnames(tmp.mat) <- snames
 #  basis.mat <- rbind(tmp.mat, control = rep(0, ncol(tmp.mat)))
 #  saveRDS(basis.mat, "../data/basis.mat.1e6.RDS") # Save it, just in case

 # Prepare files to project
 #  seltraits <- c("BASC", "EOSC", "LYMC", "MONC", "PLAV", "NEUC", "PLAC", "LEUC")
 #  pchen1 <- paste0(seltraits, "_Chen_32888493_2") %>% paste(., collapse="|")
 #  pastle <- paste0(seltraits, "_Astle") %>% paste(., collapse="|")
 #  pboth  <- paste(pchen1, pastle, sep="|")
 #  ppath  <- "../data/"
 #  mfiles <- dir(ppath, pattern=pboth)
 #  ftopr <- paste0("../data/", mfiles)
 # ## ftopr <- paste0("../data/", dir("../data/", pattern="*-ft.tsv.gz"))
 #  prds <- prepare.datasets(ftopr, snames = snames)
 # ## prds$Trait <- gsub("ERYC", "RBC", prds$Trait)
 #  fwrite(prds, "../data/prds_1e6.tsv.gz", sep ="\t")


 # Load pre-processed files, if already available
  prds <- fread("../data/prds_1e5.tsv.gz") # Just the aligned Eur Chen and Astle together in one file. 
  basis.mat <- readRDS("../data/basis.mat.1e5.RDS")


# Now we need to 
# 1. Select a number of parameters for sPCA.
# 2. Run sPCA using each of those parameters for all PCs.
# 3. Project European Chen and Astle onto each basis.
# 4. Compute the square differences between Eur_Chen and Astle for each trait, and take the mean per component.
# 5. Find the parameters for which the average mean square difference is the minimum for each PC.
#
# Select parameters

spca.param = 10^seq(-9, -5, by=0.5)
#spca.param = spca.param[args]
#param.name = formatC(spca.param, format="e", digits=2)

# The function finetune.parameters() will take care of steps 2-4
create.basis.project(prds, param=spca.param, basis.mat=basis.mat, shrinkage.DT = shrinkage.DT, destfile = "../data/projtable_8t_1e5s_9to5.tsv", snps.destfile = "../data/drivers_8t_1e5s_9to5.tsv")

# We need to also extract the driver SNPs, so let's do that too
#spca.param = 10^seq(-9, -5, by=0.5)
# spca.param = 10^seq(-9, -5, by=0.5)[c(1,4:5,9)] # Let's try with params 1e-9, 3e-8, 1e-7, 1e-5


# p.table <- fread("../data/projtable_1e5_9to5.tsv")
# d.table <- projection.stats(p.table)
# 
# #fwrite(p.table, "../data/projtable_1e5_9to5.tsv", sep="\t")
# #fwrite(d.table, "../data/difftable_1e5_9to5.tsv", sep="\t")
# #fwrite(p.table, paste0("../data/projtable_1e5_9to5_", param.name, ".tsv"), sep="\t")
# #fwrite(d.table, paste0("../data/difftable_1e5_9to5_", param.name, ".tsv"), sep="\t")
# 
