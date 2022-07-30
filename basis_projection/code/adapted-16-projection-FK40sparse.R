# This script is to apply projection of ahola.olli 1-261 datasets (reduced by adapted-15-reduction.R) and basis datasets onto
# our Ferkingstad basis (sparse, PC1-40) for a complete projection.  Significance is computed here.

# 2022-07-14
# 2022-07-19: rerun projection , as we updated project_sparse function;
# 2022-07-19: after projection, change the write.table of QC table, from ("1e-05" to "actual number"), as despite
#              convenience in reading, we need accurate overall.p for FDR computation.

library(data.table)
library(Matrix)
library(magrittr)

setDTthreads(10)

# 5519 unique SNPs


##############################################
###  projection of 1-261
##############################################

load("../Rdata/cytokine_sparseSNP_40basis-LD-ProjFun.RData")

# Create log file
mpath <- "~/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_projection/code/log/"
date <- format(Sys.time(), format="%Y%m%d")
logname <- paste0(mpath,"log_project_",date,".txt")
file.create(logname)

# Files to project
files_to_project  <- dir("../data_40_sparse", pattern = ".tsv") #261 datasets to project
nfiles <- length(files_to_project)
Trait <- rep(NA,nfiles)
nSNP <- rep(NA, nfiles)
overall_p <- rep(NA, nfiles)
mscomp <- rep(NA,nfiles)


###############################################
# Projection of 1-261 external traits		#	
##############################################

projected.table <- lapply(files_to_project, function(file){
	message("Projecting ", file)
	trait_label <- strsplit(file, "-", fixed = TRUE)[[1]][1] #trait name
	index <- which(files_to_project == file) 
	sm <- fread(paste0("../data_40_sparse/",file)) #data_40_sparse

	# Some checks
	sm <- unique(sm)
	sm[sm == ""] <- NA # Some missing data might pass as empty string. This will fix that	
	sm <- na.omit(sm, cols = c("pid", "BETA", "SE", "P")) # 
	dups <- sm$pid[duplicated(sm$pid)]

	if(length(dups) > 0){
		dupmessage= "This file has duplicated pids. I removed them prior to projection. You might want to check it."
		message(dupmessage)
		# Write to log
		write(paste0(file,". ",dupmessage), logname, append=TRUE)
		sm <- sm[!pid %in% dups] # Remove all duplicated instances, to be safe
	} # sm is 5450, less than 5519

	# A bit of QC 
	Trait[index] <<- trait_label
	nSNP[index] <<- nrow(sm)

  	projected.userdata <- tryCatch(expr = project_sparse(beta=sm$BETA, seb=sm$SE, pids=sm$pid)[,trait:=trait_label][], #project
  	                      error =function(e) {
				      failmessage <- "Projection for this file had non-zero exit status, please check. Jumping to next file..."
				      message(failmessage)
				      write(paste0(file,". ",failmessage), logname, append=TRUE)
				      return(NULL)})
  	if(is.null(projected.userdata)) {
		return(NULL)
	} else{
	projected.userdata[, proj:=NULL] # 105 * 4
  	setnames(projected.userdata, c("var.proj", "delta", "p", "trait"), c("Var.Delta", "Delta", "P", "Trait")) # Changed var.proj to var.delta
  	setcolorder(projected.userdata, c("PC", "Var.Delta", "Delta", "p.overall", "z", "P", "Trait"))

  	# More QC
  	overall_p[index] <<- projected.userdata$p.overall[1] # sprintf("%1.0e",projected.userdata$p.overall[1])
  	minc.idx <- which.min(projected.userdata$P)
  	mscomp[index] <<- sprintf("%s (%1.0e)",projected.userdata$PC[minc.idx],projected.userdata$P[minc.idx])
	}
	projected.userdata
}
) #current

projected.table[sapply(projected.table, is.null)]  <- NULL 
projected.table <- rbindlist(projected.table)
#projected.table  <- projected.table[,.(PC, Delta,Trait)]
projected.table  <- projected.table[,.(PC, Delta, Var.Delta, z, P, Trait)]  #p.overall is in QC.table

# add ratio of covered sparse SNP to QC table as a threshold
QC.table <- data.table(Trait, nSNP, overall_p, mscomp)
QC.table[, ratio_coveredSNP := nSNP / nrow(rot.pca)]


projtablename  <- paste0("Projection_cytokine_basis_sig_40sparse_", date, ".tsv") # add sig_40sparse_
qctablename  <- paste0("QC_cytokine_basis_sig_40sparse_", date, ".tsv") # add sig_40sparse_


#while(projtablename %in% dir(paste0(mpath, "cell_basis_v3_varimax/Projections"))){
#  version  <- version + 1
#  projtablename  <- paste0("Projection_cell_basis_v3_", date, "-v",version, ".tsv")
#  qctablename  <- paste0("QC_cell_basis_v3_", date, "-v",version, ".tsv")
#}

message("start writing projected.table & QC.table")

write.table(projected.table, paste0("../Projections/", projtablename), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(QC.table, paste0("../Projections/", qctablename), sep = "\t", quote = FALSE, row.names = FALSE)

message("Finished writing!")



###############################################
# Projection of basis to basis 				  #	
##############################################


# Files to project
e5  <- fread("../../basis_building/manifest/shrinkage.sparse.e5.tsv", tmpdir = "tmp") 
lists_to_project <- split(e5, e5$Trait) # A list of length 40, each element is a trait datatable to project

nfiles <- length(lists_to_project) # 40

# QC table
Trait <- unique(e5$Trait)
nSNP <- rep(NA, nfiles)
names(nSNP) <- Trait

overall_p <- rep(NA, nfiles) #current, need to figure out how to 
names(overall_p) <- Trait
mscomp <- rep(NA,nfiles)
names(mscomp) <- Trait

# Projection
projected.table <- lapply(lists_to_project, function(file){

	current.trait <- file$Trait[1]
	message("Projecting ", current.trait)

	# Some checks
	file <- unique(file)
	file[file == ""] <- NA # Some missing data might pass as empty string. This will fix that	
	file <- na.omit(file, cols = c("pid", "BETA", "SE", "P"))
	
	dups <- file$pid[duplicated(file$pid)]
	if(length(dups) > 0){
		stop("This file has duplicated pids. Which means there is duplicated pid in shrinkage table!")
	}

	# A bit of QC 
	nSNP[current.trait] <<- nrow(file)

  	projected.userdata <- tryCatch(expr = project_sparse(beta=file$BETA, seb=file$SE, pids=file$pid)[,trait:=current.trait][],
  	                      error =function(e) {
				      failmessage <- "Projection for this file had non-zero exit status, please check. Jumping to next file..."
				      message(failmessage)
				      write(paste0(current.trait,". ",failmessage), logname, append=TRUE)
				      return(NULL)})
  	if(is.null(projected.userdata)) {
		return(NULL)
	} else{
	projected.userdata[, proj:=NULL] # 105 * 4
  	setnames(projected.userdata, c("var.proj", "delta", "p", "trait"), c("Var.Delta", "Delta", "P", "Trait")) # Changed var.proj to var.delta
  	setcolorder(projected.userdata, c("PC", "Var.Delta", "Delta", "p.overall", "z", "P", "Trait"))

  	# More QC
  	overall_p[current.trait] <<- projected.userdata$p.overall[1] # sprintf("%1.0e",projected.userdata$p.overall[1])
  	minc.idx <- which.min(projected.userdata$P)
  	mscomp[current.trait] <<- sprintf("%s (%1.0e)",projected.userdata$PC[minc.idx],projected.userdata$P[minc.idx])

	}
	
	projected.userdata
}) #current

projected.table[sapply(projected.table, is.null)]  <- NULL 
projected.table <- rbindlist(projected.table)
projected.table  <- projected.table[,.(PC, Delta, Var.Delta, z, P, Trait)]  #p.overall is in QC.table

QC.table <- data.table(Trait = names(nSNP), nSNP) 
QC.table[, overall_p := overall_p[Trait]][, mscomp:=mscomp[Trait]]
QC.table[, ratio_coveredSNP := nSNP / nrow(rot.pca)]

# Write tables
projtablename  <- paste0("Projection_cytokine_basis_sig_40sparse_basisTObasis", date, ".tsv") #_40
qctablename  <- paste0("QC_cytokine_basis_sig_40sparse_basisTObasis", date, ".tsv") #_40


message("start writing projected.table & QC.table")

write.table(projected.table, paste0("../Projections/", projtablename), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(QC.table, paste0("../Projections/", qctablename), sep = "\t", quote = FALSE, row.names = FALSE)

message("Finished writing!")
