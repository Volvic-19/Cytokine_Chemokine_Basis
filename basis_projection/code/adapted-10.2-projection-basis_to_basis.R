# This script is to apply projection of Ferkingstad datasets  onto
# our Ferkingstad basis (complete, not sparse yet) for a quick visualisation test. Signficance is not computed. 

# To save the pain of reducing Ferkingstad datasets, we use BETA, SE, Pid from shrinkage table to do projection.

# Before running projection, we need to load Rdata to avoid repeated computing of : basis, 
# 2022-07-11

library(data.table)
library(Matrix)
library(magrittr)

setDTthreads(10)

# 641,079 unique SNPs


##############################################
### LOAD Rdata and run projection
##############################################

load("../Rdata/cytokine_completeSNP_40.Rdata") #_40

# Create log file
mpath <- "~/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_projection/code/log/"
date <- format(Sys.time(), format="%Y%m%d")
logname <- paste0(mpath,"log_project_",date,".txt")
file.create(logname)

# Files to project
e5  <- fread("../../basis_building/manifest/shrinkage.e5.DT.tsv.gz", tmpdir = "tmp") 
lists_to_project <- split(e5, e5$Trait) # A list of length 40, each element is a trait datatable to project

nfiles <- length(lists_to_project)

# QC table
Trait <- unique(e5$Trait)
nSNP <- rep(NA, nfiles)
names(nSNP) <- Trait


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

  	projected.userdata <- tryCatch(expr = project_nosig(beta=file$BETA, seb=file$SE, pid=file$pid)[,trait:=current.trait][],
  	                      error =function(e) {
				      failmessage <- "Projection for this file had non-zero exit status, please check. Jumping to next file..."
				      message(failmessage)
				      write(paste0(current.trait,". ",failmessage), logname, append=TRUE)
				      return(NULL)})
  	if(is.null(projected.userdata)) {
		return(NULL)
	} else{
	projected.userdata[, proj:=NULL] # 105 * 4
  	setnames(projected.userdata, c( "delta", "trait"), c("Delta", "Trait")) # Changed var.proj to var.delta
  	setcolorder(projected.userdata, c("PC", "Delta","Trait"))
	}
	
	projected.userdata
}) #current

projected.table[sapply(projected.table, is.null)]  <- NULL 
projected.table <- rbindlist(projected.table)
projected.table  <- projected.table[,.(PC, Delta,Trait)]

QC.table <- data.table(Trait = names(nSNP), nSNP)

# Write tables
projtablename  <- paste0("Projection_cytokine_basis_nosig_basisTobasis_40_", date, ".tsv") #_40
qctablename  <- paste0("QC_cytokine_basis_nosig_basisTobasis_40_", date, ".tsv") #_40


message("start writing projected.table & QC.table")

write.table(projected.table, paste0("../Projections/", projtablename), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(QC.table, paste0("../Projections/", qctablename), sep = "\t", quote = FALSE, row.names = FALSE)

message("Finished writing!")