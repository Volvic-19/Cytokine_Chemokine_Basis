# This script is to apply projection of Ohola.olli datasets (reduced by adapted-08-reduction.R) onto
# our Ferkingstad basis (complete, not sparse yet) for a quick visualisation test. Signficance is not computed. 

# Before running projection, we need to prepare Rdata to avoid repeated computing of : basis, 


library(data.table)
library(Matrix)
library(magrittr)

setDTthreads(10)

# 673,318 unique SNPs


##############################################
### LOAD Rdata and run projection
##############################################

load("../Rdata/cytokine_completeSNP.Rdata")

# Create log file
mpath <- "~/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_projection/code/log/"
date <- format(Sys.time(), format="%Y%m%d")
logname <- paste0(mpath,"log_project_",date,".txt")
file.create(logname)

# Files to project
files_to_project  <- dir("../data", pattern = ".tsv") #261 datasets to project
nfiles <- length(files_to_project)
Trait <- rep(NA,nfiles)
nSNP <- rep(NA, nfiles)

# Projection

projected.table <- lapply(files_to_project, function(file){
	message("Projecting ", file)
	trait_label <- strsplit(file, "-", fixed = TRUE)[[1]][1] #trait name
	index <- which(files_to_project == file) 
	sm <- fread(paste0("../data/",file))

	# Some checks
	sm <- unique(sm)
	sm[sm == ""] <- NA # Some missing data might pass as empty string. This will fix that	
	sm <- na.omit(sm, cols = c("pid", "BETA", "SE", "P"))
	dups <- sm$pid[duplicated(sm$pid)]

	if(length(dups) > 0){
		dupmessage= "This file has duplicated pids. I removed them prior to projection. You might want to check it."
		message(dupmessage)
		# Write to log
		write(paste0(file,". ",dupmessage), logname, append=TRUE)
		sm <- sm[!pid %in% dups] # Remove all duplicated instances, to be safe
	}

	# A bit of QC 
	Trait[index] <<- trait_label
	nSNP[index] <<- nrow(sm)

  	projected.userdata <- tryCatch(expr = project_nosig(beta=sm$BETA, seb=sm$SE, pid=sm$pid)[,trait:=trait_label][],
  	                      error =function(e) {
				      failmessage <- "Projection for this file had non-zero exit status, please check. Jumping to next file..."
				      message(failmessage)
				      write(paste0(file,". ",failmessage), logname, append=TRUE)
				      return(NULL)})
  	if(is.null(projected.userdata)) {
		return(NULL)
	} else{
	projected.userdata[, proj:=NULL] # 105 * 4
  	setnames(projected.userdata, c( "delta", "trait"), c("Delta", "Trait")) # Changed var.proj to var.delta
  	setcolorder(projected.userdata, c("PC", "Delta","Trait"))

  	# More QC
  	#overall_p[index] <<- sprintf("%1.0e",projected.userdata$p.overall[1])
  	#minc.idx <- which.min(projected.userdata$P)
  	#mscomp[index] <<- sprintf("%s (%1.0e)",projected.userdata$PC[minc.idx],projected.userdata$P[minc.idx])
	}
	projected.userdata
}
) 

projected.table[sapply(projected.table, is.null)]  <- NULL 
projected.table <- rbindlist(projected.table)
projected.table  <- projected.table[,.(PC, Delta,Trait)]

QC.table <- data.table(Trait, nSNP)

projtablename  <- paste0("Projection_cytokine_basis_nosig_", date, ".tsv")
qctablename  <- paste0("QC_cytokine_basis_nosig_", date, ".tsv")

#while(projtablename %in% dir(paste0(mpath, "cell_basis_v3_varimax/Projections"))){
#  version  <- version + 1
#  projtablename  <- paste0("Projection_cell_basis_v3_", date, "-v",version, ".tsv")
#  qctablename  <- paste0("QC_cell_basis_v3_", date, "-v",version, ".tsv")
#}

message("start writing projected.table & QC.table")

write.table(projected.table, paste0("../Projections/", projtablename), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(QC.table, paste0("../Projections/", qctablename), sep = "\t", quote = FALSE, row.names = FALSE)

message("Finished writing!")