#########################################
## DATASET FILTERING BY THE MANIFEST (e5), as part of the "quick projection"
#########################################

# Introduction: This code is meant to pre-process our big files 
#(located at 02-Liftovered/), filtering them by the SNPs in the 
# SNP manifest, to make them more manageable
# prior to project them onto the Cytokine cell basis

# This script requires "ALT_FREQ" in input data. 1-248 was generated from this script 


##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################

library(data.table)
library(magrittr)
setDTthreads(10)


#load("~/rds/rds-cew54-basis/03-Bases/cell_basis_v2/cell-basis-sparse-2.0.RData")
#setwd("/home/qz284/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_projection/code")

# Load manifest
SNP.manifest <- fread("../../basis_building/manifest/IL5_consensus_manifest_6M_e5.tsv")

######################################################
###	 Load aligner functions		 #############
######################################################

### g.complement, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
g.complement <- function (x) {
  x <- toupper(x)
  switches <- c(A = "t", T = "a", C = "g", G = "c")
  for (i in seq_along(switches)) x <- sub(names(switches)[i], switches[i], x)
  toupper(x)
}


### g.rev, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
g.rev <- function (x, sep = "/") {
  sapply(strsplit(x, sep), function(g) paste(rev(g), collapse = "/"))
}

g.class  <- function(x, y, flip_strand=TRUE){
      	diag <- rep(NA, length(x))
	diag[x == g.rev(y)]  <- "rev"
	diag[x == y] <- "nochange"
	if(flip_strand){
	diag[x == g.complement(y)]  <- "comp"
	diag[x == g.rev(g.complement(y))]  <- "revcomp"
	diag[x %in% c("A/T", "T/A", "G/C","C/G")] <- "ambig"
	}
	diag[x != y && (grepl("-", x) | x %in% c("I/D","D/I"))]  <- "indels"
	diag[is.na(diag)]  <- "impossible"
	diag
}

# Build the function
### g.align. New function I made as a wrapper to align datasets to the manifest using the previous functions. It will also flip ALT_FREQ, if available.
# Note 
g.align <- function(ds, manifest){
		
		# Copy both datasets
		d <- copy(ds)
		man <- copy(manifest)
		
		if(!all(sapply(list(d,man), is.data.table))){
			d <- data.table(d)
			man <- data.table(man)
		}
		# Sanity checks. We'll require these columns initially
		mincold <- c("CHR38", "BP38","REF","ALT", "BETA", "SE", "P","ALT_FREQ")
		if(!all(mincold %in% names(d))) stop("Minimum columns missing from dataset to be aligned, these are: ", paste(mincold, collapse=", "))
		mincolm <- mincold[1:4]
		if(!all(mincolm %in% names(man))) stop("Minimum columns missing from manifest to align dataset to be aligned to, these are: ", paste(mincolm, collapse=", "))

		# Create allele vector for manifest
		 man[,alleles:=paste(REF,ALT, sep="/")][, pid:=paste(CHR38, BP38, sep=":")]
		
#		if("ALT_FREQ" %in% names(d)) mincold <- c(mincold, "ALT_FREQ") # Include ALT_FREQ, if available
#		d <- d[, ..mincold] # Use only relevant columns
		d[,REF:=toupper(REF)][,ALT:=toupper(ALT)][,alleles:=paste(REF,ALT,sep="/")][, pid:=paste(CHR38, BP38, sep=":")]
		M <- merge.data.table(d, man[,.(pid,alleles)], by='pid', suffixes=c(".d",".m"))
        
		# Diagnose alleles
		M[, diag:=g.class(alleles.m, alleles.d)]
		if(!all(M$diag == "nochange")){
	     	    cat("Some SNPs have to be flipped. ", sum(M$diag == "rev"), " to flip, ", sum(M$diag == "comp"), " to find their complement, and ", sum(M$diag == "revcomp"), " to find their reverse complement.\n")
		    M[, alleles.f:= alleles.d]
		    M[diag == "rev", alleles.f:= unlist(g.rev(alleles.d))] # Flip reversed alleles
		    M[diag == "comp", alleles.f:= g.complement(alleles.d)] # Find complement
		    M[diag == "revcomp", alleles.f:= unlist(g.rev(g.complement(alleles.d)))] # Find rev comp
		    
		    # Remove those SNPs that couldn't be successfuly aligned
		    M[, diag2:=g.class(alleles.m, alleles.f)]
		    if(!all(M$diag2 == "nochange")){
	     	    cat("Unfortunately, ", sum(M$diag2 == "ambig"), " SNPs are ambiguous, and ", sum(M$diag2 == "impossible"), " were impossible to align. These will be removed now.\n")
			M <- M[ diag2 == "nochange"]	
		    }
		    
		    # Now update alleles, BETA, and ALT_FREQ
		    M[ diag %in% c("rev", "revcomp"), BETA:= -BETA]
		    if("ALT_FREQ" %in% names(M)){
		    	M[ diag %in% c("rev", "revcomp"), ALT_FREQ:= 1-ALT_FREQ]
		    }
		    M[, c("REF", "ALT") := tstrsplit(alleles.f, "/")] 
	            M[, c("alleles.d", "alleles.m", "alleles.f", "diag", "diag2"):=NULL]
		}else{
		    cat("all no change!")
		    M[, c("alleles.d", "alleles.m", "diag"):=NULL]
		}

	         M <- unique(M)
	         if(nrow(M) > nrow(man)) warning("Aligned file has more SNPs than the manifest. Some SNPs might be duplicated.")
	        
		 return(M)
}

# This function check.tidy 
check.tidy <- function(d){
    mincolm <- c("REF_ALT", "pid")
    
    if (!all(mincolm %in% colnames(d))){
        stop("dataset does not have minimun column REF_ALT and pid!")
    }else{
        
        # check if there is any duplicated SNP left in ds
        if(any(duplicated(d$pid))){
            
            duplicated.pid <- d$pid[which(duplicated(d$pid))]
            message(paste0("Detected duplicated pid! Removing "), length(duplicated.pid), " duplicated pid")
            d <- d[!pid %in% duplicated.pid]
        }# check if there is any odd REF_ALT (A/A, T/T, C/C, G/G) left in ds
        all.allele <- unique(d$REF_ALT)
        odd.allele <- c("A/A","T/T","C/C","G/G")
        if (any(odd.allele %in% all.allele)){
            d <- d[!REF_ALT %in% odd.allele]
            message("Detected odd REF_ALT: A/A , T/T, C/C or G/G, removed from input")
        }

		if (any(nchar(d$REF_ALT)!=3)){
			d <- d[nchar(REF_ALT)!=3]
			message("Detected odd REF_ALT: nchar > 3, removed from input")
		}
        
    }
    return(d)
}


#############################
## DATA INPUT
#############################



# Get the range of files from args
# This time we'll use a different strategy, involving array jobs.
# Since we want to control the size of the batch to be similar in all cases, we'll create another variable to store the index + number of extra jobs
file <- commandArgs(trailingOnly = TRUE)
message(paste0("Processing file ", file))

ppath <- paste0("~/rds/rds-cew54-basis/02-Processed/" , file)


#  run g.align to align, remove SNPs that are unalignable or not included in manifest
input <- fread(ppath, tmpdir = "tmp")
input <- input[, c("SNPID", "CHR38", "BP38","REF","ALT", "BETA", "SE", "P", "ALT_FREQ")]
input[,REF_ALT:=paste(REF,ALT,sep="/")][,pid:=paste(CHR38,BP38,sep=":")]

M <- g.align(input, SNP.manifest)

#  run check.tidy to remove any left duplicated pid / odd REF_ALT
M <- check.tidy(M)

# write file
rm(input)
newname <- strsplit(file, split = "-")[[1]][1]

fwrite(M, paste0("~/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_projection/data/",newname,"-ft.tsv"), sep = "\t")
cat("Done!\n")

 
