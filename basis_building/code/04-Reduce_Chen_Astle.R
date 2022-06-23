# Reducing EUR-CHEN and Astle for projecting

## Background: We need to align and filter Astle and Chen_2 (European) by the most recent manifest in order to project them and fine-tune the parameters for sparsePCA
# This script is intended to be run in the cluster.

# Guillermo Reales
# 2022-06-08

######################################################
###	 Load packages  and prepare list #############
######################################################


library(data.table)
setDTthreads(18)
library(magrittr)

arg  <- commandArgs(trailingOnly=TRUE) %>% as.numeric()

# Use selected traits only
seltraits <- c("BASC", "EOSC", "LYMC", "MONC", "PLAV", "NEUC", "PLAC", "LEUC")
#dict <- fread("../data/Chen_dict.txt") %>% .$Trait_ID_2.0
#dict <- c(dict, "ERYC") # Astle has ERYC instead of RBC

pchen1 <- paste0(seltraits, "_Chen_32888493_2") %>% paste(., collapse="|")
pastle <- paste0(seltraits, "_Astle") %>% paste(., collapse="|")
pboth  <- paste(pchen1, pastle, sep="|")
ppath  <- "~/rds/rds-cew54-basis/02-Processed/"
mfiles <- dir(ppath, pattern=pboth)
file   <- paste0(ppath, mfiles)[arg] # Select the file according to the index variable in array job
name <- paste0(strsplit(mfiles[arg], split="-")[[1]][1], "-ft.tsv.gz") # Create name to save it later

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
		mincold <- c("CHR38", "BP38","REF","ALT", "BETA", "SE", "P")
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
		}

	         M <- unique(M)
	         if(nrow(M) > nrow(man)) warning("Aligned file has more SNPs than the manifest. Some SNPs might be duplicated.")
	         M
}




#######################################################
######### Import and process ##########################
#######################################################

# Load manifest
SNP.manifest <- fread("../data/BCB4_SNP_manifest_after_Chen_8M.tsv")

# Load file
message("Working on ", file)
ds <- fread(file, tmpdir="tmp")
if("ALT_FREQ" %in% names(ds)){
ds <- ds[,.(SNPID, CHR38, BP38, REF, ALT, ALT_FREQ, BETA, SE, P)]
} else{
ds <- ds[,.(SNPID, CHR38, BP38, REF, ALT, BETA, SE, P)]

}

M <- g.align(ds, SNP.manifest)

fwrite(M, paste0("../data/", name), sep="\t")
message("Done!")




