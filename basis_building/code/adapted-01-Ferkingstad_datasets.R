# Preparing manifest-filtered (CCL8 8M) Ferkingstad datasets prior to basis creation.

# We'll now import the 40 Ferkingstad datasets of cytokine and chemokine.
# We'll filter them by the initial manifest, keep relevant columns, check alignment and that all SNPs are common among all files and put them together for the next step.

# Slurm: sbatch --array 1-40 slurm_01filter_FK

# 2022/07/07: previous mistake of REF_ALT is corrected here.
# 2022/07/10: use CCL8 initial manifest, and use 40 selected FK datasets instead of all 104.

# Load libraries
library(data.table)
setDTthreads(0)


### fk.annotate function to prepare raw Ferkingstad datasets for g.align
### According to Ferkingstad paper, not all MAF in raw data are Effect allele frequency. Therefore we need to annotate those SNPs, change to correct ALT_FREQ
### fk.annotate add ALT_FREQ, set column names, create pid.

fk.annotate <- function(fk,annotation){
        ## remove ChrX 
	fk <- fk[Chrom!="chrX"] #
	
	## add annotated ALT_FREQ by "name"
	fk.m <- merge(fk, annotation[,.(Name, effectAllele,effectAlleleFreq)], by="Name", all.x=TRUE, suffixes=c(".f",".a"))

	## sanity check: if effect allele from annotation file is all the same as raw Ferkingstad dataset
	comp <- fk.m[!is.na(effectAlleleFreq)]
	if (all(comp$effectAllele.a == comp$effectAllele.f)){
	
		message("matched effectAllele ^^")
		fk.m[, effectAllele.a:=NULL]
		fk.m[, effectAllele:=effectAllele.f]
		fk.m[, effectAllele.f:=NULL]
		fk.m[is.na(effectAlleleFreq), effectAlleleFreq:=ImpMAF] #if not in annotation file, keep original ImpMAF as effectAlleleFreq
		
	}else{
		error("mismatched effectAllele! Check fk effectAllele and annotation effectAllele!")
	}

	## change format of Chrom , create pid, create REF_ALT remove unwanted columns, convert CHR38 and BP38 to int
	fk.m[ , CHR38:=gsub("chr", "", Chrom)]
	fk.m[ , pid:=paste(CHR38, Pos, sep=":")]
	fk.m[ , BP38:=as.integer(Pos)]
	fk.m[ , CHR38:=as.integer(CHR38)]
	fk.m[ , REF_ALT:=paste(otherAllele, effectAllele, sep = "/")] 
	# Previous REF_ALT mistake corrected. 2022/07/07

	fk.m[ , c("Chrom","Pos","minus_log10_pval"):=NULL]


	## keep only 1 to 1 substitution alleles

	fk.m <- fk.m[nchar(REF_ALT)==3]

	## reset column names for g.align

	setnames(fk.m, old = c("otherAllele","effectAllele","effectAlleleFreq", "Beta","Pval"), 
		new = c("REF","ALT","ALT_FREQ","BETA","P")) 

	setcolorder(fk.m, neworder = c("pid","CHR38","BP38","REF","ALT","REF_ALT","Name","rsids","BETA",
		"P","SE","N","ImpMAF","ALT_FREQ"))

	## Remove "A/A", "T/T","C/C","G/G"
	fk.m <- fk.m[!REF_ALT %in% c("A/A","T/T","C/C","G/G")] #mistaken REF_ALT should not cause a difference here

	fk.m[order(CHR38, BP38)]
}




############################
# Helper functions   #######
############################

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
		mincold <- c("CHR38", "BP38","REF","ALT", "BETA", "SE", "P","ALT_FREQ") # mistaken REF_ALT was not required in this function
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


#########################################################

# Get file name from argument number
args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args)


# Set path and variables
ppath <- "~/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_building/raw_fk_data_40"
file <- dir(path=ppath)[i]

trait <- gsub(".txt.gz","",file)

# load manifest 
manifest <- fread("../manifest/CCL8_consensus_manifest_8M.tsv", tmpdir = "tmp")
annotation <- fread("../manifest/assocvariants.annotated.txt.gz", tmpdir = "tmp")

# process file
fk <- fread(paste0(ppath, "/", file), tmpdir = "tmp")

message("Annotating",trait)
fk.m <- fk.annotate(fk,annotation)

# Align to initial manifest
# This time will not have a lot of "no change", as our initial manifest is no longer a simple subset of Ferkingstad IL5,
# itself was aligned to Refman. So a lot of SNPs would have to be flipped.

message("Aligning", trait)
M <- g.align(fk.m, manifest)

M <- merge(M, manifest[, .(pid, ld.block)], by="pid")


# write file

fwrite(M, paste0("../manifested_data/",trait, "-ft.tsv.gz"), sep="\t")









