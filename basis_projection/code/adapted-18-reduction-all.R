########################################################################################################
## DATASET FILTERD BY THE MANIFEST (FK e5 sparse), as pre-work of the complete projection on FK-40-sparse basis #
#######################################################################################################

# This script [does not require] "ALT_FREQ" in input data. (And we don't need ALT_FREQ in our projection)

# Introduction: This code is meant to pre-process our big files 
#(located at 02-Processed/), in total [10541] external datasets, filtering them by the SNPs in the 
# final e5 sparse SNP manifest, to make them more manageable prior to project them onto the Cytokine basis

# NOTE: [ UKBB , 1836] and [ FinnGen , 2803] that are not ready to project previously are updated now, as they were pre-filtered by old manifest. 
# All traits that did not fit the coverage cutoff (80%) will be removed
# at QC step after projection.

# sbatch --array 1-10541 slurm_adapted-18-reduce [failed, as exceeds array limit?]
#  sbatch --array 1-9000 slurm_adapted-18-reduce [succeed] [stopped after 64565215_5757.err ]
# sbatch --array 5758-10541 slurm_adapted-18-reduce [failed, as exceeds array limit?]
# sbatch --array 1-4784 slurm_adapted-18-reduce [expr array_id + 5757][stop at 1021]
# sbatch --array 1022-4784 slurm_adapted-18-reduce [change to 30min time]
# sbatch --array 3557-4784 slurm_adapted-18-reduce

# 2022-07-21

##############################################
### LOAD LIBRARIES AND SET REQUIRED FUNCTIONS
##############################################

library(data.table)
library(magrittr)
setDTthreads(10)
library(IMDtools)


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

# Load manifest
SNP.manifest <- fread("../../basis_building/manifest/CCL8_consensus_manifest_6M_e5_sparse.tsv") # 5519 * 8

# This time we use g.align function from IMDtools. Which does not require ALT_FREQ, but will alter it when it is available
# test
# ds <- fread("~/rds/rds-cew54-basis/02-Processed/EGPA_Lyons_31719529_1-hg38.tsv.gz")
# ds <- ds[, c("SNPID", "CHR38", "BP38","REF","ALT", "BETA", "SE", "P")] #"ALT_FREQ" is not in output. As we don't need that for proj.
# ds[,REF_ALT:=paste(REF,ALT,sep="/")][,pid:=paste(CHR38,BP38,sep=":")]
# a <- IMDtools::g.align(ds = ds, manifest = SNP.manifest)
# a <- check.tidy(a)

# A <- fread("../data_40_sparse/EGPA_Lyons_31719529_1-ft.tsv")
#  all(a==A)
# [1] TRUE


#############################
## DATA INPUT
#############################

# Get the range of files from args

file <- commandArgs(trailingOnly = TRUE)
message(paste0("Processing file ", file))

ppath <- paste0("~/rds/rds-cew54-basis/02-Processed/" , file)


#  run g.align to align, remove SNPs that are unalignable or not included in manifest
input <- fread(ppath, tmpdir = "tmp")
input <- input[, c("SNPID", "CHR38", "BP38","REF","ALT", "BETA", "SE", "P")] #"ALT_FREQ" is not in output. As we don't need that for proj.
input[,REF_ALT:=paste(REF,ALT,sep="/")][,pid:=paste(CHR38,BP38,sep=":")]

M <- g.align(input, SNP.manifest)

#  run check.tidy to remove any left duplicated pid / odd REF_ALT
M <- check.tidy(M)

# write file 

rm(input)
newname <- strsplit(file, split = "-")[[1]][1]

fwrite(M, paste0("~/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_projection/data_40_sparse_allcollection10541/",newname,"-ft.tsv"), sep = "\t")
cat("Done!\n")



