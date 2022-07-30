
# This script is to download compete version of main GWAS table (2022-07-22, 10673 rows), save locally
# Then we need to add three extra columns : "isCyto", "Common_name" (names used in the book), and "usedinBasis"(if its one of the 40 basis). 
# This is for the convenience of labeling
# We also create a meta table for Ferkingstad, but saved separately, as FK is not included in the main GWAS table yet.

# "Common_Name": based on NCBI book to rename Ahola, Folkerson, Hugland, Traglia, Reales(meta). So that all cyto/chemokines share same short name
# "isCyto": if this trait is included in 104 basis table
# "usedinBasis": if this trait is included in 40 sparse basis

library(googlesheets4)
library(data.table)
library(readxl)
library(stringr)


setDTthreads(10)


# 2022-07-25 (one-time script)


##############################################
### Download Complete Google meta data   ###
##############################################

# Specify the user:
options(gargle_oauth_email = "qingqingzhou54@gmail.com")
gs4_deauth()

# Load meta google registry table
message("Fetching metadata...")
metadata <- read_sheet("https://docs.google.com/spreadsheets/d/16B4ANehcS4psdAFReTBQUJLYXuf5_RrpmjSPaASa2Nw/edit?usp=sharing", sheet = 1)
metadata <- metadata[!is.na(metadata$File_ID),]
metadata$Trait <- sapply(strsplit(metadata$File_ID,"-"),`[`, 1) # '[' and '[[' are extraction function in R
metadata <- data.table(metadata)
metadata[,Retrieval_date:=as.character(Retrieval_date)]

#fwrite(metadata, "../google_metadata_fullUKBBFIN", sep = "\t", na = NA)


#########################################################################
###     Add "isCyto", "Common_name", "usedinBasis"into main table     ###
#########################################################################
#metadata <- fread("../google_metadata_fullUKBBFIN")

# usedinBasis
metadata$usedinBasis <- rep(FALSE, nrow(metadata))

#  Common_name : Ahola, Folkerson
cyto.folkerson <- c("CXCL8","CD154","IL-1RA","IL-6","CCL2","TRAIL","M-CSF","SCF","IL-18",
"LIGHT","TRANCE","IL-16","CXCL6","CX3CL1","CCL20") # all common names for folkerson, not limited to 40 basis 
names(cyto.folkerson) <- c("IL8","CD40L","IL1RA","IL6","MCP1","TRAIL","MCSF","SCF","IL18",
"TNFSF14","TRANCE","IL16","CXCL6","CX3CL1","CCL20") 

cyto.ahola <- c("CCL27","CCL11","G-CSF","CXCL1","IFN-γ","IL-10",
"IL-13","IL-16","IL-17","IL-18","IL-1β","IL-1RA","IL-2","IL-4",
"IL-5","IL-6","IL-7","CXCL8","IL-9","CXCL10","CCL2","CCL7","M-CSF","CXCL9","CCL3",
"CCL5","SCF","CXCL12","TNF-α","TNF-β","TRAIL") # all common names for ahola, not limited to 40 basis . remove "IL12" and "CCL4" as they are not in FK104
names(cyto.ahola) <- c("CTACK","EOT1","GCSF","GROA","IFNG","IL10",
"IL13","IL16","IL17","IL18","IL1B","IL1RA","IL2","IL4",
"IL5","IL6","IL7","IL8","IL9","IP10","MCP1","MCP3","MCSF","MIG","MIP1A",
"RANTES","SCF","SDF1a","TNFa","TNFb","TRAIL")


#metadata[Reference == 27989323, short:=str_extract(Trait,"[^_]+")]
metadata[Reference == 27989323, Common_Name:=cyto.ahola[Trait_ID_2.0]] #ahola
#metadata[Reference == 33067605, short:=str_extract(Trait,"[^_]+")]
metadata[Reference == 33067605, Common_Name:=cyto.folkerson[Trait_ID_2.0]]#folkerson

# Common_name : Reales, Traglia
metadata[First_Author == "Reales", Common_Name:=cyto.ahola[Trait_ID_2.0]] # Reales is metaanalysis performed on Hugland & Ahola


cyto.tra <- c("CCL11", "GM-CSF", "IFN-γ", "IL-10", "IL-13", "IL-17", "IL-1α", 
"IL-1β", "IL-1RA", "IL-2", "IL-4", "IL-6", "IL-7", "CXCL8", 
"CXCL10", "CCL3", "TNF-α", "CXCL13", "CCL20", "CCL21", "CCL23", 
"CCL25", "CCL27", "CX3CL1", "CXCL5", "CXCL6", "CCL24", "CCL26", 
"CXCL1", "CCL1", "IL-16", "CXCL11", "CCL8", "CCL7", "CCL13", 
"CCL22", "CXCL9", "CCL15", "CCL19", "CXCL12")

names(cyto.tra)<- c("EOT1", "GMCSF", "IFNG", "IL10", "IL13", "IL17", "IL1A", "IL1B", 
"IL1RA", "IL2", "IL4", "IL6", "IL7", "IL8", "IP10", "MIP1A", 
"TNFa", "BCA1", "CCL20", "CCL21", "CCL23", "CCL25", "CTACK", 
"CX3CL1", "CXCL5", "CXCL6", "EOT2", "EOT3", "GROA", "I309", "IL16", 
"IP9", "MCP2", "MCP3", "MCP4", "MDC", "MIG", "MIP1D", "MIP3B", 
"SDF1AB")
metadata[Reference == 30134952, Common_Name:=cyto.tra[Trait_ID_2.0]]


# Common_name : Huglund #mannually write it down in excel, then read

cyto.hug <- c("CCL19", "CCL20", "CCL23", "CCL25", "CX3CL1", "CXCL5", "CXCL6", 
"CCL11", "Flt-3L", "CXCL1", "IFN-γ", "IL-10", "IL-13", "IL-17", 
"IL-18", "IL-1α", "IL-20", "IL-4", "IL-5", "IL-6", "IL-7", "CXCL8", 
"CXCL10", "CXCL11", "LIF", "CCL2", "CCL8", "CCL7", "CCL13", "M-CSF", 
"CXCL9", "CCL3", "OSM", "SCF", "TGF-β1", "TNF-β", "LIGHT", 
"TRAIL", "TRANCE", "TWEAK")

names(cyto.hug) <- c("CCL19", "CCL20", "CCL23", "CCL25", "CX3CL1", "CXCL5", "CXCL6", 
"EOT1", "FLT3L", "GROA", "IFNG", "IL10", "IL13", "IL17", "IL18", 
"IL1A", "IL20", "IL4", "IL5", "IL6", "IL7", "IL8", "IP10", "IP9", 
"LIF", "MCP1", "MCP2", "MCP3", "MCP4", "MCSF", "MIG", "MIP1A", 
"OSM", "SCF", "TGFb1", "TNFb", "TNFSF14", "TRAIL", "TRANCE", 
"TWEAK")
metadata[Reference == 31727947, Common_Name:=cyto.hug[Trait_ID_2.0]]

# sanity check: all "Common_Name" is within "Common_Name" from Ferkingstad
# fk.com <- read_xlsx("../../Pre-work/map_dataset_names/Ferkingstad_mapped_basis.xlsx")$Common_Name
# all(metadata[!is.na(Common_Name)]$Common_Name %in% fk.com) # TRUE

# is.Cyto
metadata[, isCyto := !is.na(Common_Name)]

# usedinBasis needs to wait untill FK.meta finished to assign value
#used.Trait <- unique(fread("../Projections/Projection_cytokine_basis_sig_40sparse_basisTObasis20220719.tsv")$Trait)
#used.Common_Name <- as.data.table(read_xlsx("../../Pre-work/map_dataset_names/Ferkingstad_mapped_basis.xlsx")[,-1])
#metadata[, usedinBasis:=Common_Name %in% used.Common_Name] # logical value. 40 "TRUE" indicates these traits are used in sparse 40 basis


##############################################
###     Make metatable for Ferkingstad     ###
##############################################
# load previously mapped fk table
name.basis <- as.data.table(read_xlsx("../../Pre-work/map_dataset_names/Ferkingstad_mapped_basis.xlsx")[,-1])[Remark!="annotation added"][,c("FileName", "Protein (short name)", "Protein (full name)", "Common_Name")]%>%
                setnames(old=c("FileName", "Protein (short name)", "Protein (full name)"), new=c("File_ID","Trait_ID_2.0","Trait_long"))

# initialize output fk metatable
cols <- colnames(metadata)
fk.metadata <- as.data.table(matrix(NA, nrow=nrow(name.basis), ncol=length(cols)))
colnames(fk.metadata) <- cols

# File_ID, Trait, Trait_ID_2.0, Trait_long, Common_Name (Common_Name already in Ferkingstad_mapped_basis.xlsx table)
fk.metadata[, File_ID:=name.basis$File_ID][, Trait_ID_2.0:= name.basis$Trait_ID_2.0][, Trait_long:=name.basis$Trait_long]
fk.metadata[, Common_Name:=name.basis$Common_Name][, Trait:=gsub(".txt.gz","",File_ID)]

# First_Author, Trait_class, Collection, Population, Year, Reference, Retrieval_date
fk.metadata$First_Author <- rep("Ferkingstad",nrow(fk.metadata))
fk.metadata$Trait_class <- rep("BMK",nrow(fk.metadata))
fk.metadata$Collection <- rep("Ferkingstad",nrow(fk.metadata))
fk.metadata$Population <- rep("European (EUR)",nrow(fk.metadata))
fk.metadata$Year <- rep(2021,nrow(fk.metadata))
fk.metadata$Reference <- rep("34857953",nrow(fk.metadata))
fk.metadata$Retrieval_date <- rep(as.character(NA),nrow(fk.metadata))

# usedinBasis, isCyto
fk.metadata$isCyto <- rep(TRUE,nrow(fk.metadata))
used.Trait <- unique(fread("../Projections/Projection_cytokine_basis_sig_40sparse_basisTObasis20220719.tsv")$Trait)
fk.metadata[, usedinBasis:=Trait %in% used.Trait] # logical value. 40 "TRUE" indicates these traits are used in sparse 40 basis

# usedinBasis: metadata 
used.common <- fk.metadata[usedinBasis==TRUE]$Common_Name
metadata[, usedinBasis:=(Common_Name %in% used.common)]

# write. Eventually, we add 3 cols: "usedinBasis" (logical), "Common_Name"(char), "isCyto"(logical)
fwrite(fk.metadata, "../google_metadata_ferkingstad", sep = "\t", na = NA)
fwrite(metadata, "../google_metadata_fullUKBBFIN", sep = "\t", na = NA) #new cols: "usedinBasis"(logical),"Common_Name"(char),"isCyto"(logical)

#all(colnames(metadata) == colnames(fk.metadata))
#[1] TRUE

#tmp <- rbind(fk.metadata, metadata) # ok