## This script is to map file names of downloaded datasets with protein name/ Uniprot ID / common protein name
## Names of the downloaded file follows the convension : "SeqId_GeneName_ProteinName.txt.gz" where SeqId is unique for each trait
## In this script and following scripts, "Common cytokine" refers to cytokine/chemokine classified in the cytokine review book.

library(dplyr)
library(reshape2)
library(readxl)
library(stringr)
library(xlsx)

setwd("~/AllMphil/internship_code/")

# read table of Ferkingstad_mannual, generated from Cytokine_chemokine_coverage.R
Ferkingstad <- read_xlsx("extract_cytokines/Ferkingstad_mannual.xlsx")[, -1]

# extract file names from url.txt file
urls <- readLines("map_dataset_names/urls.txt")
gz.urls <- urls[grep("txt.gz$",urls)] #4907, marked as "Yes" in "Included in Analysis" column

dataset_names <- gsub("https.*file\\=","",gz.urls) 


# extract unique SeqId from dataset_names
Unique_SeqID <- str_extract(dataset_names,"[^_]+\\_.*?(?=_)")

if (!all(Unique_SeqID %in% Ferkingstad$SeqId)){
    message("Check Unique_SeqID!")
}

map.df <- data.frame(SeqId = Unique_SeqID,
                     FileName = dataset_names)


# merge FileName to Ferkingstad table based on SeqId

output.all <- merge(map.df, Ferkingstad, by = "SeqId", all = TRUE)

output.basis <- subset(output.all, !is.na(output.all$FileName)) %>%
    filter(., !is.na(Common_Name)) #105 in total(100 unique Uniprot IDs), 86 unique common proteins #Hence we get 62 cytokines datasets (of  49 unique cytokines), 43 chemokine datasets (of 37 unique chemokines)


# save files

write.xlsx(output.all, "map_dataset_names/Ferkingstad_mapped_all.xlsx")
write.xlsx(output.basis, "map_dataset_names/Ferkingstad_mapped_basis.xlsx") 



#############################################################################
############# Solve duplicate problem in datasets############################
#############################################################################

basis <- read_xlsx ("map_dataset_names/Ferkingstad_mapped_basis.xlsx")[, -1]
all <- read_xlsx ("map_dataset_names/Ferkingstad_mapped_all.xlsx")[, -1]


# mannually remove "9278_9_EPYC_Epiphycan.txt.gz" as it was mislabeled as "SDF-1" with UniprotID "P48061" by the original author.
basis <- subset(basis, !basis$SeqId == "9278_9")


# Keep only one dataset when there are datasets sharing same Uniprot ID

## keep "Platelet basic protein (PBP)" as the only protein of "P02775". As PBP can be cleaved into: 
## Beta-thromboglobulin (BTG), Neutrophil-activating peptide 2 (NAP -2), Connective tissue-activating peptide III (CTAP-3)

#index.to.remove <- setdiff(which(basis$UniProt == "P02775"), which(basis$`Protein (short name)` == "PBP"))
#basis <- basis[-index.to.remove, ]

## 
#dup.uniprot <- basis$UniProt[duplicated(basis$UniProt)]
#dup.basis <- subset(basis, basis$UniProt %in% dup.uniprot)

#dup.all <- all$UniProt[duplicated(all$UniProt)]
#dup.all <- subset(all, all$UniProt %in% dup.all)

 
