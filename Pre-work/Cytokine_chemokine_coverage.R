## This script is to compare the coverage of common cytokine/chemokine across 4 papers 
## Common cytokine/chemokine is classified accroding to book: https://www.ncbi.nlm.nih.gov/books/NBK6294/
## It outputs 6 summary tables: 4 XXXX_mannual.xlsx, 1 chemokine_summary.xlsx, 1 cytokine_summary.xlsx

library(dplyr)
library(reshape2)
library(readxl)
library(stringr)
library(xlsx)

setwd("~/AllMphil/internship_code/")

## read table of desired cytokine and chemokine 
cytokine <- read_xlsx("extract_cytokines/common_cytokine.xlsx", sheet = 1)
chemokine <- read_xlsx("extract_cytokines/common_chemokine.xlsx", sheet = 1)

cytokine <- cytokine[-which(is.na(cytokine$Uniprot)),]
chemokine <- chemokine[-which(is.na(chemokine$Uniprot)),]

################################################################
############################  Ahola-Olli    ####################
################################################################

Ahola_Olli <- read_xlsx("extract_cytokines/Ahola_olli_traits.xlsx", sheet = 1, skip = 2)[-45,]

Ahola_Olli <- cbind.data.frame(Ahola_Olli[,1:2], 
                               Common_Name=NA,
                               Common_class=NA,
                               Common_Cytokine=NA,
                               Common_Chemokine=NA,
                               Ahola_Olli[,3:27])

#if "Abbreviation" is the same as cytokine/chemokine listed in the book

for (i in 1:nrow(Ahola_Olli)){
    if(Ahola_Olli$Abbreviation[i] %in% cytokine$Name){
        Ahola_Olli$Common_Name[i] <- Ahola_Olli$Abbreviation[i]
        Ahola_Olli$Common_Cytokine[i] <- T
    }else if(Ahola_Olli$Abbreviation[i] %in% chemokine$Name){
        Ahola_Olli$Common_Name[i] <- Ahola_Olli$Abbreviation[i]
        Ahola_Olli$Common_Chemokine[i] <- T
    }
}

#For those "Abbreviations" with "()" , extract the abbreviation in that bracket
index <- grep("\\(C",Ahola_Olli$Variable)

in.brackets <- unname(sapply(Ahola_Olli$Variable[index], function(x){
    in.brackets <- str_extract(x, "\\(C.*\\)")
    substring(in.brackets, 2, nchar(in.brackets)-1)
}))

Ahola_Olli$Common_Name[index] <- in.brackets 
Ahola_Olli$Common_Chemokine[index] <- T

#Additional mannually annotate
Ahola_Olli$Common_Name[8] <- "G-CSF"
Ahola_Olli$Common_Name[11] <- "IFN-γ"
Ahola_Olli$Common_Name[13] <- "IL-12"
Ahola_Olli$Common_Name[18] <- "IL-1β"
Ahola_Olli$Common_Name[19] <- "IL-1RA"
Ahola_Olli$Common_Name[31] <- "M-CSF"
Ahola_Olli$Common_Name[41] <- "TNF-α"
Ahola_Olli$Common_Name[42] <- "TNF-β"

#now all Common_Names are matched with the chemokine/cytokine table
## map Ahola_Olli traits to cytokine/chemokine table by "Name", annotate class information
for (i in 1:nrow(Ahola_Olli)){
    if (Ahola_Olli$Common_Name[i] %in% chemokine$Name){
        Ahola_Olli$Common_Chemokine[i] <- T
        Ahola_Olli$Common_class[i] <- chemokine$Class[which(chemokine$Name == Ahola_Olli$Common_Name[i])]
    }
    
    else if(Ahola_Olli$Common_Name[i] %in% cytokine$Name){
        Ahola_Olli$Common_Cytokine[i] <- T
        Ahola_Olli$Common_class[i] <- cytokine$Class[which(cytokine$Name == Ahola_Olli$Common_Name[i])]
        
    }
}


##label the availability of cytokines/chemokines in 4 papers: Ahola_Olli
cytokine <- cbind.data.frame(cytokine, Ferkingstad=NA, Ahola_Olli=NA, Folkersen=NA,Sinnott_Armstrong=NA)[,-c(4:7)]
chemokine <- cbind.data.frame(chemokine, Ferkingstad=NA, Ahola_Olli=NA, Folkersen=NA,Sinnott_Armstrong=NA)[,-c(4:7)]

for (i in 1:nrow(cytokine)){
    if (cytokine$Name[i] %in% Ahola_Olli$Common_Name)
    cytokine$Ahola_Olli[i] <- T
}

for (i in 1:nrow(chemokine)){
    if (chemokine$Name[i] %in% Ahola_Olli$Common_Name)
        chemokine$Ahola_Olli[i] <- T
}

################################################################
############################  Ferkingstad    ###################
################################################################

## read  data: 5284 row   21 col
Ferkingstad <- read_xlsx("extract_cytokines/Ferkingstad_traits.xlsx", sheet = 1,skip = 1)

##filter by: Organism == "human", Uniprot !== "NA", 4990 row  21 col
Ferkingstad <- subset(Ferkingstad, Ferkingstad$Organism == "human") %>%
    filter(., !UniProt == "NA")

##output format
Ferkingstad <- cbind.data.frame(Ferkingstad[, 1:7], 
                               Common_Name=NA,
                               Common_class=NA,
                               Common_Cytokine=NA,
                               Common_Chemokine=NA,
                               Ferkingstad[,8:21])

##map Ferkingstad traits to cytokine/chemokine table by "Uniprot ID"

for (i in 1:nrow(Ferkingstad)){
    if (Ferkingstad$UniProt[i] %in% chemokine$Uniprot){
        Ferkingstad$Common_Chemokine[i] <- T
        Ferkingstad$Common_class[i] <- chemokine$Class[which(chemokine$Uniprot == Ferkingstad$UniProt[i])[1]]
        Ferkingstad$Common_Name[i] <- chemokine$Name[which(chemokine$Uniprot == Ferkingstad$UniProt[i])[1]]
    }
    
    else if (Ferkingstad$UniProt[i] %in% cytokine$Uniprot){
        Ferkingstad$Common_Cytokine[i] <- T
        Ferkingstad$Common_class[i] <- cytokine$Class[which(cytokine$Uniprot == Ferkingstad$UniProt[i])[1]]
        Ferkingstad$Common_Name[i] <- cytokine$Name[which(cytokine$Uniprot == Ferkingstad$UniProt[i])[1]]
    }
}

##mannual: "P01563|P01567|P01562|P01570|P05014|P32881|P01571|P01568|P01569|P05013|P01566|P05015" all belong to "IFN-α"
IFN_a_ID <- unlist(str_split("P01563|P01567|P01562|P01570|P05014|P32881|P01571|P01568|P01569|P05013|P01566|P05015",
                             "\\|"))

IFN_a_index <- which(Ferkingstad$UniProt %in% IFN_a_ID)
Ferkingstad$Common_Name[IFN_a_index] <- "IFN-α"
Ferkingstad$Common_class[IFN_a_index] <- "Interferons"
Ferkingstad$Common_Cytokine[IFN_a_index] <- T


##label the availability of cytokines/chemokines in 4 papers: Ferkingstad
cytokine$Ferkingstad[which(cytokine$Name %in% Ferkingstad$Common_Name)] <- T
chemokine$Ferkingstad[which(chemokine$Name %in% Ferkingstad$Common_Name)] <- T


################################################################
############################  Folkersen    ####################
################################################################

Folkersen <- read_xlsx("extract_cytokines/Folkersen_traits.xlsx", sheet=1, skip=2)[1:92,]

Folkersen <- cbind.data.frame(Folkersen[, 1:5], 
                                Common_Name=NA,
                                Common_class=NA,
                                Common_Cytokine=NA,
                                Common_Chemokine=NA,
                                Folkersen[,6:23])


##map Folkersen traits to cytokine/chemokine table by "Uniprot ID"

for (i in 1:nrow(Folkersen)){
    if (Folkersen$UniProt[i] %in% chemokine$Uniprot){
        Folkersen$Common_Chemokine[i] <- T
        Folkersen$Common_class[i] <- chemokine$Class[which(chemokine$Uniprot == Folkersen$UniProt[i])[1]]
        Folkersen$Common_Name[i] <- chemokine$Name[which(chemokine$Uniprot == Folkersen$UniProt[i])[1]]
    }
    
    else if (Folkersen$UniProt[i] %in% cytokine$Uniprot){
        Folkersen$Common_Cytokine[i] <- T
        Folkersen$Common_class[i] <- cytokine$Class[which(cytokine$Uniprot == Folkersen$UniProt[i])[1]]
        Folkersen$Common_Name[i] <- cytokine$Name[which(cytokine$Uniprot == Folkersen$UniProt[i])[1]]
    }
}

##label the availability of cytokines/chemokines in 4 papers: Ferkingstad
cytokine$Folkersen[which(cytokine$Name %in% Folkersen$Common_Name)] <- T
chemokine$Folkersen[which(chemokine$Name %in% Folkersen$Common_Name)] <- T


################################################################
############################  Sunnott-Armstrong    ####################
################################################################


Sinnot_Armstrong <- read_xlsx("extract_cytokines/Sinnott_Armstrong_traits.xlsx", sheet = 1, skip = 2)

Sinnot_Armstrong <- cbind.data.frame(Sinnot_Armstrong[, 1], 
                              Common_Name=NA,
                              Common_class=NA,
                              Common_Cytokine=NA,
                              Common_Chemokine=NA,
                              Sinnot_Armstrong[,2:5])

# no match was found in Sinnot_Armstrong. Remains NA in cytokine/chemokine summary table.

##########################################################################
##output tidied table
write.xlsx(Ahola_Olli, file = "extract_cytokines/Ahola_olli_mannual.xlsx")
write.xlsx(Folkersen, file = "extract_cytokines/Folkersen_mannual.xlsx")
write.xlsx(Ferkingstad, file = "extract_cytokines/Ferkingstad_mannual.xlsx")
write.xlsx(Sinnot_Armstrong, file = "extract_cytokines/Sinnott_Armstrong_traits.xlsx")

write.xlsx(chemokine,file = "extract_cytokines/chemokine_summary.xlsx")
write.xlsx(cytokine,file = "extract_cytokines/cytokine_summary.xlsx")


