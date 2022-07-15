# This script is to plot the projection results of Ferkingstad(basis-40), Ahola and Folkersen together.
# Only selected top/bottom extreme traits are kept in plot.
# Traits exclusive to Ferkingstad will be hidden in plots.

# 2022-07-11

library(data.table)
library(stringr)
library(readxl)
library(ggplot2)
library(googlesheets4)
# Specify the user:
options(gargle_oauth_email = "qingqingzhou54@gmail.com")
gs4_deauth()

# Load meta google registry table

message("Fetching metadata...")
metadata <- read_sheet("https://docs.google.com/spreadsheets/d/16B4ANehcS4psdAFReTBQUJLYXuf5_RrpmjSPaASa2Nw/edit?usp=sharing", sheet = 1)
metadata <- metadata[,c("File_ID", "First_Author", "Year", "Reference","Trait_ID_2.0", "Trait_long")]
metadata <- metadata[!is.na(metadata$File_ID),]
metadata$Trait <- sapply(strsplit(metadata$File_ID,"-"),`[`, 1)
metadata <- data.table(metadata)

fwrite(metadata, "../google_metadata", sep = "\t", na = NA)

# read google_metadata
metadata <- fread("../google_metadata")

# Load projection tables and name mapping tables.
# proj.* files are projection results.
# name.* files are to map Trait name with "Common_name", which was originally from the cytokine review book. After name mapping, we can plot 
# same trait from different studies together in one plot.

proj.all <- fread("../Projections/Projection_cytokine_basis_nosig_40_20220711.tsv") #40
proj.ahola <- proj.all[grep("AholaOlli",proj.all$Trait)] # 1681 = 41 traits * 41
proj.folkerson <- proj.all[grep("Folkersen",proj.all$Trait)] # 3690 = 90 traits * 41
proj.basis <- fread("../Projections/Projection_cytokine_basis_nosig_basisTobasis_40_20220711.tsv") #40


name.ahola <- as.data.table(read_xlsx("../../Pre-work/map_dataset_names/Ahola_olli_mannual.xlsx")[,-1])
name.folkerson <- as.data.table(read_xlsx("../../Pre-work/map_dataset_names/Folkersen_mannual.xlsx")[,-1])
name.basis <- as.data.table(read_xlsx("../../Pre-work/map_dataset_names/Ferkingstad_mapped_basis.xlsx")[,-1])[Remark!="annotation added"]

name.ahola <- name.ahola[,.(Variable, Abbreviation,Common_Name)] #,Common_class,Common_Cytokine,Common_Chemokine
name.folkerson <- name.folkerson[,.(Gene, Short_annotation, Common_Name)]
name.basis <- name.basis[,c("FileName","Common_Name")][,Trait:=gsub(".txt.gz","",FileName)]

# Keep only cytokine/chemokine traits in Folkersen and Ahola (basis does not need that, because Trait is already in name.basis)
# Add "short" to projection file, for name mapping
# Add "group" to projection file, states where the data source is. 
proj.folkerson <- proj.folkerson[, short:=str_extract(Trait,"[^_]+")][, group:= "folkerson"]
proj.ahola <- proj.ahola[, short:=str_extract(Trait, "[^_]+")][, group:="ahola"]

lab <- name.basis$Common_Name
names(lab) <- name.basis$Trait
proj.basis <- proj.basis[, short:=lab[Trait]][, group:="fk"] #40 unique traits from proj.basis
proj.basis[,Common:=short]

unique(proj.folkerson$short)
unique(proj.ahola$short)

# extract cytokines in Folkerson and Ahola, rename them with "Common_name"
# Had to this step mannually

cyto.folkerson <- c("CXCL8","CD154","IL-1RA","IL-6","CCL2","TRAIL","M-CSF","SCF","IL-18",
"LIGHT","TRANCE","IL-16","CXCL6","CX3CL1","CCL20") # common names for folkerson 
names(cyto.folkerson) <- c("IL8","CD40L","IL1RA","IL6","MCP1","TRAIL","MCSF","SCF","IL18",
"TNFSF14","TRANCE","IL16","CXCL6","CX3CL1","CCL20") 

cyto.ahola <- c("CCL27","CCL11","G-CSF","CXCL1","IFN-γ","IL-10","IL-12",
"IL-13","IL-16","IL-17","IL-18","IL-1β","IL-1RA","IL-2","IL-4",
"IL-5","IL-6","IL-7","CXCL8","IL-9","CXCL10","CCL2","CCL7","M-CSF","CXCL9","CCL3","CCL4",
"CCL5","SCF","CXCL12","TNF-α","TNF-β","TRAIL") # common names for ahola
names(cyto.ahola) <- c("CTACK","EOT1","GCSF","GROA","IFNG","IL10","IL12",
"IL13","IL16","IL17","IL18","IL1B","IL1RA","IL2","IL4",
"IL5","IL6","IL7","IL8","IL9","IP10","MCP1","MCP3","MCSF","MIG","MIP1A","MIP1B",
"RANTES","SCF","SDF1a","TNFa","TNFb","TRAIL")

proj.folkerson <- proj.folkerson[, Common:=cyto.folkerson[short]]
proj.ahola <- proj.ahola[, Common:=cyto.ahola[short]]

proj.folkerson <- proj.folkerson[!is.na(Common)] #keep only cytokines in folkerson # 615= 15 traits * 41.
proj.ahola <- proj.ahola[!is.na(Common)][!Common %in% c("CCL4","IL-12")] #1271 = 31 traits * 41. #keep only cytokines in ahola (* note only CCL4 and IL-12 are not present in 104 basis, so remove it from ahola too)


# Check if all ahola, folkerson traits are in basis 40
all(unique(proj.ahola$Common) %in% proj.basis$Common) #31 traits * 41
#[1] FALSE
all(unique(proj.folkerson$Common) %in% proj.basis$Common) #15 traits * 41
#[1] FALSE

# Keep only ahola, folkerson traits that are also in basis 40
proj.folkerson <- proj.folkerson[Common %in% unique(proj.basis$Common)] # Folkersen left: 287 = 7 traits * 41
proj.ahola <- proj.ahola[Common %in% unique(proj.basis$Common)] # ahola left: 492 = 12 traits * 41


# Merged beta table
merged.data <- .rbind.data.table(proj.basis, proj.ahola, proj.folkerson)

# top.forrest function to plot top/bottom extrem basis traits
top.forrest <- function(n.top=10, n.bottom=10, merged.data, basis.group="fk", pc=1, rm.non.overlap=FALSE){
    #keep only selected PC
    merged.data <- merged.data[PC==paste0("PC",pc)]

    #keep only the top/bottom trait of basis
    basis <- merged.data[group == basis.group]
    top.traits <- unique(basis$Common[order(basis$Delta, decreasing = TRUE)][1:n.top])
    bot.traits <- unique(basis$Common[order(basis$Delta, decreasing = FALSE)][1:n.bottom])


    merged.data <- merged.data[Common %in% c(top.traits, bot.traits)]
    # remove traits from plot that are exclusive to Ferkingstad
    if(rm.non.overlap){
        shared.common <- unique(merged.data[group!=basis.group]$Common)
        merged.data <- merged.data[Common %in% shared.common]
    }

    fix.order <-  unique(merged.data[group == basis.group]$Common[order(merged.data[group == basis.group]$Delta, decreasing = FALSE)])
    merged.data$Common <- factor(merged.data$Common, levels = fix.order)

    p <- ggplot(merged.data, mapping = aes(x=Common, y=Delta, color=group))+
    geom_point()+
    coord_flip()+
    geom_hline(yintercept = 0, col="yellow", lty=2)+
    xlab(paste("top",n.top,"and", "bottom",n.bottom,"traits"))+
    labs(title = paste0("PC",pc))+
    theme_light()

    return(p)
}

# make plots for PC 1-20

for (i in 1:20){
    p <- top.forrest(15, 15, merged.data, "fk", i , rm.non.overlap = TRUE)
    ggsave(filename = paste0("../Plots/top_bot_40/rmPC",i,".png"), plot = p,width = 4, height = 4)
}

# make plots for PC 21-41

for (i in 21:41){
    p <- top.forrest(15, 15, merged.data, "fk", i , rm.non.overlap = TRUE)
    ggsave(filename = paste0("../Plots/top_bot_40/rmPC",i,".png"), plot = p,width = 4, height = 4)
}




#top.forrest(15,15,merged.data,"fk",1,rm.non.overlap = TRUE)
#top.forrest(15,15,merged.data,"fk",2,rm.non.overlap = TRUE)
#top.forrest(15,15,merged.data,"fk",3,rm.non.overlap = TRUE)
#top.forrest(15,15,merged.data,"fk",4,rm.non.overlap = TRUE)
#top.forrest(15,15,merged.data,"fk",5,rm.non.overlap = TRUE)

#ggsave("../Plots/top_bot_40/rmPC1.png", width = 4, height = 4)