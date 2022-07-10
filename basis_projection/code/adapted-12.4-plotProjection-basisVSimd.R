# This script is to plot the projection results of Ferkingstad(basis) together with IMD.
# Only selected top/bottom extreme traits are kept in plot.
# Traits exclusive to Ferkingstad will be hidden in plots.

# 2022-07-07
library(data.table)
library(stringr)
library(readxl)
library(ggplot2)
library(googlesheets4)


# read google_metadata
metadata <- fread("../google_metadata")

# Load projection tables and name mapping tables.
# proj.* files are projection results.
# name.* files are to map Trait name with "Common_name", which was originally from the cytokine review book. After name mapping, we can plot 
proj.all <- fread("../Projections/Projection_cytokine_basis_nosig_20220629.tsv")
proj.basis <- fread("../Projections/Projection_cytokine_basis_nosig_basisTobasis20220629.tsv")

name.basis <- as.data.table(read_xlsx("../../Pre-work/map_dataset_names/Ferkingstad_mapped_basis.xlsx")[,-1])[Remark!="annotation added"]
name.basis <- name.basis[,c("FileName","Common_Name")][,Trait:=gsub(".txt.gz","",FileName)]


lab <- name.basis$Common_Name
names(lab) <- name.basis$Trait
proj.basis <- proj.basis[, short:=lab[Trait]][, group:="fk"]

# extract projections of IMD from proj.all, 13 traits in total
trait.imd <- readLines("traits_to_project_early_test.txt")[249:261]
trait.imd <- gsub("-hg38.tsv.gz","",trait.imd)
proj.imd <- proj.all[Trait %in% trait.imd]
proj.imd[, short:= str_extract(Trait, "[^_]+")][, group:="imd"]

imd.lab <- c("Asthma", "Crohn's Disease", "Celiac Disease", "Eosinophilic Granulomatosis with Polyangiitis",
"ANCA negative EGPA","MPO+ EGPA", "Inflammatory Bowel Disease", "Primary Biliary Cholangitis",
"Primary Sclerosing Cholangitis", "Rheumatoid Arthritis", "Systemic Lupus Erythematosus",
"Type 1 Diabetes", "Ulcerative Colitis" )
names(imd.lab) <- c("AST", "CD","CEL","EGPA","EGPAANCAn","EGPAMPO","IBD","PBC","PSC",
"RA","SLE","T1D","UC")

proj.imd[, short:=imd.lab[short]]

# AST: Asthma
# CD: Crohn's Disease
# CEL :Celiac Disease
# EGPA: Eosinophilic Granulomatosis with Polyangiitis
# EGPAANCAn : ANCA negative EGPA
# EGPAMPO : MPO+ EGPA
# IBD : Inflammatory Bowel Disease
# PBC: Primary Biliary Cholangitis
# PSC : Primary Sclerosing Cholangitis
# RA : Rheumatoid Arthritis
# SLE : Systemic Lupus Erythematosus
# T1D : Type 1 Diabetes
# UC : Ulcerative Colitis


# Merged beta table
merged.data <- .rbind.data.table(proj.basis, proj.imd)


# top.imd.forrest function to plot top/bottom extrem basis traits together with IMDs
top.imd.forrest <- function(n.top=10, n.bottom=10, merged.data, basis.group="fk", pc=1){
    #keep only selected PC
    merged.data <- merged.data[PC==paste0("PC",pc)]

    #keep only the top/bottom trait of basis
    basis <- merged.data[group == basis.group]
    top.traits <- unique(basis$short[order(basis$Delta, decreasing = TRUE)][1:n.top])
    bot.traits <- unique(basis$short[order(basis$Delta, decreasing = FALSE)][1:n.bottom])
    imd.traits <- merged.data[group == "imd"]$short

    merged.data <- merged.data[short %in% c(top.traits, bot.traits, imd.traits)]

    fix.order.basis <-  unique(merged.data[group == basis.group]$short[order(merged.data[group == basis.group]$Delta, decreasing = FALSE)])
    fix.order.imd <- unique(merged.data[group == "imd"]$short[order(merged.data[group == "imd"]$Delta, decreasing = FALSE)])
    fix.order <- c(fix.order.imd, fix.order.basis)
    merged.data$short <- factor(merged.data$short, levels = fix.order)

    p <- ggplot(merged.data, mapping = aes(x=short, y=Delta, color=group))+
    geom_point()+
    coord_flip()+
    geom_hline(yintercept = 0, col="yellow", lty=2)+
    xlab(paste("top",n.top,"and", "bottom",n.bottom,"traits + imd"))+
    labs(title = paste0("PC",pc))+
    theme_light()

    return(p)
}

# CXCL9 and CXCL10 with imd traits
CXCL9_10.imd.forrest <- function(merged.data, pc=1){
    #keep only selected PC
    merged.data <- merged.data[PC==paste0("PC",pc)]

    #keep only "CXCL9/10" and imd traits in merged.data
    imd.traits <- merged.data[group == "imd"]$short

    merged.data <- merged.data[short %in% c("CXCL9","CXCL10", imd.traits)]

    fix.order <- unique(merged.data$short[order(merged.data$Delta, decreasing = FALSE)])
    merged.data$short <- factor(merged.data$short, levels = fix.order)

    p <- ggplot(merged.data, mapping = aes(x=short, y=Delta, color=group))+
    geom_point()+
    coord_flip()+
    geom_hline(yintercept = 0, col="yellow", lty=2)+
    xlab(paste("CXCL9, CXCL10 traits + imd"))+
    labs(title = paste0("PC",pc))+
    theme_light()

    return(p)
}

# make plots for PC 1~10
for (i in 1:10){
    top.imd.forrest(merged.data=merged.data,pc=i)
    ggsave(paste0("../Plots/imd/top_PC",i,".png"), width = 6, height = 6)
}

# make plots of CXCL9/10 for PC 1~10

for (i in 1:10){
    CXCL9_10.imd.forrest(merged.data=merged.data,pc=i)
    ggsave(paste0("../Plots/imd/CXCL9_10_PC",i,".png"), width = 6, height = 6)
}

top.imd.forrest(merged.data=merged.data,pc=5)

CXCL9_10.imd.forrest(merged.data = merged.data,pc=1)

# MS is "multiple sclerosis"
# IGAN is "IgA nephropathy"

# AST: Asthma
# CD: Crohn's Disease
# CEL :Celiac Disease
# EGPA: Eosinophilic Granulomatosis with Polyangiitis
# EGPAANCAn : ANCA negative EGPA
# EGPAMPO : MPO+ EGPA
# IBD : Inflammatory Bowel Disease
# PBC: Primary Biliary Cholangitis
# PSC : Primary Sclerosing Cholangitis
# RA : Rheumatoid Arthritis
# SLE : Systemic Lupus Erythematosus
# T1D : Type 1 Diabetes
# UC : Ulcerative Colitis