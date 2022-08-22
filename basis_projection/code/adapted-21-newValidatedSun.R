# This script is to analyze projection from Sun to our sparse Ferkingstad 40 basis ("../Rdata/cytokine_sparseSNP_40basis-LD-ProjFun.RData")
# We get more validated cyto/chemokines from Sun projection. And further look into validated cyto/chemokines 
# This script: (A) rename Sun proteins, merge Sun to Ahola/Folkersen/Hog heatmap 

# 2022-08-02

library(data.table)
library(Matrix)
library(magrittr)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(gridExtra)
library(stringr)
library(readxl)

setDTthreads(10)

#### load Sun projection, rename as Common_Name, merge Sun with  ####
load("~/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_projection/Rdata/sun-projected.RData") # datatable P

name.basis <- as.data.table(read_xlsx("../../Pre-work/map_dataset_names/Ferkingstad_mapped_basis.xlsx")[,-1])[Remark!="annotation added"][,c("FileName", "Protein (short name)", "Protein (full name)","UniProt", "Common_Name")]%>%
                setnames(old=c("FileName", "Protein (short name)", "Protein (full name)"), new=c("File_ID","Trait_ID_2.0","Trait_long"))

PT.filt <- P[fdr.overall < 0.01] # 13640 rows
PT.filt[p < 0.05, stars:= "﹡"]
PT.filt[fdr < 0.05, stars:="○"][fdr < 0.01, stars:="●"]
PT.filt[is.na(stars), stars:=""] # 275 ● , 268 ○ , 2409 ﹡

#### create "Label" for Sun projection
# Get "Common_Name" from UniProt ID
PT.name <- merge.data.table(PT.filt, name.basis[,.(UniProt, Common_Name)], by="UniProt")

# Keep Only "Common_Name"s used in 40 basis
fk.metadata <- fread("../google_metadata_ferkingstad")
used.common <- fk.metadata[usedinBasis==TRUE]$Common_Name
PT.name <- PT.name[Common_Name %in% used.common] # 1120 = 28 * 40

# make label

PT.sigviz.sun <- PT.name[,.(PC, Common_Name, protein, delta, stars)]
setnames(PT.sigviz.sun, old = "delta", new = "Delta")
# View(PT.sigviz.sun[PC=="PC1"]) # some duplicated names, need to distinguish them
PT.sigviz.sun[, Trait_class:="BMK"]
PT.sigviz.sun[, tmp:=paste0(PC, protein)]
PT.sigviz.sun <- PT.sigviz.sun[!duplicated(tmp)][, tmp:=NULL] # remove duplicated protein

make.labels.sun <- function(x){
  x[,Label:= paste0(Common_Name, " / ", "Sun")] %>% 
    .[,tmp:=1] %>%
    .[protein%in%c("CCL15.3509.1.1","CCL23.3028.36.2","CCL25.2705.5.2","CXCL12.3516.60.2"), tmp := 2]%>%
    .[,Label:=paste0(Label,"_",tmp)] %>%
    .[, tmp:=NULL]
}

make.labels.sun(PT.sigviz.sun)
PT.sigviz.sun <- PT.sigviz.sun[,.(PC, Trait_class, Label, Delta, stars)]
QC.sigviz.sun <- PT.sigviz.sun[PC=="PC1"][,.(Trait_class, Label)]

#### load previous Ahola/Folkersen/Hog plot.PT, plot.QC table ####
previous.PT <- fread("../plotdata/cytoinBasis_PT")[, .(PC, Trait_class, Label, Delta, stars)] # 720 (= 18 * 40) # previous input for pheamap function
previous.QC <- fread("../plotdata/cytoinBasis_QC")[, .(Trait_class, Label)] # 18 # previous input for pheamap function

plot.PT <- rbind(PT.sigviz.sun, previous.PT)
plot.QC <- rbind(QC.sigviz.sun, previous.QC)

#### pheatmap function ####

pheamap.plot <- function(PT.sigviz, QC.sigviz){ # PT.sigviz is the projections to show in the  Pheamap, you can define a subset

show.class <- unique(PT.sigviz$Trait_class)

hmcol <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(100)
PCorder <-  paste0("PC",1:40)

M.sig <- acast(PT.sigviz[,c("PC", "Label", "Delta")], Label ~ PC) # Melted to 543 * 40 matrix. 141 is the length of unique label
M.sig.stars <- acast(PT.sigviz[,c("PC","Label","stars")], Label ~ PC)

M.sig <- M.sig[,PCorder]
M.sig.stars <- M.sig.stars[,PCorder]

range <- max(abs(M.sig))
# [1] 0.07372398

# Let's also annotate by type of trait
Trait_class = c(BMK = "#06D6A0", IMD = "#EF476F", INF = "#FFD166", PSD = "#F9DAD0", OTH ="#26547c")
Trait_class = Trait_class[show.class]

#annviz <- data.frame(Trait_class = QC.sigviz[Trait_class %in% show.class]$Trait_class, row.names = unique(PT.sigviz[Trait_class %in% show.class]$Label)) # QC.sigviz doesn't have "Label" column, but traits are in the same order
annviz <- data.frame(Trait_class = QC.sigviz$Trait_class, row.names = QC.sigviz$Label) # QC.sigviz doesn't have "Label" column, but traits are in the same order
colannviz = list(Trait_class = Trait_class)

pheatmap(M.sig, breaks = seq(-range, range, length.out = 100), cluster_cols = FALSE, display_numbers = M.sig.stars, fontsize_row = 12, fontsize_number = 20, color = rev(hmcol), annotation_row = annviz, annotation_colors = colannviz)

}


# save plot 
external.cyto <- pheamap.plot(PT.sigviz = plot.PT, QC.sigviz = plot.QC)
ggsave("../Plots/imd_40_sparseAll10541/externalCytoHeat.png",external.cyto, width=13, height = 10)

sun.cyto <- pheamap.plot(PT.sigviz = PT.sigviz.sun, QC.sigviz = QC.sigviz.sun)
ggsave("../Plots/imd_40_sparseAll10541/sunCytoHeat.png",sun.cyto, width=10, height = 5)


# save plot data
fwrite(plot.PT, file="../plotdata/All_heatmap_cyto_PT", sep = "\t", na = NA)
fwrite(plot.QC, file="../plotdata/All_heatmap_cyto_QC", sep = "\t", na = NA)
