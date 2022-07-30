# This script is to apply projection of 10542(5894) external datasets from 02-processed/ (reduced by adapted-18-reduction.R) onto
# our Ferkingstad basis (sparse, PC1-40) for a complete projection.  Significance is computed here.

# In the plotting part, Previously projected Ferkingstad basis (from adapted-16-projection-FK40sparse) will be plotted together with
# our external datasets, but only those with : (A) QC passed 80% of coverage of sparse SNP
# (B) only traits marked as "IMD" in "Triat_class" will be shown (C) Overall significant, and sig in each PC
# will be shown in our plot

# 2022-07-26

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

setDTthreads(10)

# 5519 unique SNPs



##############################################
### LOAD sparse Rdata and run projection   ###
##############################################

load("../Rdata/cytokine_sparseSNP_40basis-LD-ProjFun.RData")

# Create log file
mpath <- "~/rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_projection/code/log/"
date <- format(Sys.time(), format="%Y%m%d")
logname <- paste0(mpath,"log_project_",date,".txt")
file.create(logname)

# Files to project
files_to_project  <- dir("../data_40_sparse_allcollection", pattern = ".tsv") # 5894 datasets to project
nfiles <- length(files_to_project)
Trait <- rep(NA,nfiles)
nSNP <- rep(NA, nfiles)
overall_p <- rep(NA, nfiles)
mscomp <- rep(NA,nfiles)


###############################################
# Projection of 1-5894 external traits		#	
##############################################

projected.table <- lapply(files_to_project, function(file){
	message("Projecting ", file)
	trait_label <- strsplit(file, "-", fixed = TRUE)[[1]][1] #trait name
	index <- which(files_to_project == file) 
	sm <- fread(paste0("../data_40_sparse_allcollection/",file)) #data_40_sparse_allcollection

	# Some checks
	sm <- unique(sm)
	sm[sm == ""] <- NA # Some missing data might pass as empty string. This will fix that	
	sm <- na.omit(sm, cols = c("pid", "BETA", "SE", "P")) # 
	dups <- sm$pid[duplicated(sm$pid)]
    ## Attention, no Inf/-Inf in Beta or SE!! #2022-07-18
    sm <- sm[!BETA %in% c(Inf, -Inf)][!SE %in% c(Inf, -Inf)]

	if(length(dups) > 0){
		dupmessage= "This file has duplicated pids. I removed them prior to projection. You might want to check it."
		message(dupmessage)
		# Write to log
		write(paste0(file,". ",dupmessage), logname, append=TRUE)
		sm <- sm[!pid %in% dups] # Remove all duplicated instances, to be safe
	} # sm is 5450, less than 5519

	# A bit of QC 
	Trait[index] <<- trait_label
	nSNP[index] <<- nrow(sm)

  	projected.userdata <- tryCatch(expr = project_sparse(beta=sm$BETA, seb=sm$SE, pids=sm$pid)[,trait:=trait_label][], #project_nosig
  	                      error =function(e) {
				      failmessage <- "Projection for this file had non-zero exit status, please check. Jumping to next file..."
				      message(failmessage)
				      write(paste0(file,". ",failmessage), logname, append=TRUE)
				      return(NULL)})
  	if(is.null(projected.userdata)) {
		return(NULL)
	} else{
	projected.userdata[, proj:=NULL] # 105 * 4
  	setnames(projected.userdata, c("var.proj", "delta", "p", "trait"), c("Var.Delta", "Delta", "P", "Trait")) # Changed var.proj to var.delta
  	setcolorder(projected.userdata, c("PC", "Var.Delta", "Delta", "p.overall", "z", "P", "Trait"))

  	# More QC
  	overall_p[index] <<- projected.userdata$p.overall[1] #sprintf("%1.0e",projected.userdata$p.overall[1])
  	minc.idx <- which.min(projected.userdata$P)
  	mscomp[index] <<- sprintf("%s (%1.0e)",projected.userdata$PC[minc.idx],projected.userdata$P[minc.idx])
	}
	projected.userdata
}
) #current

projected.table[sapply(projected.table, is.null)]  <- NULL 
projected.table <- rbindlist(projected.table)
#projected.table  <- projected.table[,.(PC, Delta,Trait)]
projected.table  <- projected.table[,.(PC, Delta, Var.Delta, z, P, Trait)]  #p.overall is in QC.table

# add ratio of covered sparse SNP to QC table as a threshold
QC.table <- data.table(Trait, nSNP, overall_p, mscomp)
QC.table[, ratio_coveredSNP := nSNP / nrow(rot.pca)]


projtablename  <- paste0("Projection_cytokine_basis_sig_40sparse_all", date, ".tsv") # add sig_40sparse_
qctablename  <- paste0("QC_cytokine_basis_sig_40sparse_all", date, ".tsv") # add sig_40sparse_


#while(projtablename %in% dir(paste0(mpath, "cell_basis_v3_varimax/Projections"))){
#  version  <- version + 1
#  projtablename  <- paste0("Projection_cell_basis_v3_", date, "-v",version, ".tsv")
#  qctablename  <- paste0("QC_cell_basis_v3_", date, "-v",version, ".tsv")
#}

message("start writing projected.table & QC.table")

write.table(projected.table, paste0("../Projections/", projtablename), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(QC.table, paste0("../Projections/", qctablename), sep = "\t", quote = FALSE, row.names = FALSE)

message("Finished writing!")


###############################################
# Apply filtering and adjust FDR		     #	
##############################################

############## inspect data, remove NA and Fingen/UKBB  ##################################################

# rename projection results, prepare them for plotting
rename_trait <- function(dt, Metatable){ # return an updated datatable
    meta <- Metatable[, .(Trait_long,Trait)]
    tmp <- merge(dt, meta, by="Trait", all.x=TRUE)
	message("anyNA in Trait_long after mege: ", anyNA(tmp$Trait_long))
	return(tmp)
}

# Load data
Metatable <- fread("../google_metadata_full")[,.(File_ID,Trait,Trait_long,First_Author, Year, Reference ,
Chip,Trait_ID_2.0, Trait_class,N0,N1,N,Population,Collection)]
projected.table <- fread("../Projections/Projection_cytokine_basis_sig_40sparse_all20220719.tsv")
QC.table <- fread("../Projections/QC_cytokine_basis_sig_40sparse_all20220719.tsv")
projected.table <- rename_trait(projected.table, Metatable)

QC.table <- merge.data.table(QC.table, Metatable, by = "Trait",all.x = TRUE)


# Inspect projection results:
anyNA(projected.table)
#[1] FALSE
anyNA(QC.table$overall_p)
#[1] TRUE

#> QC.table[which(is.na(QC.table$overall_p))]
#                       Trait nSNP overall_p mscomp ratio_coveredSNP
#1:     ANS_Cortes_23749187_1  237        NA   <NA>      0.042942562
#2:       MDD_Wray_29700475_2    6        NA   <NA>      0.001087153
#3: PEMFOL_Augusto_33991537_1  184        NA   <NA>      0.033339373
#4:    T2D_Gaulton_26551672_1   56        NA   <NA>      0.010146766
#5:         TRAIL_Reales_up_1   93        NA   <NA>      0.016850879

# These above traits failed to project, we remove them
QC.table <- QC.table[!is.na(QC.table$overall_p),] # 5889 sucessfully projected traits



c(lessthan95 = nrow(QC.table[QC.table$ratio_coveredSNP < 0.95]),
 lessthan80 = nrow(QC.table[QC.table$ratio_coveredSNP < 0.80]), 
 lessthan50 = nrow(QC.table[QC.table$ratio_coveredSNP <0.5]))/nrow(QC.table)

#lessthan95 lessthan80 lessthan50  ## 88.4% traits have less than 80% coverage of sparse SNP, because our UKBB and Fingen was filtered with IMD13 manifest.
# 0.9397181  0.8841909  0.8271353


# For now, UKBB and Fingen are not ready to project yet, so we remove them.
non.UKFIN.QC.table <- QC.table[!Reference %in% c("FinnGenR5","UKBB","PanUKBBR1")]#1250 non UKBB/Fin 
non.UKFIN.proj.table <- projected.table[Trait %in% non.UKFIN.QC.table$Trait] #1250 non UKBB/Fin

c(lessthan95 = nrow(non.UKFIN.QC.table[non.UKFIN.QC.table$ratio_coveredSNP < 0.95]),
 lessthan80 = nrow(non.UKFIN.QC.table[non.UKFIN.QC.table$ratio_coveredSNP < 0.80]), 
 lessthan50 = nrow(non.UKFIN.QC.table[non.UKFIN.QC.table$ratio_coveredSNP <0.5]))/nrow(non.UKFIN.QC.table)
#lessthan95 lessthan80 lessthan50 
# 0.7160     0.4544     0.1856  # 55% non-ukfin traits have coverage over 80%



################  Before FDR: Filtering by 80% SNP match, leaving out ImmunoChip, and removing undesired datasets  ###############
#PT.basis <- copy(non.UKFIN.proj.table)
#PT.basis[,Label:=Trait_long][, stars:=""][, c("Var.Delta", "z", "P") := NULL]

QC.table <- copy(non.UKFIN.QC.table)
proj.table <- copy(non.UKFIN.proj.table)

# Remove Immunochip datasets and those with <80% SNP match  # 
QC.filt <- QC.table[QC.table$ratio_coveredSNP >= 0.8 & !grepl("ImmunoChip|imputation",QC.table$Chip, ignore.case = TRUE), ] # 681 traits left
QC.filt <- QC.filt[ !First_Author %in% c("Inshaw","Cooper")] # 680 left
#table(QC.filt$Trait_class)
#BMK IMD INF OTH PSD 
#581  43  37  10   9 


# Create filtered projection table too
PT.filt <- merge(proj.table, QC.filt[,c("First_Author","Trait", "Trait_ID_2.0", 
"Trait_class", "Collection", "Population")], by = "Trait") # 27280 = 40 * 682 #"Trait_long",

############## Compute FDR  ##################################################
# Apply 1% FDR correction to overall p for all remaining datasets
QC.filt[, FDR.overall := p.adjust(overall_p, method = "BH"), by="Trait_class"] # 682 rows
QC.sig <- QC.filt[FDR.overall < 0.01,] # 248 rows
#table(QC.sig$Trait_class)
#BMK IMD INF PSD 
#225  16   3   2 

# Apply 5% FDR correction by trait and PC to projections, then filter by overall significant traits
PT.filt[, FDR.PC:=p.adjust(P, method = "BH"), by = c("PC","Trait_class")][, stars:=ifelse(!is.na(FDR.PC) & FDR.PC<0.05,"○","")][ FDR.PC < 0.01 , stars:="●"]

# Filter PT80 by overall significant traits
PT.sig <- PT.filt[Trait %in% QC.sig$Trait,] # 246 * 40 sig overall traits, but inside each PC, another layer of sig

# After applying FDR < 1%, only 248 datasets out of 682 were overall significant.
# From those, 18 are IMD, 3 are infectious diseases, 2 are PSD, and 228 are biomarkers.
# The cyto/chemokine basis seems to recover a mainly immune signal
# Note that many BMK in external projections are reflecting same proteins(cyto/chemokines) as our basis, except from other studies.



#######################################################################################
# QC plot of binary traits (to be done after we get complete UKBB and FinnGen)      #	
######################################################################################


#######################################################################################
#     Visualizing projections
######################################################################################

make.labels.cytochemo <- function(x){
  x[,Label:=Trait_long] %>% 
    .[, Label:=gsub(" \\(UKBB\\)", "", Label)] %>% 
    .[, Label:=gsub(" \\(FinnGen\\)", "", Label)] %>% 
    .[, Label:=gsub(" \\(FG\\)", "", Label)] %>% 
    .[, Label:=gsub("Crohn disease \\(strict definition, all UC cases excluded\\)", "Crohn's, strict", Label)] %>% 
    .[, Label:=gsub("Eosinophilic Granulomatosis with Polyangiitis", "EGPA", Label)] %>% 
    .[ First_Author == "FinnGen", Label:=paste0(str_trunc(Label, width = 50), " / ", First_Author)] %>% 
    .[ First_Author != "FinnGen", Label:=paste0(str_trunc(paste0(Label, " (", Population,")"), width = 50), " / ", First_Author)] %>% 
    .[, Label:=gsub("Neale", "UKBB", Label)] %>% 
    .[, Label:=gsub("PanUKBB", "UKBB", Label)] # Addition to accommodate PanUKBB
  
}


# Remove unwanted trait : As "HBC_Chen_32888493_6" has extreme Delta for some reason.
# Remove chips!!! from our collection.

unwanted <- c("HBC_Chen_32888493_6")
PT.sigviz <- PT.sig[!grepl(unwanted, Trait),] #245 traits

PT.sigviz <- make.labels.cytochemo(PT.sigviz)
QC.sigviz <- QC.sig[Trait %in% unique(PT.sigviz$Trait),]
nrow(QC.sigviz)
#245


# Check duplicated labels, as in , different traits sharing same trait_long + first author
# such as ""
#which(M.sig[,"PC1"]>1)
#PT.sigviz[Label == "Monocyte Count (Trans-ethnic (EUR, AFR, EAS, SAS, HA)) / Chen"][PC == "PC1"]
#PT.sigviz[Label == "C-X-C motif chemokine 5 (CXCL5, ENA78) levels (European (... / Hoglund"][PC == "PC1"]
#duplicate.labels <- PT.sigviz$Label[duplicated(PT.sigviz[PC=="PC2"]$Label)]
PT.sigviz[,Label:=paste0(Label,"_",substr(Trait,nchar(Trait), nchar(Trait)))]



#####################################     Plot  Complete Pheamap   ############################################################
hmcol <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(100)
PCorder <-  paste0("PC",1:40)

M.sig <- acast(PT.sigviz[,c("PC", "Label", "Delta")], Label ~ PC) # Melted to 141 * 40 matrix. 141 is the length of unique label

M.sig.stars <- acast(PT.sigviz[,c("PC","Label","stars")], Label ~ PC)

M.sig <- M.sig[,PCorder]
M.sig.stars <- M.sig.stars[,PCorder]

range <- max(abs(M.sig))
# [1] 27686165095 ?!
# [1] 0.07372398

#View(PT.sigviz)
# found some very extreme Delta in Chen datasets, why?
# exclue "HBC_Chen_32888493_6"!

# Let's also annotate by type of trait
annviz <- data.frame(Trait_class = QC.sigviz$Trait_class, row.names = unique(PT.sigviz$Label)) # QC.sigviz doesn't have "Label" column, but traits are in the same order
colannviz = list(Trait_class = c(BMK = "#06D6A0", IMD = "#EF476F", INF = "#FFD166", PSD = "#F9DAD0"))

png(file="../Plots/imd_40_sparseAll/pheamap.png", width = 1000, height = 5000)
pheatmap(M.sig, breaks = seq(-range, range, length.out = 100), cluster_cols = FALSE, display_numbers = M.sig.stars, fontsize_row = 12, fontsize_number = 20, color = rev(hmcol), annotation_row = annviz, annotation_colors = colannviz)
dev.off()
#fontsize_number: change circle diameter

##########################################  make pheamap of IMD only  ####################################################

hmcol <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(100)
PCorder <-  paste0("PC",1:40)

M.sig.imd <- acast(PT.sigviz[Trait_class =="IMD"][,c("PC", "Label", "Delta")], Label ~ PC) # Melted to 141 * 40 matrix. 141 is the length of unique label

M.sig.stars.imd <- acast(PT.sigviz[Trait_class =="IMD"][,c("PC","Label","stars")], Label ~ PC)

M.sig.imd <- M.sig.imd[,PCorder]
M.sig.stars.imd <- M.sig.stars.imd[,PCorder]

range <- max(abs(M.sig.imd))
# [1] 0.01066709 max abs delta

# Let's also annotate by type of trait
annviz <- data.frame(Trait_class = QC.sigviz[Trait_class =="IMD"]$Trait_class, row.names = unique(PT.sigviz[Trait_class =="IMD"]$Label)) # QC.sigviz doesn't have "Label" column, but traits are in the same order
colannviz = list(Trait_class = c(IMD = "#EF476F"))


png(file="../Plots/imd_40_sparseAll/pheamap_imd.png", width = 1000, height = 500)
pheatmap(M.sig.imd, breaks = seq(-range, range, length.out = 100), cluster_cols = FALSE, display_numbers = M.sig.stars.imd, fontsize_row = 12, fontsize_number = 20, color = rev(hmcol), annotation_row = annviz, annotation_colors = colannviz)
dev.off()


##########################################  make pheamap of Non-BMK projections only  ####################################################

hmcol <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(100)
PCorder <-  paste0("PC",1:40)

M.sig.nonBMK <- acast(PT.sigviz[Trait_class !="BMK"][,c("PC", "Label", "Delta")], Label ~ PC) # Melted to 141 * 40 matrix. 141 is the length of unique label

M.sig.stars.nonBMK <- acast(PT.sigviz[Trait_class !="BMK"][,c("PC","Label","stars")], Label ~ PC)

M.sig.nonBMK <- M.sig.nonBMK[,PCorder]
M.sig.stars.nonBMK <- M.sig.stars.nonBMK[,PCorder]

range <- max(abs(M.sig.nonBMK))
# [1] 0.04303982 max abs delta

#View(PT.sigviz)
# found some very extreme Delta in Chen datasets, why?
# exclue "HBC_Chen_32888493_6"!

# Let's also annotate by type of trait
annviz <- data.frame(Trait_class = QC.sigviz[Trait_class !="BMK"]$Trait_class, row.names = unique(PT.sigviz[Trait_class !="BMK"]$Label)) # QC.sigviz doesn't have "Label" column, but traits are in the same order
colannviz = list(Trait_class = c(IMD = "#EF476F", INF = "#FFD166", PSD = "#F9DAD0"))

png(file="../Plots/imd_40_sparseAll/pheamap_nonBMK.png", width = 1000, height = 600)
pheatmap(M.sig.nonBMK, breaks = seq(-range, range, length.out = 100), cluster_cols = FALSE, display_numbers = M.sig.stars.nonBMK, fontsize_row = 12, fontsize_number = 20, color = rev(hmcol), annotation_row = annviz, annotation_colors = colannviz)
dev.off()

##########################################  make pheamap of external Cyto/chemokine projections only  ####################################################

hmcol <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(100)
PCorder <-  paste0("PC",1:40)

PT.sigviz.bmk <- PT.sigviz[Trait_class == "BMK"]
QC.sigviz.bmk <- QC.sigviz[Trait_class == "BMK"]

# Find out the cyto/chemokine traits_long
ahola.traits <- c("CTACK","EOT1","GCSF","GROA","IFNG","IL10","IL12",
"IL13","IL16","IL17","IL18","IL1B","IL1RA","IL2","IL4",
"IL5","IL6","IL7","IL8","IL9","IP10","MCP1","MCP3","MCSF","MIG","MIP1A","MIP1B",
"RANTES","SCF","SDF1a","TNFa","TNFb","TRAIL")

folkersen.traits <- c("IL8","CD40L","IL1RA","IL6","MCP1","TRAIL","MCSF","SCF","IL18",
"TNFSF14","TRANCE","IL16","CXCL6","CX3CL1","CCL20")

hoglund.traits <- unique(PT.sigviz.bmk[First_Author == "Hoglund"]$Trait_ID_2.0)
hoglund.traits <- hoglund.traits[!hoglund.traits %in% c("4EBP1","ADA","AXIN1","bNGF","CASP8","CD244","CD40","CD5","CD6","CDCP1","CST5", "DNER","ENRAGE",
"FGF19","FGF21","FGF23","FGF5","FLT3L","GDNF","HGF","IL10RB","IL15RA","IL18R1","IL20RA","IL22RA1","IL2RB",
"MMP1","MMP10","NRTN","NT3","OPG","PDL1","SIRT2","SLAMF1","STAMBP","SULT1A1","UPA","VEGFA")]

safe.ID2 <- union(union(ahola.traits, folkersen.traits), hoglund.traits)

PT.sigviz.cyto <- PT.sigviz.bmk[First_Author %in% c("Hoglund","Folkersen","AholaOlli")][Trait_ID_2.0 %in% safe.ID2]
QC.sigviz.cyto <- QC.sigviz.bmk[First_Author %in% c("Hoglund","Folkersen","AholaOlli")][Trait_ID_2.0 %in% safe.ID2]

# plot external cyto/chemokines

M.sig.cyto <- acast(PT.sigviz.cyto[,c("PC", "Label", "Delta")], Label ~ PC) # Melted to 141 * 40 matrix. 141 is the length of unique label

M.sig.stars.cyto <- acast(PT.sigviz.cyto[,c("PC","Label","stars")], Label ~ PC)

M.sig.cyto <- M.sig.cyto[,PCorder]
M.sig.stars.cyto <- M.sig.stars.cyto[,PCorder]

range <- max(abs(M.sig.cyto))
# [1] 0.03464276 max abs delta

# Let's also annotate by type of trait
annviz <- data.frame(Trait_class = QC.sigviz.cyto$Trait_class, row.names = unique(PT.sigviz.cyto$Label)) # QC.sigviz doesn't have "Label" column, but traits are in the same order
colannviz = list(Trait_class = c(BMK = "#06D6A0"))


png(file="../Plots/imd_40_sparseAll/pheamap_cyto.png", width = 1000, height = 600)
pheatmap(M.sig.cyto, breaks = seq(-range, range, length.out = 100), cluster_cols = FALSE, display_numbers = M.sig.stars.cyto, fontsize_row = 12, fontsize_number = 20, color = rev(hmcol), annotation_row = annviz, annotation_colors = colannviz)
dev.off()


###############################################
# Pie chart to show the proportion change before/after FDR adjust #	
##############################################

#BMK IMD INF OTH PSD  (after 80% SNP search)
#581  43  37  10   9 

#BMK IMD INF PSD (after  overall p value adjusted for FDR)
#225  16   3   2

pie <- data.frame(value=c(581,43,37,9,10), type=c("BMK","IMD","INF","PSD","OTH"))
pie2 <- data.frame(value=c(225,16,3,2), type=c("BMK","IMD","INF","PSD"))

plot.pie <- ggplot(pie, aes(x = "", y = value, fill = type)) +
  geom_col()+
  geom_col(color = "black") +
    geom_text(aes(label = value),
            position = position_stack(vjust = 0.3)) +
    coord_polar(theta = "y") +
	scale_fill_brewer(palette="Set3") +
	labs(title = "Count of traits of different classes - Before FDR adjust")+
	ylab(paste0("Number of projected traits ( imd counts for ", round(100* pie$value[which(pie$type == "IMD")] / sum(pie$value),3),"% )"))

plot.pie2 <- ggplot(pie2, aes(x = "", y = value, fill = type)) +
  geom_col()+
  geom_col(color = "black") +
    geom_text(aes(label = value),
            position = position_stack(vjust = 0.3)) +
    coord_polar(theta = "y") +
	scale_fill_brewer(palette="Set3") +
	labs(title = "Count of traits of different classes - After FDR filter")+
	ylab(paste0("Number of projected traits ( imd counts for ", round(100* pie2$value[which(pie2$type == "IMD")] / sum(pie2$value),3),"% )"))

# save pie charts

ggsave("../Plots/imd_40_sparseAll/traitPie_beforeFDR.png",plot.pie, width = 8,height = 8)
ggsave("../Plots/imd_40_sparseAll/traitPie_afterFDR.png",plot.pie2, width = 8, height = 8)




###############################################
# Forrest plot:  Load previous projection of basis to basis,				  #	
##############################################
PT.delta <- PT.sigviz[stars != "",] #3746 significant PC traits to plot

# Read basis-basis projection
PT.basis <- fread("../Projections/Projection_cytokine_basis_sig_40sparse_basisTObasis20220719.tsv")

#
PT.basis[, Trait_class:="Ferkingstad"]
PT.delta <- rbind(PT.basis, PT.delta, fill = TRUE)
PT.delta[,Var.Delta:=as.numeric(Var.Delta)]
PT.delta[Var.Delta != 0, ci:=sqrt(Var.Delta) * 1.96][is.na(ci),ci:=0][is.na(Label), Label:=Trait] #label 
PT.delta[Trait_class=="Ferkingstad", ci:=0]
PCs <- unique(PT.delta$PC)


# prepare dt (could be a subset of entire projection, e.g. "dt <- all.dt[Trait_class %in% c("IMD","Ferkingstad")]")

forest.plot <- function(x="PC1", PT.delta, show.class=c("IMD","Ferkingstad")){

all.dt <- PT.delta[PC == x,][order(Delta, decreasing = FALSE),][,colours:=ifelse(ci == 0, "red", ifelse(Trait_class == "IMD", "#26547c", "black"))]
dt <- all.dt[Trait_class %in% show.class]

fplot <- ggplot(dt, aes(x = reorder(Label, -Delta), y = Delta, ymin=Delta-ci, ymax=Delta+ci, colour = colours))+
  geom_pointrange()+
  scale_colour_manual(values = c("red" = "red", "black" = "black", "#26547c" = "#26547c"))+
  geom_hline(yintercept = 0, col="red", lty=2)+
  coord_flip()+
  xlab("Traits")+
  ylab("Delta")+
  ggtitle(paste("Delta Plot ", x, sep = ""))+
  theme_light()+
  theme(axis.text.y = element_text(colour = rev(dt$colours)), legend.position = "none")
  fplot
}

# forest.plot(PT.delta=PT.delta, x="PC5")


# save forrest plots

for (i in 1:length(PCs)){
  tmp <- forest.plot(x=PCs[i], PT.delta = PT.delta, show.class=c("IMD","Ferkingstad"))
  ggsave(filename = paste0("../Plots/imd_40_sparseAll/","IMD_FK_Forest_PC",i,".png"),plot = tmp)
}


# Plots to make:
#  1) forrest plot : Ferkingstad + IMDs 
#  2) forrest plot : Ferkingstad + nonBMKs
#  3) forrest plot : Ferkingstad + cytochemos
#  4) heatmap : Ferkingstad + IMDs ; IMDs
#  5) heatmap : Ferkingstad + nonBMKs ; nonBMKs
#  6) heatmap : Ferkingstad + cytochemos : cytochemos
