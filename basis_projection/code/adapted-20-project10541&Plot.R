# This script is to apply projection of 10541 external datasets from 02-processed/ (reduced by adapted-18-reduction.R) onto
# our Ferkingstad basis (sparse, PC1-40) for a complete projection.  Significance is computed here.

# In the plotting part, Previously projected Ferkingstad basis (from adapted-16-projection-FK40sparse) will be plotted together with
# our external datasets, but only those with : (A) QC passed 80% of coverage of sparse SNP
# (B) Overall significant, and sig in each PC
# will be shown in our plot

# 2022-07-27

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
files_to_project  <- dir("../data_40_sparse_allcollection10541", pattern = ".tsv") # 10541 datasets to project
nfiles <- length(files_to_project)
Trait <- rep(NA,nfiles)
nSNP <- rep(NA, nfiles)
overall_p <- rep(NA, nfiles)
mscomp <- rep(NA,nfiles)


# ############## test projection of "ANS_Cortes_23749187_1-ft.tsv"#############
# file <- "ANS_Cortes_23749187_1-ft.tsv" #237 SNPs
# file <- "MDD_Wray_29700475_2-ft.tsv" #6 SNPs
# file <- "PEMFOL_Augusto_33991537_1-ft.tsv" #185 SNPs
# file <- "T2D_Gaulton_26551672_1-ft.tsv" #56 SNPs
# file <- "O15_BREAST_LACT_OTHER_DIS_FinnGen_FinnGenR5_1-ft.tsv" #222 SNPs, but works in project_sparse

# trait_label <- strsplit(file, "-", fixed = TRUE)[[1]][1] #trait name
# sm <- fread(paste0("../data_40_sparse_allcollection10541/",file)) #data_40_sparse_allcollection10541

# 	sm <- unique(sm)
# 	sm[sm == ""] <- NA # Some missing data might pass as empty string. This will fix that	
# 	sm <- na.omit(sm, cols = c("pid", "BETA", "SE", "P")) # 
# 	dups <- sm$pid[duplicated(sm$pid)]
# 	sm <- sm[!BETA %in% c(Inf, -Inf)][!SE %in% c(Inf, -Inf)]


# x <- project_sparse(beta=sm$BETA, seb=sm$SE, pids=sm$pid)[,trait:=trait_label]

# beta=sm$BETA
# seb=sm$SE
# pids=sm$pid

# project_sparse <- function (beta=sm$BETA, seb=sm$SE, pids=sm$pid) {
#     if (length(beta) != length(seb) || length(beta) != length(pids) || 
#         !length(beta)) 
#         stop("arguments must be equal length vectors > 0")
#     if (!all(pids %in% SNP.manifest$pid)) 
#         stop("all pids must be members of sparse basis (SNP.manifest$pid)")
#     if (length(pids) < 0.95 * nrow(rot.pca)) 
#         warning("more than 5% sparse basis snps missing")
#     b <- beta * shrinkage[pids] - beta.centers[pids]
#     proj <- b %*% rot.pca[pids, ]
#     # ATTENTION: Before running this function, remove Inf/-Inf from beta and SE!!
#     # ATTENTION: seb needs to be NA free! otherwise will introduce NA into the matrix computation
#     v <- seb * shrinkage[pids] * rot.pca[pids, ] # v is not variance, it is an intermediate term, like a projection of seb, parallel to b and proj
#     var.proj <- t(v) %*% LD[pids, pids] %*% v # why NA?
#     ctl <- (-beta.centers[pids]) %*% rot.pca[pids, ]
#     delta <- (proj - ctl)[1, ]
#     chi2 <- (t(delta) %*% solve(var.proj) %*% delta)[1, 1]
#     ret <- data.table::data.table(PC = colnames(proj), proj = proj[1, 
#         ], var.proj = Matrix::diag(var.proj), delta = delta, 
#         p.overall = stats::pchisq(chi2, df = 40, lower.tail = FALSE)) # NOTE: remember to update df as number of PCs!
#     ret$z = ret$delta/sqrt(ret$var.proj)
#     ret$p = stats::pnorm(abs(ret$z), lower.tail = FALSE) * 2
#     copy(ret)
# }




###############################################
# Projection of 1-10541 external traits		#	
##############################################

projected.table <- lapply(files_to_project, function(file){
	message("Projecting ", file)
	trait_label <- strsplit(file, "-", fixed = TRUE)[[1]][1] #trait name
	index <- which(files_to_project == file) 
	sm <- fread(paste0("../data_40_sparse_allcollection10541/",file)) #data_40_sparse_allcollection10541

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


projtablename  <- paste0("Projection_cytokine_basis_sig_40sparse_all10541", date, ".tsv") # add sig_40sparse_
qctablename  <- paste0("QC_cytokine_basis_sig_40sparse_all10541", date, ".tsv") # add sig_40sparse_


message("start writing projected.table & QC.table")

write.table(projected.table, paste0("../Projections/", projtablename), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(QC.table, paste0("../Projections/", qctablename), sep = "\t", quote = FALSE, row.names = FALSE)

message("Finished writing!")



###############################################
# Apply filtering and adjust FDR		     #	
##############################################

############## inspect data, remove NA and Fingen/UKBB  ##################################################

# Load data
Metatable <- fread("../google_metadata_fullUKBBFIN")[,.(File_ID,Trait,Trait_long,First_Author, Year, Reference ,
Chip,Trait_ID_2.0, Trait_class,N0,N1,N,Population,Collection,usedinBasis,Common_Name,isCyto)]  # add "usedinBasis","Common_Name","isCyto"
projected.table <- fread("../Projections/Projection_cytokine_basis_sig_40sparse_all1054120220727.tsv")
QC.table <- fread("../Projections/QC_cytokine_basis_sig_40sparse_all1054120220727.tsv")

projected.table <- merge.data.table(projected.table, Metatable, by = "Trait",all.x = TRUE)
QC.table <- merge.data.table(QC.table, Metatable, by = "Trait",all.x = TRUE)


# Inspect projection results:
anyNA(projected.table$Delta)
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
QC.table <- QC.table[!is.na(QC.table$overall_p),] # 10536 sucessfully projected traits



c(lessthan95 = nrow(QC.table[QC.table$ratio_coveredSNP < 0.95]),
 lessthan80 = nrow(QC.table[QC.table$ratio_coveredSNP < 0.80]), 
 lessthan50 = nrow(QC.table[QC.table$ratio_coveredSNP <0.5]))/nrow(QC.table)

#lessthan95 lessthan80 lessthan50  ## 46.73% traits have less than 80% coverage of sparse SNP, because some of our old UKBB and Fingen was filtered with IMD13 manifest.
# 0.4984814  0.4673500  0.4354594


# UKBBR5 is the old UKBB datasets which were pre-filtered with IMD-13 manifest, so we remove them
non.UKFIN.QC.table <- QC.table[!Reference %in% c("FinnGenR5","PanUKBBR1")]#6180 non UKBB/Fin #PanUKBBR1 is old version, now we have PanUKBBR2; Similarly, we have FinnGenR7 instead of FinnGenR5
non.UKFIN.proj.table <- projected.table[Trait %in% non.UKFIN.QC.table$Trait] #6180 non UKBB/Fin

c(lessthan95 = nrow(non.UKFIN.QC.table[non.UKFIN.QC.table$ratio_coveredSNP < 0.95]),
 lessthan80 = nrow(non.UKFIN.QC.table[non.UKFIN.QC.table$ratio_coveredSNP < 0.80]), 
 lessthan50 = nrow(non.UKFIN.QC.table[non.UKFIN.QC.table$ratio_coveredSNP <0.5]))/nrow(non.UKFIN.QC.table)
#lessthan95 lessthan80 lessthan50 
#0.14498382 0.09190939 0.03754045 # 90.8% non-ukfin traits have coverage over 80%




################  Before FDR: Filtering by 80% SNP match, leaving out ImmunoChip, and removing undesired datasets  ###############
#PT.basis <- copy(non.UKFIN.proj.table)
#PT.basis[,Label:=Trait_long][, stars:=""][, c("Var.Delta", "z", "P") := NULL]

QC.table <- copy(non.UKFIN.QC.table)
proj.table <- copy(non.UKFIN.proj.table)

# Remove Immunochip datasets and those with <80% SNP match  # 
QC.filt <- QC.table[QC.table$ratio_coveredSNP >= 0.8 & !grepl("ImmunoChip|imputation",QC.table$Chip, ignore.case = TRUE), ] # 5611 traits left
QC.filt <- QC.filt[ !First_Author %in% c("Inshaw","Cooper")] # 5610 left
#table(QC.filt$Trait_class)
#BMK  CAN  IMD  INF  OTH  PSD 
#646  598  481  373 3147  365 


# Create filtered projection table too
PT.filt <- proj.table[Trait %in% QC.filt$Trait] # 224400 = 40 * 5610 #"Trait_long" #merge(proj.table, QC.filt[,c("First_Author","Trait", "Trait_ID_2.0", 
#"Trait_class", "Collection", "Population")], by = "Trait")


############## Compute FDR  ##################################################

# Apply 1% FDR correction to overall p for all remaining datasets
QC.filt[, FDR.overall := p.adjust(overall_p, method = "BH"), by="Trait_class"] # 5610 rows
QC.sig <- QC.filt[FDR.overall < 0.01,] # 544 rows
#table(QC.sig$Trait_class)
#BMK IMD INF OTH PSD 
#284  79  16 155  10  

# Apply 5% FDR correction by trait and PC to projections, then filter by overall significant traits
PT.filt[, FDR.PC:=p.adjust(P, method = "BH"), by = c("PC","Trait_class")][P < 0.05, stars:="﹡"]

PT.filt[FDR.PC<0.05, stars:="○"][ FDR.PC < 0.01 , stars:="●"] # P value, <0.05, less strict

PT.filt[is.na(stars), stars:=""] # IF no "*" : PT.filt[, FDR.PC:=p.adjust(P, method = "BH"), by = c("PC","Trait_class")][, stars:=ifelse(!is.na(FDR.PC) & FDR.PC<0.05,"○","")][ FDR.PC < 0.01 , stars:="●"]


# Filter PT80 by overall significant traits
PT.sig <- PT.filt[Trait %in% QC.sig$Trait,] # 544 * 40 sig overall traits, but inside each PC, another layer of sig

# After applying FDR < 1%, only 544 datasets out of 5610 were overall significant.
# From those, 79 are IMD, 16 are infectious diseases, 10 are PSD, 155 are OTH(mostly removed after filtering), and 284 are biomarkers.
# The cyto/chemokine basis seems to recover a mainly immune signal
# Note that many BMK in external projections are reflecting same proteins(cyto/chemokines) as our basis, except from other studies.


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
    .[, Label:=gsub("PanUKBB", "UKBB", Label)] %>% # Addition to accommodate PanUKBB
    .[isCyto==TRUE,Label:= paste0(Common_Name, " / ", First_Author )] %>% # for cyto/cheomos, label them as "Common_Name + Author"
    .[,Label:=paste0(Label,"_",substr(Trait,nchar(Trait), nchar(Trait)))] %>% 
    .[Trait=="I9_TIA_FinnGen_FinnGenR7_1", Label := paste0(Label,"(I9)")] %>% #duplicated label
    .[Trait=="K11_UC_STRICT2_FinnGen_FinnGenR7_1", Label := paste0(Label, "(2)")] %>% #duplicated label
    .[Trait=="I9_DISVEINLYMPH_FinnGen_FinnGenR7_1", Label := paste0(Label, "(I9)")]
}


# Metatable

# Remove unwanted trait : As "HBC_Chen_32888493_6" has extreme Delta for some reason.
# Remove chips!!! from our collection.

unwanted <- c("HBC_Chen_32888493_6")
PT.sigviz <- PT.sig[!grepl(unwanted, Trait),] # 543 traits

PT.sigviz <- make.labels.cytochemo(PT.sigviz)
QC.sigviz <- QC.sig[Trait %in% unique(PT.sigviz$Trait),]
nrow(QC.sigviz)
#543

# Check duplicated labels, as in , different traits sharing same trait_long + first author
any(duplicated(PT.sigviz[PC=="PC1"]$Label))
# FALSE

# add "Label" to QC.sigviz, for convenience of phemap plotting
QC.sigviz <- merge.data.table(QC.sigviz, PT.sigviz[PC=="PC1"][,.(Trait, Label)], by="Trait", all.x = TRUE)

# Plots to make:
#  1) forrest plot : Ferkingstad + IMDs 
#  2) forrest plot : Ferkingstad + nonBMKs
#  3) forrest plot : Ferkingstad + cytochemos
#  4) heatmap : Ferkingstad + IMDs ; IMDs
#  5) heatmap : Ferkingstad + nonBMKs ; nonBMKs
#  6) heatmap : Ferkingstad + cytochemos : cytochemos


###############################
#####      Phemaps         ####
###############################

#####################################     Plot  Complete Pheamap   ############################################################

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


# make pheamap  : IMDs
plot.PT <- PT.sigviz[Trait_class=="IMD"]
plot.QC <- QC.sigviz[Trait_class=="IMD"]

#png(file="../Plots/imd_40_sparseAll10541/pheatmap_IMD.png", width = 10, height = 15)
#pheamap.plot(PT.sigviz = plot.PT, QC.sigviz = plot.QC)
#dev.off()

IMDheat <- pheamap.plot(PT.sigviz = plot.PT, QC.sigviz = plot.QC)
ggsave("../Plots/imd_40_sparseAll10541/pheatmap_IMD.png",IMDheat, width=13, height = 15)


# make pheamap  : cytochemos
plot.PT <- PT.sigviz[isCyto==TRUE]
plot.QC <- QC.sigviz[isCyto==TRUE]

png(file="../Plots/imd_40_sparseAll10541/pheatmap_CYTO.png", width = 1000, height = 800)
pheamap.plot(PT.sigviz = plot.PT, QC.sigviz = plot.QC)
dev.off()


# make pheamap  : no Reales,
plot.PT <- PT.sigviz[isCyto==TRUE][First_Author!="Reales"]
plot.QC <- QC.sigviz[isCyto==TRUE][First_Author!="Reales"]

png(file="../Plots/imd_40_sparseAll10541/pheatmap_CYTO_noReales_omi.png", width = 1000, height = 600)
pheamap.plot(PT.sigviz = plot.PT, QC.sigviz = plot.QC)
dev.off()


# make pheamap  : no Reales, only showing cyto/chemokines used in basis
plot.PT <- PT.sigviz[isCyto==TRUE][First_Author!="Reales"][usedinBasis==TRUE]
plot.QC <- QC.sigviz[isCyto==TRUE][First_Author!="Reales"][usedinBasis==TRUE]

        # save plotdata
fwrite(plot.PT, file="../plotdata/cytoinBasis_PT", sep = "\t", na = NA)
fwrite(plot.QC, file="../plotdata/cytoinBasis_QC", sep = "\t", na = NA)


png(file="../Plots/imd_40_sparseAll10541/pheatmap_CYTO_noReales_omi_onlybasis.png", width = 1000, height = 600)
pheamap.plot(PT.sigviz = plot.PT, QC.sigviz = plot.QC)
dev.off()

# make pheamap  : no Reales, no Hoglund

plot.PT <- PT.sigviz[isCyto==TRUE][!First_Author %in% c("Reales","Hoglund")]
plot.QC <- QC.sigviz[isCyto==TRUE][!First_Author %in% c("Reales","Hoglund")]

png(file="../Plots/imd_40_sparseAll10541/pheatmap_CYTO_noReales_noHug.png", width = 1000, height = 400)
pheamap.plot(PT.sigviz = plot.PT, QC.sigviz = plot.QC)
dev.off()

# make pheamap  : non BMKs
plot.PT <- PT.sigviz[Trait_class!="BMK"]
plot.QC <- QC.sigviz[Trait_class!="BMK"]

png(file="../Plots/imd_40_sparseAll10541/pheatmap_nonBMK.png", width = 1000, height = 3300)
pheamap.plot(PT.sigviz = plot.PT, QC.sigviz = plot.QC)
dev.off()



###########################################################################
###      Forrest plot:  Load previous projection of basis to basis 		 ##	
###########################################################################

PT.delta <- PT.sigviz[stars !="",] #3746 significant PC traits to plot

# Read basis-basis projection
PT.basis <- fread("../Projections/Projection_cytokine_basis_sig_40sparse_basisTObasis20220719.tsv")
meta.fk <- fread("../google_metadata_ferkingstad")

PT.basis <- merge.data.table(PT.basis, meta.fk, by="Trait", all.x=TRUE) %>%
.[, Label:=paste0(Common_Name, " / ", "Ferkingstad")] %>%
.[, Trait_class:="Ferkingstad"]

PT.delta <- rbind(PT.basis, PT.delta, fill = TRUE)
PT.delta[,Var.Delta:=as.numeric(Var.Delta)]
PT.delta[Var.Delta != 0, ci:=sqrt(Var.Delta) * 1.96][is.na(ci),ci:=0][is.na(Label), Label:=Trait] #label 
PT.delta[Trait_class=="Ferkingstad", ci:=0]
PCs <- unique(PT.delta$PC)


# prepare dt (could be a subset of entire projection, e.g. "dt <- all.dt[Trait_class %in% c("IMD","Ferkingstad")]")

forest.plot <- function(x="PC1", PT.delta, show.class=c("IMD","Ferkingstad")){

dt <- PT.delta[PC == x,][Trait_class %in% show.class][order(Delta, decreasing = FALSE),][,colours:=ifelse(ci == 0, "red", ifelse(Trait_class == "IMD", "#26547c", "black"))]
dt$Label <- factor(dt$Label, levels = unique(dt$Label))
#reorder(Label, -Delta)

mycol <- copy(dt[-which(duplicated(dt$Label))]$colours)

fplot <- ggplot(dt, aes(x = Label, y = Delta, ymin=Delta-ci, ymax=Delta+ci, colour = colours))+
  geom_pointrange()+
  scale_colour_manual(values = c("red" = "red", "black" = "black", "#26547c" = "#26547c"))+
  geom_hline(yintercept = 0, col="red", lty=2)+
  coord_flip()+
  xlab("Label")+
  ylab("Delta")+
  ggtitle(paste("Delta Plot ", x, sep = ""))+
  theme_light()+
  theme(axis.text.y = element_text(colour = mycol), legend.position = "none")
  fplot
}

# forest.plot(PT.delta=PT.delta, x="PC11")


# forest plots : IMDs + Ferkingstad

for (i in 1:length(PCs)){
  tmp <- forest.plot(x=PCs[i], PT.delta = PT.delta, show.class=c("IMD","Ferkingstad"))
  ggsave(filename = paste0("../Plots/imd_40_sparseAll10541/","IMD_FK_Forest_PC",i,".png"), width=12, height=10,plot = tmp)
}

i=11
tmp <- forest.plot(x=PCs[i], PT.delta = PT.delta, show.class=c("IMD","Ferkingstad"))
  ggsave(filename = paste0("../Plots/imd_40_sparseAll10541/","IMD_FK_Forest_PC",i,".png"),width=12, height=10, plot = tmp)

# forest plots : External cytokines + Ferkingstad
for (i in 1:length(PCs)){
tmp <- forest.plot(x=PCs[i], PT.delta = PT.delta[usedinBasis==TRUE], show.class = c("BMK","Ferkingstad"))
  ggsave(filename = paste0("../Plots/imd_40_sparseAll10541/","CYTO_FK_Forest_PC",i,".png"), width=12, height=10,plot = tmp)
}




###########################################################################
###      Pie plot: Before and After p.adjust 		 ##	
###########################################################################

#BMK  CAN  IMD  INF  OTH  PSD (after 80% SNP search)
#646  598  481  373 3147  365 

#BMK IMD INF OTH PSD (after  overall p value adjusted for FDR)
#284  79  16 155  10  


pie <- data.frame(value=c(481,646,373,3147,365,598), Type=c("IMD","BMK","INF","OTH","PSD","CAN"))
pie2 <- data.frame(value=c(79,284,16,155,10,0), Type=c("IMD","BMK","INF","OTH","PSD","CAN"))

plot.pie <- ggplot(pie, aes(x = "", y = value, fill = Type)) +
  geom_col()+
  geom_col(color = "black") +
    geom_text(aes(label = value),
            position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
	scale_fill_brewer(palette="Set3") +
	labs(title = "Count of traits of different classes - Before FDR adjust")+
	ylab(paste0("Number of projected traits ( imd counts for ", round(100* pie$value[which(pie$Type == "IMD")] / sum(pie$value),3),"% )"))

plot.pie2 <- ggplot(pie2, aes(x = "", y = value, fill = Type)) +
  geom_col()+
  geom_col(color = "black") +
    geom_text(aes(x=1.2, label = value),
            position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
	scale_fill_brewer(palette="Set3") +
	labs(title = "Count of traits of different classes - After FDR filter")+
	ylab(paste0("Number of projected traits ( imd counts for ", round(100* pie2$value[which(pie2$Type == "IMD")] / sum(pie2$value),3),"% )"))

# save pie charts

ggsave("../Plots/imd_40_sparseAll10541/traitPie_beforeFDR.png",plot.pie, width = 8,height = 8)
ggsave("../Plots/imd_40_sparseAll10541/traitPie_afterFDR.png",plot.pie2, width = 8, height = 8)


# print diseases
imd7 <- unique(PT.delta[PC %in% c("PC7")][Trait_class=="IMD"]$Label)


