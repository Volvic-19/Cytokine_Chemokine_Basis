### Visualising projection results across parameters

# Guillermo Reales
# 2022-06-20

# In this script we'll use basis + projection of Chen (Eur) and Astle to evaluate the different parameters

library(data.table)
setDTthreads(8)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(reshape2)

projection.stats <- function(p.table){

	message("Creating summary statistics from projection tables...")
	var_snps <- unique(p.table[,.(PC, Param, Var.exp, SNP.PC, SNP.total)])
	d1 <- p.table[, .(sqdiff = diff(Delta)^2), by=c("Trait2", "PC", "Param")][, .(meansqdiff.Delta = mean(sqdiff)), by=c("PC", "Param")]
	d2 <- p.table[, .(mean.Delta = mean(abs(Delta)), sd.Delta = sd(abs(Delta))), by=c("PC","Param")]
	diff.table <- merge(d1, d2, by=c("PC", "Param"))
	diff.table <- merge(diff.table, var_snps, by=c("PC", "Param"))

}

p.table <- fread("../data/projtable_8t_1e5s_9to5.tsv")
d.table <- projection.stats(p.table) 

summary(d.table) 
  
m=melt(d.table, id.vars=c("PC","Param"))
head(m)
m[,scaled_value:=(value - min(value))/(max(value) - min(value)), by=c("PC","variable")][, PC:= factor(PC, levels=paste0("PC",1:15))]
  
pt <- ggplot(m[variable!="nSNPs"], aes(x=Param, y=scaled_value, col=variable))  +
    geom_path() + geom_point() +
    facet_wrap(~PC,nrow=3) +
    theme_cowplot(10)+
    scale_x_log10()+
    background_grid(major="xy", minor="x")+
    theme(legend.title=element_blank())
pt

ggsave("../figures/finetuning_1e5_9to5_8t.png", pt, height=1500, width=2200, units="px",  bg="white")


# Heatmaps

ph <- p.table[,.(PC, Param, Trait, Delta)]

php <- split(ph, by = "Param")
phm <- lapply(php, function(x){ acast(x[,.(PC, Trait, Delta)], PC ~ Trait) })


phhm <- lapply(1:9, function(x){ 
		       ds <- phm[[x]]
		       name <- names(phm)[x]
		       range <- max(abs(ds))
		       pheatmap(ds,  breaks = seq(-range, range, length.out = 100), cluster_cols = TRUE, cluster_rows=FALSE, main=name )
})

phhm[[1]]
phhm[[2]]
phhm[[3]]
phhm[[4]]
phhm[[5]]
phhm[[6]]
phhm[[7]]
phhm[[8]]
phhm[[9]]

