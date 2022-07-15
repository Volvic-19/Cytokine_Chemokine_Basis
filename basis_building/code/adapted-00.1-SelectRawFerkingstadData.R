# This script is to select and save raw Ferkingstad data into a directory for further usage.
# Selected traits meet the flollowing requirement : have a peak within the region of 100kb upstream + gene length + 100kb downstream
# In total we select 40 out of 104 traits, or 38 out of 100 proteins.
# Table was generated from basis_building/code/get_Ferkingstad_gene_regions.R

library(data.table)

# load ps table
ps <- fread("../../Pre-work/Selected_traits/Ferkingstad_peaks.tsv")
selected.traits <- unique(ps$FileName)
writeLines(selected.traits, "../../Pre-work/Selected_traits/selected_Ferkingstad_list.tsv")

# copy selected datasets from raw_fk_data to raw_fk_data_40
lapply(selected.traits, function(x){
    file.copy(from = paste0("../raw_fk_data/", x),
    to = "../raw_fk_data_40", overwrite = TRUE)
})
