# Decide the GWAS summary statistics datasets to use

1. Compare 4 papers (Ferkingstad, Ahola-Olli, Lasse Folkersen, Sinnott-Armstrong) based on their coverage of common cytokine/chemokine.

   * Generate table of common cytokine/chemokines from book  [<Cytokines, Chemokines and Their Receptors>](https://www.ncbi.nlm.nih.gov/books/NBK6294/)
   * Map the availabe traits from 4 papers with common cytokine/chemokine table based on Uniprot ID. Summary of coverage is in [Summary of coverage](https://docs.google.com/spreadsheets/d/1lXT4dmS9ogFrl2nv5Qr2KHEoDk5g2vf2kOv9OzG8El4/edit#gid=1947277513)
   * Ferkingstad datasets show the best coverage of cytokine/chemokine: 50 (out of 52) common cytokine, 37 (out of 40) common chemokine, covering all subclasses listed in the book table.

  * Relevant Code: [Cytokine_chemokine_coverage.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/Pre-work/Cytokine_chemokine_coverage.R)

# Download desired cytokine/chemokine datasets from Ferkingstad

1. Map downloaded files with protein names via unique SeqId
   * Datasets of 4,907 traits were downloaded from [decode](https://www.decode.com/summarydata/)
   * Map file names with protein names via unique SeqId
   
2. Select desired datasets to use
   * Filter datasets for common cytokines/chemokines
   * Remove duplicated datasets (i.e. two or more files sharing the same Uniprot ID): Manhattan plots were made to compare duplicated datasets. We intended to keep only one for each protein trait, but as quality of GWAS stats are close and there is no clear ground to eliminate duplicates, we keep all of them (in total 104 traits) in our current basis.
   * Copy desired datasets (.txt.gz files) and save them under /rds/project/rds-YRWZsjDGyaU/Projects/Cytokine_Chemokine_Basis/basis_building
  
  * Relevant Code: [map_dataset_names.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/Pre-work/map_dataset_names.R)
 
 
# Create Manifest: list of SNPs to use in our cytokine/chemokine basis

1. Intersect of **UKBB Asthma**, **PanUKBB Asthma**, **FinnGen Asthma**
  * Download 3 datasets, liftover two UKBB files from hg19 to hg38 build.[ongoing]
  * 
