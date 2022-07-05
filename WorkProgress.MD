
# 1. Pre-work

## Decide the GWAS summary statistics datasets to use

1. Compare 4 papers (Ferkingstad, Ahola-Olli, Lasse Folkersen, Sinnott-Armstrong) based on their coverage of common cytokine/chemokine.

   * Generate table of common cytokine/chemokines from book  [<Cytokines, Chemokines and Their Receptors>](https://www.ncbi.nlm.nih.gov/books/NBK6294/)
   * Map the availabe traits from 4 papers with common cytokine/chemokine table based on Uniprot ID. Summary of coverage is in [Summary of coverage](https://docs.google.com/spreadsheets/d/1lXT4dmS9ogFrl2nv5Qr2KHEoDk5g2vf2kOv9OzG8El4/edit#gid=1947277513)
   * Ferkingstad datasets show the best coverage of cytokine/chemokine: 50 (out of 52) common cytokine, 37 (out of 40) common chemokine, covering all subclasses listed in the book table.

  * Relevant Code: [Cytokine_chemokine_coverage.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/Pre-work/Cytokine_chemokine_coverage.R)

## Download desired cytokine/chemokine datasets from Ferkingstad

1. Map downloaded files with protein names via unique SeqId
   * Datasets of 4,907 traits were downloaded from [decode](https://www.decode.com/summarydata/)
   * Map file names with protein names via unique SeqId
   
2. Select desired datasets to use
   * Filter datasets for common cytokines/chemokines
   * Remove duplicated datasets (i.e. two or more files sharing the same Uniprot ID): Manhattan plots were made to compare duplicated datasets. We intended to keep only one for each protein trait, but as quality of GWAS stats are close and there is no clear ground to eliminate duplicates, we keep all of them (in total 104 traits) in our current basis.
   * Copy desired datasets (in total 104 .txt.gz files) and save them under rds/rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_building/raw_fk_data
  
  * Relevant Code: [map_dataset_names.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/Pre-work/map_dataset_names/map_dataset_names.R)
 
 
 # 2. Basis-building

## Create Manifest

1. Intersect of **UKBB Asthma**, **PanUKBB Asthma**, **FinnGen Asthma** and **IL5 (a randomly selected Ferkingstad trait)** to create the initial manifest
    * Input : 11071_1_IL5_IL_5.txt.gz, Asthma_PanUKBB_processed-hg38.tsv.gz, Asthma_Neale_processed-hg38.tsv.gz, Asthma_FinnGen.tsv.gz
    * Relevant Code : [adapted-00-create_IL5_consensus_SNP_list.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_building/code/adapted-00-create_IL5_consensus_SNP_list.R)
    * Output : rds-cew54-basis/Projects/Cytokine_Chemokine_Basis/basis_building/manifest/IL5_consensus_manifest_8M.tsv

2. 





