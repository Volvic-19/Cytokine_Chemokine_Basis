
# 1. Pre-work

## Decide the GWAS summary statistics datasets to use

1. Compare 4 papers (Ferkingstad, Ahola-Olli, Lasse Folkersen, Sinnott-Armstrong) based on their coverage of common cytokines/chemokines.

   * Generate the table of common cytokines/chemokines from the review book by Mark J. Cameron and David J. Kelvin.[<Cytokines, Chemokines and Their Receptors>](https://www.ncbi.nlm.nih.gov/books/NBK6294/)
   * Map the availabe traits from 4 papers with common cytokine/chemokine table based on Uniprot ID. Comparison of coverage is summarised in [this google spreadsheet](https://docs.google.com/spreadsheets/d/1lXT4dmS9ogFrl2nv5Qr2KHEoDk5g2vf2kOv9OzG8El4/edit#gid=1947277513)
   * Ferkingstad datasets show the best coverage of cytokine/chemokine: 50 (out of 52) common cytokine, 37 (out of 40) common chemokine, covering all subclasses of cytokines/chemokines listed in the book table.

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

## Create Manifest and process Ferkingstad datasets

1. Take intersect SNPs of **UKBB Asthma**, **PanUKBB Asthma**, **FinnGen Asthma** and **IL5 (a randomly selected Ferkingstad trait)** to create the initial manifest, add LD block. We align "REF" and "ALT" in our initial manifest with "Reference_1000GP_manifest.tsv.gz". In our initial manifest, there is only one allele in one pid position. But in Reference_1000GP_manifest.tsv.gz, same pid can have different alleles. Be aware of that.
    * Input : `11071_1_IL5_IL_5.txt.gz`, `Asthma_PanUKBB_processed-hg38.tsv.gz`, `Asthma_Neale_processed-hg38.tsv.gz`, `Asthma_FinnGen.tsv.gz`,`Reference_1000GP_manifest.tsv.gz`
    * Relevant Code : [adapted-00-create_IL5_consensus_SNP_list.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_building/code/adapted-00-create_IL5_consensus_SNP_list.R)
    * Output : `basis_building/manifest/IL5_consensus_manifest_8M.tsv`

2. Note that not all ImpMAF in raw Ferkingstad data refers to effect allele frequency. Therefore we first mapped and created `ALT_FREQ` using `assocvariants.annotated.txt.gz`, which is an annotation file provided by Ferkingstad to address this issue. Next we aligned all Ferkingstad datasets with initial manifest (IL5_consensus_manifest_8M.tsv), and eventually got 104 "manifested_data", in which only SNPs shared with initial manifest were kept.
    * Input : `assocvariants.annotated.txt.gz`, `IL5_consensus_manifest_8M.tsv`, `basis_building/raw_fk_data/`
    * Relevant Code : [adapted-01-Ferkingstad_datasets.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_building/code/adapted-01-Ferkingstad_datasets.R) ` Previous mistake of "REF_ALT" has been corrected in this new run. `
    * Output : `basis_building/manifested_data/`

3.  Merge 104 Ferkingstad datasets into one table. This is an intermediate step for slurm convenience.
    * Input : `basis_building/manifested_data/`
    * Relevant Code : [adapted-02-mergetables.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_building/code/adapted-02-mergetables.R)
    * Output : `basis_building/manifest/gwas.DT.tsv.gz`

4.  Remove SNPs with ALT_FREQ < 1% or >99%. This is an intermediate step to help size down `gwas.DT.tsv.gz`
    * Input : `basis_building/manifest/gwas.DT.tsv.gz`
    * Relevant Code : [adapted-03.1-FilterAlt_Freq.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_building/code/adapted-03.1-FilterAlt_Freq.R)
    * Output : `basis_building/manifest/filtered.1.99.gwas.DT.tsv.gz`

5.  Keep only SNPs present in all 104 Ferkingstad traits. This is an intermediate step. Non-essential columns were also removed to save computing memory, as otherwise Slurm job could easily be terminated by HPC.
    * Input : `basis_building/manifest/filtered.1.99.gwas.DT.tsv.gz`
    * Relevant Code : [adapted-03.2-overlap-removeCol.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_building/code/adapted-03.2-overlap-removeCol.R)
    * Output : `basis_building/manifest/overlap.gwas.DT.tsv.gz`


## Compute Shrinkage table and PCA

1.  Compute shrinkage, output complete shrinkage table, shrinkage table with threshold set on 1e-5, 1e-6
    * Input : `basis_building/manifest/overlap.gwas.DT.tsv.gz`
    * Relevant Code : [adapted-03.3-Compute_shrinkage.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_building/code/adapted-03.3-Compute_shrinkage.R)
    * Output : `basis_building/manifest/shrinkage.DT.tsv.gz`, `shrinkage.e5.DT.tsv.gz`, `shrinkage.e6.DT.tsv.gz`, 

2.  Update manifest. This is an intermediate step to keep only SNPs shared by all 104 Ferkingstad traits in our manifest.
    * Input : `basis_building/manifest/IL5_consensus_manifest_8M.tsv`
    * Relevant Code :[adapted-04-update_manifest.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_building/code/adapted-04-update_manifest.R)
    * Output : `basis_building/manifest/IL5_consensus_manifest_6M.tsv`


3.  Apply threshold on shrinkage table. This is an intermediate step to size down shrinkage table and SNPs to use.
    * Input : `basis_building/manifest/shrinkage.e5.DT.tsv.gz`
    * Relevant Code :[adapted-05-Threshold_shrinkagetable.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_building/code/adapted-05-Threshold_shrinkagetable.R)
    * Output : `basis_building/manifest/shrinkage.e4.DT.tsv.gz`, `shrinkage.e3.DT.tsv.gz`, `shrinkage.e2.DT.tsv.gz`

4.  PCA. We take shrinkage table with threshold at 1e-5 as input, and save the computed basis.
    * Input : `basis_building/manifest/shrinkage.e5.DT.tsv.gz`
    * Relevant Code :[adapted-06-PCA.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_building/code/adapted-06-PCA.R)
    * Output : `basis_building/PCA/cytokine_basis.RDS`, `cytokine_basis_matrix.RDS`

5.  Update manifest. This is an intermediate step, after update, manifest will only keep SNPs whose shrinkage value  is > 1e-5 in shrinkage table. (i.e. This is final version of manifest, where all SNPs are the same as our basis. In projection process, all external datasets will be reduced and aligned with this final manifest before projection.)
    * Input : `basis_building/manifest/IL5_consensus_manifest_6M.tsv`
    * Relevant Code :[adapted-07-update_manifest_e5.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_building/code/adapted-07-update_manifest_e5.R)
    * Output : `IL5_consensus_manifest_6M_e5.tsv`



 # 3. Basis-projection
 
1.  Reduction: Before projection, reduce and align external datasets with final manifest (e5). List of external datasets to project is [here](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_projection/code/traits_to_project_early_test.txt) We projected in total 261 external datasets, the first no. 1-248 datasets were processed through `adapted-08.1-reduction.R`, the rest no. 249-261 datasets were processed with `adapted-08.2-reduction.R` due to format difference. I noticed that many datasets have less than 95% coverage of basis SNPs, possibly because we are not using sparse SNPs and the total size of SNP is still quite large.
    * Input : `IL5_consensus_manifest_6M_e5.tsv`, `~/rds/rds-cew54-basis/02-Processed/`
    * Relevant Code : [adapted-08.1-reduction.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_projection/code/adapted-08.1-reduction.R)  [adapted-08.2-reduction.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_projection/code/adapted-08.2-reduction.R)
    * Output : `basis_projection/data/`

2.  Pack projection function and basis into one Rdata, to save future repeatitive computation. 
    * Input :`cytokine_basis.RDS`, `cytokine_basis_matrix.RDS`, `shrinkage.e5.DT.tsv.gz`, `IL5_consensus_manifest_6M_e5.tsv`
    * Relevant Code: [adapted-09-saveBasisRdata.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_projection/code/adapted-09-saveBasisRdata.R)
    * Output : `/Rdata/cytokine_completeSNP.Rdata`

3.  Projection: 261 external datasets 
    * Input : `cytokine_completeSNP.Rdata`, `basis_projection/data/`
    * Relevant Code: [adapted-10.1-projection-1-261.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_projection/code/adapted-10.1-projection-1-261.R)
    * Output : `basis_projection/Projections`

4.  Projection: Ferkingstad datasets (We don't need to reduce them, instead we can use shrinkage table which has been aligned to manifest)
    * Input : `cytokine_completeSNP.Rdata`, `shrinkage.e5.DT.tsv.gz`
    * Relevant Code: [adapted-10.2-projection-basis_to_basis.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_projection/code/adapted-10.2-projection-basis_to_basis.R)
    * Output : `basis_projection/Projections`

5.  Plot projections, we plot Ferkingstad projection together with Ahola / Folkerson projections. In order to unify names used in different studies, we mapped each trait name with "Common_name", which are the names used in [<Cytokines, Chemokines and Their Receptors>] book. As the plots look a little odd, I have checked the "odd data points" in projection table, and can confirm that the odd results are genuinely because of computed delta value, not any error introduced in plotting.
    * Input : `basis_projection/Projections`
    * Relevant Code: [adapted-11-plotprojection-intersectAholaFolkersen.R](https://github.com/Volvic-19/Cytokine_Chemokine_Basis/blob/main/basis_projection/code/adapted-11-plotprojection-intersectAholaFolkersen.R)
    * Output :`basis_projection/Plots/`






