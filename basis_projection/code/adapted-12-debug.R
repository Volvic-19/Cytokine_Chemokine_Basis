# This script is to debug the inconsistency in some traits during the plotting of Ahola and Ferkingstad projection results
# To check: "IL1beta","IL6","IL7"

library(data.table)
setDTthreads(10)

# Load rotation matrix, get the top SNPs in PC1
load("../Rdata/cytokine_completeSNP.Rdata")

str(rot.pca) #num [1:673318, 1:105]
rotation.DT <- as.data.table(rot.pca)

top.snp <- dimnames(rot.pca)[[1]][order(abs(rotation.DT$PC1), decreasing = TRUE)][1:3] #[1] "2:3592552"  "2:3594771"  "6:96614917" "17:35954653" "3:56815721"
top.snp.value <- rotation.DT$PC1[order(abs(rotation.DT$PC1), decreasing = TRUE)][1:3] #[1] -0.96476711 -0.12001866  0.06993009 -0.04595821  0.04285648

# Check top SNPs in  Ferkingstad data,  Ahola data, for traits: "IL1beta","IL6","IL7"
IL1b_fk_raw <- fread("../../basis_building/raw_fk_data/3037_62_IL1B_IL_1b.txt.gz") #30848421       12
fk_reduced <- fread("../../basis_building/manifest/shrinkage.e5.DT.tsv.gz", tmpdir = "tmp")

IL1b_ah_raw <- fread("~/rds/rds-cew54-basis/02-Processed/IL1B_AholaOlli_27989323_1-hg38.tsv.gz",tmpdir = "tmp")
IL1b_ah_reduced <- fread("../data/IL1B_AholaOlli_27989323_1-ft.tsv") #637042 11


IL1b_fk_raw <- IL1b_fk_raw[, pid:=paste0(gsub("chr","",Chrom), ":", Pos)]

IL1b_fk_raw[pid %in% top.snp]
#> IL1b_fk_raw[pid %in% top.snp]
#   Chrom      Pos              Name     rsids effectAllele otherAllele    Beta
#1:  chr2  3592552  chr2:3592552:T:C rs6542680            T           C -0.3681
#2:  chr2  3594771  chr2:3594771:C:T rs3820897            C           T -0.3607
#3:  chr6 96614917 chr6:96614917:G:T  rs211161            G           T  0.1139
#         Pval minus_log10_pval       SE     N  ImpMAF        pid #annotation
#1: 4.148e-273        272.38216 0.010425 35332 0.18362  2:3592552   0.81638
#2: 4.683e-265        264.32948 0.010371 35333 0.18633  2:3594771   0.81367
#3:  4.118e-04          3.38531 0.032244 35343 0.01526 6:96614917   0.98474

fk_reduced[pid %in% top.snp][Trait == "3037_62_IL1B_IL_1b"]
#          pid REF_ALT    BETA          P       SE     N ALT_FREQ ld.block
#1:  2:3592552     T/C -0.3681 4.148e-273 0.010425 35332  0.81638      107
#2:  2:3594771     C/T -0.3607 4.683e-265 0.010371 35333  0.81367      107
#3: 6:96614917     G/T  0.1139  4.118e-04 0.032244 35343  0.98474      553
#                Trait ss_emp_maf_se      sdY          ppi        wj     ws_ppi
#1: 3037_62_IL1B_IL_1b      2.753808 1.083423 1.000000e+00 1.0000000 0.42441299
#2: 3037_62_IL1B_IL_1b      2.731972 1.083423 9.557728e-09 1.0000000 0.05337502
#3: 3037_62_IL1B_IL_1b      8.381990 1.083423 6.718616e-03 0.1466297 0.26577102
#    shrinkage       metric
#1: 0.15411859 -0.056731054
#2: 0.01953717 -0.007047059
#3: 0.03170739  0.003611472


IL1b_ah_raw[SNPID %in% c("rs6542680","rs3820897","rs211161")]
#          pid     SNPID CHR38     BP38 CHR       BP REF ALT        BETA
#1:  2:3640142 rs6542680     2  3592552   2  3640142   C   T 0.004214766
#2:  2:3642361 rs3820897     2  3594771   2  3642361   T   C 0.002897652
#3: 6:97062793  rs211161     6 96614917   6 97062793   T   G 0.011854030
#           SE Direction      P HetPVal orig.alleles ALT_FREQ ref.alleles
#1: 0.02528860      +--+ 0.9544  0.7770          T/C 0.651515         C/T
#2: 0.02515689      -++- 0.9956  0.8174          T/C 0.651515         T/C
#3: 0.06467032      -++- 0.7099  0.3794          T/G 0.969697         T/G
#   orig.BETA orig.SE
#1:   -0.0032  0.0192
#2:    0.0022  0.0191
#3:    0.0090  0.0491

IL1b_ah_reduced[pid %in% top.snp]
#> IL1b_ah_reduced[pid %in% top.snp]
#          pid     SNPID CHR38     BP38 REF ALT        BETA         SE      P
#1:  2:3592552 rs6542680     2  3592552   C   T 0.004214766 0.02528860 0.9544
#2:  2:3594771 rs3820897     2  3594771   T   C 0.002897652 0.02515689 0.9956
#3: 6:96614917  rs211161     6 96614917   T   G 0.011854030 0.06467032 0.7099
#   ALT_FREQ REF_ALT
#1: 0.651515     C/T
#2: 0.651515     T/C
#3: 0.969697     T/G



# "IL6"
IL6_fk_raw <- fread("../../basis_building/raw_fk_data/4673_13_IL6_IL_6.txt.gz")
IL6_ah_raw <- fread("~/rds/rds-cew54-basis/02-Processed/IL6_AholaOlli_27989323_1-hg38.tsv.gz",tmpdir = "tmp")
IL6_ah_reduced <- fread("../data/IL6_AholaOlli_27989323_1-ft.tsv")
IL6_folkersen_raw <- fread("~/rds/rds-cew54-basis/02-Processed/IL6_Folkersen_33067605_1-hg38.tsv.gz")
IL6_folkersen_reduced <- fread("../data/IL6_Folkersen_33067605_1-ft.tsv")


IL6_fk_raw <- IL6_fk_raw[, pid:=paste0(gsub("chr","",Chrom), ":", Pos)]

IL6_fk_raw[pid %in% top.snp]
#   Chrom      Pos              Name     rsids effectAllele otherAllele    Beta
#1:  chr2  3592552  chr2:3592552:T:C rs6542680            T           C -0.0214
#2:  chr2  3594771  chr2:3594771:C:T rs3820897            C           T -0.0207
#3:  chr6 96614917 chr6:96614917:G:T  rs211161            G           T  0.1395
#         Pval minus_log10_pval       SE     N  ImpMAF        pid
#1: 3.3551e-02          1.47429 0.010069 35332 0.18362  2:3592552
#2: 3.8283e-02          1.41699 0.009991 35333 0.18633  2:3594771
#3: 5.7380e-06          5.24124 0.030755 35343 0.01526 6:96614917

fk_reduced[pid %in% top.snp][Trait == "4673_13_IL6_IL_6"]
#         pid REF_ALT    BETA          P       SE     N ALT_FREQ ld.block
#1:  2:3592552     T/C -0.0214 3.3551e-02 0.010069 35332  0.81638      107
#2:  2:3594771     C/T -0.0207 3.8283e-02 0.009991 35333  0.81367      107
#3: 6:96614917     G/T  0.1395 5.7380e-06 0.030755 35343  0.98474      553
#              Trait ss_emp_maf_se      sdY          ppi        wj     ws_ppi
#1: 4673_13_IL6_IL_6      2.753808 1.032703 5.374745e-05 0.1257900 0.42441299
#2: 4673_13_IL6_IL_6      2.731972 1.032703 4.770087e-05 0.1257900 0.05337502
#3: 4673_13_IL6_IL_6      8.381990 1.032703 2.311915e-01 0.4024496 0.26577102
#    shrinkage        metric
#1: 0.15411859 -0.0032981379
#2: 0.01953717 -0.0004044195
#3: 0.03170739  0.0044231810

IL6_ah_raw[SNPID %in% c("rs6542680","rs3820897","rs211161")]
#          pid     SNPID CHR38     BP38 CHR       BP REF ALT        BETA
#1:  2:3640142 rs6542680     2  3592552   2  3640142   C   T 0.010219025
#2:  2:3642361 rs3820897     2  3594771   2  3642361   T   C 0.009722956
#3: 6:97062793  rs211161     6 96614917   6 97062793   T   G 0.040082389
#           SE Direction      P HetPVal orig.alleles ALT_FREQ ref.alleles
#1: 0.01637028       +-- 0.5437  0.5510          T/C 0.651515         C/T
#2: 0.01627107       -++ 0.5604  0.5720          T/C 0.651515         T/C
#3: 0.04166981       +++ 0.3389  0.9237          T/G 0.969697         T/G
#   orig.BETA orig.SE
#1:   -0.0103  0.0165
#2:    0.0098  0.0164
#3:    0.0404  0.0420


IL6_ah_reduced[pid %in% top.snp]
#          pid     SNPID CHR38     BP38 REF ALT        BETA         SE      P
#1:  2:3592552 rs6542680     2  3592552   C   T 0.010219025 0.01637028 0.5437
#2:  2:3594771 rs3820897     2  3594771   T   C 0.009722956 0.01627107 0.5604
#3: 6:96614917  rs211161     6 96614917   T   G 0.040082389 0.04166981 0.3389
#   ALT_FREQ REF_ALT
#1: 0.651515     C/T
#2: 0.651515     T/C
#3: 0.969697     T/G

IL6_folkersen_raw[pid %in% top.snp]


IL6_folkersen_reduced[pid %in% top.snp]
#          pid          SNPID CHR38     BP38 REF ALT        BETA         SE
#1:  2:3592552  2_3640142_C_T     2  3592552   C   T 0.007551039 0.01009225
#2:  2:3594771  2_3642361_C_T     2  3594771   T   C 0.007914069 0.01009225
#3: 6:96614917 6_97062793_G_T     6 96614917   T   G 0.011544376 0.03942513
#        P ALT_FREQ REF_ALT
#1: 0.4543   0.7462     C/T
#2: 0.4308   0.7446     C/T
#3: 0.7692   0.9803     G/T