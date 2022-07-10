# This script is to test if run standarization on raw Ferkingstad dataset will make a difference on beta value
# sdY.correction function is adapted from https://github.com/GRealesM/IMDtools/blob/master/R/gr-wrenches.R#L45

library(data.table)
setDTthreads(10)

##### Function to standarize beta

sdY.correction <- function(beta, se, maf, n) { # This standarization process is often carried on right after downloading raw datasets
    if(length(beta) != length(maf) | length(se) != length(maf)) stop("Beta, SE and/or MAF are not the same length")
    oneover <- 1/se^2
    nvx <- 2 * n * maf * (1-maf) #imf does not matter
    m <- lm(nvx ~ oneover - 1)
    cf <- coef(m)[['oneover']]
    if(cf < 0)
        stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
    message("Estimated sdY is ", sqrt(cf))
    BETA  <- beta/sqrt(cf)
    SE  <- se/sqrt(cf)
    list(BETA,SE)
}


##### load Ferkingstad raw dataset.

IL1b_fk_raw <- fread("../../basis_building/raw_fk_data/3037_62_IL1B_IL_1b.txt.gz") #30848421       12
# Chrom, Pos, Name, rsids, effectAllele, otherAllele, Beta, Pval, minus_log10_pval, SE,N, ImpMAF

a <- sdY.correction(IL1b_fk_raw$Beta, IL1b_fk_raw$SE, IL1b_fk_raw$ImpMAF, IL1b_fk_raw$N)
IL1b_fk_raw[, corrected.Beta := unlist(a[[1]])][,corected.SE := unlist(a[[2]])]

cor(IL1b_fk_raw$Beta, IL1b_fk_raw$corrected.Beta)
# [1] 1

#IL1b_fk_raw
#          Chrom       Pos                                 Name        rsids
#       1:  chr1     23597                       chr1:23597:G:A         <NA>
#       2:  chr1     92875                       chr1:92875:C:T  rs193157612
#       3:  chr1    107682                      chr1:107682:G:C  rs879827054
#       4:  chr1    113324                     chr1:113324:TA:T rs1235985416
#       5:  chr1    113394                      chr1:113394:G:T rs1366027486
#      ---                                                                  
#30848417:  chrX 156029849                   chrX:156029849:G:A         <NA>
#30848418:  chrX 156029857                   chrX:156029857:G:C rs1384946096
#30848419:  chrX 156029880 chrX:156029880:CATGTGC:CATGTGCTTAGGG rs1440518544
#30848420:  chrX 156029926                   chrX:156029926:C:G rs1378190828
#30848421:  chrX 156029936                   chrX:156029936:A:T rs1434518991
#          effectAllele   otherAllele    Beta     Pval minus_log10_pval       SE
#       1:            G             A -0.0198 0.967768          0.01423 0.490004
#       2:            C             T  0.2501 0.595796          0.22490 0.471481
#       3:            G             C  0.7104 0.263016          0.58002 0.634687
#       4:           TA             T  0.0774 0.165432          0.78138 0.055803
#       5:            G             T  0.0040 0.949865          0.02234 0.063617
#      ---                                                                      
#30848417:            G             A  0.0786 0.367708          0.43450 0.087258
##30848418:            G             C -0.2455 0.051602          1.28733 0.126128
#30848419:      CATGTGC CATGTGCTTAGGG  0.0677 0.237086          0.62509 0.057261
#30848420:            C             G -0.0753 0.166988          0.77731 0.054488
#30848421:            A             T -0.0037 0.865258          0.06285 0.021805
#              N  ImpMAF corrected.Beta corected.SE
#       1: 35287 0.00011   -0.018315468  0.45326528
#       2: 35287 0.00012    0.231348412  0.43613107
#       3: 35287 0.00011    0.657136792  0.58710048
#       4: 35287 0.00539    0.071596829  0.05161909
#       5: 35287 0.00461    0.003700095  0.05884723
#      ---                                         
#30848417: 30587 0.00195    0.072706858  0.08071571
#30848418: 30587 0.00094   -0.227093303  0.11667138
#30848419: 30587 0.00365    0.062624100  0.05296778
#30848420: 30587 0.00485   -0.069654280  0.05040269
#30848421: 30587 0.03115   -0.003422587  0.02017014

# After trying to standarize beta in Ferkingstad raw data, no dramatic shrifting of data. "corrected.Beta" very close to
# original "Beta"
