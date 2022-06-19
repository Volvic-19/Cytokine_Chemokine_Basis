# Compute shrinkage for gwas.DT

# Once we have filtered and aligned Chen datasets, we'll create a single Chen file with all 8 datasets.

# 2022-06-08
# Guillermo Reales

# Load libraries
library(data.table)
setDTthreads(15)

#########################################################

############################
# Helper functions   #######
############################

# These functions are similar to those in Blood cell basis v2

maf_se_estimate_sample_size_quant <- function(N,se.beta,f){
  se_maf_ss <- sqrt(2 * N) * se.beta
  idx <- which(is.infinite(se_maf_ss) | is.nan(se_maf_ss) | se_maf_ss>100)
  se_maf_ss[idx] <- maf_se_estimate(f[idx])
  se_maf_ss
}

wakefield_pp_quant <- function(beta, se.beta, N, sdY, sd.prior, pi_i=1e-4) { # Adapted
  # compute V
  V <- se.beta^2
  # Compute z too
  z <- beta/se.beta
  # Multiply prior by sdY
  prior <- sd.prior*sdY
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- prior^2 / (prior^2 + V)
  ## Approximate BF
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ## tABF - to add one we create another element at the end of zero for which pi_i is 1
  tABF <- c(lABF,0)
  vpi_i<-c(rep(pi_i,length(lABF)),1)
  sBF <- logsum(tABF + log(vpi_i))
  exp(lABF+log(pi_i)-sBF)
}

compute_shrinkage_metrics_quant <- function(DT, pi_i = 1e-04,sd.prior = 0.15){
  # Create a copy of the DT to play with it
  tmp.shrinkage.DT <- copy(DT)
  # Calculate ss_emp_maf_se - Variance in the population per SNP
  tmp.shrinkage.DT[, ss_emp_maf_se:=maf_se_estimate_sample_size_quant(n, se, maf), by = pid]
  tmp.shrinkage.DT[, ss_emp_maf_se:= mean(abs(ss_emp_maf_se)), by = pid]
  # Calculate sdY - Trait standard deviation
  tmp.shrinkage.DT[, sdY:=sdY.est(se^2, maf, n), by = trait]
  # Calculate ppi (PPti in the paper) - Posterior probability of each SNP in each region to be causal
  tmp.shrinkage.DT[, ppi:=wakefield_pp_quant(beta, se, n, sdY, sd.prior, pi_i), by = c('trait', 'ld.block')]
  # Calculate wj (Vtr in the paper) - Sum of all ppi per region, to get the probability that the region contains a causal SNP
  tmp.shrinkage.DT[, wj:=sum(ppi), by = c('trait', 'ld.block')]
  # Calculate ws_ppi (Wi in the paper) - SNP weight, a weighted average of PPti
  tmp.shrinkage.DT[, ws_ppi:=sum(ppi * wj)/sum(wj), by = c('pid', 'ld.block')]
  # Finally, we divide the result by the variance in the population per SNP in order to get the final shrinkage metric
  tmp.shrinkage.DT[, shrinkage := ws_ppi/ss_emp_maf_se, by = 'pid']
  # And multiply by beta to get our metric for each SNP
  tmp.shrinkage.DT[, metric := shrinkage * beta]
  return(tmp.shrinkage.DT)
}

### Unchanged from original

maf_se_estimate <- function(f){
  #1/sqrt(f * (1-f))
  sqrt(1/f + 1/(1-f)) * 2
}

sdY.est <- function(vbeta, maf, n) {
  warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  cf <- coef(m)[['oneover']]
  if(cf < 0)
    stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
  return(sqrt(cf))
}

logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form)
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}


#########################################################
###### Load file and apply shrinkage    #################
#########################################################


gwas.DT <- fread("../data/gwas.DT.tsv.gz", tmpdir="tmp")

# Change columns names momentarily
setnames(gwas.DT, c("BETA", "SE", "ALT_FREQ", "N", "Trait"), c("beta", "se", "maf", "n", "trait"))

# Apply shrinkage!!
shrinkage.DT <- compute_shrinkage_metrics_quant(gwas.DT)
summary(shrinkage.DT)


# Change names back
setnames(shrinkage.DT, c("beta", "se", "maf", "n", "trait"), c("BETA", "SE", "ALT_FREQ", "N", "Trait") )
shrinkage.DT

# Filter shrinkage table by columns 
shr.wk4 <- shrinkage.DT[ shrinkage > 1e-4 ,.(pid, Trait, shrinkage, metric)]
shr.wk5 <- shrinkage.DT[ shrinkage > 1e-5 ,.(pid, Trait, shrinkage, metric)]
shr.wk6 <- shrinkage.DT[ shrinkage > 1e-6 ,.(pid, Trait, shrinkage, metric)]
shr.wk2 <- shrinkage.DT[ shrinkage > 1e-2 ,.(pid, Trait, shrinkage, metric)]


# Save shrinkage metrics
fwrite(shrinkage.DT, "../data/shrinkage.DT.tsv.gz", sep="\t")
fwrite(shr.wk4, "../data/shrinkage.DT.1e4.tsv.gz", sep="\t")
fwrite(shr.wk5, "../data/shrinkage.DT.1e5.tsv.gz", sep="\t")
fwrite(shr.wk6, "../data/shrinkage.DT.1e6.tsv.gz", sep="\t")
fwrite(shr.wk2, "../data/shrinkage.DT.1e2.tsv.gz", sep="\t")
