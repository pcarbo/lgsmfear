# This script maps QTLs using genotype and phenotype data from the F2
# and F34 crosses.

# SCRIPT PARAMETERS
# -----------------
qtl.method   <- "hk"       # Which QTL mapping method to use in qtl.
relatedness  <- "markers"  # How to estimate relatedness.
map.function <- "Haldane"  # Map function to use for interval mapping.
num.perm     <- 1000       # Number of replicates for permutation test.
jitter.amt   <- 1e-6       # Amount by which marker positions are adjusted.
deps         <- 0.01       # Small number added to diagonal entries of
                           # variance components in combined analysis
                           # using pedigree data.

# Save the results to this file.
resultsfile <- "gwscan.rr.RData"

# Separately map QTLs for these phenotypes, and use these covariates
# in QTL mapping for all phenotypes.
phenotypes <- c("pretrainfreeze","freezetocontext","freezetocue",
                "distperiphery","distcenter","propcenter")
covariates <- c("sex","age","albino","agouti")

# Initialize the random number generator.
set.seed(7)

# Load packages and function definitions.
library(abind)
library(qtl)
capture.output(library(QTLRel))
source("misc.R")
source("read.data.R")
source("data.manip.R")
source("mapping.tools.R")

# LOAD DATA
# ---------
# Load the genotype, phenotype and marker data for the combined
# cohort. I remove the first two columns of the genotype data frame
# containing the ID and generation, and set the row names for the
# genotype data to the mouse IDs. I adjust the map positions
# (i.e. genetic distances) slightly so that no two markers have the
# same position. I also reorder the rows of the phenotype and genotype
# data so that F2 mice are in the top rows and the F34 mice are in the
# bottom rows.
cat("Loading data.\n")
pheno           <- read.phenotypes("../data/pheno.csv")
pheno           <- pheno[order(pheno$generation),]
geno            <- read.genotypes("../data/geno.csv")
geno            <- geno[order(geno$generation),]
row.names(geno) <- geno$id
geno            <- geno[-(1:2)]
map             <- read.map("../data/map.csv",chromosomes)
map             <- jitter.gendist(map,jitter.amt)

# Read in the pedigree data. First I read in the pedigree table, and
# remove the first two rows of the table, which correspond to the
# inbred founders (generation F0). Since QTLRel does not allow the
# founders to be inbred, I artificially create inbred founds by
# including several generations of self-mating to produce the male and
# female founders. In practice, we find that identity and kinship
# coefficients for the F0, F1 and F2 now match approximately what you
# would get if the founders were inbred (with ~1% error).
if (relatedness == "pedigree") {
  cat("Loading pedigree.\n")
  ped <- read.pedigree("../data/ped.csv")
  ped <- ped[-(1:2),]
  ped <- rbind(read.pedigree("../data/inbred.ped.csv"),ped)
  row.names(ped) <- NULL
}

# DISCARD SAMPLES
# ---------------
# I remove three mice from the F2 cohort because a large proportion of
# their genotypes are not available.
rows  <- which(!is.element(pheno$id,c("648A","738A","811A")))
pheno <- pheno[rows,]
geno  <- geno[rows,]

# GET F2 and F34 CROSSES
# ----------------------
# Get the rows of the genotype and phenotype matrices corresponding to
# the mice in the F2 and F34 crosses.
F2.rows  <- which(pheno$generation == "F2")
F34.rows <- which(pheno$generation == "F34")
n2       <- length(F2.rows)
n34      <- length(F34.rows)

# DISCARD MARKERS
# ---------------
# For now, I drop the X chromosome from the analysis. I also remove
# all markers that are *not* genotyped in the F34 cross.
markers <- which(map$chr != "X" & !all.missing.col(geno[F34.rows,]))
geno    <- geno[,markers]
map     <- transform(map[markers,],chr = droplevels(chr))

# Get the markers genotyped in the F2 cross. Since I'm only retaining
# markers that are genotyped in the F34 cross, any marker genotyped in
# the F2 cross is also genotyped in the F34 cross.
F2.markers <- which(!all.missing.col(geno[F2.rows,]))

# TRANSFORM TRAITS
# ----------------
# I transform pretrainfreeze, freezetocontext, freezetocue and
# propcenter to the log-odds scale using the logit(x) function (in
# base 10). I transform sex, age, and create traits albino and agouti
# from the coat colour data.
cols        <- c("pretrainfreeze","freezetocontext",
                 "freezetocue","propcenter")
flogit      <- function (x) logit10(project.onto.interval(x,0.01,0.99))
pheno[cols] <- lapply(pheno[cols],flogit)
pheno       <- transform(pheno,
                         sex    = factor2integer(sex) - 1,
                         age    = age - mean(age),
                         albino = as.integer(coatcolor == "W"),
                         agouti = as.integer(coatcolor == "A"))

# COMPUTE GENETIC MATRICES
# ------------------------
# Calculate the genetic matrices from the identity coeffcients for
# both the F2 and F34 cohorts.
if (relatedness == "pedigree") {
  cat("Calculating genetic matrices from identity coefficients.\n")
  capture.output(
    F2.idcf <- cic(ped,ids = pheno$id[F2.rows],df = 100,
                   ask = FALSE,verbose = FALSE))
  load("../data/F34.idcf.RData")
  F2.gm  <- genMatrix(F2.idcf)
  F34.gm <- genMatrix(F34.idcf)

  # Combine the F2 and F34 genetic matrices. I add a small scalar to
  # the diagonal entries of the variance components to ensure that the
  # covariance matrix for the polygenic effects remains symmetric
  # positive definite.
  n   <- n2 + n34
  ids <- c(rownames(F2.gm$AA),rownames(F34.gm$AA))
  gm  <- list(AA = rbind(cbind(F2.gm$AA,ones(n2,n34)),
                        cbind(ones(n34,n2),F34.gm$AA)),
              DD = rbind(cbind(F2.gm$DD,ones(n2,n34)/4),
                         cbind(ones(n34,n2)/4,F34.gm$DD)))
  gm$AA <- gm$AA + deps * eye(n)
  gm$DD <- gm$DD + deps * eye(n)
  dimnames(gm$AA) <- list(ids,ids)
  dimnames(gm$DD) <- list(ids,ids)
}

# COMPUTE GENOTYPE PROBABILITIES
# ------------------------------
# Compute the conditional genotype probabilities using QTLRel. To
# accomplish this, we need to replace the genotypes with allele counts
# (AA, AB, BB become 1, 2, 3, respectively), and we replace any 
# missing values with zeros.
cat("Calculating probabilities of missing genotypes.\n")
G  <- genotypes2counts(geno)
G  <- zero.na(G)
gp <- combine.genoprob(genoProb(G[F2.rows,],map,step = Inf,
                                method = map.function,gr = 2),
                       genoProb(G[F34.rows,],map,step = Inf,
                                method = map.function,gr = 34))

# Initialize storage for QTL mapping results (gwscan), parameter
# estimates of variance components from QTLRel analysis (vcparams),
# parameter estimates of additive QTL effects for individual markers
# (additive) and dominance QTL effects (dominance), and permutation
# tests (perms).
gwscan    <- list(F2.qtl   = NULL,
                  F2.rel   = empty.scanone(map[F2.markers,]),
                  F34      = empty.scanone(map),
                  combined = empty.scanone(map))
r         <- list(F2       = empty.scanone(map[F2.markers,]),
                  F34      = empty.scanone(map),
                  combined = empty.scanone(map))
additive  <- r
dominance <- r
pve       <- r
vcparams  <- list(F2 = NULL,F34 = NULL,combined = NULL)
perms     <- list(F2 = NULL,F34 = NULL,combined = NULL)

# ANALYSIS WITH F2 CROSS USING QTL
# --------------------------------
cat("QTL MAPPING WITH F2 CROSS USING qtl\n")
out <- map.cross.qtl(pheno,geno,phenotypes,covariates,gp,
                     F2.rows,F2.markers,qtl.method,num.perm,
                     verbose = TRUE)
gwscan$F2.qtl <- out$gwscan
perms$F2      <- out$perms

# PERMUTATION TESTS FOR F34 AND COMBINED CROSSES
# ----------------------------------------------
# Perform permutation tests to calculate thresholds for significance
# for F34 cross, ignoring relatedness between mice.
cat("CALCULATING SIGNIFICANCE THRESHOLDS FOR F34 CROSS USING qtl\n")
perms$F34 <- map.cross.qtl(pheno,geno,phenotypes,covariates,gp,
                           F34.rows,NULL,qtl.method,num.perm,
                           verbose = FALSE)$perms

# Perform a permutation test to calculate thresholds for significance
# for combined F2 and F34 cohort, ignoring relatedness between mice.
cat("CALCULATING SIGNIFICANCE THRESHOLDS FOR COMBINED COHORT USING qtl\n")
perms$combined <- map.cross.qtl(pheno,geno,phenotypes,covariates,gp,
                                qtl.method = qtl.method,num.perm = num.perm,
                                verbose = FALSE)$perms

# Repeat for each phenotype.
for (phenotype in phenotypes) {
  cat("PHENOTYPE =",toupper(phenotype),"\n")
  cols <- c(phenotype,covariates)
  
  # ANALYSIS WITH F2 CROSS USING QTLRel
  # -----------------------------------
  # Only analyze samples (i.e. rows of the genotype and phenotype
  # matrices) for which the phenotype and all the covariates are
  # observed.
  cat("(a) QTL mapping with F2 cross using QTLRel\n")
  F2.rows <- which(pheno$generation == "F2" &
                   none.missing.row(pheno[cols]))
  if (relatedness == "pedigree")
    out <- map.cross.ped(pheno,geno,map,phenotype,covariates,F2.gm,gp,
                         F2.rows,F2.markers,verbose = TRUE)
  else if (relatedness == "markers")
    out <- map.cross.rr(pheno,geno,map,phenotype,covariates,gp, 
                        F2.rows,F2.markers,verbose = TRUE)
  
  # Get the parameter estimates and QTL mapping results.
  gwscan$F2.rel[[phenotype]] <- out$gwscan$lod
  additive$F2[[phenotype]]   <- out$gwscan$additive
  dominance$F2[[phenotype]]  <- out$gwscan$dominance
  pve$F2[[phenotype]]        <- out$gwscan$pve
  vcparams$F2                <- rbind(vcparams$F2,out$vcparams)

  # ANALYSIS WITH F34 CROSS USING QTLRel
  # ------------------------------------
  # Only analyze samples (i.e. rows of the genotype and phenotype
  # matrices) for which the phenotype and all the covariates are
  # observed.
  cat("(b) QTL mapping with F34 cross using QTLRel\n")
  F34.rows <- which(pheno$generation == "F34" &
                    none.missing.row(pheno[cols]))
  if (relatedness == "pedigree")
    out <- map.cross.ped(pheno,geno,map,phenotype,covariates,F34.gm,gp,
                         rows = F34.rows,verbose = TRUE)
  else if (relatedness == "markers")
    out <- map.cross.rr(pheno,geno,map,phenotype,covariates,gp,
                        rows = F34.rows,verbose = TRUE)

  # Get the parameter estimates and QTL mapping results.
  gwscan$F34[[phenotype]]    <- out$gwscan$lod
  additive$F34[[phenotype]]  <- out$gwscan$additive
  dominance$F34[[phenotype]] <- out$gwscan$dominance
  pve$F34[[phenotype]]       <- out$gwscan$pve
  vcparams$F34               <- rbind(vcparams$F34,out$vcparams)

  # ANALYSIS WITH F2 AND F34 SAMPLES USING QTLRel
  # ---------------------------------------------
  cat("(c) QTL mapping with combined F2 and F34 cohorts using QTLRel\n")
  rows <- c(F2.rows,F34.rows)
  if (relatedness == "pedigree")
    out <- map.cross.ped(pheno,geno,map,phenotype,covariates,
                         gm,gp,rows = rows,verbose = TRUE)
  else if (relatedness == "markers")
    out <- map.cross.rr(pheno,geno,map,phenotype,covariates,gp, 
                        rows = rows,verbose = TRUE)

  # Get the parameter estimates and QTL mapping results.
  gwscan$combined[[phenotype]]    <- out$gwscan$lod
  additive$combined[[phenotype]]  <- out$gwscan$additive
  dominance$combined[[phenotype]] <- out$gwscan$dominance
  pve$combined[[phenotype]]       <- out$gwscan$pve
  vcparams$combined               <- rbind(vcparams$combined,out$vcparams)
}

# Add labels to the parameter estimates.
vcparams <- within(vcparams,{
  F2       <- data.frame(F2,check.names = FALSE,row.names = phenotypes)
  F34      <- data.frame(F34,check.names = FALSE,row.names = phenotypes)
  combined <- data.frame(combined,check.names = FALSE,row.names = phenotypes)
})

# SAVE RESULTS TO FILE
# --------------------
cat("Saving results to file.\n")
save(file = resultsfile,
     list = c("relatedness","map.function","num.perm","qtl.method",
              "covariates","phenotypes","gwscan","additive",
              "dominance","pve","vcparams","perms"))
