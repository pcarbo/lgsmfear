# This script maps QTLs on the X chromosome using genotype and
# phenotype data from the combined F2 and F34 crosses. Marker-based
# estimates of pairwise relatedness are used here, so there is no need
# to load the pedigree data.

# SCRIPT PARAMETERS
# -----------------
qtl.method   <- "hk"       # Which QTL mapping method to use in qtl.
map.function <- "Haldane"  # Map function to use for interval mapping.
jitter.amt   <- 1e-6       # Amount by which marker positions are adjusted.

# Save the results to this file.
resultsfile <- "gwscan.X.RData"

# Separately map QTLs for these phenotypes, and use these covariates
# in QTL mapping for all phenotypes. Note that sex is treated as an
# interactive covariate.
phenotypes <- c("pretrainfreeze","freezetocontext","freezetocue",
                "distperiphery","distcenter","propcenter")
covariates <- c("age","albino","agouti")

# Initialize the random number generator.
set.seed(7)

# Load packages and function definitions.
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
# data so that the rows of the tables list the F2 females, F2 males,
# F34 females, and F34 males in that order, from top of the table to
# bottom.
cat("Loading data.\n")
pheno           <- read.phenotypes("../data/pheno.csv")
geno            <- read.genotypes("../data/geno.csv")
rows            <- order(pheno$generation,pheno$sex)
pheno           <- pheno[rows,]
geno            <- geno[rows,]
row.names(geno) <- geno$id
geno            <- geno[-(1:2)]
map             <- read.map("../data/map.csv",chromosomes)
map             <- jitter.gendist(map,jitter.amt)

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
# the male and female mice in the F2 and F34 crosses.
sex         <- pheno$sex
F2.females  <- which(pheno$generation == "F2"  & sex == "F")
F2.males    <- which(pheno$generation == "F2"  & sex == "M")
F34.females <- which(pheno$generation == "F34" & sex == "F")
F34.males   <- which(pheno$generation == "F34" & sex == "M")
F2.rows     <- c(F2.females,F2.males)
F34.rows    <- c(F34.females,F34.males)

# DISCARD MARKERS
# ---------------
# I also remove markers that are *not* genotyped in the F34 cross.
markers <- which(!all.missing.col(geno[F34.rows,]))
geno    <- geno[,markers]
map     <- map[markers,]

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
                         age    = age - mean(age),
                         albino = as.integer(coatcolor == "W"),
                         agouti = as.integer(coatcolor == "A"))

# COMPUTE GENOTYPE PROBABILITIES
# ------------------------------
# Compute the conditional genotype probabilities separately for males
# and females in F2 and F34 crosses. To accomplish this, we need to
# replace the genotypes with allele counts (AA, AB, BB become 1, 2, 3,
# respectively), and we replace any missing values with zeros.
cat("Calculating probabilities of missing genotypes.\n")
G  <- genotypes2counts(geno)
G  <- zero.na(G)
gp <- combine.genoprob(genoProb(G[F2.females,],map,step = Inf,
                                method = map.function,gr = 2),
                       genoProb(G[F2.males,],map,step = Inf,
                                method = map.function,gr = 2))
gp <- combine.genoprob(gp,genoProb(G[F34.females,],map,step = Inf,
                                method = map.function,gr = 34))
gp <- combine.genoprob(gp,genoProb(G[F34.males,],map,step = Inf,
                                method = map.function,gr = 34))

# Initialize storage for QTL mapping results (gwscan) and parameter
# estimates of variance components from QTLRel analysis (vcparams).
gwscan   <- list()
vcparams <- list()

# QTL MAPPING WITH F2 AND F34 SAMPLES USING QTLRel
# ------------------------------------------------
cat("QTL mapping with combined F2 and F34 cohorts using QTLRel.\n")
for (phenotype in phenotypes) {
  cat("PHENOTYPE =",toupper(phenotype),"\n")
  cols <- c(phenotype,covariates)
  rows <- which(none.missing.row(pheno[cols]))
  out  <- map.cross.rr.chr(pheno[rows,],geno[rows,],map,phenotype,
                           covariates,"X",subset.genoprob(gp,rows),
                           int.covariates = "sex",verbose = FALSE)

  # Get the parameter estimates and QTL mapping results.
  gwscan[[phenotype]]   <- out$gwscan
  vcparams[[phenotype]] <- out$vcparams
}

# Turn the list of parameter estimates into a data frame.
vcparams <- do.call(rbind,vcparams)

# SAVE RESULTS TO FILE
# --------------------
cat("Saving results to file.\n")
save(file = resultsfile,
     list = c("map.function","qtl.method","covariates",
              "phenotypes","gwscan","vcparams"))
