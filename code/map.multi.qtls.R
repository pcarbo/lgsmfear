# This script assesses support for multiple QTLs using the combined
# genotype and phenotype data from the F2 and F34 crosses. This script
# requires a list of SNPs that show the strongest support for being
# QTLs (qtls.csv). Here we use marker-based estimates of relatedness,
# so there is no need to load the pedigree data.

# SCRIPT PARAMETERS
# -----------------
qtl.method   <- "hk"       # Which QTL mapping method to use in qtl.
map.function <- "Haldane"  # Map function to use for interval mapping.
jitter.amt   <- 1e-6       # Amount by which marker positions are adjusted.

# Save the results to this file.
resultsfile <- "gwscan.multi.RData"

# Use these covariates in QTL mapping for all phenotypes. Here I am
# conditioning on two additional covariates: the allele count of the
# QTL identified in a prior step of the analysis (additive effect);
# and a binary variable that is 1 if the individual is heterozygous at
# the marker (dominance effect).
covariates <- c("sex","age","albino","agouti","qtla","qtld")

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

# Load the candidate QTLs we will condition on in these analyses.
qtls <- read.csv("../data/qtls.csv",comment.char="#",check.names = FALSE,
                 stringsAsFactors = FALSE)

# DISCARD SAMPLES
# ---------------
# I remove three mice from the F2 cohort because a large proportion of
# their genotypes are not available.
rows  <- which(!is.element(pheno$id,c("648A","738A","811A")))
pheno <- pheno[rows,]
geno  <- geno[rows,]

# Get the rows of the genotype and phenotype matrices corresponding to
# the mice in the F2 and F34 crosses.
F2.rows  <- which(pheno$generation == "F2")
F34.rows <- which(pheno$generation == "F34")

# DISCARD MARKERS
# ---------------
# For now, I drop the X chromosome from the analysis. I also remove
# all markers that are not genotyped in the F34 cross.
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

# Initialize storage for the QTL mapping results.
gwscan   <- list()
vcparams <- list()

# Repeat for each QTL identified in the initial genome-wide scan.
numqtls <- nrow(qtls)
cat("Assessing support for",numqtls,"multi-QTL models.\n")
for (i in 1:numqtls) {

  # Get the phenotype and chromosome (k) to analyze, and the index of
  # the marker that is included as a covariate in the QTL mapping.
  phenotype <- qtls$phenotype[i]
  chr       <- qtls$chr[i]
  j         <- which(map$refSNP == qtls$refSNP[i])
    
  cat(i,". QTL mapping on chromosome ",chr," for phenotype '",
      qtls$phenotype[i],"',\n",sep="")
  cat("conditioned on genotype of SNP",map$refSNP[j],"at",
      format(map$pos[j]/1e6,digits=5),"Mb.\n")

  # Add two phenotype columns, one with the mean allele count of the
  # QTL polymorphism, and another with the probability that the
  # individual is heterozygous at the marker.
  pAB        <- gp$pr[,2,j]
  pBB        <- gp$pr[,3,j]
  pheno$qtla <- pAB + 2*pBB
  pheno$qtld <- pAB
  
  # Map QTLs on chromsome with combined F2 and F34 cohorts using QTLRel.
  cols <- c(phenotype,covariates)
  rows <- which(none.missing.row(pheno[cols]))
  out  <- map.cross.rr.chr(pheno[rows,],geno[rows,],map,phenotype,covariates,
                           chr,subset.genoprob(gp,rows),verbose = FALSE)
  
  # Get the parameter estimates and QTL mapping results.
  gwscan[[i]]   <- out$gwscan
  vcparams[[i]] <- out$vcparams
}

# SAVE RESULTS TO FILE
# --------------------
cat("Saving results to file.\n")
save(file = resultsfile,
     list = c("map.function","qtl.method","covariates","qtls",
              "gwscan","vcparams"))
