# SUMMARY
# -------
# Some functions for manipulating QTL experiment data. Here is an
# overview of the functions defined in this file:
#
#   factor2integer(x)
#   genotypes2counts(d)
#   counts2genotypes(d)
#   rel2qtl(pheno,geno,map)
#   qtl2rel(cross)
#   combine.genoprob(gp1,gp2)
#   rel2qtl.genoprob(cross,gp) 
#   qtl2rel.genoprob(cross)
#   jitter.gendist(map,j)
#   empty.scanone(map,label)
#   subset.genoprob(gp,samples,markers)
#
library(abind)
library(plyr)
library(qtl)
library(QTLRel)

# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Convert a factor to a vector of integers.
factor2integer <- function (x)
  match(x,levels(x))

# ----------------------------------------------------------------------
# Convert the genotypes to a matrix of allele counts.
genotypes2counts <- function (d)
  as.matrix(data.frame(lapply(d,factor2integer),
                       row.names = rownames(d),
                       check.names = FALSE))

# ----------------------------------------------------------------------
# Convert the allele counts to genotypes (stored in a data frame).
counts2genotypes <- function (d) {
  f <- function (x) factor(x,labels = c("AA","AB","BB"))
  return(data.frame(lapply(as.data.frame(d),f),
                    row.names = dimnames(d)[[1]],
                    check.names = FALSE))
}

# ----------------------------------------------------------------------
# Convert the QTL experiment data from the format used in the QTLRel
# library to the format used in qtl library. The return value is a
# 'cross' object that keeps track of all the data in a single QTL
# experiment; for more details, see the help for function read.cross
# in the qtl library.
#
# Here I assume that the alleles are A and B and the genotypes are
# labeled as AA, AB and BB. qtl requires a paternal grandmother
# ("pgm") phenotype, so a column in the table for this phenotype is
# included if it is missing, and the entries of this column are set to
# zero.
#
# Note that I currently do not properly handle the X chromosome.
rel2qtl <- function (pheno, geno, map) {

  # Get the set of chromosomes.
  chromosomes <- levels(map$chr)
  
  # Convert the genotypes to a matrix of integers.
  geno <- genotypes2counts(geno)
  
  # Add the paternal grandmother ("pgm") phenotype if it does not
  # already exist.
  if (!is.list.element("pgm",pheno))
    pheno <- cbind(pheno,pgm = 0)
  
  # Initialize the cross object, setting the genotypes to an empty list.
  cross <- list(geno = list(),pheno = pheno)
  class(cross) <- c("f2","cross")

  # Set the alleles.
  attributes(cross)$alleles <- c("A","B")
  
  # Split the markers by chromosome.
  for (i in chromosomes) {

    # Get the markers on the chromosome.
    markers <- which(map$chr == i)
    
    # Get the genotype and map data corresponding to these markers.
    # Note that I need to coerce the genotype data to be a matrix in
    # the exceptional case when there is only a single marker
    # genotyped on the chromosome.
    m                <- map[markers,]
    d                <- list(data = as.matrix(geno[,markers]),map = m$dist)
    colnames(d$data) <- m$snp
    names(d$map)     <- m$snp

    # I need to treat the X chromosome a little differently.
    if (i == "X")
      class(d) <- "X"
    else
      class(d) <- "A"

    # Store the genotypes for markers on the chromosome.
    cross$geno[[i]] <- d
  }
  
  return(cross)
}

# ----------------------------------------------------------------------
# Convert the QTL experiment data from the format used in the qtl
# library to the format used in QTLRel library. The return value is a
# list with two elements, 'pheno' and 'geno', containing the phenotype
# and genotype data, respectively.
#
# Note that I currently do not properly handle the X chromosome.
qtl2rel <- function (cross) {

  # Get the set of chromosomes.
  chromosomes <- names(cross$geno)
  
  # Merge the genotypes into a single matrix.
  geno <- lapply(cross$geno,function (x) x$data)
  geno <- do.call(cbind,geno)
  
  # Convert the genotypes to factors.
  geno <- counts2genotypes(geno)

  # Return a list with two elements: one containing the phenotype
  # data, and another containing the genotype data.
  return(list(pheno = cross$pheno,geno = geno))
}

# ----------------------------------------------------------------------
# Combine the genotype probabilities from two collections of samples
# (individuals) into a single, larger collection of genotype
# probabilities. This involves combining over the arrays of genotype
# probabilities (the 'pr' list element) along their rows; each row
# corresponds to a sample.
combine.genoprob <- function (gp1, gp2) {
  gp    <- gp1
  gp$pr <- abind(gp1$pr,gp2$pr,along = 1)
  return(gp)
}

# ----------------------------------------------------------------------
# Convert the genotype probabilities from the format used by QTLRel to
# the format used by qtl, and store the genotype probabilities in the
# cross object. Note that I currently do not properly handle the X
# chromosome.
rel2qtl.genoprob <- function (cross, gp) {

  # Transpose the array.
  prob <- aperm(gp$pr,c(1,3,2))

  # Get the set of chromosomes.
  chromosomes <- unique(gp$chr)

  # Split the genotype probabilities by chromosome.
  for (i in chromosomes) {

    # Get the data for all markers on this chromosome.
    geno <- cross$geno[[i]]
    
    # Get the markers on the chromosome.
    markers <- which(gp$chr == i)

    # Get the genotype probabilities for these markers.
    geno$prob                 <- prob[,markers,]
    attributes(geno$prob)$map <- geno$map
    dimnames(geno$prob)       <- list(NULL,
                                      names(geno$map),
                                      c("AA","AB","BB"))

    # Store the new genotype data.
    cross$geno[[i]] <- geno
  }
  
  return(cross)
}

# ----------------------------------------------------------------------
# Convert the genotype probabilities from the format used by qtl to
# the format used by QTLRel. Note that I currently do not properly
# handle the X chromosome.
qtl2rel.genoprob <- function (cross) {

  # This is the return value.
  gp <- list()

  # Get the SNP labels and genetic distances.
  gp$dist        <- do.call(c,lapply(cross$geno,function (x) x$map))
  gp$snp         <- do.call(c,lapply(cross$geno,function (x) names(x$map)))
  names(gp$dist) <- NULL

  # Get the set of chromosomes.
  chromosomes <- names(cross$geno)

  # Get the chromosome numbers.
  chr <- list()
  for (i in chromosomes) {
    numsnps  <- length(cross$geno[[i]]$map)
    chr[[i]] <- rep(i,times = numsnps)
  }
  gp$chr <- do.call(c,chr)
  names(gp$chr) <- NULL

  # Get the genotype probabilities.
  gp$pr <- aperm(abind(lapply(cross$geno,function (x) x$prob),along = 2),
                 c(1,3,2))
  
  return(gp)
}

# ----------------------------------------------------------------------
# Adjust the map positions (i.e. genetic distances) slightly (by the
# amount j) so that no two markers have the same position.
jitter.gendist <- function (map, j) {

  # Get the set of chromosomes.
  chromosomes <- levels(map$chr)

  # Repeat for each chromosome.
  for (i in chromosomes) {

    # Get the markers on the chromosome.
    markers  <- which(map$chr == i)
    n        <- length(markers)

    # Adjust the genetic distances slightly according to scalar j.
    map[markers,"dist"] <- map[markers,"dist"] + cumsum(rep(j,times = n))
  }

  return(map)
}

# ----------------------------------------------------------------------
# This function creates a "scanone" object that will be used to store
# QTL mapping results based on the SNP data provided as input. (For
# further info, see the 'scanone' function in the 'qtl' package.) 
empty.scanone <- function (map) {
  d            <- map[c("chr","dist")]
  row.names(d) <- map[["snp"]]
  names(d)     <- c("chr","pos")
  class(d)     <- c("scanone","data.frame")
  return(d) 
}

# ----------------------------------------------------------------------
# Get the genotype probabilities for the specified samples and markers.
subset.genoprob <- function (gp, samples = NULL, markers = NULL) {

  # If the set of markers is the null object, use all the markers.
  if (is.null(markers)) {
    n       <- length(gp$chr)
    markers <- 1:n
  }

  # If the set of samples is the null object, use all the samples.
  if (is.null(samples)) {
    n       <- nrow(gp$pr)
    samples <- 1:n
  }
  
  class(gp) <- c("Pr","list")
  return(within(gp,{
    pr   <- pr[samples,,markers]
    chr  <- chr[markers]
    dist <- dist[markers]
    snp  <- snp[markers]
  }))
}
