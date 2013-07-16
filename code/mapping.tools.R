# SUMMARY
# -------
# Some functions to analyze the QTL experiment data. Here is an
# overview of the functions defined in this file:
# 
#   mean.geno(gp)
#   center.vec(a)
#   center.cols(A)
#   rr.matrix(gp)
#   map.cross.rel(pheno,geno,map,gp,vc,phenotype,
#                 covariates,int.covariates)
#   map.cross.qtl(pheno,geno,phenotypes,covariates,gp,rows,
#                 markers,qtl.method,num.perm,verbose)
#   map.cross.ped(pheno,geno,map,phenotype,covariates,
#                 gm,gp,rows,markers,verbose)
#   map.cross.rr.chr(pheno,geno,map,phenotype,covariates,
#                    chr,gp,int.covariates,verbose)
#   map.cross.rr(pheno,geno,map,phenotype,covariates,gp, 
#                rows,markers,verbose)
# 
library(qtl)
library(QTLRel)

# Return the matrix of expected allele counts given the genotype
# probabilities (the output of function genoProb in the 'QTLRel'
# package). 
mean.geno <- function (gp) {
  p1 <- gp$pr[,2,]   # Probability of AB genotype.
  p2 <- gp$pr[,3,]   # Probability of BB genotype.
  return(p1 + 2*p2)  # Calculate mean genotype.
}

# ----------------------------------------------------------------------
# Subtract the mean from all the entries of the vector 'a' so that
# mean(a) = 0.
center.vec <- function (a)
    return(a - mean(a))

# ----------------------------------------------------------------------
# Subtract the mean of each column of matrix A so that that the
# entries of each column add up to zero.
center.cols <- function (A)
    return(apply(A,2,center.vec))

# ----------------------------------------------------------------------
# Use all the markers to estimate the n x n pairwise relatedness
# matrix, which I define to be 2 times the matrix of kinship
# coefficients (for the definition of kinship coefficients, see Lange,
# 2002). To allow for uncertainty in the genotypes at the markers, I
# compute the expected value of this matrix using the genotype
# probabilities provided as output from the 'genoProb' function in
# QTLRel.
rr.matrix <- function (gp) { 

  # Get the number of samples (n) and the number of markers (p).
  d <- dim(gp$pr)
  n <- d[1]
  p <- d[3]
  
  # Get the genotype probabilities.
  pAA <- gp$pr[,1,]
  pAB <- gp$pr[,2,]
  pBB <- gp$pr[,3,]

  # Get the probability of the genotype AB averaged over all the markers.
  mAB <- matrix(rep(rowMeans(pAB),times = n),n,n)

  # Return the expected value of the pairwise relatedness matrix.
  return(2*(matrix.square(pAA)/p + matrix.square(pBB)/p) 
         + mAB + t(mAB) - matrix.square(pAB)/p)
}

# ----------------------------------------------------------------------
# Analyze the QTL experiment data using 'scanOne' from the 'QTLRel'
# library. The function returns a data frame containing the QTL
# mapping results, specifically the LOD scores and the estimated
# additive and dominance effects corresponding to individual markers.
map.cross.rel <- function (pheno, geno, map, gp, vc, phenotype,
                           covariates, int.covariates = NULL) {
  if (!is.null(int.covariates))
    int.covariates <- pheno[,int.covariates]
  out    <- scanOne(pheno[,phenotype],pheno[,covariates],geno,gp,
                    vc,intcovar = int.covariates,test = "None")
  gwscan <- empty.scanone(map)

  # Get the LOD scores and the proportion of variance explained by each SNP. 
  gwscan$lod <- out$p/(2*log(10))
  gwscan$pve <- out$v/100
  
  # Get the additive and dominance effects of all the SNPs.
  gwscan$additive  <- sapply(out$parameters,"[[","a")
  gwscan$dominance <- sapply(out$parameters,"[[","d")
  
  return(gwscan)
}

# ----------------------------------------------------------------------
# Analyze the QTL experiment data using 'scanone' from the 'qtl'
# library, but with the phenotype and genotype data in the format used
# by QTLRel. If the rows and markers are not specified, all the rows
# and markers are used. If the conditional genotype probabilities are
# not provided, they are computed using the qtl library. The return
# value is a list with two elements: a data frame containing the QTL
# mapping results, specifically the LOD scores (gwscan), and the
# permutation tests (perms).
map.cross.qtl <- function (pheno, geno, phenotypes, covariates, gp = NULL,
                           rows = NULL, markers = NULL, qtl.method = "em",
                           num.perm = 1, verbose = TRUE) {

  # If the rows and markers are not specified, use all the rows and
  # markers.
  if (is.null(rows))
    rows <- 1:nrow(pheno)
  if (is.null(markers))
    markers <- 1:ncol(geno)
  
  # Convert the F2 cross data to the format used by qtl.
  data <- rel2qtl(pheno[rows,],geno[rows,markers],map[markers,])

  # Get the conditional genotype probabilities. If not specified,
  # compute them using qtl.
  if (is.null(gp))
    data <- calc.genoprob(data,step = 0)
  else {
    gp   <- subset.genoprob(gp,rows,markers)
    data <- rel2qtl.genoprob(data,gp) 
  }
  
  # Map QTLs for all phenotypes using qtl.
  if (verbose)
    cat("- Mapping QTLs.\n")
  suppressWarnings(
    gwscan <- scanone(data,pheno.col = phenotypes,
                      addcovar = data$pheno[covariates],
                      model = "normal",method = qtl.method,
                      use = "all.obs"))

  # Perform a permutation test to calculate thresholds for significance.
  if (verbose)
    cat("- Calculating significance thresholds.\n")
  suppressWarnings(
    perms <- scanone(data,pheno.col = phenotypes,
                     addcovar = data$pheno[covariates],
                     model = "normal",method = qtl.method,
                     use = "all.obs",n.perm = num.perm,
                     verbose = FALSE))

  # Return the QTL mapping results and the permutation tests.
  return(list(gwscan = gwscan,perms = perms))
}

# ----------------------------------------------------------------------
# Analyze the QTL experiment data using 'scanOne' from the 'QTLRel'
# library, in which the relatedness of the individuals is estimated
# using the pedigree. If the rows and markers are not specified, all
# the rows and markers are used. The return value is a list containing
# two elements: a data frame containing the QTL mapping results
# (gwscan), and the parameter estimates from the model fitting of the
# variance components (vcparams).
map.cross.ped <- function (pheno, geno, map, phenotype, covariates,
                           gm, gp, rows = NULL, markers = NULL,
                           verbose = TRUE) {
  
  # If the rows and markers are not specified, use all the rows and
  # markers.
  if (is.null(rows))
    rows <- 1:nrow(pheno)
  if (is.null(markers))
    markers <- 1:ncol(geno)
  n <- length(rows)
  
  # Get the QTL experiment data for the selected samples and markers.
  if (verbose)
    cat("- Acquiring",n,"samples and",length(markers),"markers.\n")
  pheno <- pheno[rows,]
  geno  <- geno[rows,markers]
  gp    <- subset.genoprob(gp,rows,markers)
  map   <- map[markers,]
  
  # Estimate the variance components for mice in F2 cross.
  if (verbose)
    cat("- Estimating variance components.\n")
    
  # Get the rows of the genetic matrix corresponding to the samples
  # used in the analysis.
  k <- match(pheno$id,row.names(gm$AA))
    
  # Use the relatedness matrix estimated from the pedigree data to
  # estimate the variance components.
  r <- estVC(pheno[,phenotype],pheno[,covariates],
             v = list(AA = gm$AA[k,k],DD = gm$DD[k,k],
             AD = NULL,HH = NULL,MH = NULL,EE = diag(n)))

  # Once we have the variance components estimated, build a matrix
  # from the variance components.
  vc <- r$par["AA"] * r$v[["AA"]] +
        r$par["DD"] * r$v[["DD"]] +
        r$par["EE"] * r$v[["EE"]]

  # Map QTLs for given phenotype accounting for relatedness among mice
  # in the pedigree.
  if (verbose)
    cat("- Mapping QTLs.\n")
  gwscan <- map.cross.rel(pheno,geno,map,gp,vc,phenotype,covariates)
  
  # Return the QTL mapping results and the parameter estimates
  # corresponding to these variance components.
  return(list(gwscan = gwscan,vcparams = r$par))
}

# ----------------------------------------------------------------------
# Map QTLs for all markers on a single chromosome using 'scanOne' from
# the 'QTLRel' library, in which pairwise relatedness is estimated
# using all markers except the markers on the same chromosome as the
# one being analyzed. The return value is a list containing two
# elements: a data frame containing the QTL mapping results (gwscan),
# and the parameter estimates from the model fitting of the variance
# components (vcparams).
map.cross.rr.chr <- function (pheno, geno, map, phenotype, covariates,
                              chr, gp, int.covariates = NULL,
                              verbose = TRUE) {

  # Get the markers on the chromosome.
  markers <- which(map$chr == chr)
  if (verbose)
    cat("    * Chromosome ",chr," (",length(markers)," markers)\n",sep = "")
  
  # Compute the (expected) relatedness matrix using all markers except
  # the markers on the current chromosome.
  R <- rr.matrix(subset.genoprob(gp,markers = which(map$chr != chr)))
  dimnames(R) <- list(pheno$id,pheno$id)
  
  # Use the relatedness matrix estimated from the marker data to
  # estimate the variance components.
  n <- nrow(pheno)
  r <- estVC(pheno[,phenotype],pheno[,covariates],
             v = list(AA = R,DD = NULL,AD = NULL,HH = NULL,
                      MH = NULL,EE = diag(n)))
    
  # Once we have the variance components estimated, build a matrix
  # from the variance components. Note that the "AA" component here is
  # the n x n relatedness matrix estimated from the marker data (R),
  # and the "EE" component is the n x n identity matrix.
  vc <- r$par["AA"] * r$v[["AA"]] +
        r$par["EE"] * r$v[["EE"]]
    
  # Map QTLs in the chromosome.
  gwscan <- map.cross.rel(pheno,geno[markers],map[markers,],
                          subset.genoprob(gp,markers = markers),
                          vc,phenotype,covariates,
                          int.covariates = int.covariates)

  # Return a list containing the QTL mapping results (gwscan) and the
  # parameter estimates of the variance components (vcparams).
  return(list(gwscan = gwscan,vcparams = r$par))
}

# ----------------------------------------------------------------------
# Analyze the QTL experiment data using 'scanOne' from the 'QTLRel'
# library, in which pairwise relatedness is estimated using all
# selected markers except the markers on the same chromosome as the
# one being analyzed. If rows and markers are not specified, all the
# rows and markers are used in QTL mapping and estimating pairwise
# relatedness. The return value is a list containing two elements: a
# data frame containing the QTL mapping results (gwscan), and the
# parameter estimates from the model fitting of the variance
# components (vcparams).
map.cross.rr <- function (pheno, geno, map, phenotype, covariates, gp, 
                          rows = NULL, markers = NULL, verbose = TRUE) {

  # If the rows and markers are not specified, use all the rows and
  # markers.
  if (is.null(rows))
    rows <- 1:nrows(pheno)
  if (is.null(markers))
    markers <- 1:ncol(geno)

  # Get the set of chromosomes.
  chromosomes <- levels(map$chr)
  
  # Get the QTL experiment data for the selected samples and markers.
  if (verbose)
    cat("- Acquiring",length(rows),"samples and",length(markers),"markers.\n")
  pheno <- pheno[rows,]
  geno  <- geno[rows,markers]
  gp    <- subset.genoprob(gp,rows,markers)
  map   <- map[markers,]
  
  # Compute the (expected) relatedness matrix using all markers, and
  # estimate the variance components using this relatedness
  # matrix. This step is taken so that we can record the fitted
  # parameter estimates for the variance components using the whole
  # genome.
  if (verbose) {
    cat("- Estimating variance components and mapping QTLs.\n")
    cat("    * Whole genome\n")
  }
  R           <- rr.matrix(gp)
  dimnames(R) <- list(pheno$id,pheno$id)
    
  # Use the relatedness matrix estimated from the marker data to
  # estimate the variance components, and store the parameters
  # corresponding to these variance components.
  n        <- nrow(pheno)
  vcparams <- estVC(pheno[,phenotype],pheno[,covariates],
                    v = list(AA = R,DD = NULL,AD = NULL,HH = NULL,
                             MH = NULL,EE = diag(n)))$par

  # Map QTLs separately for each chromosome.
  gwscan <- list()
  for (i in chromosomes)
    gwscan[[i]] <- map.cross.rr.chr(pheno,geno,map,phenotype,covariates,
                                    i,gp,verbose = verbose)$gwscan

  # Merge the QTLRel mapping results.
  r           <- do.call(rbind,gwscan)
  rownames(r) <- do.call(c,lapply(gwscan,rownames))

  # Return the QTL mapping results and parameter estimates
  # corresponding to these variance components.
  return(list(gwscan = r,vcparams = vcparams))
}
