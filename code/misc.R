# SUMMARY
# -------
# This file contains a number of useful symbol and function
# definitions that do not fit under any one category. Here is an
# overview of the functions defined in this file:
#
#   is.list.element(x,lst)
#   none.missing.row(d)
#   all.missing.col(d)
#   logit10(x)
#   project.onto.interval(x,a,b)
#   zero.na(A)
#   eye(n)
#   ones(n,m)
#   blockdiag(A,B)
#   offdiag(A)
#   matrix.square(A)
# 

# SYMBOL DEFINITIONS
# ----------------------------------------------------------------------
# Shorthand for machine epsilon.
eps <- .Machine$double.eps

# The chromosomes that make up the mouse genome.
chromosomes <- c(as.character(1:19),"X")

# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Returns true if x is an element of the list.
is.list.element <- function (x, lst)
  is.element(x,names(lst))

# ----------------------------------------------------------------------
# For each row of the matrix or data frame, returns true if all the
# entries in the row are provided (not missing).
none.missing.row <- function (d)
  rowSums(is.na(d)) == 0

# ----------------------------------------------------------------------
# Returns a logical vector such that an entry is TRUE if all the
# entries in the corresponding list element (or column of the table)
# are missing---that is, all the entries are set to NA.
all.missing.col <- function (d)
  sapply(d,function(x) all(is.na(x)))

# ----------------------------------------------------------------------
# The logit function.
logit10 <- function (x)
  log10((x + eps)/(1 - x + eps))

# ----------------------------------------------------------------------
# Projects points to the interval [a,b].
project.onto.interval <- function (x, a, b)
  pmin(b,pmax(a,x))

# ----------------------------------------------------------------------
# Set the missing entries (NA) of a matrix to zero.
zero.na <- function (A) {
  A[is.na(A)] <- 0
  return(A)
}

# ----------------------------------------------------------------------
# Return the n x n identity matrix.
eye <- function (n)
  diag(rep(1,n))

# ----------------------------------------------------------------------
# Create an n x m matrix in which every matrix element is 1.
ones <- function (n, m)
  matrix(1,n,m)

# ----------------------------------------------------------------------
# Creates the block diagonal matrix
#
#   | A 0 |
#   | 0 B |
# 
# from two square matrices, A and B.
blockdiag <- function (A, B) {

  # Create a matrix of zeros to store the result.
  na <- nrow(A)
  nb <- nrow(B)
  n  <- na + nb
  D  <- matrix(0,n,n)

  # Copy the entries from blocks A and B.
  D[1:na,1:na] <- A
  D[1:nb,1:nb] <- B

  return(D)
}

# ----------------------------------------------------------------------
# Return a vector containing the upper triangular entries of A, or the
# off-diagonal entries of symmetric matrix A.
offdiag <- function (A)
  return(A[upper.tri(A)])
  
# ----------------------------------------------------------------------
# Returns the matrix product A*A'.
matrix.square <- function (A)
  return(A %*% t(A))
