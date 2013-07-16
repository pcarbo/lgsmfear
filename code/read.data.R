# SUMMARY
# -------
# This file contains several functions for reading the QTL experiment
# data from files in comma-delimited ("csv") format. Here is an overview 
# of the functions defined in this file:
#
#   read.pedigree(file)
#   read.phenotypes(file)
#   read.genotypes(file)
#   read.map(file,chromosomes)
#   read.F2.data(file)
#
library(qtl)

# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Returns a data frame containing the pedigree data stored in a CSV file.
read.pedigree <- function (file)
  read.csv(file,comment.char="#",as.is = c(1,4,5),check.names = FALSE)

# ----------------------------------------------------------------------
# Returns a data frame containing the phenotype data stored in a CSV file.
read.phenotypes <- function (file)
  read.csv(file,comment.char="#",as.is = 1,check.names = FALSE)

# ----------------------------------------------------------------------
# Returns a data frame containing the genotype data stored in a CSV file.
read.genotypes <- function (file)
  read.csv(file,comment.char="#",as.is = 1,check.names = FALSE)

# ----------------------------------------------------------------------
# Returns a data frame containing the marker data stored in a CSV file.
# Here I convert the chromosomes to factors manually (as opposed to
# letting the read.csv function handle this) to make sure that the
# chromosomes are ordered properly.
read.map <- function (file, chromosomes)
  within(read.csv(file,comment.char="#",check.names = FALSE,
                  stringsAsFactors = FALSE),
         chr <- factor(chr,chromosomes))

# ----------------------------------------------------------------------
# Read and output the QTL experiment data for the F2 mice.
read.F2.data <- function (file)
  read.cross(format = "csv",file = file,comment.char = "#",
             na.strings = "NA",genotypes = c("AA","AB","BB"))
