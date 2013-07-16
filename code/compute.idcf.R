# This script generates the identity coefficients for the F2 and F34
# crosses.
library(QTLRel)
source("misc.R")
source("read.data.R")

# Load the pedigree and phenotype data.
cat("Loading data.\n")
pheno <- read.phenotypes("../data/pheno.csv")
ped   <- read.pedigree("../data/ped.csv")
ped   <- ped[-(1:2),]
ped   <- rbind(read.pedigree("../data/inbred.ped.csv"),ped)

# Get the ids of the F2 and F34 crosses.
F2.ids  <- subset(pheno,generation == "F2")$id
F34.ids <- subset(pheno,generation == "F34")$id

# Compute the identity coefficients for the F2 cross.
cat("Computing identity coefficients for the F2 cross.\n");
F2.idcf <- cic(ped,ids = F2.ids,df = 100,ask = FALSE,verbose = FALSE)
save(list = c("F2.idcf"),file = "F2.idcf.RData")

# Compute the identity coefficients for the F34 cross.
cat("Computing identity coefficients for the F34 cross.\n");
F34.idcf <- cic(ped,ids = F34.ids,df = 3,ask = FALSE,verbose = TRUE)
save(list = c("F34.idcf"),file = "F34.idcf.RData")

