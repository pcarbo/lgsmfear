# QTL mapping of fear conditioning traits in LG/J x SM/J AIL mouse study

###Objectives

Overview of what is contained in repository goes here.

###License

Copyright (c) 2013, Peter Carbonetto

The lgsmfear project repository is free software: you can redistribute
it and/or modify it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html) as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
**without any warranty**; without even the implied warranty of
**merchantability** or **fitness for a particular purpose**. See
[LICENSE](LICENSE) for more details.

###Getting started

Things to explain in this section: (1) how to download the code and
data to your computer; (2) which R libraries to install (QTLRel is one
of them, of course); (3) how to run the script to compute QTL mapping
results for all phenotypes; (4) how to display the results.

###Overview of data files

Here is a brief summary of the files in the [data](data) directory:

+ [pheno.csv](data/pheno.csv) Phenotype data from 3-day fear
conditioning study for 490 mice from the F2 cross, and 687 mice from
the F34 cross. Includes other information such as gender, age and coat
colour.

+ [geno.csv](data/geno.csv) A large table giving the genotypes sampled
at 4608 markers (single nucleotide polymorphisms, or SNPs) in the F2
and F34 crosses. Most of the genotypes are marked as missing in the F2
mice because a subset of only 162 SNPs were genotyped in these mice.

+ [map.csv](data/map.csv) Information about SNPs genotyped in mouse
advanced intercross line. Information for each SNP includes chromosome
number, base pair position on chromosome, refSNP identifier (if
available), and genetic distance estimate.

+ [ped.csv](data/ped.csv) Pedigree data for the mouse advanced
  intercross line.

+ [inbred.ped.csv](data/inbred.ped.csv) Pedigree fragment used to
define inbred founders. QTLRel assumes that the alleles of founders
in the pedigree are not identical by descent (IBD). This pedigree
fragment is added to the pedigree to circumvent this restriction.

+ [qtls.csv](data/qtls.csv) SNPs with the strongest support for being
QTLs based on the initial QTL mapping.

+ [F34.idcf.RData](data/F34.idcf.RData) Identity coefficients for F34
  cross calculated using the **cic** function in QTLRel.

###Overview of R source code files

Here is a brief summary of the files in the [code](code) directory

The code subdirectory contains several files of interest. Here are the
main ones of interest:

+ [map.qtls.R](code/map.qtls.R) Description of this file goes here.

+ [plot.gwscan.R](code/plot.gwscan.R) Description of this file goes here.

+ [compute.idcf.R](code/compute.idcf.R) Description of this file goes here.

+ [map.X.qtls.R](code/map.X.qtls.R) Description of this file goes here.

+ [map.X.perms.R](code/map.X.perms.R) Description of this file goes here.

+ [map.multi.qtls.R](code/map.multi.qtls.R) Description of this file goes here.

+ [read.data.R](code/read.data.R) Description of this file goes here.

+ [data.manip.R](code/data.manip.R) Description of this file goes here.

+ [mapping.tools.R](code/mapping.tools.R) Description of this file goes here.

+ [misc.R](code/misc.R) Description of this file goes here.

###Credits

The R code implementing the analysis procedures was developed by:<br>
[Peter Carbonetto]((http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics<br>
University of Chicago<br> 
July 2013
