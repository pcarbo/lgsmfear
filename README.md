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

this program is distributed in the hope that it will be useful, but
**without any warranty**; without even the implied warranty of
**merchantability** or **fitness for a particular purpose**. See
[LICENSE](LICENSE) for more details.

###Getting started

Things to explain in this section: (1) how to download the code and
data to your computer; (2) which R libraries to install (QTLRel is one
of them, of course); (3) how to run the script to compute QTL mapping
results for all phenotypes; (4) how to display the results.

###Overview of data files

Here is a brief summary of the files in the data directory:

+ **pheno.csv** Phenotype data from 3-day fear conditioning study for
490 mice from the F2 cross, and 687 mice from the F34 cross. Includes
other information such as gender, age and coat colour.

+ **geno.csv** A large table giving the genotypes sampled at 4608
markers (single nucleotide polymorphisms, or SNPs) in the F2 and F34
crosses. Most of the genotypes are marked as missing in the F2 mice
because a subset of only 162 SNPs were genotyped in these mice.

+ **map.csv** 

###Overview of R source code files

The code subdirectory contains several files of interest. Here are the
main ones of interest:

+ **map.qtls.R** Description of this file goes here.

+ **plot.gwscan.R** Description of this file goes here.

###Who

The R code implementing the analysis procedure for this mouse study 
was developed by:<br>
[Peter Carbonetto]((http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics<br>
University of Chicago<br> 
July 2013
