# load libraries
chooseCRANmirror(ind=73)
# Load packages
Pckg.Lst <-c("vroom","readxl","tidyverse", "bigmemory", "biganalytics", "parallel", "MASS")
package.check <- lapply(
  Pckg.Lst,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)}})

# load gapit functions to convert to numeric format.
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

#
myG = vroom("../VCFfilter/4TradStep6.hmp.txt", col_names=F)

#load phenotype file
pheno <- read.csv("GWA.Phenotype.csv", head=T)
# convert phenotype to dataframe
pheno <- as.data.frame(pheno)
colnames(pheno)[1] <- "taxa"
# subset to traits of interest
pheno <- pheno[,c(1,5,7:10)]

# run GAPIT model
myGAPIT <- GAPIT(
  Y=pheno,
  G=myG,
  PCA.total=1, # previous analyses indicate first PCA explains 10% and the rest < 2%
  model=c("MLM"))