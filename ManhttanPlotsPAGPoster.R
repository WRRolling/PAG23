# Create manhattan plot
# Author: William Rolling
#####################################################################################
# Section 1.0: Prepare R environment
# 1.1 Clear R environment
rm(list=ls()) 
# 1.2 Set Cran Mirror if any packages need to be installed. 
chooseCRANmirror(ind=73)
# 1.3 list all packages to load or install
Pckg.Lst <-c("vroom",
             "tidyverse",
             "qqman")
#1.4 Install packages. Acknowledgement to https://vbaliga.github.io for install function
package.check <- lapply( 
  Pckg.Lst, 
  FUN = function(x) { 
    if (!require(x, character.only = TRUE)) { 
      install.packages(x, dependencies = TRUE) 
      library(x, character.only = TRUE)}})
#1.5 Create/Set Working Directory - project specific
#####################################################################################
# Section 2.0 Prepare data for qqman
#2.1 Find File of interest
File.Name <- "GAPIT.MLM.Alpha.Beta.GWAS.Results.csv"
Location <- list.files("~", File.Name, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
#2.2 Read in Genome-wide analyses output
Results <- read_csv(Location[1], col_names = T)
#2.3 Format Input data for qqman
  # Only need four columns of data including SNP name, physical position and significance
qqdat <- subset(Results, select = c(SNP, Chromosome, Position, P.value ))
# 2.4 Rename columns for qqman package
colnames(qqdat) <- c("SNP", "CHR", "BP", "P")
#####################################################################################
# Section 3: Highlight regions putatively invovled in Carotenoid accumulation
# 3.1 Identify SNPs at Or locus
Chr03 <- which(qqdat$CHR == 3)
Chr03 <- qqdat[Chr03,]
Ordf <- Chr03[which(Chr03$BP > 5000000 & Chr03$BP < 5200000),]
# 3.2 Identify SNPs at Y locus
Chr05 <- which(qqdat$CHR == 5)
Chr05 <- qqdat[Chr05,]
Ydf <- Chr05[which(Chr05$BP > 29935615 & Chr05$BP < 30065947),]   
# 3.3 Identify SNPs at Y2 Locus
Chr07 <- which(qqdat$CHR == 7)
Chr07 <- qqdat[Chr07,]
#3.4 Combine & Format SNPs at known loci
Y2df <- Chr07[which(Chr07$BP > 38656562 & Chr07$BP < 39560030),]   
OrSNPs <- as.character(Ordf$SNP)
YSNPs <- as.character(Ydf$SNP)
Y2SNPs <- as.character(Y2df$SNP)
KnownLoci <- c(OrSNPs, YSNPs, Y2SNPs)    
#####################################################################################
# Section 5: Plot Results of GWA with GBS data
# 5.1 Set a custom color to use in plot - USDA Green
Color1 <- rgb(red = 0, green = 0.349, blue = 0.255)
# 5.2 Set a custom color to use in plot -  USDA Blue
Color2 <- rgb(red = 0, green = 0.321, blue = 0.616)
# 5.3 fix(manhattan) can highlight SNPs a different color but you have to update script
HighlightUpdate <- readline(prompt="Do you want bright green SNPs of interest?")
if (HighlightUpdate == "Yes") {
    print("I Like Color Too!")
  } else {
    print ("Stop & perform fix(manhattan) function")
  }

# 5.4 Draw Manhattan plot
png(filename="GBS.Base.png", width =800, height = 400)
manhattan(qqdat, main="Genotype by Sequencing", ylim=c(0,12), cex=0.8, 
          cex.axis=0.9, col=c(Color1, Color2), suggestiveline=F, genomewideline=6.230954, 
          chrlabs=c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09"),
          highlight = KnownLoci)
dev.off()
  # set y axis limits with ylim
  # set size of manhattan "dots" with cex
  # Set your colors - I used custom colors I found in a USDA logo
  # Add significance threshold. Mine came from simpleM - Modified Boneferroni 
  # Label you chromosomes as you want
  # highlight SNPs at known loci
#####################################################################################
# Section 6.0 Find results of GWAS with whole genome resequencing data
# 6.1 Find File of interest
File.Name <- "GAPIT.Association.GWAS_Results.MLM.Ratio.alpha_betaRSQ.csv"
Location <- list.files("~", File.Name, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
# 6.2 Read in Genome-wide analyses output
Results <- read_csv(Location[1], col_names = T)
# 6.3 Format Input data for qqman
# Only need four columns of data including SNP name, physical position and significance
qqdat <- subset(Results, select = c(SNP, Chromosome, Position, P.value ))
# 6.4 Rename columns for qqman package
colnames(qqdat) <- c("SNP", "CHR", "BP", "P")
#####################################################################################
# Section 7: New Genotype file means we have to reidentify SNPS at Known Loci 
# Or locus
Chr03 <- which(qqdat$CHR == 3)
Chr03 <- qqdat[Chr03,]
Ordf <- Chr03[which(Chr03$BP > 5000000 & Chr03$BP < 5200000),]
# Y locus
Chr05 <- which(qqdat$CHR == 5)
Chr05 <- qqdat[Chr05,]
Ydf <- Chr05[which(Chr05$BP > 29935615 & Chr05$BP < 30065947),]   
# Y2 Locus
Chr07 <- which(qqdat$CHR == 7)
Chr07 <- qqdat[Chr07,]
Y2df <- Chr07[which(Chr07$BP > 38656562 & Chr07$BP < 39560030),]   
OrSNPs <- as.character(Ordf$SNP)
YSNPs <- as.character(Ydf$SNP)
Y2SNPs <- as.character(Y2df$SNP)
KnownLoci <- c(OrSNPs, YSNPs, Y2SNPs)    

#####################################################################################
# Section 8: Draw Manhattan Plot
png(filename="RSeq.png", width =800, height = 400)
manhattan(qqdat, main="Whole Genome Resequencing", ylim=c(0,12), cex=0.8, 
          cex.axis=0.9, col=c(Color1, Color2), suggestiveline=F, genomewideline=7.863069706, 
          chrlabs=c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09"),
          highlight = KnownLoci)
dev.off()
#####################################################################################
# Switch to ploting results when genotype files and phenotype files are switched 
#####################################################################################
#####################################################################################
# Section 9.0 Find data with 2018 GBS data and phenotypic data from 2019
  # Same Genotype as used in step 5 above
# 9.0 Set a custom colors to use in plot - Carrot Orange & USDA Blue
Color1 <- rgb(red = 0.929, green = 0.49, blue = 0.192)
Color2 <- rgb(red = 0, green = 0.321, blue = 0.616)
#9.1 Find File of interest
File.Name <- "2018GBS.2019Pheno"
Location <- list.files("~", File.Name, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
#9.2 Read in Genome-wide analyses output
Results <- read_csv(Location[1], col_names = T)
#9.3 Format Input data for qqman
# Only need four columns of data including SNP name, physical position and significance
qqdat <- subset(Results, select = c(SNP, Chromosome, Position, P.value ))
#9.4 Rename columns for qqman package
colnames(qqdat) <- c("SNP", "CHR", "BP", "P")
#####################################################################################
# Section 9: Highlight regions putatively invovled in Carotenoid accumulation
# Or locus
Chr03 <- which(qqdat$CHR == 3)
Chr03 <- qqdat[Chr03,]
Ordf <- Chr03[which(Chr03$BP > 5000000 & Chr03$BP < 5200000),]
# Y locus
Chr05 <- which(qqdat$CHR == 5)
Chr05 <- qqdat[Chr05,]
Ydf <- Chr05[which(Chr05$BP > 29935615 & Chr05$BP < 30065947),]   
# Y2 Locus
Chr07 <- which(qqdat$CHR == 7)
Chr07 <- qqdat[Chr07,]
Y2df <- Chr07[which(Chr07$BP > 38656562 & Chr07$BP < 39560030),]   
OrSNPs <- as.character(Ordf$SNP)
YSNPs <- as.character(Ydf$SNP)
Y2SNPs <- as.character(Y2df$SNP)
KnownLoci <- c(OrSNPs, YSNPs, Y2SNPs)    

#####################################################################################
# Section 9: Plot Manhattan Plot
png(filename="2018GBS.2019Pheno.png", width =600, height = 500 )
manhattan(qqdat, main="2019 Phenotype - 2018 Genotype by Sequencing", ylim=c(0,12), cex=0.8, 
          cex.axis=0.9, col=c(Color1, Color2), suggestiveline=F, genomewideline=6.230954, 
          chrlabs=c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09"),
          highlight = KnownLoci)
dev.off()

#####################################################################################
# Section  10.0 repeat the above ploting for using whole genome resequencing w/ 2019 phenotype
#10.1 Find File of interest
File.Name <- "RSQGeno.2019Pheno"
Location <- list.files("~", File.Name, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
#10.2 Read in Genome-wide analyses output
Results <- read_csv(Location[1], col_names = T)
#10.3 Format Input data for qqman
# Only need four columns of data including SNP name, physical position and significance
qqdat <- subset(Results, select = c(SNP, Chromosome, Position, P.value ))
# 10.4 Rename columns for qqman package
colnames(qqdat) <- c("SNP", "CHR", "BP", "P")
# 10.5 reidentify SNPs in known loci
Chr03 <- which(qqdat$CHR == 3)
Chr03 <- qqdat[Chr03,]
Ordf <- Chr03[which(Chr03$BP > 5000000 & Chr03$BP < 5200000),]
# Y locus
Chr05 <- which(qqdat$CHR == 5)
Chr05 <- qqdat[Chr05,]
Ydf <- Chr05[which(Chr05$BP > 29935615 & Chr05$BP < 30065947),]   
# Y2 Locus
Chr07 <- which(qqdat$CHR == 7)
Chr07 <- qqdat[Chr07,]
Y2df <- Chr07[which(Chr07$BP > 38656562 & Chr07$BP < 39560030),]   
OrSNPs <- as.character(Ordf$SNP)
YSNPs <- as.character(Ydf$SNP)
Y2SNPs <- as.character(Y2df$SNP)
KnownLoci <- c(OrSNPs, YSNPs, Y2SNPs)  
# 10.6 Plot Manhattan Figure
png(filename = "2019Pheno.2018RSQ.png", width =600, height = 500)
manhattan(qqdat, main="2019 Phenotype - 2018 Whole Genome Resequencing", ylim=c(0,12), cex=0.8, 
          cex.axis=0.9, col=c(Color1, Color2), suggestiveline=F, genomewideline=7.863069706, 
          chrlabs=c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09"),
          highlight = KnownLoci)
dev.off()
#####################################################################################
# Section 11.0 Final Figure of using 2019 GBS data with 2018 phenotype
# 11.1 find file of interest
File.Name <- "2019Geno.2018Pheno"
Location <- list.files("~", File.Name, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
#11.2 Read in Genome-wide analyses output
Results <- read_csv(Location[1], col_names = T)
#11.3 Format Input data for qqman
# Only need four columns of data including SNP name, physical position and significance
qqdat <- subset(Results, select = c(SNP, Chromosome, Position, P.value ))
# 11.4 Rename columns for qqman package
colnames(qqdat) <- c("SNP", "CHR", "BP", "P")
# 11.5 Reidentify SNPs at Known Loci
Chr03 <- which(qqdat$CHR == 3)
Chr03 <- qqdat[Chr03,]
Ordf <- Chr03[which(Chr03$BP > 5000000 & Chr03$BP < 5200000),]
# Y locus
Chr05 <- which(qqdat$CHR == 5)
Chr05 <- qqdat[Chr05,]
Ydf <- Chr05[which(Chr05$BP > 29935615 & Chr05$BP < 30065947),]   
# Y2 Locus
Chr07 <- which(qqdat$CHR == 7)
Chr07 <- qqdat[Chr07,]
Y2df <- Chr07[which(Chr07$BP > 38656562 & Chr07$BP < 39560030),]   
OrSNPs <- as.character(Ordf$SNP)
YSNPs <- as.character(Ydf$SNP)
Y2SNPs <- as.character(Y2df$SNP)
KnownLoci <- c(OrSNPs, YSNPs, Y2SNPs)  
# 11.6 Create Manhattan Plot
png(filename="2018Pheno.2019Geno.png",  width =600, height = 500)
manhattan(qqdat, main="2018 Pheno - 2019 Genotype by Sequencing", ylim=c(0,12), cex=0.8, 
          cex.axis=0.9, col=c(Color1, Color2), suggestiveline=F, genomewideline=6.2, 
          chrlabs=c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09"),
          highlight = KnownLoci)
dev.off()
# The end :)
quit()

