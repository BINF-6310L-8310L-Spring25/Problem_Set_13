# Problem_Set_13

This week we will be conducting a GWAS using rice data. The tutorial is adapted from https://whussain2.github.io/Materials/Teaching/GWAS_R_2.html

## Step 1 - Install packages

We will need the following packages 
```
library(rrBLUP)
library(BGLR)
library(DT)
library(SNPRelate)
library(dplyr)
library(qqman)
library(poolR)
```

To install SNPRelate you will need to use BiocManger

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SNPRelate")
```

More info here: https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html 

To install poolR you will need to pull it from github

```
install.packages("remotes")
remotes::install_github("ozancinar/poolR")
```

More info here: https://rdrr.io/github/ozancinar/poolR/ 

## Step 2 - Read in the data

We have 3 genotype files 
sativas412.fam (names of acessions/genotypes)
sativas412.map (chromosome, position, and marker names)
sativas412.bed (allele data)

The allele data is encoded as follows
1: Heterozygous
0: Homozygous for minor allele
3: Homozygous for major allele
2: Missing data 

### Read in genotype data 

- read in the .ped file using ```read_ped()```
- read in the .fam and .map files as tables


