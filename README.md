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
library(poolr)
```

**Note** in the tutorial they use pool**R** which is the precurser to poolr. I ran into an issue with the lockdir - this was solved by using ```options("install.lock"=FALSE)``` before install

To install SNPRelate you will need to use BiocManger

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SNPRelate")
```

More info here: https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html 

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

### Extract and edit marker data

We need to recode the genome data and turn it into a matrix. If you get stuck see the original tutorial 

- The ped object contains 3 things: the dimenions of the dataset (p and n) and a single list of the genotypes (x)
- Save the genotype data (x) into a new object
- Recode the genotype data such that
- - 2 becomes NA
- - 3 becomes 2
- Create a genotype matrix with a ```nrow=p, ncol=n, byrow=TRUE```
- Transpose the genotype matrix using ```t()```

### Read in the Phenotype Data

We have 1 phenotype file RiceDiversity_44K_Phenotypes_34traits_PLINK.txt
The NFSTV ID is the name of the acessions
There are 36 different traits

-Read in the phenotype data but be sure to use ```stringsAsFactors = False```

### Answer the following questions 

- How many individuals were genotyped?
- How many NAs are in our genotype data?
- How many NAs are in our phenotype data?



## Step 2 Prepare data for GWAS 

We are going to conduct our GWAS on plant height stored in ```Plant.height```

We need to process and prep our data first

- Subsettting: Subset the **phenotype and genotype data** to only include the species that are NOT NA in Plant.height. We are going to assume that the Genotype data is in the **same order** as the phenotype data
- Genotype NA: We still have NA values in our genotype. We want to replace the NAs with the column means. You should implement a loop that loops through all of the columns and replaces the NA with the column mean 
- Rare Alleles: We need to compute the minor allele frequency and then remove the ones that have a maf < 0.05 (aka maf >= 0.05) from the genotype file 
- Subset the Map data to exclude markers with MAF < 0.05 *assume they are in the same order as the gentoype data above*

We have subset the Genotye, Phenotype and Map data

Calculating the maf
```
p <- colSums(Geno)/(2 * nrow(Geno))
maf <- ifelse(p > 0.5, 1 - p, p)
```

### Answer the following questions
- How many individuals are left for analysis
- How many SNPs are left for analysis 


## Step 3 Analyze population structure with a PCA

We are goin gto use the SNPRelate package for this step. Conduct this analysis on our **filtered data**

- Create a gds file using ```snpgdsCreateGeno()``` This requires the options below

1. ```genmat``` - our genotype matrix
2. ```sample.id``` - this is stored in V2 of the Family data
3. ```snp.id``` - this is stored in the V2 of the Map data
4. ```snp.chromosome``` - this is stored in the V1 of the Map data
5. ```snp.position``` - this is stored in the V4 of the Map data
6. ```snpfirstdim``` - set this equal to false 

- Read in the saved gds file
- Summarize the data

A sample of this analysis is below

```
# create gds formate file with marker and sample ids and save it as 44k.gds
snpgdsCreateGeno("44k.gds", genmat = Geno1, sample.id = sample, snp.id = snp.id, 
    snp.chromosome = MAP1$V1, snp.position = MAP1$V4, snpfirstdim = FALSE)
# Now open the 44k.gds file
geno_44k <- snpgdsOpen("44k.gds")
snpgdsSummary("44k.gds")
```

