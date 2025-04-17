# Problem_Set_13

This week we will be conducting a GWAS using rice data. The tutorial is adapted from https://whussain2.github.io/Materials/Teaching/GWAS_R_2.html

All files you need are located in the Github repository OR the Canvas Module

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

&nbsp;
&nbsp;

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

&nbsp;

### Read in genotype data 

- read in the .ped file using ```read_ped()```
- read in the .fam and .map files as tables

&nbsp;

### Extract and edit marker data

We need to recode the genome data and turn it into a matrix. If you get stuck, see the original tutorial. 

- The ped object contains 3 things: the dimensions of the dataset (p and n) and a single list of the genotypes (x)
- Save the genotype data (x) into a new object
- Recode the genotype data such that
- - 2 becomes NA
- - 3 becomes 2
- Create a genotype matrix with a ```nrow=p, ncol=n, byrow=TRUE```
- Transpose the genotype matrix using ```t()```

&nbsp;

### Read in the Phenotype & Population Data

We have 1 phenotype file RiceDiversity_44K_Phenotypes_34traits_PLINK.txt
The NFSTV ID is the name of the acessions
There are 36 different traits

-Read in the phenotype data but be sure to use ```stringsAsFactors = False```

We also have population data saved in RiceDiversity.44K.germplasm.csv

- Read in the population data - this is **comma** delimited

&nbsp;

# Question 2A
- How many individuals were genotyped?

# Question 2B
- How many NAs are in our genotype data?

# Question 2C
- How many NAs are in our phenotype data?

&nbsp;
&nbsp;

## Step 3 Prepare data for GWAS 

We are going to conduct our GWAS on plant height stored in ```Plant.height```

We need to process and prep our data first

- Subsetting: Subset the **phenotype, fam, and genotype data** to only include the species that are NOT NA in Plant.height. We are going to assume that the Genotype and FAM data is in the **same order** as the phenotype data
- Genotype NA: We still have NA values in our genotype. We want to replace the NAs with the column means. You should implement a loop that loops through all of the columns and replaces the NA with the column mean 
- Rare Alleles: We need to compute the minor allele frequency and then remove the ones that have a maf < 0.05 (aka maf >= 0.05) from the genotype file 
- Subset the Map data to exclude markers with MAF < 0.05 *assume they are in the same order as the genotype data above*

We have subset the Genotye, Phenotype, Fam, and Map data

Calculating the maf
```
p <- colSums(Geno)/(2 * nrow(Geno))
maf <- ifelse(p > 0.5, 1 - p, p)
```

# Question 3A
- How many individuals are left for analysis

# Question 3B
- How many SNPs are left for analysis 

&nbsp;

## Step 4 Analyze population structure with a PCA

We are going to use the SNPRelate package for this step. Conduct this analysis on our **filtered data**

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
- Conduct a PCA using ```snpgsdPCA()``` example: ```pca <- snpgdsPCA(geno_44k, snp.id = colnames(Geno1))```
- Save your PCA data in a dataframe and visualize it using code similar to below

```
pca <- data.frame(sample.id = edited_fam_data$V2, EV1 = pca$eigenvect[, 1], EV2 = pca$eigenvect[, 
    2], EV3 = pca$eigenvect[, 3], EV3 = pca$eigenvect[, 4], stringsAsFactors = FALSE)
# Plot the PCA
plot(pca$EV2, pca$EV1, xlab = "eigenvector 3", ylab = "eigenvector 4")
```

- Now we are going to add the population data to the data frame containing our PCA data. The samples *are not in the same order* this time so we will use ```full_join()``` from the dplyr package. See below for example

```pca<-left_join(pca, Pop, by=c("sample.id"="NSFTV.ID"))```

- Finally, visualize the PCA with the population colored using the code below or ggplot. 

```
plot(pca$EV1, pca$EV2, xlab = "PC1", ylab = "PC2", col = c(1:6)[factor(pca$Sub.population)])
```

&nbsp;

# Question 4A

Which population has the most variation on the PCA? 

# Question 4B

Which two populations are distinguished on the second PC/Eigenvector? 

&nbsp;
&nbsp;

## Step 5 - Conduct the GWAS

The GWAS function requires two data frames that are in the format below:

```
pheno	
Data frame where the first column is the line name (gid). The remaining columns can be either a phenotype or the levels of a fixed effect. Any column not designated as a fixed effect is assumed to be a phenotype.

geno	
Data frame with the marker names in the first column. The second and third columns contain the chromosome and map position (either bp or cM), respectively, which are used only when plot=TRUE to make Manhattan plots. If the markers are unmapped, just use a placeholder for those two columns. Columns 4 and higher contain the marker scores for each line, coded as {-1,0,1} = {aa,Aa,AA}. Fractional (imputed) and missing (NA) values are allowed. The column names must match the line names in the "pheno" data frame.
```

- Format a final phenotype and genotype data table. I have included an example below 


```
ph_geno.final<-data.frame(marker = edited_map[, 2], chrom = edited_map[, 1], pos = edited_map[,4], t(edited_geno - 1), check.names = FALSE) 

colnames(ph_geno.final)<-c("marker","chrom","pos",edited_family$V2)

ph_pheno.final<-as.data.frame(edited_phenotype[,c("NSFTVID","Plant.height")])

```

- Run the GWAS! This may take a second to run

```GWAS <- GWAS(ph_pheno.final, ph_geno.final, min.MAF = 0.05, P3D = TRUE, plot = TRUE)```


## Step 6 - Multiple Test Correction & Visualization 

The multiple-test correction method takes a while to compute. So, instead, we will use the MTC cutoff calculated in the tutorial. 

- Use the following commands to access the GWAS results ```GWAS_1 <- GWAS %>% filter(Plant.height != "0")```
- Access the list of significant snps using ```GWAS_1 %>% filter(Plant.height < 1e-04)```
- Visualize the GWAS using the command below

```
manhattan(x = GWAS_1, chr = "chrom", bp = "pos", p = "Plant.height", snp = "marker", col = c("blue4", 
    "orange3"), suggestiveline = -log10(1e-04), logp = TRUE)
```
    
Answer the questions below

# Question 6A
- How many loci pass our MTC and are significantly associated with Plant height?

# Question 6B
- What gene is associated with the variant found on chromosome 8? Use the Rice Genome Browser to find the gene (http://rice.uga.edu/cgi-bin/gbrowse/rice/) that contains the associated variant found on Chromosome 8 (note the format for searching a landmark region MUST follow a strict format. If the associated locus was on chromosome 2 and position 145 you would search for "Chr2:145..146". If you do not have a variant on chromosome 8, report any variant and include the variant position. 



