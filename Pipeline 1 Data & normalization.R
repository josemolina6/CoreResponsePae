# *** Pipeline 1: DATA SOURCE AND NORMALIZATION ***

# Work: A first perturbome of Pseudomonas aeruginosa: Identification of core genes related to multiple perturbations by a machine learning approach
# Developer: Jose A. Molina-Mora, Ph.D.

#Packages: Installation
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#library(BiocManager) 
#BiocManager::install("GEOquery")
#BiocManager::install("affy")

#Packages
library(GEOquery)
library(affy)

# ------------------------------------------------------------------------------
# OPTION 1. DOWNLOAD RAW FILES AND USE "RMA" TO NORMALIZE
# This option was used in the study. Check other options for other applications. 

# Select the working directory, for example:
#setwd("C:/...")

# Download and unzip all files by serie number. Directory is called "data". 
#If download was not possible, a manual step was done to complete the assay.
for (i in c("GSE2430", "GSE3090", "GSE4152", "GSE5443", "GSE7402", "GSE10604", "GSE12738", "GSE13252", "GSE14253" ,"GSE36753")){
  print(i)
  getGEOSuppFiles(i, fetch_files = TRUE, makeDirectory = FALSE)
  untar(paste(i,"_RAW.tar", sep=""), exdir = "data")
}

#List of ".CEL" files and unzip them
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep = "/"), gunzip)
cels = list.files("data/", pattern = "CEL")

# Join data and plot raw values
raw.data =  ReadAffy(verbose = FALSE, filenames = paste("data", cels, sep = "/"))
boxplot(raw.data)

# perform RMA normalization (log2)
data.rma.norm = rma(raw.data)

# Get the expression estimates for each array
rma = exprs(data.rma.norm)

# Take a look at the result (5 rows and 3 columes)
rma[1100:1105, 1:3]

write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t")

# Plot normalizated values
boxplot(rma, col= c("pink"))

# Arrangement: Samples were manually verified and curated if necessary. 
# Metadata were incorporated, including class (control or perturbation), 
# GEO-IDs (serie & sample), and the name of the perturbation.  
# The final archive "Complete_dataset.csv" (with all the information) is available with this pipeline 1.
# See Pipelines 2 and 3 for subsequent steps. 

# ------------------------------------------------------------------------------
# OPTION 2. DOWNLOAD NORMALIZED DATA:

# It is possible to download and join normalized SOFT data (but 
#normalization is independent among series), for example: 

# First, produce an data.array with the first case
gse<-getGEO("GSE2430")
e1<-as.data.frame(exprs(gse$GSE2430_series_matrix.txt.gz))
e1$ID<-rownames(e1)

# For the other assays:
for (i in c("GSE3090", "GSE4152", "GSE5443", "GSE7402", "GSE10604", "GSE12738", "GSE13252", "GSE14253" ,"GSE36753")){
  gse<-getGEO(i)
  e2<-as.data.frame(exprs(gse[[1]]))
  e2$ID<-rownames(e2)  
  e1<-merge(e1, e2, by = "ID")
}

# If we plot them, we saw an independent normalization by series... 
boxplot(log(e1[,-1]), col="gold")
# Thus, RMA was preferred. 

# ------------------------------------------------------------------------------
# OTHER OPTIONS

# Get data by platform
getGEO("GPL84") #plataforma

# Get data by serie
getGEO("GSE5443") #Serie

# ------------------------------------------------------------------------------
# ---------------------THE-END--------------------------------------------------

