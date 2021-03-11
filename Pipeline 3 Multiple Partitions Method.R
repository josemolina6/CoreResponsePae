# *** Pipeline 3: MULTIPLE PARTITIONS (MP) METHOD ***

# Work: A first perturbome of Pseudomonas aeruginosa: Identification of core genes related to multiple perturbations by a machine learning approach
# Developer: Jose A. Molina-Mora, Ph.D.

# For this pipeline, the file "Complete_dataset.csv" and Script "Classif_Independent.R" are required. It must be in the same directory. 
# Set the working directory with the files. 

#setwd("C:/...")

# NOTE: By defaul, the KNN is the selected classifier. 
#See lines 46-48 to change to SVM or RF algorithms. 

# Packages: Installation
#install.packages("ggplot2")
#install.packages("caret")
#install.packages("mlbench")
#install.packages("dplyr")
library(ggplot2)
library(caret)
library(mlbench)
library(dplyr)

#--------------------------------------------------------------
# PART 1. Datasets 

#Load data
preDataset <- read.csv("Complete_Dataset.csv") # It includes metadata
Dataset<-preDataset[,c(8:107,1)] #c(8:5555,1) #Select expression data & class only
dim(Dataset)[2]

# The multiple random partitions are executed using the Script "Classif_Independent.R"
# See the Script for details of the data partitioning.

#Number of replicates to try. In each replicate, a random data partitioning is run. 
kreplicas=100 # In the study we used 100 replicates, to test you can use 10. 

# TopK genes to consider (see Methods, k=56 in this case). 
ktop=56

#Metrics to consider
metrics<-c("Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull", "AccuracyPValue","McnemarPValue","Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision", "Recall","F1","Prevalence" ,"Detection Rate","Detection Prevalence","Balanced Accuracy")

#Classification algorithms: Select one of the algorithms
set.seed(10)
algorithm="knn"
#algorithm="svmRadial"
#algorithm="rf"

# At this point, you can open the Script "Classif_Independent.R" to study the step within. 
# If you just want to run the pipeline, continue with PART 2.

#--------------------------------------------------------------
# PART 2. TRAINING USING MACHINE-LEARNING AND RANKING OF GENES

# The Script "Classif_Independent.R" is used here. The script makes the random partition,
# the training, ranking of genes, selection of ranked topk  genes and assessment 
# using testing dataset for each kreplicates or iteration. 

# Call the function
source("Classif_Independent.R")

# Run the kreplicates for the topk genes with multiple partitions
MPtable <- as.matrix(as.data.frame.list(replicate(kreplicas, Classif_Independent(Dataset,algorithm,ktop))))
rownames(MPtable)<-c(metrics,sprintf("Gen%d",1:ktop),"Class")
colnames(MPtable)<-sprintf("Run%d",1:kreplicas)
#MPtable

# Explore accuracy and kappa values
#hist(as.double(MPtable[1,]), main=paste("Accuracy by",algorithm,"with",kreplicas,"replicates and top genes:", ktop), xlab = "Accuracy", ylab="Frequency")
#hist(as.double(MPtable[2,]), main=paste("Kappa by",algorithm,"with",kreplicas,"replicates" ), xlab = "Valor kappa", ylab="Frecuencia")
#boxplot(as.double(MPtable[1,]), main=paste("Accuracy by",algorithm,"with",kreplicas,"replicates and top genes:", ktop))

# Other exploratory plots
par(mfrow = c(3,6))
for (i in 1:18){
  boxplot(as.double(MPtable[i,]),main=metrics[i], ylim=c(0,1))}

# Select gene names and count them among the replicates
GenesTable<-MPtable[-c(1:18),]
GeneCount<-as.data.frame(sort(table(GenesTable), decreasing = TRUE))
GeneCount[1:10,] # The class is a internal control that must be equal to kreplicates

colores<-c("green","yellow","","","","","","","","","","red","orange","blue")

# Assess metrics among replicates
par(mfrow = c(1,5))
for (i in c(1,2,12,13,14)){
boxplot(as.double(MPtable[i,]),main=metrics[i], ylim=c(0,1), col=colores[i])
}
title(paste("Metrics for", algorithm, "with",kreplicas,"replicates and top genes:", ktop), outer = TRUE)


#--------------------------------------------------------------
# PART 3. DEFINITIVE LISTS AND PLOTS 

#Tables and final list
Genes<-as.data.frame.array(GeneCount)

# List of the consensus topK genes by the classifier after kreplicates
TopGenes<-as.data.frame(Genes[2:{ktop+1},])

write.csv(TopGenes, file = paste("List_genes",algorithm,"-",ktop,"replicates.csv"))

# Reorder the data
colnames(TopGenes)<-c("x","y")
TopGenes=TopGenes %>%
  arrange(y) %>%
  mutate(x=factor(x,x))

# Plot
ggplot(TopGenes, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y), color="black", size=0.3 ) +
  geom_point( color=rgb(43, 107, 225, max=255), size= 1.3 ) +
  theme_classic() +
  coord_flip() +
  theme(
    legend.position="none",
    axis.ticks.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  ) +
  xlab("") +
  ylab("Frequency")

# ------------------------------------------------------------------------------
# ---------------------THE-END--------------------------------------------------



