# *** Function: Classif_Independent function ***

# Work: A first perturbome of Pseudomonas aeruginosa: Identification of core genes related to multiple perturbations by a machine learning approach
# Developer: Jose A. Molina-Mora, Ph.D.

# This Script is called using Classif_Independent(FData,Falgorithm,Fktop), in which
# The Dataset (FData), the classifier (Falgorithm) and the top K genes (Fktop) must be provided. 
# In the Dataset, the column name with the classes must be called "Class".  

Classif_Independent <- function(FData,Falgorithm,Fktop) { #Initial "F" to identify the variables within the function  

# START HERE if you are studying the pipeline step by step:
# You can run this function with the parameters of Part1 of Pipeline3. Run
# the FData,Falgorithm, and Fktop values.  
   
#Falgorithm=algorithm  
#Fktop=ktop
#FData=Dataset
nc<-dim(FData)[2] # number of columns

# RANDOM data partitioning: 80% training dataset and 20% testing dataset. 
index<-createDataPartition(FData$Class, p=0.80, list = FALSE)

# Create dataset using the index: this will rank all the genes 
Dtraining <- FData[index, ]
Dtesting<- FData[-index,]

# Parameters of the training scheme
control <- trainControl(method="cv", number=10)#, classProbs = TRUE)

# train the model
model <- train(Class~., data=Dtraining, method=Falgorithm, preProcess="scale", trControl=control)

# estimate variable importance
importance <- varImp(model, scale=FALSE)

#Create the index of ranked genes
if(Falgorithm=="rf") {
  IndexRank <-data.frame(sort(importance$importance$Overall, index.return = TRUE, decreasing = TRUE)[2])
  # For RF, it is required to use importance$importance$Overall
  } else {
  IndexRank <-data.frame(sort(importance$importance$Control, index.return = TRUE, decreasing = TRUE)[2])
}

# Index of the topK ranked genes
Indexk<-t(IndexRank[1:Fktop,]) 
#IndexRank[1:4,]
#Indexk[1,]

#Re-build dataset using the ranked topk genes
DtrainingRanked<-Dtraining[,c(Indexk[1,],nc)]
DtestingRanked<-Dtesting[,c(Indexk[1,],nc)]

# Run algorithm using 10-fold cross validation using the topK genes
control <- trainControl(method="cv", number=10)#, classProbs = TRUE)
metric <- "Accuracy"
modelR <- train(Class~., data=DtrainingRanked, method=Falgorithm, metric=metric, trControl=control)

# Assess the model on the testing dataset
predictions <- predict(modelR, DtestingRanked)
confusion<-confusionMatrix(predictions, DtestingRanked$Class)
#confusion$overall
#confusion$byClass

# Metrics of the model with topK genes
Prediction_metrics<-t(as.data.frame(c(confusion$overall,confusion$byClass)))

# List of ranked topK genes
Genesk<-as.data.frame(colnames(DtrainingRanked))
GenesRanked<-t(as.data.frame(Genesk))
#colnames(GenesRanked)<-c(sprintf("Gen%d",1:Fktop),"Class")

# Join all results in a single vector: the 18 metric and the list of topk genes
output<-as.data.frame(c(Prediction_metrics,GenesRanked))

return(output)
}


