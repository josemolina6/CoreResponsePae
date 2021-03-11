# *** Pipeline 2: SINGLE PARTITION (SP) METHOD ***

# Work: A first perturbome of Pseudomonas aeruginosa: Identification of core genes related to multiple perturbations by a machine learning approach
# Developer: Jose A. Molina-Mora, Ph.D.

# For this pipeline, the file "Complete_dataset.csv" is required. It must be in the same directory. 
# Set the working directory with the Dataset and the Pipeline 2. 

#setwd("C:/...")

# NOTE: By defaul, the KNN is the selected classifier. 
# See lines 39-41 to change to SVM or RF algorithms. 

# Packages: Installation
#install.packages("caret")
#install.packages("ROCR")
#install.packages("rpart")
#install.packages("rattle")

# Packages
library(caret)
library(ROCR)
library(rpart)
library(rattle)

#--------------------------------------------------------------
# PART 1. Datasets and partitions 

#Load data
Dataset <- read.csv("Complete_Dataset.csv")

# Make a single partition, using the the pre-asigned value in "Set" column. 
# Sets includes only the expression values and the class (control or perturbation)
Dtraining<-Dataset[Dataset$Set=="Training",c(8:5555,1)]
Dtesting<-Dataset[Dataset$Set=="Testing",c(8:5555,1)]

# Distribution of samples by class and dataset
x<-prop.table(table(Dataset$Class, Dataset$Set),2)
barplot(x,col = c("lightblue","gold"),legend = TRUE, ylim = c(0,1.7), args.legend = list(bty = "n", x = "topright", ncol = 1))

#--------------------------------------------------------------
# PART 2. TRAINING USING MACHINE-LEARNING AND RANKING OF GENES

# Parameters of the training scheme: Run algorithms using 10-fold cross validation
control <- trainControl(method="cv", number=10, classProbs=TRUE)
metric <- "Accuracy"

#Classification algorithms: Select one of the algorithms
set.seed(10)
algorithm="knn"
#algorithm="svmRadial"
#algorithm="rf"

# Training step and results 
model <- train(Class~., data=Dtraining, method=algorithm, metric=metric, trControl=control)
model$method
model$results

# Genes ranking by "importance"
importance <- varImp(model, scale=FALSE)
head(importance)

#Create index to sort genes by importance:
# For SMV and KNN use importance$importance$Control, for RF use importance$importance$Overall

if(algorithm=="rf") {
  IndexRank <-t(data.frame(sort(importance$importance$Overall, index.return = TRUE, decreasing = TRUE)[2]))
} else {
  IndexRank <-t(data.frame(sort(importance$importance$Control, index.return = TRUE, decreasing = TRUE)[2]))
}

# Apply the index among the dataset (to sort by importance)
DtrainingRanked<-Dtraining[,c(IndexRank[1,],5549)]
DtestingRanked<-Dtesting[,c(IndexRank[1,],5549)]

#Verify the order
colnames(DtrainingRanked)[1:10]
colnames(DtestingRanked)[1:10]

# Plot of the top 6 genes: box and whisker plots for each attribute
featurePlot(x=DtrainingRanked[,c(4,5,6,1,2,3)], y=DtrainingRanked[,5549], plot="box", main = paste("Top 6 genes with algorithm", algorithm), col = c("red","green"))

#--------------------------------------------------------------
# PART 3. ASSESSING 200 MODELS USING THE TOP 200 RANKED GENES 

# Using the ranked genes, evaluate the performance of the 200 top genes
kvalue=200 # Number of top genes to consider

kEvaluation<-matrix(,nrow=kvalue, ncol=19) # 18 metrics + one for the K value
colnames(kEvaluation)<-c("Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull", "AccuracyPValue","McnemarPValue","Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision", "Recall","F1","Prevalence" ,"Detection Rate","Detection Prevalence","Balanced Accuracy", "K")

# Training 200 models, each time adding one gene (from top1 to top200)
for (k in 1:kvalue){ 
  DtrainK<-DtrainingRanked[,c(1:k,5549)]
  DtestK<-DtestingRanked[,c(1:k,5549)]
  fit.algorK <- train(Class~., data=DtrainK, method=algorithm, metric=metric, trControl=control)
  predictionsK <- predict(fit.algorK, DtestK)
  StatisticsK<-confusionMatrix(predictionsK, DtestK$Class)
  kEvaluation[k,1:18]<-t(as.data.frame(c(StatisticsK$overall,StatisticsK$byClass)))
  kEvaluation[k,19]<-k
}

#Statistics
StatisticsK$byClass

# Transform output to dataframe
EvalK<-as.data.frame(kEvaluation)

# Plot the metrics (Accuracy, kappa and F1 for all the 200 models)
par(mfrow = c(1,1))
plot(EvalK$Accuracy,type="o",pch=1,col="green",main = paste("Evaluation by ranked genes with algorithm",algorithm,sed=""),xlab = "Top genes",ylab = "Metrics", ylim=c(0,1))
lines(EvalK$F1,type = "o", pch=2, col="blue")
lines(EvalK$Kappa,type = "o", pch=3, col="red")
legend("bottomright",legend=c("Accuracy","F1","Kappa"),pch=c(1,2,3),col=c("green","blue","red"))

#--------------------------------------------------------------
# PART 4. EVALUATION OF THE SELECTED TOP-K GENES 

# After the K value was selected (consensus by all the three classfication algorithms), run again the model for the K top genes. See Methods. 
ks=56 #Selected K top genes 

# Run the model for the K top genes
# Subsets
DtrainKS<-DtrainingRanked[,c(1:ks,5549)]
DtestKS<-DtestingRanked[,c(1:ks,5549)]
# Classification
fit.algorKS <- train(Class~., data=DtrainKS, method=algorithm, metric=metric, trControl=control)
importance <- varImp(fit.algorKS, scale=FALSE) 
plot(importance, main = paste('Top ', ks, "genes with algorithm", algorithm))

# Predict performance using the testing dataset
predictionsKS <- predict(fit.algorKS, DtestKS)
StatisticsKS<-confusionMatrix(predictionsKS, DtestKS$Class)
StatisticsKS
StatisticsKS$overall
StatisticsKS$byClass

# Prediction for ROC curve
predictionsKSroc<- predict(fit.algorKS, DtestKS,type = "prob")[,2] #prob. class=yes
predict.rocr  <- prediction (predictionsKSroc,DtestKS$Class)
perf.rocr     <- performance(predict.rocr,"tpr","fpr") #True and False postivie.rate

# Plot the ROC curve
auc <- as.numeric(performance(predict.rocr ,"auc")@y.values)
par(mfrow = c(1,1))
plot(perf.rocr,type='o', col = "red",main = paste('Area Under Curve =',round(auc,2),"with algorithm", algorithm), ylim=c(0,1.01),xlim=c(0,1))  
abline(a=0, b=1)

# List of predicted vrs actual classes in the Testing dataset:
comparison<-as.data.frame(cbind(predictionsKS,DtestKS$Class))
colnames(comparison)<-c("Prediction","Actual_class")

# ------------------------------------------------------------------------------
# ---------------------THE-END--------------------------------------------------

