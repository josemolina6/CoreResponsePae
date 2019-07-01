#CLASIFICACION DE PSEUDOMONAS AERUGINOSA BAJO PERTURBACIONES

#setwd("C:/Users/INSPIRON 14 7000/Desktop")


#install.packages("caret")
#install.packages("mlbench")
library(caret)
library(mlbench)
 
#preDatos <- read.csv("Datos pert_5600g replicas.csv")
#Datos<-preDatos[,c(1:100,5549)]
Datos <- read.csv("05 Datos completos oct18.csv")

dim(Datos)[2]

#numero de intentos
kreplicas=100
ktop=56
#metrics<-c("Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull", "AccuracyPValue","McnemarPValue","Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value","Prevalence" ,"Detection Rate","Detection Prevalence","Balanced Accuracy")
metrics<-c("Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull", "AccuracyPValue","McnemarPValue","Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision", "Recall","F1","Prevalence" ,"Detection Rate","Detection Prevalence","Balanced Accuracy")


#------------------------------------------------------------------------------------------
algorithm="svmRadial"
#algorithm="rf"  

Result1 <- as.matrix(as.data.frame.list(replicate(kreplicas, ClasifPAIndep1(Datos,algorithm,ktop))))


rownames(Result1)<-c(metrics,sprintf("Gen%d",1:ktop),"Clase")
colnames(Result1)<-sprintf("Run%d",1:kreplicas)
#Result

R_accuracy1<-as.double(Result1[1,])
#R_accuracy
R_kappa1<-as.double(Result1[2,])
#dim(Result)

#Grafica
par(mfrow = c(1,1))
hist(R_accuracy1, main=paste("Evaluation by",algorithm,"with",kreplicas,"replicates and top", ktop, "genes"), xlab = "Accuracy", ylab="Frequency")
hist(R_kappa1, main=paste("Kappa by",algorithm,"with",kreplicas,"replicates" ), xlab = "Kappa value", ylab="Frecuency")
row
GenesTable1<-Result1[-c(1:18),]

GeneCount1<-as.data.frame(sort(table(GenesTable1), decreasing = TRUE))
GeneCount1[1:10,]

boxplot(R_accuracy1,main=paste("Accuracy by SVM with",kreplicas,"replicates and top genes:", ktop))

par(mfrow = c(3,6))
for (i in 1:18){
  R_prueba1<-as.double(Result1[i,])
  boxplot(R_prueba1,main=metrics[i], ylim=c(0,1))
}


par(mfrow = c(1,5))
for (i in c(1,2,12,13,14)){
  R_prueba1<-as.double(Result1[i,])
  boxplot(R_prueba1,main=metrics[i], ylim=c(0,1))
}
title(paste("Metric for SVM"), outer = TRUE)



#------------------------------------------------------------------------------------------
#algorithm="svmRadial"
algorithm="knn"  

Result2 <- as.matrix(as.data.frame.list(replicate(kreplicas, ClasifPAIndep1(Datos,algorithm,ktop))))


rownames(Result2)<-c(metrics,sprintf("Gen%d",1:ktop),"Clase")
colnames(Result2)<-sprintf("Run%d",1:kreplicas)
#Result

R_accuracy2<-as.double(Result2[1,])
#R_accuracy
R_kappa2<-as.double(Result2[2,])
#dim(Result)

#Grafica
par(mfrow = c(1,1))
hist(R_accuracy2, main=paste("Accuracy by",algorithm,"with",kreplicas,"replicates and top genes:", ktop), xlab = "Exactitud en la clasificaci?n", ylab="Frecuencia")
hist(R_kappa2, main=paste("Kappa by",algorithm,"with",kreplicas,"replicates" ), xlab = "Valor kappa", ylab="Frecuencia")

GenesTable2<-Result2[-c(1:18),]

GeneCount2<-as.data.frame(sort(table(GenesTable2), decreasing = TRUE))
GeneCount2[1:10,]

boxplot(R_accuracy2, main=paste("Accuracy by KNN with",kreplicas,"replicates and top genes:", ktop))
boxplot(R_accuracy1, R_kappa1)

par(mfrow = c(3,6))
for (i in 1:18){
  R_prueba2<-as.double(Result2[i,])
  boxplot(R_prueba2,main=metrics[i], ylim=c(0,1))
}


par(mfrow = c(1,5))
for (i in c(1,2,12,13,14)){
  R_prueba2<-as.double(Result2[i,])
  boxplot(R_prueba2,main=metrics[i], ylim=c(0,1))
}
title(paste("Metric for KNN"), outer = TRUE)


#------------------------------------------------------------------------------------------
#algorithm="svmRadial"
algorithm="rf"  

Result3 <- as.matrix(as.data.frame.list(replicate(kreplicas, ClasifPAIndep2(Datos,algorithm,ktop))))


rownames(Result3)<-c(metrics,sprintf("Gen%d",1:ktop),"Clase")
colnames(Result3)<-sprintf("Run%d",1:kreplicas)
#Result

R_accuracy3<-as.double(Result3[1,])
#R_accuracy
R_kappa3<-as.double(Result3[2,])
#dim(Result)

#Grafica
par(mfrow = c(1,1))
hist(R_accuracy3, main=paste("Accuracy by",algorithm,"with",kreplicas,"replicates and top genes:", ktop), xlab = "Exactitud en la clasificaci?n", ylab="Frecuencia")
hist(R_kappa3, main=paste("Kappa by",algorithm,"with",kreplicas,"replicates" ), xlab = "Valor kappa", ylab="Frecuencia")

GenesTable3<-Result3[-c(1:18),]

GeneCount3<-as.data.frame(sort(table(GenesTable3), decreasing = TRUE))
GeneCount3[1:10,]

boxplot(R_accuracy3, main=paste("Accuracy by RF with",kreplicas,"replicates and top genes:", ktop))


par(mfrow = c(3,6))
for (i in 1:18){
  R_prueba3<-as.double(Result3[i,])
  boxplot(R_prueba3,main=metrics[i], ylim=c(0,1))
}


par(mfrow = c(1,5))
for (i in c(1,2,12,13,14)){
  R_prueba3<-as.double(Result3[i,])
  boxplot(R_prueba3,main=metrics[i], ylim=c(0,1))
}
title(paste("Metric for RF"), outer = TRUE)



#--------------------------------------------------------------------------------------------------------
  
par(mfrow = c(1,1))

boxplot(R_accuracy1,R_kappa1, xlab="Algorithms",ylab="Accuracy", col=c("red","yellow"))

boxplot(R_accuracy1,R_accuracy2,R_accuracy3, col=c("blue","red","green"),xlim=c(0.5,3.8),main="Accuracy", xlab="Algorithms",ylab="Accuracy")
legend("bottomright",legend=c("SVM","KNN","RF"),pch=c(15,15,15),col=c("blue","red","green"))


#-------------------------------------------------------------------------------------------------------
  
#TABLAS
GenesSVM<-as.data.frame.array(GeneCount1)
GenesKNN<-as.data.frame.array(GeneCount2)
GenesRF<-as.data.frame.array(GeneCount3)

ListaSVM<-as.data.frame(GenesSVM[2:57,1])
dim(ListaSVM)
ListaKNN<-as.data.frame(GenesKNN[2:57,1])
dim(ListaKNN)
ListaRF<-as.data.frame(GenesRF[2:57,1])
dim(ListaRF)

Listas<-c(ListaSVM,ListaKNN, ListaRF)

write.csv(Listas, file = "RandomListaGenes.csv")
