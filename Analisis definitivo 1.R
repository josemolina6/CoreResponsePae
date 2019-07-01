#ANALISIS DE ALGORITMOS DE RANKEO Y SELECCION DE GENES
#Desarrollado por Jose Molina Mora

library(caret)
library(ROCR)    # para curva ROC
#library(rpart)   # para arbol decision en caso de que se use RF de la parte final
#library(rattle)  # para data set, y arbol decision


#Load data
Dtraining <- read.csv("05 Datos entrenamiento oct18.csv")
Dvalidation <- read.csv("05 Datos validacion oct18.csv")
#Dtraining <- read.csv("05 Datos completos oct18.csv")

dim(Dtraining)
dim(Dvalidation)

# Distribucion de la clase en trainig
percentage <- prop.table(table(Dtraining$Class)) * 100
cbind(freq=table(Dtraining$Clase), percentage=percentage)
#summary(Dtraining)

#Distribucion de la clase en validación
percentage2 <- prop.table(table(Dvalidation$Class))*100
cbind(freq=table(Dvalidation$Clase),percentage=percentage2)
#summary(Dvalidation)

# split input and output
x <- Dtraining[,1:5548]
y <- Dtraining[,5549]
sx <- Dtraining[,c(4,5,6,1,2,3)]#Estaba 1:6, pero es para que se ordene en el gráfico

#boxplot
par(mfrow=c(1,6))
for(i in 1:6) {
  boxplot(x[,i], main=names(Dtraining)[i])
}


#plot y
plot(y)

#install.packages("ellipse")

# scatterplot matrix
featurePlot(x=sx, y=y, plot="ellipse")

# box and whisker plots for each attribute
featurePlot(x=sx, y=y, plot="box")

# density plots for each attribute by class value
scales <- list(x=list(relation="free"), y=list(relation="free"))
featurePlot(x=sx, y=y, plot="density", scales=scales)

#ALGORITMOS

# Run algorithms using 10-fold cross validation
control <- trainControl(method="cv", number=10, classProbs=TRUE)
metric <- "Accuracy"

#Clasification algorithms
# a) linear algorithms
#set.seed(7)DA ERROR CON LDA EN ESTE CASO, hat colinealidad
#fit.lda <- train(Class~., data=Dtraining, method="lda", metric=metric, trControl=control)
# b) nonlinear algorithms
# CART
#fit.cart <- train(Class~., data=Dtraining, method="rpart", metric=metric, trControl=control)
# kNN
fit.knn <- train(Class~., data=Dtraining, method="knn", metric=metric, trControl=control)
# c) advanced algorithms
# SVM
#fit.svm <- train(Class~., data=Dtraining, method="svmRadial", metric=metric, trControl=control)
# Random Forest
#set.seed(7)
#fit.rf <- train(Class~., data=Dtraining, method="rf", metric=metric, trControl=control)

# summarize accuracy of models
results <- resamples(list(knn=fit.knn, svm=fit.svm, rf=fit.rf, lda=fit.lda, cart=fit.cart))
summary(results)


#el mejor estaría dando con RF pero al evaluar los datos no da bien, por lo que se prueban los otros
# y da muy bien con lda y svm

# compare accuracy of models
dotplot(results)

#RANKING PARA ALGORITMO

#ALGORITMO SELECCIONADO
model=fit.knn
kalgoritmo="knn"
#model=fit.knn
#kalgoritmo="knn"
#model=fit.rf
#kalgoritmo="rf"

importance <- varImp(model, scale=FALSE)
#head(importance)
#INDEX
IndexRank <-data.frame(sort(importance$importance$Control, index.return = TRUE, decreasing = TRUE)[2])
#Para SMV y KNN importance$importance$Control, para RF/cart usar $Overall (RF solo con Overall)
Ranking5548<-t(IndexRank)

DtrainingRanked<-Dtraining[,c(Ranking5548[1,],5549)]
DvalidationRanked<-Dvalidation[,c(Ranking5548[1,],5549)]

#GRAFICOS DE LOS TOP

# split input and output
xR <- DtrainingRanked[,1:5548]
yR <- DtrainingRanked[,5549]
sxR <- DtrainingRanked[,c(4,5,6,1,2,3)]

#boxplot
par(mfrow=c(1,6))

for(i in 1:6) {
  boxplot(xR[,i], main=names(DtrainingRanked)[i])
}
title(paste("Top 6 genes with algorithm", kalgoritmo), outer=TRUE) 


# scatterplot matrix
featurePlot(x=sxR, y=yR, plot="ellipse",main = paste("Top 6 genes with algorithm", kalgoritmo))

# box and whisker plots for each attribute
featurePlot(x=sxR, y=yR, plot="box", main = paste("Top 6 genes with algorithm", kalgoritmo))

# density plots for each attribute by class value
scales <- list(x=list(relation="free"), y=list(relation="free"))
featurePlot(x=sxR, y=yR, plot="density", scales=scales, main = paste("Top 6 genes with algorithm", kalgoritmo))



#Ver los primeros 30 datos
k=30
Dtrain30<-DtrainingRanked[,c(1:k,5549)]
Dval30<-DvalidationRanked[,c(1:k,5549)]

fit.algor30 <- train(Class~., data=Dtrain30, method=kalgoritmo, metric=metric, trControl=control)
importance <- varImp(fit.algor30, scale=FALSE) 
plot(importance, main = paste('Top ', k, "genes with algorithm", kalgoritmo))

predictions30 <- predict(fit.algor30, Dval30)
Statistics30<-confusionMatrix(predictions30, Dval30$Class)
#Statistics30
#Statistics30$overall
#Statistics30$byClass

#EVALUACION DE DIF SUBSET DE GENES TOP
#Numero de genes top a evaluar
kvalue=200

  
kEvaluacion<-matrix(,nrow=kvalue, ncol=19) #Son 18 parametros y uno más para calculos que uno quiera hacr en una columna
colnames(kEvaluacion)<-c("Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull", "AccuracyPValue","McnemarPValue","Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Precision", "Recall","F1","Prevalence" ,"Detection Rate","Detection Prevalence","Balanced Accuracy", "K")

for (k in 1:kvalue){  #c(30,35,100)
  #k=32
  DtrainK<-DtrainingRanked[,c(1:k,5549)]
  DvalK<-DvalidationRanked[,c(1:k,5549)]
  
  fit.algorK <- train(Class~., data=DtrainK, method=kalgoritmo, metric=metric, trControl=control)
  predictionsK <- predict(fit.algorK, DvalK)
  StatisticsK<-confusionMatrix(predictionsK, DvalK$Class)
  kEvaluacion[k,1:18]<-t(as.data.frame(c(StatisticsK$overall,StatisticsK$byClass)))
  #Si no quiere agregar extras en columnas , puede usar kEvaluacion[k,]<-t(as.data.frame(c(StatisticsK$overall,StatisticsK$byClass)))
}

StatisticsK$byClass


#Para agregar otro parámetro a dedo, por ejemplo el valor F1 calculado o el valor K (SOLO DE PRUEBA)
for (k in 1:kvalue){
  kEvaluacion[k,19]<-k
#PARA F1 calculado  kEvaluacion[k,19]<-(2*kEvaluacion[k,12]*kEvaluacion[k,13]/(kEvaluacion[k,12]+kEvaluacion[k,13]))
}

EvalK<-as.data.frame(kEvaluacion)



par(mfrow = c(1,1))
plot(EvalK$Accuracy,type="o",pch=1,col="green",main = paste("Evaluation by ranked genes with algorithm",kalgoritmo,sed=""),xlab = "Top genes",ylab = "Metrics", ylim=c(0,1))
lines(EvalK$F1,type = "o", pch=2, col="blue")
lines(EvalK$Kappa,type = "o", pch=3, col="red")
legend("bottomright",legend=c("Accuracy","F1","Kappa"),pch=c(1,2,3),col=c("green","blue","red"))



#Para el valor K escogido:

ks=67
DtrainKS<-DtrainingRanked[,c(1:ks,5549)]
DvalKS<-DvalidationRanked[,c(1:ks,5549)]

fit.algorKS <- train(Class~., data=DtrainKS, method=kalgoritmo, metric=metric, trControl=control)
importance <- varImp(fit.algorKS, scale=FALSE) 
plot(importance, main = paste('Top ', ks, "genes with algorithm", kalgoritmo))

predictionsKS <- predict(fit.algorKS, DvalKS)
StatisticsKS<-confusionMatrix(predictionsKS, DvalKS$Class)
StatisticsKS
StatisticsKS$overall
StatisticsKS$byClass

#PREDICCION

#control <- trainControl(method="cv", number=10, classProbs=TRUE) SE AGREGO CLASSPRO TRUE PARA SVM
fit.algorKS <- train(Class~., data=DtrainKS, method=kalgoritmo, metric=metric, trControl=control)
# ------------------------------------------------------------------------------
predictionsKSroc<- predict(fit.algorKS, DvalKS,type = "prob")[,2] #prob. clase=yes
predict.rocr  <- prediction (predictionsKSroc,DvalKS$Class)
perf.rocr     <- performance(predict.rocr,"tpr","fpr") #True y False postivie.rate

# GRAFICO CURVA ROC
# ------------------------------------------------------------------------------
auc <- as.numeric(performance(predict.rocr ,"auc")@y.values)
par(mfrow = c(1,1))
plot(perf.rocr,type='o', col = "red",main = paste('Area Bajo la Curva =',round(auc,2),"with algorithm", kalgoritmo), ylim=c(0,1.01),xlim=c(0,1))  
abline(a=0, b=1)

# GRAFICO ARBOL DECISION
# ------------------------------------------------------------------------------
#fancyRpartPlot(fit.algorKS) 



#GRAFICA DE RESULTADOS TOTALES POST APERTURA DE ARCHIVOS GUARDADOS
results <- resamples(list(knn=fit.knn, svm=fit.svm, rf=fit.rf, lda=fit.lda, cart=fit.cart))
dotplot(results)

par(mfrow=c(1,6))
for(i in 1:6) {
   boxplot(xR[,i], main=names(DtrainingRanked)[i])
 }
title(paste("Top 6 genes with algorithm", kalgoritmo), outer=TRUE)

featurePlot(x=sxR, y=yR, plot="density", scales=scales, main = paste("Top 6 genes with algorithm", kalgoritmo))

par(mfrow = c(1,1))
plot(EvalK$Accuracy,type="o",pch=1,col="green",main = paste("Evaluation by ranked genes with algorithm",kalgoritmo,sed=""),xlab = "Top genes",ylab = "Metrics", ylim=c(0,1))
lines(EvalK$F1,type = "o", pch=2, col="blue")
lines(EvalK$Kappa,type = "o", pch=3, col="red")
legend("bottomright",legend=c("Accuracy","F1","Kappa"),pch=c(1,2,3),col=c("green","blue","red"))

plot(importance, main = paste('Top ', ks, "genes with algorithm", kalgoritmo))
plot(perf.rocr,type='o', col = "red",main = paste('Area Bajo la Curva =',round(auc,2),"with algorithm", kalgoritmo), ylim=c(0,1.01),xlim=c(0,1))  
abline(a=0, b=1)


#COMPARACION CASO A CASO DEL SET DE VALIDACION
PRED<-as.data.frame(c(predictionsKS))
PRED2<-as.data.frame(c(DvalKS$Class))
juntos<-as.data.frame(c(PRED,PRED2))
View(juntos)



#Extracción de listas
ListaGenes<-as.data.frame(rownames(importance$importance))
write.csv(ListaGenes, file="ListaGenesRF.csv")
