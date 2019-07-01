#FUNCI?N PARA LA PARTICI?N DE DATOS (80/20) y hacer selecci?n de atributos UNICAMENTE
#con los datos de entrenamiento (sin incluir los de validaci?n).
#Con la llamada de la funci?n ClasifPAIndep, la partici?n se hace al azar en cada corrida, y la
#selecci?n de atributos se hace cada vez, dejando el top 56 (1%) de los genes rankeados.
#Con esos 56 atributos, se entrenan algoritmos de clasificaci?n (con validaci?n cruzada de 10x)
# y se hace la evaluaci?n con el conjunto de validaci?n. 

#setwd("C:/Users/INSPIRON 14 7000/Desktop")

ClasifPAIndep2 <- function(Datos,algoritmo,valork) {

algoritmo=algorithm  
valork=ktop
nc<-dim(Datos)[2]

# Partir los datos en 80% training y 20% validaci?n ALEATORIAMENTE CADA TIRO!
validation_index<-createDataPartition(Datos$Class, p=0.80, list = FALSE)

#Con preDTraining se har? la selecci?n de atributos, y luego el index se pasa por preDValidat 
preDtraining <- Datos[validation_index, ]
preDvalidation<- Datos[-validation_index,]

#RANKING DE ATRIBUTOS
# prepare training scheme

#??????????????????????????????????????????
#control <- trainControl(method="repeatedcv", number=10, repeats=3)
control <- trainControl(method="cv", number=10)#, classProbs = TRUE)

# train the model
model <- train(Class~., data=preDtraining, method=algoritmo, preProcess="scale", trControl=control)

# estimate variable importance
importance <- varImp(model, scale=FALSE)

#Crear el index de los datos basado en el ranking



IndexRank <-data.frame(sort(importance$importance$Overall, index.return = TRUE, decreasing = TRUE)[2])
#Usar importance$importance$Overall en el caso de rf, $Control para SVM y KNN
Featk<-t(IndexRank[1:valork,])
#IndexRank[1:4,]
#Feat56[1,]
#IndexRank

#Reconstruir los datos con los genes seleccionados (para eso se usa el index)
Dtraining<-preDtraining[,c(Featk[1,],nc)]
Dvalidation<-preDvalidation[,c(Featk[1,],nc)]

Genesk<-as.data.frame(colnames(Dtraining))
out1<-Genesk

# Run algorithms using 10-fold cross validation
control <- trainControl(method="cv", number=10)#, classProbs = TRUE)
metric <- "Accuracy"

#CORRER ALGORITMOS

# a) linear algorithms
#set.seed(7)
#fit.lda <- train(Clase~., data=Dtraining, method="lda", metric=metric, trControl=control)
# b) nonlinear algorithms
# CART
#set.seed(7)
#fit.cart <- train(Clase~., data=Dtraining, method="rpart", metric=metric, trControl=control)
# kNN
#set.seed(7)
#fit.knn <- train(Clase~., data=Dtraining, method="knn", metric=metric, trControl=control)
# c) advanced algorithms
# SVM
#set.seed(7)
#fit.svm <- train(Clase~., data=Dtraining, method="svmRadial", metric=metric, trControl=control)
# Random Forest
#set.seed(7)
fit <- train(Class~., data=Dtraining, method=algoritmo, metric=metric, trControl=control)

#?summarize accuracy of models
#results <- resamples(list(lda=fit.lda,cart=fit.cart, knn=fit.knn, svm=fit.svm, rf=fit.rf))
#results <- resamples(list(cart=fit.cart, knn=fit.knn, svm=fit.svm, rf=fit.rf))

#summary(results)

# compare?accuracy of models
#dotplot(results)

# estimate skill of model on the validation dataset
predictions <- predict(fit, Dvalidation)
VAL<-confusionMatrix(predictions, Dvalidation$Class)
#VAL
#VAL$overall
#VAL$byClass

output1<-t(as.data.frame(c(VAL$overall,VAL$byClass)))
output1
output2<-t(as.data.frame(Genesk))
output2
#colnames(output2)<-c(sprintf("Gen%d",1:56),"Clase")

#prueba 
#out3<-list(output1,output2)

out3<-as.data.frame(c(output1,output2))

return(out3)
}


