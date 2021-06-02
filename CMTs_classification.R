# Biomarkers of canine mammary tumors in metabolomics data analysis 
# Classification using logistic regression, random forest, support vector machine, and linear discriminant analysis
# Cases are labeled with 1 and controls are 0

require(Epi)
require(randomForest)
require(MASS)
require(caret)
require(ROCR)
require(e1071)
require(pROC)
require(grid)

training<-read.csv("/CMTs_Training.csv", header=TRUE, row.names=1)
training$Label<-as.factor(training$Label)
test<-read.csv("/CMTs_Test.csv", header=TRUE, row.names=1)
test$Label<-as.factor(test$Label)

###########################
###########################
#Logistic Regression
set.seed(1234)

#ctrl<-trainControl(method="repeatedcv", number=5, savePredictions = TRUE) #5-fold CV with iteration 5
mod_fit<-train(Label~., 
               data=training, method="glm", family="binomial")
mod_fit

mod_pred_r<-predict(mod_fit, test) #Accuracy (Cutoff:0.5)
mod_pred_r
confusionMatrix(mod_pred_r, test$Label, positive="1")

mod_pred<-predict(mod_fit, test, type="prob") #AUC
mod_pred

###########################
###########################
#Random Forest (RF)
set.seed(1234)

#Multiple comparisons of options
#Choose the options having the best performance in training datset
########
metric<-"Accuracy"
customRF<-list(type="Classification", library="randomForest", loop=NULL)
customRF$parameters<-data.frame(parameter=c("mtry", "ntree"), class=rep("nuberic", 2), label=c("mtry", "ntree"))
customRF$grid<-function(x,y,len=NULL, search="grid"){}
customRF$fit<-function(x,y,wts, param, lev, last, weights, classProbs, ...){
  randomForest(x, y, mtry=param$mtry, ntree=param$ntree, ...)
}

customRF$predict<-function(modelFit, newdata, preProc=NULL, submodels=NULL)
  predict(modelFit, newdata)
customRF$prob<-function(modelFit, newdata, preProc=NULL, submodels=NULL)
  predict(modelFit, newdata, type="prob")
customRF$sort<-function(x) x[order(x[,1]),]
customRF$levels<-function(x) x$classes
control<-trainControl(method="boot", search="grid")

tunegrid<-expand.grid(.mtry=c(1:3), .ntree=c(100, 200, 300, 400))
custom<-train(Label~sn.Glycero.3.phosphocholine+Dimethylamine+Acetone+X2.Hydroxybutyrate+Succinate+Carnitine, 
              data=training, method=customRF, metric=metric, tuneGrid=tunegrid, trControl=control)
custom
plot(custom)
########

set.seed(1234)
rf_met<-randomForest(Label~.,
                     data=training, mtry=floor(sqrt(3)), ntree=400, proximity=TRUE, importance=TRUE)
print(rf_met)
plot(rf_met)
legend(locator(1), legend=c("Benign", "Average", "Normal"), lty=c(2,1,2), col=c(2,1,3))
confusionMatrix(rf_met, training$Label, positive="1")

#Result of variable importances with mean decrease accuracy score
importance(rf_met)
varImpPlot(rf_met)

metPred<-predict(rf_met, newdata=test)
table(metPred, test$Label)
confusionMatrix(metPred, test$Label, positive="1")
modPred<-predict(rf_met, newdata=test, type="prob")
modPred

###########################
###########################
#Suppor Vector Machine (SVM)
set.seed(1234)

svm_train_x<-training[,-1]
svm_train_y<-training[,1]

svm_test_x<-test[,-1]
svm_test_y<-test[,1]

tune.svm(Label~., 
         data=training, gamma=2^(-5:5), cost=2^(-5:5))
svm_model_after_tune<-svm(Label~.,
                          data=training, kernel="radial", gamma=0.03125, cost=2, probability=TRUE)
summary(svm_model_after_tune)
svm_model_after_tune

pred_test<-predict(svm_model_after_tune, test, type="prob", probability = TRUE)

table(pred_test, svm_test_y)
pred_test
confusionMatrix(pred_test, svm_test_y, positive="1")
roc_test<-ROC(test=attr(pred_test, "probabilities")[,1], stat=test[,1], plot="ROC")

###########################
###########################
#LDA
set.seed(1234)
ld<-lda(formula=Label~., data=training)
ld_pred<-predict(ld, test)
summary(ld_pred$class)
xtab<-table(ld_pred$class, test$Label)
confusionMatrix(xtab, positive="1")
ld_pred

###########################
###########################
#ROC curves
#Performance of training set using each machine learning method
#Run both lr_train and lr_pr_train for check the performance in training dataset

lr_train<-predict(mod_fit, type="prob")[,2] #LR
lr_pr_train<-prediction(lr_train, training$Label) #LR

#lr_train<-predict(rf_met, type="prob") #RF
#lr_pr_train<-prediction(lr_train[,2], training$Label) #RF

#lr_train<-predict(svm_model_after_tune, svm_train_x, probability=TRUE) #SVM
#lr_pr_train<-prediction(attr(lr_train, "probabilities")[,1], svm_train_y)#SVM

lr_pr_train
trainingAcc<-performance(lr_pr_train, "acc")
max(trainingAcc@y.values[[1]])
performance(lr_pr_train, measure="auc")@y.values[[1]]#AUC

#training set sens. spec.
lr_tf_train<-performance(lr_pr_train, "sens", "spec")
ix<-which.min(abs(lr_tf_train@alpha.values[[1]]-0.5))
sensitivity<-lr_tf_train@y.values[[1]][ix]
specificity<-lr_tf_train@x.values[[1]][ix]
sensitivity
specificity
lr_tf_train@alpha.values[[1]][ix]

#test set AUC
pred<-prediction(mod_pred[,2], test$Label) #LR
#pred<-prediction(modPred[,2], test$Label) #RF
#pred<-prediction(attr(pred_test, "probabilities")[,1], svm_test_y) #SVM
#pred<-prediction(ld_pred$posterior[,2], test$Label) #LDA
pred

performance(pred, "auc")@y.values[[1]] #AUC = y.values
pref<-performance(pred, "prec", "rec")
pref@x.values[[1]] #Recall value
performance(pred, "acc", "cutoff")

plot(performance(pred, "tpr", "fpr")) #AUC
plot(performance(pred, "prec", "rec")) #AUPR

ex<-performance(pred, "sens", "spec")
ix1<-which.min(abs(ex@alpha.values[[1]]-0.5))
sensitivity<-ex@y.values[[1]][ix1]
specificity<-ex@x.values[[1]][ix1]
sensitivity
specificity
ex@alpha.values[[1]][ix1]

testAcc<-performance(pred, "acc")
max(testAcc@y.values[[1]])
testAcc@y.values[[1]][ix1]

lr_f1_test<-performance(pred, "f") #f1 score
ixi<-which.min(abs(lr_f1_test@x.values[[1]]-0.5))
f1score<-lr_f1_test@y.values[[1]][ixi]
f1score

#Drawing multiple ROC with ROCR
#Input the ROC data into the variables
lr_roc<-roc(test$Label, mod_pred$`1`)
rf_roc<-roc(test$Label, modPred[,2])
svm_roc<-roc(test$Label, attr(pred_test, "probabilities")[,1])
lda_roc<-roc(test$Label, ld_pred$posterior[,2])
roc4_list<-list("Logistic regression"=lr_roc, "Random forest"=rf_roc, "Support vector machine"=svm_roc, "Linear discriminant analysis"=lda_roc)
roc.test(lr_roc, rf_roc, svm_roc)

#The plot of ROC curves
ggroc(roc4_list, legacy.axes=TRUE, aes=c("colour"))+
  scale_color_discrete(labels=c("Random forest (AUC=0.72)", "Linear discriminant analysis (AUC=0.738)", "Logistic regression (AUC=0.74)", "Support vector machine (AUC=0.762)"))+
  labs(x="1-Specificity", y="Sensitivity")+
  geom_abline(linetype="dashed")+
  theme_bw()+
  theme(legend.position = c(0.7, 0.25), legend.title=element_text(size=13), legend.background=element_rect(color="black"))+
  labs(color='Classification methods')+
  geom_line(size=1)

