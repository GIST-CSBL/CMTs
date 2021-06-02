# Biomarkers of canine mammary tumors in metabolomics data analysis 
# Multivariate analysis for selecting metabolic biomarker candidates.
# Consist of Partial Least Squared-Discriminate Analysis (PLS-DA) and calculation of variable importance in projection (VIP) scores.

#Import packages
require(ropls)

#Import data
#DataM: Metabolite dataset. All values in this file were imputed and log2 transformed.
#SampleM: Sample information.
#VariableM: Metabolite information.
set.seed(1234)
dataM<-read.csv("/Datasets/CMTs_DataMatrix.csv", header=TRUE, row.names=1)
sampleM<-read.csv("/Datasets/CMTs_SampleMetadata.csv", header=TRUE, row.names=1)
variableM<-read.csv("/Datasets/CMTs_VariableMetadata.csv", header=TRUE, row.names=1)

listM<-list(dataM, sampleM, variableM)
names(listM)<-c("dataMatrix", "sampleMetadata", "variableMetadata")

#PCA
#set.seed(1234)
#listM.pca<-opls(dataM)
#StateFc<-sampleM[, "State"]
#plot(listM.pca, typeVc="x-score", parAsColFcVn=StateFc, parEllipsesL=TRUE)

#PLS-DA
listM.plsda<-opls(dataM, StateFc, predI=2, permI=200)
plot(listM.plsda, typeVc="x-score", parAsColFcVn=StateFc, parEllipsesL=TRUE)

#Calcuate the VIP scores
VIPs<-getVipVn(listM.plsda)

