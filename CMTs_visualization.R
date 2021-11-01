# Visualization of preprocessed metabolomics dataset using principal component analysis (PCA)
# Input data should not have missing values for avoiding errors

#Import packages
require(ggplot2)
require(ggfortify)

#Import data
#dataM: Preprocessed metabolite dataset. All values in this file were imputed using proper method.
set.seed(1234)
dataM<-read.csv("/CMTs_DataMatrix.csv", header=TRUE, row.names=1)
dataM$group<-gsub("1", "Benign", dataM$group)
dataM$group<-gsub("2", "Malignant", dataM$group)
dataM$group<-gsub("3", "Normal", dataM$group)

#Log2 transformation
#For the PCA, data was transformed using log2 
dataM[,2:56]<-log(dataM[,2:56]+1, base=2)

#PCA
pca_dataM<-prcomp(dataM[,-c(1,2)], scale.=TRUE)
autoplot(pca_dataM)
autoplot(pca_dataM, data=dataM, colour='breeds', frame.type='norm', frame.colour="group", shape='group', frame.alpha=0.05)+
  theme_bw()+
  theme(legend.key = element_rect(fill = "white"))+
  scale_shape_manual(values=c(17, 16, 8))

