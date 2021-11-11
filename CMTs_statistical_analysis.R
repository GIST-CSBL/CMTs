# Visualization of preprocessed metabolomics dataset using principal component analysis (PCA)
# Input data should not have missing values for avoiding errors

#Import packages
require(ggplot2)
require(ggfortify)
require(factoextra)

#Import data
#dataM: Preprocessed metabolite dataset. All values in this file were imputed using proper method.
set.seed(1234)
dataM<-read.csv("/CMTs_DataMatrix.csv", header=TRUE, row.names=1)
dataM$group<-gsub("1", "Benign", dataM$group)
dataM$group<-gsub("2", "Malignant", dataM$group)
dataM$group<-gsub("3", "Normal", dataM$group)

#If you need to exclude the outlier, please run the below code
#x,y: the row number of outlier samples
#dataM<-dataM[-c(x,y),]

#Log2 transformation
#For the PCA, data was transformed using log2 
#Please ignore if you want to visualize the raw data
log_dataM<-dataM
log_dataM[,3:57]<-log(log_dataM[,3:57]+1, base=2)

#PCA for the raw data
pca_dataM<-prcomp(dataM[,-c(1,2)], scale.=TRUE) #remove groups and breeds features
fviz_pca_ind(pca_dataM,
             col.ind=dataM$group,
             addEllipses = FALSE,
             repel=TRUE)

#PCA for the log2 transformed data
pca_log_dataM<-prcomp(log_cmts[,-c(1,2)], scale.=TRUE) #remove groups and breeds features
fviz_pca_ind(pca_log_dataM,
             col.ind=log_dataM$breeds,
             addEllipses = FALSE,
             repel=TRUE,
             labelsize=3)

#Export the plot for TIFF
tiff("PLS-DA.tiff", units="mm", width=190, height=180, res=300)
fviz_pca_ind(pca_log_dataM,
             col.ind=log_dataM$breeds,
             addEllipses = FALSE,
             repel=TRUE,
             labelsize=3)
dev.off()
