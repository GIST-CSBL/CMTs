# Outlier identification and visualization of preprocessed metabolomics dataset using z-score and principal component analysis (PCA), respectively
# Input data should not have missing values for avoiding errors

#Import packages
require(ggplot2)
require(factoextra)

#Import data
#dataM: Preprocessed metabolite dataset. All values in this file were imputed using proper method.
set.seed(1234)
dataM<-read.csv("/CMTs_DataMatrix.csv", header=TRUE, row.names=1)
dataM$group<-gsub("1", "Benign", dataM$group)
dataM$group<-gsub("2", "Malignant", dataM$group)
dataM$group<-gsub("3", "Normal", dataM$group)

#### Outlier detection ####
#If you need to exclude the outlier, please run the below code
dataM_z<-dataM[,-c(1,2)]
z_scores <- as.data.frame(sapply(dataM_z, function(dataM_z) (abs(dataM_z-mean(dataM_z))/sd(dataM_z))))
rownames(z_scores)<-rownames(dataM_z)
z_scores$counts<-rowSums(z_scores>3.29) #Outlier criteria: 3.29
#dataM<-dataM[-c(x,y),] #x,y: the row number of outlier samples

#### PCA plot ####
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
tiff("PCA.tiff", units="mm", width=190, height=180, res=300)
fviz_pca_ind(pca_log_dataM,
             col.ind=log_dataM$breeds,
             addEllipses = FALSE,
             repel=TRUE,
             labelsize=3)
dev.off()
