# Outlier identification and visualization of preprocessed metabolomics dataset using principal component analysis (PCA)
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

#Log2 transformation
#For the PCA, data was transformed using log2 
#Please ignore if you want to visualize the raw data
log_dataM<-dataM
log_dataM[,3:57]<-log(log_dataM[,3:57]+1, base=2)

#### PCA plot ####
#PCA for the raw data
pca_dataM<-prcomp(dataM[,-c(1,2)], scale.=TRUE) #remove groups and breeds features
fviz_pca_ind(pca_dataM,
             col.ind=dataM$group,
             addEllipses = FALSE,
             repel=TRUE,
            labelsize=3,
            legend.title="Groups")

#PCA for the log2 transformed data
#To eliminate the outlier samples, please rid out them using the index
#Colored by breeds
log_dataM_NoOutliers<-log_dataM[-c(2, 5),]
pca_log_dataM<-prcomp(log_dataM_NoOutliers[,-c(1,2)], scale.=TRUE) #remove groups and breeds features
fviz_pca_ind(pca_log_dataM,
             col.ind=log_dataM$breeds,
             addEllipses = FALSE,
             repel=TRUE,
             labelsize=3,
             legend.title="Breeds")

#Colored by ages
fviz_pca_ind(pca_age,
             col.ind=log_cmts_age$age,
             addEllipses = FALSE,
             repel=TRUE,
             labelsize=3,
             gradient.cols=c("#97BDD6", "#7692B5", "black"),
             legend.title="Age(month)")

#Export the plot for TIFF
tiff("PCA.tiff", units="mm", width=190, height=180, res=300) #draw PCA plot in 300dpi TIFF file with 190x180 (two-columns figure)
fviz_pca_ind(pca_log_dataM,
             col.ind=log_dataM$breeds,
             addEllipses = FALSE,
             repel=TRUE,
             labelsize=3)
dev.off()
