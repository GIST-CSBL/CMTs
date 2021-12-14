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
log_dataM[,4:58]<-log(log_dataM[,4:58]+1, base=2)

#### PCA plot ####
#PCA for the raw data
pca_dataM<-prcomp(dataM[,-c(1,3)], scale.=TRUE) #remove groups and breeds features
fviz_pca_ind(pca_dataM,
             col.ind=dataM$group,
             addEllipses = FALSE,
             repel=TRUE,
             labelsize=5,
             legend.title="Groups")+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        axis.text=element_text(size=14),
        text = element_text(size = 14))

#PCA for the log2 transformed data
#To eliminate the outlier samples, please rid out them using the index (i.e., the number of rows)
log_dataM_NoOutliers<-log_dataM[-c(2, 5),]
pca_log_dataM<-prcomp(log_dataM_NoOutliers[,-c(1,3)], scale.=TRUE) #remove groups, breeds, and age features

#Colored by groups
fviz_pca_ind(pca_log_dataM,
             col.ind=log_dataM_NoOutliers$group,
             addEllipses = TRUE,
             ellipse.level=0.95,
             label="none",
             pointsize=2,
             mean.point=FALSE,
             ellipse.alpha=0,
             legend.title="Groups")+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        axis.text=element_text(size=14),
        text = element_text(size = 14))

#Colored by breeds
fviz_pca_ind(pca_log_dataM,
             col.ind=log_dataM_NoOutliers$breeds,
             addEllipses = FALSE,
             ellipse.level=0.95,
             label="none",
             pointsize=3,
             mean.point=FALSE,
             legend.title="Breeds")+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        axis.text=element_text(size=14),
        text = element_text(size = 14))

#Colored by ages
fviz_pca_ind(pca_log_dataM,
             col.ind=log_dataM_NoOutliers$age,
             addEllipses = FALSE,
             ellipse.level=0.95,
             label="none",
             pointsize=3,
             mean.point=FALSE,
             gradient.cols=c("#97BDD6", "#7692B5", "black"),
             legend.title="Age(month)")+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        axis.text=element_text(size=14),
        text = element_text(size = 14))

#Export the plot for TIFF
#For exporting the plot, please copy and paste the code "fviz_~" between tiff() and dev.off()
tiff("PCA.tiff", units="mm", width=190, height=180, res=300) #draw PCA plot in 300dpi TIFF file with 190x180 (two-columns figure)
fviz_pca_ind(pca_log_dataM,
             col.ind=log_dataM_NoOutliers$breeds,
             addEllipses = FALSE,
             repel=TRUE,
             labelsize=3,
             legend.title="Breeds")+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        axis.text=element_text(size=14),
        text = element_text(size = 14))
dev.off()
