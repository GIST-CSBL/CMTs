# Outlier identification and visualization of preprocessed metabolomics dataset using z-score and principal component analysis (PCA), respectively
# Input data should not have missing values for avoiding errors

#Import packages
require(ggplot2)
require(factoextra)
require(dendextend)
library(colorspace)

#Import data
#dataM: Preprocessed metabolite dataset. All values in this file were imputed using proper method.
set.seed(1234)
dataM<-read.csv("/CMTs_DataMatrix.csv", header=TRUE, row.names=1)
dataM$group<-gsub("1", "Benign", dataM$group)
dataM$group<-gsub("2", "Malignant", dataM$group)
dataM$group<-gsub("3", "Normal", dataM$group)

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
pca_log_dataM<-prcomp(log_dataM[,-c(1,2)], scale.=TRUE) #remove groups and breeds features
fviz_pca_ind(pca_log_dataM,
             col.ind=log_dataM$breeds,
             addEllipses = FALSE,
             repel=TRUE,
             labelsize=3)

#Export the plot for TIFF
tiff("PCA.tiff", units="mm", width=190, height=180, res=300) #draw PCA plot in 300dpi TIFF file with 190x180 (two-columns figure)
fviz_pca_ind(pca_log_dataM,
             col.ind=log_dataM$breeds,
             addEllipses = FALSE,
             repel=TRUE,
             labelsize=3)
dev.off()

#### Dendrogram ####
hc<-hclust(dist(log_dataM[,-c(1,2)]), method="average")
dend<-as.dendrogram(hc)

labels_colors(dend)<-
  rainbow_hcl(3)[
    as.numeric(log_cmts[,1])[order.dendrogram(dend)]]

dend<-set(dend, "labels_cex", 0.7)
plot(dend,
     main="Clustered CMTs data set",
     horiz=FALSE, nodePar=list(cex=.05))
legend("topright", inset=.07, legend=c("Benign", "Malignant", "Normal"), fill=rainbow_hcl(3))

tiff("Hierarchical_300_groups.tiff", units="mm", width=380, height=180, res=300) #1-column: 90mm, 2-column:190mm
plot(dend,
     main="Clustered CMTs data set",
     horiz=FALSE, nodePar=list(cex=.05))
legend("topright", inset=.07, legend=c("Benign", "Malignant", "Normal"), fill=rainbow_hcl(3))

dev.off()
