# Outlier identification and visualization of preprocessed metabolomics dataset using hierarchical clustering
# Input data should not have missing values for avoiding errors

require(dendextend)
require(colorspace)

#Import data
#dataM: Preprocessed metabolite dataset. All values in this file were imputed using proper method.
set.seed(1234)
dataM<-read.csv("/CMTs_DataMatrix.csv", header=TRUE, row.names=1)

#Log2 transformation
#For the PCA, data was transformed using log2 
#Please ignore if you want to visualize the raw data
log_dataM<-dataM
log_dataM[,3:57]<-log(log_dataM[,3:57]+1, base=2)

#### Dendrogram ####
hc<-hclust(dist(log_dataM[,-c(1,2)]), method="average")
dend<-as.dendrogram(hc)

labels_colors(dend)<-
  rainbow_hcl(3)[
    as.numeric(log_dataM[,1])[order.dendrogram(dend)]]

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
