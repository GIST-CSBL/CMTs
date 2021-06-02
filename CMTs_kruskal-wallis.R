# Biomarkers of canine mammary tumors in metabolomics data analysis 
# Check the distribution of clinical features in training and test dataset
# Import packages
require(dunn.test)
require(agricolae)
require(MASS)

# Mann-Whitney U test for numerical features
# The table has rows(samples) and columns(features)
clinical<-read.csv("/Datasets/CMTs_clinical_training_test_mwU.csv", header=TRUE, row.names=1)
wilcox.test(clinical[1:3,3], clinical[4:6,3], correct=FALSE)

# Chi-squared test for categorical features
# The input table has rows(Training, Test) and columns (categories) with the number of corresponding samples
# The input table contains the number of samples in 
table<-read.csv("/Datasets/CMTs_clinical_training_test_chisquare.csv", header=TRUE, row.names=1)
table<-as.matrix(table)
chisq.test(table, correct=FALSE)
