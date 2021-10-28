# 1H NMR based urinary metabolites profiling dataset of canine mammary tumors 
R scripts for the preprocessing of urinary metabolite profiling data from the canine mammary tumors (CMTs) \
R version: 3.6.3 \
IDE: RStudio 

## Note
Before run the codes, missing values in the raw metabolite data should be imputed with proper methods. \
In this study, the data was preprocessed using Quantile Regression Imputation of Left-Censored data (QRILC) method. \
You can perform the imputation from https://metabolomics.cc.hawaii.edu/software/MetImp/

## Workflow
The toy examples(.csv) for the codes are contained in Datasets folder. 

### Prepare the preprocessed metabolite data 
If you use your own metabolite data, you should get a copmlete dataset that imputed and filtered the metabolites.
In this case, metabolites that has the missing value over 20% and mismatched with KEGG ID are excluded from the following analysis.

### CMTs_multivariate analysis.R
To find the metabolites that contribute to classificaiton of each group, we performed the Partial Least Squares-Dicriminant Analysis (PLS-DA).
The importance of each variable are calculated with Variable Importance in Projection (VIP) values.
