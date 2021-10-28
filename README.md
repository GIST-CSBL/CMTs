# 1H NMR based urinary metabolites profiling dataset of canine mammary tumors 
R scripts for the identification of urinary metabolic biomarker candidates in canine mammary tumors (CMTs) \
R version: 3.6.3 \
IDE: RStudio 

## Note
Before run the codes, missing values in the raw metabolite data should be imputed with proper methods. \
In this study, the data was preprocessed using Quantile Regression Imputation of Left-Censored data (QRILC) method. \
You can perform the imputation from https://metabolomics.cc.hawaii.edu/software/MetImp/

## Workflow
The toy examples(.csv) for the codes are contained in Datasets folder. 

### 0. Prepare the preprocessed metabolite data 
If you use your own metabolite data, you should get a copmlete dataset that imputed and filtered the metabolites.
In this case, metabolites that has the missing value over 20% and mismatched with KEGG ID are excluded from the following analysis.

### 1. CMTs_univariate analysis.R 
For the selection of biomarker candidates, metabolites that has the p-value lower than 0.01 were remained.\
There are additional codes for false discovery rate and fold change value.

### 2. CMTs_multivariate analysis.R
To find the metabolites that contribute to classificaiton of each group, we performed the Partial Least Squares-Dicriminant Analysis (PLS-DA).
The importance of each variable are calculated with Variable Importance in Projection (VIP) values.

### 3. CMTs_Mann-Whitney_Chi-squared.R
After choose the biomarker candidates according to the p-value and VIP thresholds (p-value < 0.01, VIP > 1), samples were separated into training and test dataset. \
For the pair comparison between control and case, we checked the distribution of samples using Mann-Whitney U test and Chi-squared test.

### 4. CMTs_classification.R
To find the most optimal combination of biomarker candidates, various classification methods are applied. \
In this study, we applied logistic regression, linear discriminant anlaysis for linear classification and random forest, support vector machine for nonlinear classification.

### 5. CMTs_plots.R
Finally, we visualized the patterns of metabolites from healthy control to benign and malignant states. 
