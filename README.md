# 1H NMR based urinary metabolites profiling dataset of canine mammary tumors 
R scripts for the preprocessing of urinary metabolite profiling data from the canine mammary tumors (CMTs) \
R version: 3.6.3 \
IDE: RStudio 

## Abstract
Identifying efficient and sensitive biomarkers for non-invasive tests is one of the most important problems in cancer diagnosis. To deal with this problem, metabolomics is widely applied for biomarker research that detects abnormal changes in cancer patients. Canine mammary tumors exhibit physiological characteristics identical to those in human breast cancer and serve as a useful animal model to conduct breast cancer research. Here, we aimed to provide a reliable large-scale metabolite dataset collected from dogs with mammary tumors using proton nuclear magnetic resonance spectroscopy. We identified 55 metabolites in urine samples from 22 benign, 84 malignant, and 49 healthy control subjects. This dataset provides the potential biomarker candidates for mammary tumors in dogs and insights into cancer-specific metabolic alterations that share similar molecular characteristics.

## Workflow
The toy examples(.csv) for the codes are contained in Datasets folder. 

### Missing value imputation using MetImp
Missing values in the raw metabolite profiling data should be imputed with proper methods. \
In this study, the data was preprocessed using Quantile Regression Imputation of Left-Censored data (QRILC) method. \
We perform the imputation from https://metabolomics.cc.hawaii.edu/software/MetImp/ with minor modification. 

### Prepare the preprocessed metabolite data 
If you use your own metabolite data, you should get a copmlete dataset that imputed and filtered the metabolites.
In this case, metabolites that has the missing value over 20% and mismatched with KEGG ID are excluded from the following analysis.

### CMTs_multivariate analysis.R
To visualize the samples' distribution and check the outliers according to the factors, we used principal component analysis (PCA).
