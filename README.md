# 1H NMR based urinary metabolites profiling dataset of canine mammary tumors 
R scripts for the preprocessing of urinary metabolite profiling data from the canine mammary tumors (CMTs) \
R version: 3.6.3 \
IDE: RStudio 

## Abstract
The identification of efficient and sensitive biomarkers for non-invasive tests is one of the major challenges in cancer diagnosis. To address this challenge, metabolomics is widely applied for identifying biomarkers that detects abnormal changes in cancer patients. Canine mammary tumors exhibit physiological characteristics identical to those in human breast cancer and serve as a useful animal model to conduct breast cancer research. Here, we aimed to provide a reliable large-scale metabolite dataset collected from dogs with mammary tumors, using proton nuclear magnetic resonance spectroscopy. We identified 55 metabolites in urine samples from 20 benign, 87 malignant, and 49 healthy control subjects. This dataset provides details of mammary tumor-specific metabolites in dogs and insights into cancer-specific metabolic alterations that share similar molecular characteristics.

## Workflow

### Missing value imputation using MetImp
Missing values in the raw metabolite profiling data should be imputed with proper methods. \
In this study, the data were preprocessed using the Quantile Regression Imputation of Left-Censored data (QRILC) method. \
We perform the imputation from https://metabolomics.cc.hawaii.edu/software/MetImp/ with minor modification of the code in the https://github.com/WandeRum/MVI-evaluation. 

### Prepare the preprocessed metabolite data 
If you use your own metabolite data, you should get a copmlete dataset that imputed and filtered the metabolites.
In this case, metabolites that has the missing value over 20% and mismatched with KEGG ID are excluded from the following analysis.

### CMTs_Dendrogram.R
To identify outliers, we used hierarchical clustering. The group value should be numerical value for coloring by index.

### CMTs_PCA.R
To identify outliers and visualize the samples' distribution, we used principal component analysis (PCA). This code contains not only basic PCA but also coloring by breeds and age.
