# Biomarkers of canine mammary tumors in metabolomics data analysis 

# Univariate analysis for selecting metabolic biomarker candidates
# Consist of t-test, Benjamini-Hochberg procedure, and log2 fold change calculation

# Read the pre-processed metabolite data as .csv file
# Normalize the values using log2 transformation if the data is raw concentration
table<-read.csv("/Datasets/CMTs_preprocessed_metabolic_data.csv",header=TRUE, row.names=1)
#table<-log(table+1, base=2)

# Calculate the p-value between case and control using Student t-test
# Set the range of case and control
table["p_value"]<-NA
for(i in (1:length(table[,1])))
{
  table[i,]$p_value=t.test(table[i,1:3],table[i,4:6], alternative="two.sided",var.equal=TRUE)[3]$p.value
}
x<-vapply(table$p_value,length,1L)
table<-table[rep(rownames(table),x),]
table$p_value<-unlist(table$p_value,use.names=FALSE)

# Calculate the q-value as correction of p-value using Benjamini-Hochberg procedure (BH)
table["FDR"]<-NA
table["FDR"]<-p.adjust(table$p_value,method="BH")

# Calculate the log2 fold change value.
# Since all values in the table are log2 transformed, fold change values can be obtained by subtraction.
# Consider the exceptions using if-else.
table["LogFC"]<-NA 
for(i in (1:length(table[,1])))
{
  average_tumor=mean(as.numeric(table[i,1:3]))
  average_normal=mean(as.numeric(table[i,4:6]))
  
  if(average_normal==0)
  {
    if(average_tumor==0)
    {
      table[i,]$LogFC<-"NAN"
    }
    else
    {
      table[i,]$LogFC<-Inf
    }
  }
  else if(average_tumor==0)
  {
    table[i,]$LogFC<--Inf
  }
  else
  {
    table[i,]$LogFC<-average_tumor-average_normal
  }
}
