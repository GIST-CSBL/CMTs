# Biomarkers of canine mammary tumors in metabolomics data analysis 
# Boxplots with jitter for the metabolic biomarker candidates showing the statistical significance between groups

require(GGally)
require(ggplot2)
require(rcartocolor)
require(grDevices)
require(colorspace)
require(reshape2)
require(ggsignif)
require(ggpubr)
require(tidyverse)
require(ggpol)


#Boxplot with jitter
##########
candidates<-read.csv("/Datasets/CMTs_Biomarker_candidates_list_with_labeled_samples.csv", header=TRUE, row.names=1)
candidates$Label<-as.factor(candidates$Label)
clinical_l.m<-melt(candidates, id.var="Label")
clinical_l.m<-transform(clinical_l.m, Label=factor(Label, levels=c("Normal", "Benign", "Malignant")))
ggplot(data = clinical_l.m, aes(x=Label, y=value, fill=Label)) + 
  ylab(bquote(Log[2]~'Concentration'))+
  facet_wrap( ~ variable, scales="free") +
  theme_bw()+
  theme(legend.position="none", axis.title.x=element_blank())+
  scale_color_discrete_sequential(palette = "Blues", nmax = 10, order = c(3,5,8))+
  scale_fill_discrete_sequential(palette = "Blues", nmax = 10, order = c(3,5,8))+
  geom_boxjitter(aes(x=Label, y=value, fill=Label), jitter.shape=21, jitter.color=NA, outlier.color=NA, errorbar.draw=TRUE)+
  geom_signif(comparisons = list(c("Normal", "Benign"), c("Normal", "Malignant"), c("Benign", "Malignant")),
              map_signif_level=TRUE,
              step_increase = 0.1,
              vjust=0.1,
              test="t.test",
              textsize=3)
