# this script is to draw ROC curve
# made by hirota 2022/11/22

# activate package to draw ROC curve
library(Epi)

# make directory
setwd("C:/Rdata")
dir.create("20221122_ROC_curve_of_pvalue_after_cutoff")

# inport result of correlation analysis between expression level of 410 RBPs and 377 miRNAs
# this result is located at "https://github.com/Ryosuke-Hirota/20221122_ROC_curve_of_pvalue_after_cutoff"
setwd("C:/Rdata/20221108_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier")
c.result <-read.table("summary_of_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier.txt",sep="\t",header = T,stringsAsFactors = F)

# edit list and convert p.value (because I couldn't draw correct ROC curve without converting p.value)
c.result[,1] <-paste0(c.result[,1],"_vs_",c.result[,2])
c.result[,4] <-log10(c.result[,4])*-1

# inport list of physical interaction
# this list is located at "https://github.com/Ryosuke-Hirota/20221017_ROC_curve_with_list_of_functional_or_physical_interactions"
setwd("C:/Rdata/20221017_ROC_curve_for_cutoff_with_functional_interactions")
phy.list <-read.table("list_of_treiber_physical_interaction_between_RBP_and_miRNA.txt",sep="\t",header = T,stringsAsFactors = F)

# cutoff singnal > 3 and edit list
phy.list <-phy.list[phy.list[,5]>3,]
phy.list[,1] <-paste0(phy.list[,4],"_vs_",phy.list[,2])
phy.list <-phy.list[,c(1,5)]


# set value of cutoff
cutoff <-seq(50,900,50)

setwd("C:/Rdata/20221122_ROC_curve_of_pvalue_after_cutoff")

# draw ROC curve
for (i in 1:length(cutoff)) {
  # remove NA rows
  result <-c.result
  result <-subset(result,!is.na(result[,3]))
  
  # annotate postive and negative by physical interaction
  m <-match(phy.list[,1],result[,1])
  m <-na.omit(m)
  result[m,6] <-1
  result[-m,6] <-0
  colnames(result)[6] <-"physical_interaction"
  
  # cutoff number of cell line
  result <-result[result[,5]>=cutoff[i],]
  
  # draw ROc curve
  pdf(paste0("ROC_curve_pvalue_cutoff_",cutoff[i],".pdf"))
  ROC(test=result$p.value, stat=result$physical_interaction, plot="ROC")
  dev.off()
}
