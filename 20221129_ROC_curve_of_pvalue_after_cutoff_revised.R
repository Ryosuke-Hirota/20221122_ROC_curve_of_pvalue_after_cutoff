# this script is to draw ROC curve with number of combination with/without physical interaction
# 2022/11/29 made
# 2022/11/30 revised

setwd("C:/Rdata")
dir.create("20221129_ROC_curve_of_pvalue_after_cutoff_revised")

# import result of correlation analysis between expression level of 410 RBPs and 377 miRNAs
# this result is located at "https://github.com/Ryosuke-Hirota/20221122_ROC_curve_of_pvalue_after_cutoff"
setwd("C:/Rdata/20221108_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier")
c.result <-read.table("summary_of_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier.txt",sep="\t",header = T,stringsAsFactors = F)

# edit list and convert p.value (because I couldn't draw correct ROC curve without converting p.value)
c.result[,1] <-paste0(c.result[,1],"_vs_",c.result[,2])
c.result[,4] <-log10(c.result[,4])*-1
c.result <-c.result[,-2]
colnames(c.result)[1] <-"combination"

# Import list of physical interaction
# This list is located at "https://github.com/Ryosuke-Hirota/20221122_ROC_curve_of_pvalue_after_cutoff"
# Attention : this list is made by pri-miRNA. Thus, when you correspond this list to the correlation analysis using mature miRNA, duplicated rows appear. 
setwd("C:/Rdata/20221129_ROC_curve_of_pvalue_after_cutoff_revised")
phy.list <-read.table("list_of_treiber_physical_interaction_with_not_considering_primary_transcript.txt",sep="\t",header = T,stringsAsFactors = F)

# edit list of physical interaction 
phy.list[,5] <-paste0(phy.list[,3],"_vs_",phy.list[,2])
colnames(phy.list)[5] <-"combination"
phy.list <-phy.list[,c(5,4,1)]

# merge correlation analysis result and list of physical interaction
merge.df <-merge(c.result,phy.list,by="combination",all=T)
merge.df <-subset(merge.df,!is.na(merge.df[,2]))

# annotate whether treiber's combinations match with CCLE correlation analysis result
merge.df[!is.na(merge.df[,5]),7] <-"match"
merge.df[is.na(merge.df[,5]),7] <-"not_match"
colnames(merge.df)[7] <-"match_with_CCLE"
merge.df[,5] <-ifelse(is.na(merge.df[,5]),0,merge.df[,5])

# check number of combinations with/without physical interaction
# 454(unique:431) treiber's combinations with physical interaction(score>3) are matched with CCLE correlation analysis result
# 5760(unique:4420) treiber's combinations without physical interaction(score<=3) are matched with CCLE correlation analysis result
phy <-merge.df[merge.df[,5]>3&merge.df[,7]=="match",1]
phy <-phy[-1]
length(unique(phy))

no.phy <-merge.df[merge.df[,5]<=3&merge.df[,7]=='match',1]
no.phy <-no.phy[-1]
length(unique(no.phy))

# set value of cutoff
cutoff <-seq(50,900,50)

setwd("C:/Rdata/20221129_ROC_curve_of_pvalue_after_cutoff_revised")

# draw ROC curve
for (i in 1:length(cutoff)) {
  
  result <-merge.df
  
  # annotate positive and negative by physical interaction
  result[result[,5]>3,8] <-1
  result[result[,5]<=3,8] <-0
  colnames(result)[8] <-"physical_interaction"
  
  # cutoff number of cell line
  result <-result[result[,4]>=cutoff[i],]
  
  # count number of combinations with/without physical interaction
  # caution : this number is duplicated. If you wanna know non-duplicated number, write additional script.
  count <-as.data.frame(table(result[,7],result[,8]))
  p <-count[3,3]
  np <-count[1,3]
  
  # draw ROC curve
  pdf(paste0("ROC_curve_pvalue_cutoff_",cutoff[i],".pdf"))
  ROC(test=result$p.value, stat=result$physical_interaction, plot="ROC")
  mtext(text = paste0("physical interaction : ",p," , no physical interaction : ",np))
  dev.off()
  }
