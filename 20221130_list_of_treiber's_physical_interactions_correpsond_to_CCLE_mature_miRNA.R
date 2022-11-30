# this script is to generate list of treiber's physical interactions correpsond to CCLE mature miRNA
# 2022/11/30 made

# import table of CCLE mature miRNA correspond pri-miRNA
# this table is located at ""
setwd("C:/Rdata/CCLE_data")
ccle.miRNA <-read.table("table_of_pri-miRNA_of_CCLE.txt",sep="\t",header = T,stringsAsFactors = F)

# import list of treiber's physical interactions
# this list is located at ""
setwd("C:/Rdata")
treiber.score <-read.table("treiber_heatmap_score.txt",sep="\t",header = T,stringsAsFactors = F)
colnames(treiber.score)[2] <-"primary"

# merge table and list
ccle.treiber <-merge(treiber.score,ccle.miRNA)
ccle.treiber <-ccle.treiber[,c(1,5,2,3)]

# output 
setwd("C:/Rdata/20221129_ROC_curve_of_pvalue_after_cutoff_revised")
write.table(ccle.treiber,"list_of_treiber_physical_interaction_with_not_considering_primary_transcript.txt",sep="\t",row.names = F,quote = F)
