# this script is to generate talbe about outliers in correlation analysis between expression level of 410 RBP and 377 miRNA
# 2022/12/07 made

# import table about expression levels of RBPs and miRNAs
# these tables are located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\CCLE_Data"
setwd("C:/Rdata/CCLE_data")
RBP.df <-read.table("RNAseq_of_RBP.txt",sep="\t",header = T,stringsAsFactors = F)
RBP.df <-RBP.df[,c(2:411,1,412)]

miRNA.df <-read.table("CCLE_miRNAseq.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)
miRNA.df <-miRNA.df[,c(2:378,1,379)]

# set function for caluculating minimum value of expression level
find.min <-function(x,y){
  for (i in 1:y) {
    m <-x[x[,i]!=0,i]
    mv <-min(m)
    if(i==1){
      mvs <-mv
    }else{
      mvs <-append(mvs,mv)
    }}
  return(min(mvs))
}

# investigate minimum value of each expression level
r.min <-find.min(RBP.df,410)
m.min <-find.min(miRNA.df,377)

# set function for caluculating outlier
ex.outlier <-function(x,y){
  q <-as.numeric(quantile(x[,y]))
  iqr <-IQR(x[,y])
  outlier1 <-q[2]-iqr*1.5
  outlier2 <-q[4]+iqr*1.5
  outliers <-append(outlier1,outlier2)
  return(outliers)
}

sm <-as.data.frame(matrix(nrow = 410*377,ncol = 6))
colnames(sm) <-c("RBP","miRNA","RBP_outlier_small","RBP_outlier_large","miRNA_outlier_small","miRNA_outlier_small")

# 
for (i in 1:410) {
  for (k in 1:377) {
    # extract expression levels of RBP and miRNA
    RBP <-RBP.df[,c(i,411,412)]
    miRNA <-miRNA.df[,c(k,378,379)]
    RBP.miRNA <-merge(RBP,miRNA,by=c("CELL","Site_Primary"))
    RBP.miRNA <-RBP.miRNA[,c(3,4,1,2)]
    
    # count number of cell line without RBP expression
    zero.cell <-nrow(RBP.miRNA[RBP.miRNA[,1]==0,])
    
    # if number of cell line without RBP expression is higher than 100, add minimum RBP expression level
    # if number of cell line without RBP expression is lower than 100, exclude cell line without RBP expression
    if(zero.cell>100){
      RBP.miRNA[,1] <-ifelse(RBP.miRNA[,1]==0,r.min,RBP.miRNA[,1])
    }else{
      RBP.miRNA <-RBP.miRNA[RBP.miRNA[,1]!=0,]
    }
    
    # log2 each expression level 
    RBP.miRNA[,1] <-log2(RBP.miRNA[,1])
    RBP.miRNA[,2] <-log2(RBP.miRNA[,2])
    
    # investigate each outlier 
    r.outlier <-ex.outlier(RBP.miRNA,1)
    m.outlier <-ex.outlier(RBP.miRNA,2)

    # write summary
    sm[(i-1)*377+k,1] <-colnames(RBP.miRNA)[1]
    sm[(i-1)*377+k,2] <-colnames(RBP.miRNA)[2]
    sm[(i-1)*377+k,3] <-r.outlier[1]
    sm[(i-1)*377+k,4] <-r.outlier[2]
    sm[(i-1)*377+k,5] <-m.outlier[1]
    sm[(i-1)*377+k,6] <-m.outlier[2]
    }}

# output summary
setwd("C:/Rdata/20221108_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier")
write.table(sm,"table_about_outliers_in_correlaiton_analysis_between_RBP_level_and_miRNA_level.txt",sep="\t",row.names = F,quote = F)
