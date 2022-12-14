# this script is to draw plot of correlation analysis between expression level of 410 RBP and 377 miRNA with excluding outlier
# 2022/12/06 made

library(ggplot2)

# import table about expression levels of RBPs and miRNAs
# these tables are located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\CCLE_Data"
setwd("C:/Rdata/CCLE_data")
RBP.df <-read.table("RNAseq_of_RBP.txt",sep="\t",header = T,stringsAsFactors = F)
RBP.df <-RBP.df[,c(2:411,1,412)]

miRNA.df <-read.table("CCLE_miRNAseq.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)
miRNA.df <-miRNA.df[,c(2:378,1,379)]

# set function for calculating minimum value of expression level
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

# set function for calculating outlier
ex.outlier <-function(x,y){
  q <-as.numeric(quantile(x[,y]))
  iqr <-IQR(x[,y])
  outlier1 <-q[2]-iqr*1.5
  outlier2 <-q[4]+iqr*1.5
  outliers <-append(outlier1,outlier2)
  return(outliers)
}

sm <-as.data.frame(matrix(nrow = 410*377,ncol = 6))
colnames(sm) <-c("RBP","miRNA","r","p.value","number_of_cell_line","number_of_cell_line_without_expression")

# draw plot
for (i in 1:410) {
  for (k in 1:377) {
    # extract expression levels of RBP and miRNA
    RBP <-RBP.df[,c(i,411,412)]
    miRNA <-miRNA.df[,c(k,378,379)]
    RBP.miRNA <-merge(RBP,miRNA,by=c("CELL","Site_Primary"))
    RBP.miRNA <-RBP.miRNA[,c(3,4,1,2)]
    
    #  minimum expression level of miRNA
    RBP.miRNA[,2] <-ifelse(RBP.miRNA[,2]==0,m.min,RBP.miRNA[,2])
    
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
    
    # exclude outlier of RBP expression level
    RBP.miRNA <-RBP.miRNA[RBP.miRNA[,1]>=r.outlier[1],]
    RBP.miRNA <-RBP.miRNA[RBP.miRNA[,1]<=r.outlier[2],]
    
    # exclude outlier of miRNA expression level
    RBP.miRNA <-RBP.miRNA[RBP.miRNA[,2]>=m.outlier[1],]
    RBP.miRNA <-RBP.miRNA[RBP.miRNA[,2]<=m.outlier[2],]
    
    # column of Site_Primary as factor 
    RBP.miRNA[,4] <-as.factor(RBP.miRNA[,4])
    
    # calculate correlation coefficient
    r <-cor.test(RBP.miRNA[,1],RBP.miRNA[,2],method = "pearson")

    # write summary
    sm[(i-1)*377+k,1] <-colnames(RBP.miRNA)[1]
    sm[(i-1)*377+k,2] <-colnames(RBP.miRNA)[2]
    sm[(i-1)*377+k,3] <-signif(r$estimate,3)
    sm[(i-1)*377+k,4] <-signif(r$p.value,3)
    sm[(i-1)*377+k,5] <-nrow(RBP.miRNA)
    sm[(i-1)*377+k,6] <-zero.cell
    }}

# add column to order by p.value
sm[,7] <-log10(sm[,4])*-1
sm <-sm[order(sm[,3],sm[,7],decreasing = T),]
sm <-sm[,-7]

# output summary
setwd("C:/Rdata/20221108_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier")
write.table(sm,"summary_of_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier_and_considering_no_expression.txt",sep="\t",row.names = F,quote = F)

#### script for drawing plot
dir.create(colnames(RBP.miRNA)[1])
setwd(paste0("C:/Rdata/20221108_correlaiton_between_RBP_level_and_miRNA_level_with_excluding_outlier/",colnames(RBP.miRNA)[1]))

cp <-ggplot(data=RBP.miRNA,aes(x=RBP.miRNA[,1],y=RBP.miRNA[,2]))+
  geom_point(aes(color=Site_Primary))+ #draw scatter plot
  stat_smooth(method = "lm", se = FALSE, colour = "black", size = 0.5)+ # draw regression line
  scale_color_manual(values =unique(RBP.miRNA[,4]))+ # color by site primary
  labs(title=paste0("r =",signif(r$estimate,3),", p = ",signif(r$p.value,3),", n = ",nrow(RBP.miRNA)),x=colnames(RBP.miRNA)[1]
       ,y=colnames(RBP.miRNA)[2])+ #write main title and la title
  guides(color = guide_legend(ncol = 1))+ #make legend 1line
  theme(legend.text = element_text(size=5)) + #change size 
  theme_bw()+
  theme(legend.background = element_rect(fill = "white", colour = "black"))

cp <-cp+
  theme(plot.title=element_text(hjust = 0.5),axis.title = element_text(size = 15),axis.text = element_text(size = 13),
        title = element_text(size = 15),legend.title = element_text(hjust = 0.5))


# draw plot
ggsave(paste0(colnames(RBP.miRNA)[1],"_vs_",colnames(RBP.miRNA)[2],".pdf"),width =10 ,height =7.5 )