library(EBSeq)
library(gplots)
library(RColorBrewer)
library(zoo)
source('~/code/cb/LoadFunctions.R')
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
datapath <- "~/code/data/cb/IncTempComparison/"
setwd(datapath)

data <- read.csv2(paste0(datapath,"genes.no_mt.ec.tab"), header=T, sep="\t",row.names=1)
data <- data[,-ncol(data)]
data <- as.matrix(data)
storage.mode(data) <- "numeric"
g = barplot(colSums(data),main = "Total Number of Reads")
text(g,colSums(data)-100000 ,colSums(data),cex = .5) 
data <- GetNormalizedMat(data,MedianNorm(data))
colnames(data) <- gsub("X","",colnames(data))
gs=rownames(data)[grepl(pattern = "HOX",rownames(data))]
gs=gs[!grepl('P|R|S',gs)]
gs=gs[order(as.numeric(gsub("HOX[A-Z]","",gs)))]
gs=gs[order(gsub("[0-9]\1{,2}","",gs))]
heatmap.2(data[gs,],trace = "none",Rowv = F,Colv = F,col = cols,srtCol = 45, main = "Hox Genes \n Log2 Normalized ECs")
neural_list <- c(neural_list,c("DLK1","LIX1","NR2F1","GLI1","GLI2","HES5"))
neural_list <- c(neural_list,gs)

heatmap.2(cor(data[,order(conditions[2,])],method = 's'),col = cols,trace='none',Rowv = F,Colv = F)

conditions <- sapply(colnames(data),function(x) unlist(strsplit(x,'_')))
conditions[3,] <- as.integer(gsub('d','',conditions[3,]))
datasets <- lapply(unique(conditions[2,]), function(x){
  data[,conditions[2,]==x]
} )
names(datasets) <- unique(conditions[2,])
tps <- lapply(unique(conditions[2,]), function(x){
  q=conditions[3,]
  q[conditions[2,]==x]
} )
names(tps) <- unique(conditions[2,])

write.pdf({
  mypar(3,4)
  for(gn in neural_list[neural_list%in% rownames(datasets[[1]])]){
    zooz <-  lapply(names(datasets), function(x){zoo(datasets[[x]][gn,],as.numeric(tps[[x]]))})
    #print(merge(zooz[[1]],zooz[[2]]))
    plot.zoo(Reduce( merge, zooz),plot.type = "single",col = 1:5, lty = 1,main=paste(gn),xlab = "days",ylab = "Normalized Counts" )
    legend("topleft", unique(conditions[2,]), col = 1:5, lty = 1)
  }}
  ,filename = "~/Desktop/IncTempComparisonNeurals.pdf")

write.pdf({
  mypar(3,4)
  for(gn in cc_list[cc_list%in% rownames(datasets[[1]])]){
    zooz <-  lapply(names(datasets), function(x){zoo(datasets[[x]][gn,],as.numeric(tps[[x]]))})
    #print(merge(zooz[[1]],zooz[[2]]))
    plot.zoo(Reduce( merge, zooz),plot.type = "single",col = 1:5, lty = 1,main=paste(gn),xlab = "days",ylab = "Normalized Counts" )
    legend("topleft", unique(conditions[2,]), col = 1:5, lty = 1)
  }}
  ,filename = "~/Desktop/IncTempComparisonCC.pdf")





write.table(neural_list,file = '~/Desktop/20161111NeuralSelections.tsv',sep = '\t',quote = F,row.names = F)
