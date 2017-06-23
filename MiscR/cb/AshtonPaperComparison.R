library(EBSeq)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)

timeconverth <- function(tp){
  tp[grepl('d', tp)] <- as.integer(gsub('d','',tp[grepl('d', tp)]))*24
  tp[grepl('h', tp)] <- as.integer(gsub('h','',tp[grepl('h', tp)]))
  as.integer(tp)
}

humandatapath <- "~/code/data/cb/AshtonPaperComparison/Human/"
setwd(humandatapath)
humandata <- read.csv2(paste0(humandatapath,"genes.no_mt.ec.tab"), header=T, sep="\t",row.names=1)
humandata <- humandata[,-ncol(humandata)]
humandata <- as.matrix(humandata)
storage.mode(humandata) <- "numeric"
g = barplot(colSums(humandata),main = "Total Number of Reads")
text(g,colSums(humandata)-100000 ,colSums(humandata),cex = .5) 
humandata <- GetNormalizedMat(humandata,MedianNorm(humandata))
colnames(humandata) <- gsub("X","",colnames(humandata))
gs=rownames(humandata)[grepl(pattern = "HOX",rownames(humandata))]
gs=gs[!grepl('P|R|S',gs)]
gs=gs[order(as.numeric(gsub("HOX[A-Z]","",gs)))]
gs=gs[order(gsub("[0-9]\1{,2}","",gs))]
heatmap.2(humandata[gs,],trace = "none",Rowv = F,Colv = F,col = cols,srtCol = 45, main = "Hox Genes \n Log2 Normalized ECs")
curTP <- colnames(humandata)
curTP <- timeconverth(curTP)
humanTP <- curTP
curMat <- humandata

plotlist <- lapply( neural_selections,function(gn){
  if(gn %in% rownames(curMat)){
    gg <-  ggplot(as.data.frame(cbind("tme"=curTP, "value"=curMat[gn,])), aes(x=tme, y=value)) + 
      geom_line() +
      geom_point(size=.8) +
      theme_bw(base_size = 6) +
      theme(legend.position="none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      #scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
      labs(x = "Hour", y = "Log2Count", 
           title = gn,family="arial",size=6)
    gg
  }
})
pdf(paste0("~/Desktop/","HumanCHIR",".pdf"), onefile = T)
marrangeGrob(grobs=Filter(Negate(is.null), plotlist) , nrow=4, ncol=4)
dev.off()



mousedatapath <- "~/code/data/cb/AshtonPaperComparison/Mouse/"
setwd(mousedatapath)

mousedata <- read.csv2(paste0(mousedatapath,"genes.no_mt.ec.tab"), header=T, sep="\t")
mousedata <- as.matrix(mousedata)
rownames(mousedata) <- mousedata[,1]
mousedata <- mousedata[!duplicated(mousedata[,1]),]
mousedata <- mousedata[,-c(1,ncol(mousedata))]
mousedata <- matrix(as.numeric(mousedata) , nrow = nrow(mousedata), ncol = ncol(mousedata), dimnames =dimnames(mousedata))
rownames(mousedata) <- toupper(rownames(mousedata))
g = barplot(colSums(mousedata),main = "Total Number of Reads")
text(g,colSums(mousedata)-100000 ,colSums(mousedata),cex = .5) 
mousedata <- GetNormalizedMat(mousedata,MedianNorm(mousedata))
colnames(mousedata) <- gsub("X","",colnames(mousedata))
gs=rownames(mousedata)[grepl(pattern = "HOX",rownames(mousedata))]
gs=gs[!grepl('P|R|S',gs)]
gs=gs[order(as.numeric(gsub("HOX[A-Z]","",gs)))]
gs=gs[order(gsub("[0-9]\1{,2}","",gs))]
heatmap.2(mousedata[gs,],trace = "none",Rowv = F,Colv = F,col = cols,srtCol = 45, main = "Hox Genes \n Log2 Normalized ECs")
curTP <- gsub("mEpi8_","",colnames(mousedata))
curTP <- timeconverth(curTP)
mouseTP <- curTP
curMat <- mousedata

plotlist <- lapply( neural_selections,function(gn){
  if(gn %in% rownames(curMat)){
    gg <-  ggplot(as.data.frame(cbind("tme"=curTP, "value"=curMat[gn,])), aes(x=tme, y=value)) + 
      geom_line() +
      geom_point(size=.8) +
      theme_bw(base_size = 6) +
      theme(legend.position="none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      #scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
      labs(x = "Hour", y = "Log2Count", 
           title = gn,family="arial",size=6)
    gg
  }
})
pdf(paste0("~/Desktop/","MouseCHIR",".pdf"), onefile = T)
marrangeGrob(grobs=Filter(Negate(is.null), plotlist) , nrow=4, ncol=4)
dev.off()

std.heatmap(cor.compare(humandata,mousedata,method="spearman"))

mouseScores <- peakScore(mousedata,mouseTP)
humanScores <- peakScore(humandata,humanTP)

sort(mouseScores,T)
sort(humanScores,T)

gnz <- union(names(humanScores),names(mouseScores))

pdf(file = "~/Desktop/humanpeaksCHIR.pdf")
mypar(2,1)
for(g in names(sort(humanScores,T))[1:48]){
  plot(humanTP,humandata[g,],xlab="hrs",ylab="normalized counts",main=paste("Human",g,humanScores[g]))
}
dev.off()

pdf(file = "~/Desktop/mousepeaksCHIR.pdf")
mypar(2,1)
for(g in names(sort(mouseScores,T))[1:48]){
  plot(mouseTP,mousedata[g,],xlab="hrs",ylab="normalized counts",main=paste("Mouse",g,mouseScores[g]))
}
dev.off()


pdf(file = "~/Desktop/sharedpeaksCHIR.pdf")
mypar(2,1)
for(g in names(sort(apply(cbind(humanScores[gnz],mouseScores[gnz]),1,min),T))[1:48]){
  plot(mouseTP,mousedata[g,],xlab="hrs",ylab="normalized counts",main=paste("Mouse",g,mouseScores[g]))
  plot(humanTP,humandata[g,],xlab="hrs",ylab="normalized counts",main=paste("Human",g,humanScores[g]))
}
dev.off()
