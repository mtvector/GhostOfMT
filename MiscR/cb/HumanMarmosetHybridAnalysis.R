source("~/code/cb/LoadFunctions.R")
source("http://bioconductor.org/biocLite.R")
library(zoo)
library(EBSeq)
library(biomaRt)
library(reshape2)
neural_list <- toupper(sort(as.character(read.csv("~/code/data/cb/markers/FullNeuralGeneList.csv",header = T)[,1])))
neural_related_genes <- make.names(toupper(as.character(read.csv("~/code/data/cb/markers/neuralrelatedgenes.csv",header = T)[,1])))
pluripotent_genes <- toupper(sort(as.character(read.csv("~/code/data/general/marker_lists/ES_Germ/ES.txt",header = T)[,1])))
tpmdata <- read.csv2("~/code/data/cb/H1vCJvHybridDiff/1_TPM_all_samples.txt", header=T, sep="\t",row.names=1)
data <- read.csv2("~/code/data/cb/H1vCJvHybridDiff/1_ExpectCounts_all_samples 2.txt", header=T, sep="\t",row.names=1)
data <- as.matrix(data[,which(colnames(data)!='Transcript_ID')])
data <- matrix(as.numeric(data) , nrow = nrow(data), ncol = ncol(data), dimnames =dimnames(data))
data <- data[,colSums(data)>1e5]
tpmdata <- as.matrix(tpmdata[,which(colnames(tpmdata)!='Transcript_ID')])
tpmdata <- matrix(as.numeric(tpmdata) , nrow = nrow(tpmdata), ncol = ncol(tpmdata), dimnames =dimnames(tpmdata))
tpmdata <- tpmdata[,colSums(tpmdata)>1e5]

condits <- colnames(data)

#condits <- gsub(pattern = "_[0-9]+","",condits)
withtp <- grepl(pattern = "_d[0-9]+",condits)
conditswithtp <- condits[withtp]
tp <- withtp
tp[!withtp]<-0
tp[withtp] <- sapply(strsplit(condits,"_d+")[withtp],"[[",2)
tp <- as.numeric(tp)

condits <- sapply(strsplit(condits,"_d"),"[[",1)

humanData <- data[grepl("hg19",rownames(data)),]
marmData <- data[grepl("CJ367",rownames(data)),]
rownames(humanData) <- make.names(gsub("\\|hg19","",rownames(humanData)))
rownames(marmData) <-make.names(gsub("\\|CJ367","",rownames(marmData)))

humantpmData <- tpmdata[grepl("hg19",rownames(tpmdata)),]
marmtpmData <- tpmdata[grepl("CJ367",rownames(tpmdata)),]
rownames(humantpmData) <- make.names(gsub("\\|hg19","",rownames(humantpmData)))
rownames(marmtpmData) <- make.names(gsub("\\|CJ367","",rownames(marmtpmData)))

rawHumanData <- humanData
rawMarmData <- marmData

humanData <- GetNormalizedMat(humanData,MedianNorm(data))
marmData <- GetNormalizedMat(marmData,MedianNorm(data))

std.heatmap(cor(humanData[,tp==0]))
std.heatmap(cor(marmData[,tp==0]))
std.heatmap(cor(marmData))

neural_selections <- read.csv2(file.path("~/code/data/cb/NBNogData/MixingTest/CustomGenes.txt"),sep = "\t",header = T)
neural_selections <- toupper(as.character(neural_selections[,1]))

figcondit <- unique(condits[withtp])
figcondit <- figcondit[-c(2,3)]

barplot2(rbind(colSums(rawMarmData),colSums(rawHumanData)),las=2,cex.names = 1.1,ylab = "Reads",col = 1:2)
legend("topleft", c("Marmoset","Human"),fill = 1:2)

pdf(file = "~/Desktop/FusionNeurals.pdf")
mypar(1,1)
for(g in c(neural_list,pluripotent_genes)){
  if(g%in% rownames(marmtpmData) & g %in% rownames(humantpmData)){
    marmg <- lapply(1:length(figcondit),function(x){
      zoo(marmData[g,condits==figcondit[x]],tp[condits==figcondit[x]] ) 
    })
    humg <- lapply(1:length(figcondit),function(x){
      zoo(humanData[g,condits==figcondit[x]],tp[condits==figcondit[x]] ) 
    })
    #plot(na.approx(Reduce(merge.zoo, humg)),plot.type = "single",col = 1:10,type = "l",ylab = "TPM", xlab = "Day",main = paste("Human",g))
    #legend("topleft", figcondit , col = 1:10, lty = 1)
    #plot(na.approx(Reduce(merge.zoo, marmg)),plot.type = "single",col = 1:10,type = "l",ylab = "TPM", xlab = "Day",main = paste("Marmoset",g))
    #legend("topleft", figcondit , col = 1:10, lty = 1)
    plot(na.approx(Reduce(merge.zoo, c(marmg,humg ))),plot.type = "single",col = rep(1:5,2),lty = c(rep(1,5),rep(2,5)),type = "l",ylab = "NormalizedCounts", xlab = "Day",main = paste(g))
    legend("topleft", paste(c(rep("marm",5), rep("human",5)) ,rep(figcondit,2)) ,lty =c(rep(1,5),rep(2,5))  , col = rep(1:5,2))
  }
}
dev.off()

humanData["PAX6",condits=="Human.4"]


mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "cjacchus_gene_ensembl")
MarmResults1 <- getBM(attributes = c( 'hgnc_symbol','chromosome_name'), filters = 'hgnc_symbol',
                 values = rownames(marmData)[!grepl("ENSCJAG",rownames(marmData))] , mart = mart)

MarmResults2 <- getBM(attributes = c( 'ensembl_gene_id','chromosome_name'), filters = 'ensembl_gene_id',
                      values = rownames(marmData)[grepl("ENSCJAG",rownames(marmData))] , mart = mart)

MarmResults <- rbind(MarmResults1,setNames(MarmResults2, names(MarmResults1)))

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")
HumanResults <- getBM(attributes = c('hgnc_symbol', 'chromosome_name'), filters = 'hgnc_symbol',
                 values = rownames(humanData) , mart = mart)

HumanResults<- HumanResults[HumanResults$hgnc_symbol%in% rownames(rawHumanData),]

table(HumanResults$chromosome_name)[nchar(names(table(HumanResults$chromosome_name)))<3]
table(MarmResults$chromosome_name)[nchar(names(table(MarmResults$chromosome_name)))<3]

#If hybrid5_1 were to be missing a chromosome 6

dat <- rawMarmData
res <- MarmResults
containcondit <- (grepl("CJ",condits) | grepl("Hybrid",condits))

matt1 <-  sapply(names(table(res$chromosome_name))[nchar(names(table(res$chromosome_name)))<3], function(n){
  #mat <- colSums(dat[res[res$chromosome_name==n,]$hgnc_symbol, containcondit])/colSums(dat[,containcondit])
  mat <- log(colSums(dat[res[res$chromosome_name==n,]$hgnc_symbol, containcondit])/colSums(dat[,containcondit]),2)
  #std.heatmap(outer(mat,mat,FUN = "-"),main=n)
})
matt1 <- t(matt1)
heatmap.2((matt1-rowMins(matt1)),col=cols,trace="none", main="Marmoset Chromosomal Expression\nlog2(Fraction of reads)\nScaled as difference from rowMin")
#heatmap.2((matt1-rowMeans(matt1))/rowSds(matt1),col=cols,trace="none")

dat <- rawHumanData
res <- HumanResults
containcondit <- (grepl("Human",condits) | grepl("Hybrid",condits))

#dat[res[res$chromosome_name==6,]$hgnc_symbol, "Hybrid5_1"] <- dat[res[res$chromosome_name==6,]$hgnc_symbol,"Hybrid5_1"]/2

matt2 <-  sapply(names(table(res$chromosome_name))[nchar(names(table(res$chromosome_name)))<3], function(n){
  #mat <- colSums(dat[res[res$chromosome_name==n,]$hgnc_symbol, containcondit])/colSums(dat[,containcondit])
  mat <- log(colSums(dat[res[res$chromosome_name==n,]$hgnc_symbol, containcondit])/colSums(dat[,containcondit]),2)
  #std.heatmap(outer(mat,mat,FUN = "-"),main=n)
})
matt2 <- t(matt2)
heatmap.2((matt2-rowMins(matt2)),col=cols,trace="none", main="Human Chromosomal Expression\nlog2(Fraction of reads)\nScaled as difference from rowMin")
#heatmap.2((matt2-rowMeans(matt2))/rowSds(matt2),col=cols,trace="none")
#heatmap.2(rbind((matt1-rowMins(matt1))/rowSds(matt1) ,(matt2-rowMins(matt2))/rowSds(matt2)),col=cols,trace="none")

colnames(matt2) <- gsub("_[0-9]+","",condits[containcondit])
chrom <- factor(rownames(matt))
spec <- factor(colnames(matt)[1:6])
rowttests(matt[,1:6],spec[1:6])$p.value*nrow(matt)
mmatt <- melt(matt)
l <- lm(value~Var1*Var2 ,data = as.data.frame(mmatt))
anova(l)

marmData[MarmResults[MarmResults$chromosome_name=="Y",]$hgnc_symbol,]
humanData[HumanResults[HumanResults$chromosome_name=="Y",]$hgnc_symbol,]


marmTC <- rawMarmData[,4:12]
marmTP <- tp[4:12]
marmTC <- marmTC[,order(marmTP)]
marmTP <- sort(marmTP)
marmTC <- GetNormalizedMat(marmTC,MedianNorm(marmTC))

plot(marmTP,marmTC["PAX6",])

colnames(marmTC) <- gsub("X","",colnames(marmTC))
gs=rownames(marmTC)[grepl(pattern = "HOX",rownames(marmTC))]
gs=gs[!grepl('P|R|S',gs)]
gs=gs[order(as.numeric(gsub("HOX[A-Z]","",gs)))]
gs=gs[order(gsub("[0-9]\1{,2}","",gs))]
heatmap.2(marmTC[gs,],trace = "none",Rowv = F,Colv = F,col = cols,srtCol = 45, main = "Hox Genes \n Log2 Normalized ECs")

gs=neural_selections[neural_selections%in%rownames(marmTC)]
heatmap.2(log(marmTC[gs,]+1,2),trace = "none",Rowv = F,Colv = F,col = cols,srtCol = 45, main = "Hox Genes \n Log2 Normalized ECs")

gs=neural_selections[neural_selections%in%rownames(datasets[[11]])]
heatmap.2(log(datasets[[11]][gs,1:33]+1,2),trace = "none",Rowv = F,Colv = F,col = cols,srtCol = 45, main = "Hox Genes \n Log2 Normalized ECs")
