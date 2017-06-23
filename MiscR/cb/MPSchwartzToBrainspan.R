source("http://bioconductor.org/biocLite.R")
install.packages("~/code/MT_BIOC/Enrichment-test/pkgs/EACI_0.0.1.tar.gz",repos=NULL,type = "source")
library(EACI)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(reshape2)
library(ExpressionAtlas)
library(biomaRt)
source("~/code/cb/LoadFunctions.R")
source("~/code/cb/LoadData.R")

aggregateRowCols <- function(q){
  m <-  t(apply(q,1, function(x){
    a <- aggregate(x,list(colnames(q)),mean)[,2]
    names(a) <- aggregate(x,list(colnames(q)),mean)[,1]
    a
  }))
  m <- m[,unique(colnames(q))]
  
  m <-  t(apply(m,2, function(x){
    a <- aggregate(x,list(rownames(m)),mean)[,2]
    names(a) <- aggregate(x,list(rownames(m)),mean)[,1]
    a
  }))
  m <- m[unique(rownames(m)),]
  t(m)
}

if(F){
BrainSpan <- read.csv2(file = "~/code/data/cb/NBNogData/brainspan_devcourse/expression_matrix.csv",sep = ",",header = F,row.names = 1)
BrainSpanCols <-  read.csv2(file = "~/code/data/cb/NBNogData/brainspan_devcourse/columns_metadata.csv",sep = ",",header = T,row.names = 1)
BrainSpanRows <-  read.csv2(file = "~/code/data/cb/NBNogData/brainspan_devcourse/rows_metadata.csv",sep = ",",header = T,row.names = 1)
BrainSpan <- as.matrix(BrainSpan)
BrainSpan <- matrix(as.numeric(BrainSpan) , nrow = nrow(BrainSpan), ncol = ncol(BrainSpan), dimnames =dimnames(BrainSpan))
rownames(BrainSpan) <-  make.names(BrainSpanRows$gene_symbol,unique = T)
colnames(BrainSpan) <- paste(BrainSpanCols$structure_acronym,BrainSpanCols$age)
BrainSpanFPKM <- BrainSpan
BrainSpan <- (sweep(BrainSpan, 2, colSums(BrainSpan), "/") * 10^6) 

rownames(BrainSpan) <- rownames(BrainSpanFPKM) <-toupper(rownames(BrainSpan))
BrainSpanCols$structure_name <- as.character(BrainSpanCols$structure_name)
BrainEA <- ExpressionAtlas::getAtlasExperiment("E-MTAB-4840")
BrainEA <- BrainEA[[names(BrainEA)]]
}

studerData <- read.csv2("~/code/data/cb/NBNogData/GSE56796.results.hg19/genes.no_mt.ec.tab", header=T, sep="\t",row.names=1)
studerData <- studerData[,-ncol(studerData)]
studerData <- as.matrix(studerData)
storage.mode(studerData) <- "numeric"

colnames(studerData) <- gsub("X[0-9]+_","",colnames(studerData))
colnames(studerData) <- gsub("d","",colnames(studerData))
colnames(studerData) <- gsub("Day.","",colnames(studerData))
colnames(studerData) <- gsub(".[A-Z+]","",colnames(studerData))
colnames(studerData) <- paste("CCON",colnames(studerData))
studerDataLo <- studerData[,25:ncol(studerData)]
studerData <- studerData[,1:24]

studerDataTPM <- read.csv2("~/code/data/cb/NBNogData/GSE56796.results.hg19/genes.no_mt.tpm.renorm.tab", header=T, sep="\t",row.names=1)
studerDataTPM <- studerDataTPM[,-ncol(studerDataTPM)]
studerDataTPM <- as.matrix(studerDataTPM)
storage.mode(studerDataTPM) <- "numeric"

colnames(studerDataTPM) <- gsub("X[0-9]+_","",colnames(studerDataTPM))
colnames(studerDataTPM) <- gsub("d","",colnames(studerDataTPM))
colnames(studerDataTPM) <- gsub("Day.","",colnames(studerDataTPM))
colnames(studerDataTPM) <- gsub(".[A-Z+]","",colnames(studerDataTPM))
colnames(studerDataTPM) <- paste("CCON",colnames(studerDataTPM))
studerDataTPMLo <- studerDataTPM[,25:ncol(studerDataTPM)]
studerDataTPM <- studerDataTPM[,1:24]


#? Mikes data
msData <- read.csv2("~/code/data/cb/NBNogData/MPSchwartzAnalysis/mstl069.TPMtxt.txt",sep = "\t",row.names = 1,header = T)
msData <- msData[,-ncol(msData)]
msData <- as.matrix(msData)
storage.mode(msData) <- "numeric"
rownames(msData) <- toupper(rownames(msData))
colnames(msData) <- gsub("X","",colnames(msData))
colnames(msData) <- gsub(".CB","",colnames(msData))
msData <- msData[,!grepl("SBNOG",colnames(msData))]
msData[,order(msCondits)]

msDataEC <- read.csv2("~/code/data/cb/NBNogData/MPSchwartzAnalysis/mstl069EC.txt",sep = "\t",row.names = 1,header = T)
msDataEC <- as.matrix(msDataEC)
storage.mode(msDataEC) <- "numeric"
msDataEC <- msDataEC[,!grepl("SBNOG",colnames(msDataEC))]
rownames(msDataEC) <- toupper(rownames(msDataEC))
colnames(msDataEC) <-colnames(msData)
msCondits <-  sapply(colnames(msDataEC), strsplit,"[.]")
msCondits <- apply(sapply(msCondits,"[",c(2,3,1)),2,function(x)Reduce(paste,x) )
msCondits <- gsub("SHPEG","XSHPEG",msCondits)
colnames(msDataEC) <- msCondits
msDataEC <- msDataEC[,order(msCondits)]

write.table(msDataEC[,order(msCondits)],file = "~/Desktop/EC.txt",sep = "\t" ,quote = F)

BEA <- BrainEAassay <- assay(BrainEA)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
rnSymbol <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),filters = c("ensembl_gene_id"),values =rownames(assay(BrainEA)) ,mart = mart) 
inds <- rnSymbol$hgnc_symbol%in%union(rownames(datasets[[1]]),msData)
BEA <-  BEA[rnSymbol$ensembl_gene_id,]
BEA <-  BEA[inds,]
colnames(BEA) <- paste(colData(BrainEA)$organism_part,colData(BrainEA)$developmental_stage)
BEA <- BEA[,order(colnames(BEA))]
rownames(BEA) <- rnSymbol$hgnc_symbol[inds]
BEA.sub <- BEA[,grepl("Carnegie Stage",colnames(BEA))]
colnames(BEA.sub) <- gsub("Carnegie Stage ","CS",colnames(BEA.sub))
BEA.sub <- BEA.sub[,!grepl("fragment",colnames(BEA.sub))]


msExtraDataEC <- read.csv2("~/code/data/cb/NBNogData/MPSchwartzAnalysis/3D_Neural_Construct_samples_EC.txt",sep = "\t",row.names = 1,header = T)
msExtraDataEC <- as.matrix(msExtraDataEC)
storage.mode(msExtraDataEC) <- "numeric"
rownames(msExtraDataEC) <- toupper(rownames(msExtraDataEC))
msExtraConditsMat <-  sapply(colnames(msExtraDataEC), strsplit,"_")
msExtraCondits <- toupper(apply(sapply(msExtraConditsMat,"[",1:5),2,function(x)Reduce(paste,x) ))
msExtraConditsMat <- sapply(msExtraConditsMat,"[",1:6)
colnames(msExtraDataEC) <- msExtraCondits

msExtraData2EC <- read.csv2("~/code/data/cb/NBNogData/MPSchwartzAnalysis/d21multicomponentneuralconstructEC.txt",sep = "\t",row.names = 1,header = T)
msExtraData2EC <- as.matrix(msExtraData2EC)
storage.mode(msExtraData2EC) <- "numeric"
rownames(msExtraData2EC) <- toupper(rownames(msExtraData2EC))
msExtraConditsMat2 <-  sapply(colnames(msExtraData2EC), strsplit,"_")
msExtraCondits2 <- toupper(apply(sapply(msExtraConditsMat2,"[",1:5),2,function(x)Reduce(paste,x) ))
colnames(msExtraData2EC) <- gsub("\ S[1-9]|NA|E[1-9]","",msExtraCondits2)

organoidDataRPKM <- read.csv2("~/code/data/cb/NBNogData/MPSchwartzAnalysis/GSE80264_RNA-Seq_results_Ensembl.gene.rpkm.txt",sep = "\t",row.names = 1,header = T)
organoidDataRPKM <- as.matrix(organoidDataRPKM)
storage.mode(organoidDataRPKM) <- "numeric"
rownames(organoidDataRPKM) <- toupper(rownames(organoidDataRPKM))

organoidDataEC <- read.csv2("~/code/data/cb/NBNogData/MPSchwartzAnalysis/GSE80264.results.hg19/genes.no_mt.ec.tab",sep = "\t",row.names = 1,header = F)
organoidDataEC <- as.matrix(organoidDataEC )
organoidDataEC <- organoidDataEC[,-ncol(organoidDataEC)]
storage.mode(organoidDataEC ) <- "numeric"
rownames(organoidDataEC ) <- toupper(rownames(organoidDataEC ))
colnames(organoidDataEC) <- colnames(organoidDataRPKM)


#cbData is from first tera control, published
cbData <- datasets[[1]]
colnames(cbData) <- tpDatasets[[1]]

cbDataEC <- read.csv2("~/code/data/cb/NBNogData/MPSchwartzAnalysis/GSE90053_H1_EC.txt",sep = "\t",row.names = 1,header = T)
#cbDataEC <- cbDataEC[,-ncol(cbDataEC)]
cbDataEC <- as.matrix(cbDataEC)
storage.mode(cbDataEC) <- "numeric"
rownames(cbDataEC) <- toupper(rownames(cbDataEC))
colnames(cbDataEC) <- gsub("X","",colnames(cbDataEC))
colnames(cbDataEC) <- tpDatasets[[1]]

cbDataEC <- cbDataEC[,1:23]
cbData <- cbData[,1:23]
cbDsub <- cbDataEC[,c("0","2","4","6","8","10","14","20","24","30","34","36")]

BSS <-  md(BrainSpan[,order(BrainSpanCols$structure_name)],rownames(msData))
BSS <-BSS[,grepl( "pcw",colnames(BSS))]
BSS <-BSS[,grepl( "\ 8|\ 9|10|11|12|13|14|15",colnames(BSS))]

pdf("~/Desktop/MPS_EDA.pdf")
std.heatmap(cor.compare(cbData,msData,min = 1,method="spe"))
std.heatmap(cor.compare(studerData,msData,min = 1,method="spe"))
std.heatmap(cor.compare(BSS,BEA.sub,min = 0,method="spe"),cexRow=.5,cexCol=.5)
std.heatmap(cor.compare(BEA.sub,msDataEC,min = 1),cexRow=.5)
std.heatmap(cor.compare(BSS,cbData,min = 1,method="spe"),cexRow=.5)
std.heatmap(cor.compare(BEA.sub,cbDataEC,min = 1),cexRow=.5)
dev.off()

sort(rowMeans(cor.compare(BEA.sub,msData,min = 0,method="spe")),T)
sort(apply(cor.compare(BSS,msData,min = 0,method="spe"),1,max),T)

sort(rowMeans(cor.compare(BEA.sub,cbData,min = 0,method="spe")),T)
sort(apply(cor.compare(BSS,cbData,min = 0,method="spe"),1,max),T)

rnmBSS <- rn.merge(t(t(sort(apply(cor.compare(BSS,cbData,min = 1,method="spe"),1,max),T))), t(t(sort(apply(cor.compare(BSS,msData,min = 1,method="spe"),1,max),T))))
rnmBEA <- rn.merge(t(t(sort(apply(cor.compare(BEA.sub,cbDataEC,min = 1,method="spe"),1,max),T))), t(t(sort(apply(cor.compare(BEA.sub,msDataEC,min = 1,method="spe"),1,max),T))))
colnames(rnmBEA) <- c("CB","MPS")
colnames(rnmBSS) <- c("CB","MPS")

write.table(rnmBEA[order(rnmBEA[,1]),],"~/Desktop/HDBRMaxs.txt",sep="\t",col.names =T,row.names = T,quote = F)
write.table(rnmBSS[order(rnmBSS[,2]),],"~/Desktop/BrainspanMaxs.txt",sep="\t",col.names =T,row.names = T,quote = F)


BSSCors <- rn.merge(cor.compare(BSS,cbData,min = 1,method="spe"),cor.compare(BSS,msData,min = 1,method="spe"))
iii <- Reduce(intersect,list(rownames(cbDataEC),rownames(BEA.sub),rownames(msDataEC)))
#HDBRCors <-  rn.merge(cor.compare(BEA.sub[iii,],cbDataEC[iii,],min = .00000001,method="spe"), cor.compare(BEA.sub[iii,],msDataEC[iii,],min = 1,method="spe"))
HDBRCors <-  cor.compare(BEA.sub[iii,],cbind( cbDataEC[iii,],msDataEC[iii,]),min = .00000001,method="spe")

iii <- Reduce(intersect,list(rownames(msExtraDataEC),rownames(BEA.sub),rownames(msDataEC)))
#HDBRCors <-  rn.merge(cor.compare(BEA.sub[iii,],cbDataEC[iii,],min = .00000001,method="spe"), cor.compare(BEA.sub[iii,],msDataEC[iii,],min = 1,method="spe"))
HDBRXCors <-  cor.compare(BEA.sub[iii,],cbind( msExtraDataEC[iii,],msDataEC[iii,]),min = .00000001,method="spe")

iii <- Reduce(intersect,list(rownames(studerData),rownames(BEA.sub),rownames(msDataEC),rownames(cbDsub)))
HDBRSCors <-  cor.compare(BEA.sub[iii,],cbind( cbDsub[iii,],studerData[iii,],msDataEC[iii,]),min = .00000001,method="spe")

iii <- Reduce(intersect,list(rownames(studerDataTPM),rownames(BSS),rownames(msData),rownames(cbDsub)))
BSSSCors <-  cor.compare(BSS[iii,],cbind( cbData[iii,],studerDataTPM[iii,],msData[iii,]),min = .00000001,method="spe")

iii <- Reduce(intersect,list(rownames(organoidDataEC),rownames(BEA.sub),rownames(msDataEC),rownames(studerDataTPM),rownames(msExtraData2EC),rownames(cbDataEC)))
HDBROCors <- cor.compare(BEA.sub[iii,],cbind(cbDsub[iii,],studerData[iii,], organoidDataEC[iii,3:6],msExtraData2EC[iii,],msDataEC[iii,]),min = .000000001,method="spe") 
std.heatmap(aggregateRowCols(HDBROCors),cexCol=.3,cexRow=.3)


std.heatmap(HDBROCors,cexCol=.3,cexRow=.3)

f <- factor(msExtraConditsMat[4,])
cdit <- rbind(gsub("d|D","",msExtraConditsMat[2,]),
gsub("NPC|new","",msExtraConditsMat[3,]),
sapply(f, function(x) which(levels(f) == x)),
as.integer(!is.na(msExtraConditsMat[5,]=="EC")) )
colnames(cdit) <- 1:ncol(cdit)

ggplot(data = melt(cdit))+
  geom_text(aes(x=Var2*30,y=Var1*30,label=value),size=1)+
  theme_void()





msExtraConditsMat[4,]

write.table(HDBRCors,"~/Desktop/HDBRCors.txt",sep="\t",col.names =T,row.names = T,quote = F)
write.table(BSSCors,"~/Desktop/BrainspanCors.txt",sep="\t",col.names =T,row.names = T,quote = F)
write.table(cor(BEA.sub,BEA.sub,min = .0001,method="spe"),"~/Desktop/HDBRvsSelf.txt",sep="\t",col.names =T,row.names = T,quote = F)


iii <- Reduce(intersect,list(rownames(cbDataEC),rownames(BEA),rownames(msDataEC)))
#HDBRCors <-  rn.merge(cor.compare(BEA.sub[iii,],cbDataEC[iii,],min = .00000001,method="spe"), cor.compare(BEA.sub[iii,],msDataEC[iii,],min = 1,method="spe"))
HDBRCors <-  cor.compare(BEA[iii,],cbind( cbDataEC[iii,],msDataEC[iii,]),min = .0000000001,method="spe")

meandHDBRCors <- aggregateRowCols(HDBRCors)
meandHDBRCors <- t(meandHDBRCors[, unique(colnames(BEA))])
meandHDBRCorsF <- t(meandHDBRCors[grepl("forebrain",rownames(meandHDBRCors)),])
colnames(meandHDBRCorsF) <- gsub("forebrain\ ","", colnames(meandHDBRCorsF))
meandHDBRCorsF <- meandHDBRCorsF[,order(colnames(meandHDBRCorsF))]
meandHDBRCorsF <- meandHDBRCorsF[-seq(1,18,by=2),]
std.heatmap(meandHDBRCorsF[c(1:14,16,15,17,18,20,19,21,22),],cexCol=.5,cexRow=.4,main="Forebrain Spearman Correlations")


iii<- Reduce(intersect,list(rownames(BEA),rownames(msDataEC)),rownames(cbDataEC))
BEA.fore <- BEA[iii,grepl("forebrain",colnames(BEA))&!grepl("midbrain|fragment",colnames(BEA))]
BEAf.rank <- apply(BEA.fore,2,rank)
colnames(BEAf.rank)[19]
sort(BEAf.rank[,11]-rank(msData[iii,1]),decreasing = T)
sort((BEAf.rank[,19]-rank(msData[iii,1])) - (BEAf.rank[,11]-rank(msData[iii,1])))
BEA.fore["EOMES",]
msData["EOMES",]
cbData["EOMES",]
ea13 <- eacitest(BEAf.rank[,11]-rank(msData[iii,1]),lib = "org.Hs.eg",idtype = "SYMBOL")
ea19 <- eacitest(BEAf.rank[,19]-rank(msData[iii,1]),lib = "org.Hs.eg",idtype = "SYMBOL")
ea1319 <- eacitest(BEAf.rank[,19]-BEAf.rank[,11],lib = "org.Hs.eg",idtype = "SYMBOL")
eaDiff <- eacitest((BEAf.rank[,19]-rank(msDataEC[iii,1])) - (BEAf.rank[,11]-rank(msDataEC[iii,1])),lib = "org.Hs.eg",idtype = "SYMBOL")
ea2DDiff <- eacitest((BEAf.rank[,19]-rank(cbDataEC[iii,23])) - (BEAf.rank[,19]-rank(msDataEC[iii,1])),lib = "org.Hs.eg",idtype = "SYMBOL")
ea2DDiv <- eacitest(1/signLog((BEAf.rank[,19]-rank(cbDataEC[iii,23]))) - signLog(BEAf.rank[,19]-rank(msDataEC[iii,1])) ,lib = "org.Hs.eg",idtype = "SYMBOL")

signLog <- function(x,base=10){ log(abs(x)+1,base)*sign(x)}

ea13$setscores[order(ea13$setscores$set.mean,decreasing = F),"Term"]
ea19$setscores[order(ea19$setscores$set.mean,decreasing = T),"Term"]
ea1319$setscores[order(ea1319$setscores$set.mean,decreasing = T),"Term"]

ea2DDiff$setscores[order(ea2DDiff$setscores$set.mean,decreasing = T)[1:20],"Term"]
ea2DDiv$setscores[order(ea2DDiv$setscores$set.mean,decreasing = F)[1:20],"Term"]


rnmBEA[order(rnmBEA[,1]),]
rnmBSS[order(rnmBSS[,2]),]








AutOrganoid <-  getGEO("GSE61476",file = "~/Desktop/GSE61476_RAW/GSM1505861_Sample_1123-01_23_4_31_004.counts.txt.gz")


AutOrganoid <-  getGEO("GSE57595",AnnotGPL = T)
AO <-  exprs(AutOrganoid[[1]])
a <-  featureData(AutOrganoid[[1]])
a <- pData(a)
rownames(AO) <- a$`Gene symbol`
AO <- AO[nchar(rownames(AO))>0,]
