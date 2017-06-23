#Axolotl-DanioRerio Divergence 429.6 MYA
#Axo-Human Divergence 355.7 MYA
path <- "~/code/data/axo"
loadLibraries()


#load blastema files
blastemaFile <- file.path(path, "tanaka_norrna_rsem_human_gene_expressionForUnixJustJuvTCECs.txt")
ECBAll <- read.table(blastemaFile,header=TRUE)
rownames(ECBAll) <-ECBAll$symbol
ECB<- as.matrix(ECBAll[,9:ncol(ECBAll)])
sizesB <- MedianNorm(ECB[rowMeans(ECB)>5,])
NormECB <- GetNormalizedMat(ECB[rowMeans(ECB)>5,], sizesB)


#Load Human Alignment Data File
stagesH <- c("Stage_1","Stage_1","Stage_2","Stage_2","Stage_2", "Stage_3","Stage_3","Stage_3","Stage_4","Stage_4","Stage_4","Stage_5","Stage_5","Stage_5","Stage_6",	"Stage_6",	"Stage_6",	"Stage_7",	"Stage_7",	"Stage_7","Stage_8",	"Stage_8",	"Stage_8","Stage_9",	"Stage_9",	"Stage_9",	"Stage_10",	"Stage_10",	"Stage_10",	"Stage_11",	"Stage_11",	"Stage_12","Stage_12","Stage_12","Stage_14","Stage_14","Stage_14","Stage_16","Stage_16","Stage_16","Stage_19","Stage_19","Stage_19","Stage_24","Stage_24","Stage_24","Stage_40","Stage_40", "Stage_40")
stagesH <- factor(stagesH, levels = c("Stage_1","Stage_1","Stage_2","Stage_2","Stage_2", "Stage_3","Stage_3","Stage_3","Stage_4","Stage_4","Stage_4","Stage_5","Stage_5","Stage_5","Stage_6",	"Stage_6",	"Stage_6",	"Stage_7",	"Stage_7",	"Stage_7","Stage_8",	"Stage_8",	"Stage_8","Stage_9",	"Stage_9",	"Stage_9",	"Stage_10",	"Stage_10",	"Stage_10",	"Stage_11",	"Stage_11",	"Stage_12","Stage_12","Stage_12","Stage_14","Stage_14","Stage_14","Stage_16","Stage_16","Stage_16","Stage_19","Stage_19","Stage_19","Stage_24","Stage_24","Stage_24","Stage_40","Stage_40", "Stage_40"))


filenameH = file.path(path,"Human_Annotation/ExpectCounts_all_samples.txt")
#Expected Counts All Samples
ECHAll <- read.table(filenameH,header=TRUE)
#Assign rownames as gene names
rownames(ECHAll) <- ECHAll[,1]
ECH <- as.matrix(ECHAll[-1:-2])
ECHnonavg <- ECH[rowMeans(ECH)>5,]
nonavgsizes <- MedianNorm(ECHnonavg)
normNonAvgH <-  GetNormalizedMat(ECHnonavg, nonavgsizes)

NormMatchH <- matrix(nrow=0,ncol=ncol(ECH))
NormMatchB <- matrix(nrow=0,ncol=ncol(ECB))
for(i in 1:nrow(ECH)){
  A <- ECH[i,]
  B <- ECB[which(rownames(ECB)==ECHAll[i,1]),]
  if(length(B) > 0 && mean(A)>5 && mean(B)>5){
    NormMatchH <-rbind(NormMatchH,A)
    NormMatchB <-rbind(NormMatchB,B)
    rownames(NormMatchH)[nrow(NormMatchH)] <-rownames(ECH)[i] 
    rownames(NormMatchB)[nrow(NormMatchB)] <-rownames(ECH)[i]
  }
}
MNH <- MedianNorm(NormMatchH)
MNB <- MedianNorm(NormMatchB)
NormMatchH <- GetNormalizedMat(NormMatchH,MNH)
NormMatchB <- GetNormalizedMat(NormMatchB,MNB)
NormMatchH <- sapply(1:nrow(NormMatchH), function(i){
  tapply(NormMatchH[i,], stagesH, mean, simplify=TRUE)
})
NormMatchH <- t(NormMatchH)
NormMatchH <-  NormMatchH[, colSums(is.na(NormMatchH)) != nrow(NormMatchH)]
rownames(NormMatchH) <- rownames(NormMatchB)
colnames(NormMatchB) <- c("0h","3h", "6h","12h","1d","3d","5d","7d","10d", "14d","21d", "28d")

#Get the means of comparable columns
ECH <- sapply(1:nrow(ECH), function(i){
  tapply(ECH[i,], stagesH, mean, simplify=TRUE)
})
ECH <- t(ECH)
ECH <-  ECH[, colSums(is.na(ECH)) != nrow(ECH)]
rownames(ECH) <- rownames(ECHAll)
conditH <- factor(unique(stagesH))

dim(matchH)
#Create unNormalized matrices of the genes expressed in both Blastema and H aligned embryo
matchH <- matrix(nrow=0,ncol=ncol(ECH))
matchB <- matrix(nrow=0,ncol=ncol(ECB))
for(i in 1:nrow(ECH)){
  A <- ECH[i,]
  B <- ECB[which(rownames(ECB)==ECHAll[i,1]),]
  if(length(B) > 0 && mean(A)>5 && mean(B)>5){
    matchH <-rbind(matchH,A)
    matchB <-rbind(matchB,B)
    rownames(matchH)[nrow(matchH)] <-rownames(ECH)[i] 
    rownames(matchB)[nrow(matchB)] <-rownames(ECH)[i]
  }
}
colnames(matchB) <- c("0h","3h", "6h","12h","1d","3d","5d","7d","10d", "14d","21d", "28d")

timepointsH <- c("1:2","2:3","3:4","4:5","5:6","6:7","7:8","8:9","9:10","10:11","11:12","12:14","14:16","16:19","19:24","24:40")
timepointsB <- c("0h:3h","3h:6h", "6h:12h","12h:1d","1d:3d","3d:5d","5d:7d","7d:10d","10d:14d", "14d:21d","21d:28d")
timepointsB <- c("0h","3h", "6h","2h","1d","3d","5d","7d","10d", "14d","21d", "28d")
timepointsH <- c("1","2","3","4","5","6","7","8","9","10","11","12","14","16","19","24","40")
#timepointsH <- c("1:2","2:3","3:4","4:5","5:6","6:7","7:8","8:9","9:10","10:11","11:12","12:14","14:16","16:19","19:24","24:40")

TFs <- c("ZIC2", "USP44" , "ATF3", "EGR1", "ETS2", "FOS", "FOXO1","JUN", "JUND", "KLF4", "KLF6", "MYC", "ZFP36", "SOX2", "POU5F1", "NANOG")
HMG <- c("HMGA1","SOX11", "SOX4", "HMGA2")
Limbs <- c("PRRX1", "HOXD10", "PRDM1", "SALL1", "TBX18", "SALL4", "HOXD11", "GLI3", "SALL3", "TGFB1", "TNC", "SHH", "EED", "APC", "SMAD4", "JARID2", "ZIC2", "HMGA2", "SUFU")
oncogenes <- read.table( file.path(path,"MSKCCOncogenes.txt"), header = TRUE, sep = "\t")
oncogenes <- as.character(oncogenes$Gene.Symbol)

###############DIFFERENTIAL EXPRESSION###############

#Subset to only 2 conditions, and test consecutive subsets
ECHsubFrame <- as.data.frame(ECH)
EBOutECH <- sapply(1:(ncol(ECH)-1), function(i){
  condit <- factor(colnames(ECHsubFrame[(i):(i+1)]))  
  subECH <- ECH[,(i):(i+1)]
  Sizes=MedianNorm(subECH)
  EBOut <- EBTest(subECH, Conditions = condit,sizeFactors = Sizes, maxround=5)
  return(EBOut)
})
EBOutResultsH <- sapply(1:(ncol(ECH)-1), function(i){
  EBDERes=GetDEResults(EBOutECH[,i], FDR=0.05)
})
mypar(1,1)
#The number of genes differentially expressed between consecutive time points i and i+1
barplot2(sapply(1:(ncol(ECH)-1), function(i) length(EBOutResultsH[,i]$DEfound)), ylab="# of DE Genes",xlab="Timepoint Stage", names.arg = timepointsH , main= "Embryo Consecutive Step DE")

####DOING THE SAME FOR BLASTEMA
#Subset to only 2 conditions, and test consecutive subsets
# ECBsubFrame <- as.data.frame(ECB)
# EBOutECB <- sapply(1:(ncol(ECB)-1), function(i){
#   condit <- factor(colnames(ECBsubFrame[(i):(i+1)]))  
#   subECB <- ECB[,(i):(i+1)]
#   Sizes=MedianNorm(subECB)
#   EBOut <- EBTest(subECB, Conditions = condit,sizeFactors = Sizes, maxround=5)
#   return(EBOut)
# })
# EBOutResultsB <- sapply(1:(ncol(ECB)-1), function(i){
#   EBDERes=GetDEResults(EBOutECB[,i], FDR=0.05)
# })

#Subset to only 2 conditions, and test consecutive subsets
ECBsubFrame <- as.data.frame(ECB)
EBOutECB <- sapply(seq(2,12), function(i){
  condit <- factor(colnames(ECBsubFrame[c((1),(i))]))  
  subECB <- ECB[,c((1),(i))]
  Sizes=MedianNorm(subECB)
  EBOut <- EBTest(subECB, Conditions = condit,sizeFactors = Sizes, maxround=5)
  return(EBOut)
})
EBOutResultsB <- sapply(1:(ncol(ECB)-1), function(i){
  EBDERes=GetDEResults(EBOutECB[,i], FDR=0.05)
})




mypar(1,1)
#The number of genes differentially expressed between consecutive time points i and i+1
barplot2(sapply(1:(ncol(ECB)-1), function(i) length(EBOutResultsB[,i]$DEfound)), ylab="# of DE Genes",xlab="Timepoint Index", names.arg =timepointsB , main= "Blastema Consecutive Step DE")


foldChangesH <-sapply(1:ncol(EBOutECH), function(i){
  PostFC(EBOutECH[,i])
})

foldChangesB <-sapply(1:ncol(EBOutECB), function(i){
  PostFC(EBOutECB[,i])
})  
####################DESEQ2#######################

ECHi <- matrix(as.integer(as.matrix(ECHAll[-1:-2])), nrow=nrow(ECHAll[-1:-2]),ncol=ncol(ECHAll[-1:-2]))
ECHi
stagesH <- c("Stage_1","Stage_1","Stage_2","Stage_2","Stage_2", "Stage_3","Stage_3","Stage_3","Stage_4","Stage_4","Stage_4","Stage_5","Stage_5","Stage_5","Stage_6",	"Stage_6",	"Stage_6",	"Stage_7",	"Stage_7",	"Stage_7","Stage_8",	"Stage_8",	"Stage_8","Stage_9",	"Stage_9",	"Stage_9",	"Stage_10",	"Stage_10",	"Stage_10",	"Stage_11",	"Stage_11",	"Stage_12","Stage_12","Stage_12","Stage_14","Stage_14","Stage_14","Stage_16","Stage_16","Stage_16","Stage_19","Stage_19","Stage_19","Stage_24","Stage_24","Stage_24","Stage_40","Stage_40", "Stage_40")
conditHi <- factor(stagesH)
rownames(ECHi) <- rownames(ECHAll)
colnames(ECHi) <- colnames(ECHAll[-1:-2])
ddsH <- DESeqDataSetFromMatrix(ECHi, DataFrame(conditHi), ~ conditHi)
# standard analysis
ddsH <- DESeq(ddsH)
resH<- sapply(1:(length(conditH)-1), function(i){
  results(ddsH, contrast = c("conditHi", as.character(conditH[i]),  as.character(conditH[i+1])), independentFiltering = T)
})

ECBi <- matrix(as.integer(ECB), nrow=nrow(ECB),ncol=ncol(ECB))
rownames(ECBi) <- rownames(ECB)
colnames(ECBi) <- colnames(ECB)
conditB <- factor(colnames(matchB))
ddsB <- DESeqDataSetFromMatrix(ECBi, DataFrame(conditB), ~ conditB)
# standard analysis
ddsB <- DESeq(ddsB)
resB<- sapply(1:(ncol(ECBi)-1), function(i){
   results(ddsB, contrast = c("conditB", as.character(conditB[i]),  as.character(conditB[i+1])))
})
?"DESeq"
results(ddsB)

DEUpHi <- sapply(1:(length(conditH)-1), function(i) {-resH[[i]]$log2FoldChange > 1.4})
DEDownHi <- sapply(1:(length(conditH)-1), function(i) {-resH[[i]]$log2FoldChange < -1.4})
DEUpBi <- sapply(1:(ncol(ECBi)-1), function(i) {-resB[i][[1]]$log2FoldChange > 1.4})
DEDownBi <- sapply(1:(ncol(ECBi)-1), function(i) {-resB[i][[1]]$log2FoldChange < -1.4})

DEUpHi
DEUpBi
sapply(1:(ncol(ECBi)-1), function(i) {rownames(ECBi)[which(DEUpBi[,i])]})
DEUpHi <- sapply(1:(length(conditH)-1), function(i) {resH[[i]]$pval < .05})
DEDownHi <- sapply(1:(length(conditH)-1), function(i) {resH[i][[1]]$pval < .05})
DEUpBi <- sapply(1:(ncol(ECBi)-1), function(i) {resB[i][[1]]$pval < .05})
DEDownBi <- sapply(1:(ncol(ECBi)-1), function(i) {resB[i][[1]]$pval < .05 })

DEUpHi <- sapply(1:(length(conditH)-1), function(i) {rownames(ECH)[which(DEUpHi[,i])]})
DEDownHi <- sapply(1:(length(conditH)-1), function(i) {rownames(ECH)[which(DEDownHi[,i])]})
DEUpBi <- sapply(1:(ncol(ECBi)-1), function(i) {rownames(ECBi)[which(DEUpBi[,i])]})
DEDownBi <- sapply(1:(ncol(ECBi)-1), function(i) {rownames(ECBi)[which(DEDownBi[,i])]})

DEUpHi
mypar(1,1)
barplot2(sapply(1:(ncol(ECBi)-1), function(i) {sum(-resB[i][[1]]$log2FoldChange >1.3, na.rm = T)}) , ylab="# of DE Genes",xlab="Timepoint Index", names.arg =timepointsB , main= "Blastema Consecutive Step (+) DE")
barplot2(sapply(1:(ncol(ECBi)-1), function(i) {sum(-resB[i][[1]]$log2FoldChange < -1.3, na.rm = T)}) , ylab="# of DE Genes",xlab="Timepoint Index", names.arg =timepointsB , main= "Blastema Consecutive Step (-) DE")
barplot2(sapply(1:(length(conditH)-1), function(i) {sum(resH[[i]]$log2FoldChange < -1.3, na.rm = T)}) , ylab="# of DE Genes",xlab="Timepoint Index", names.arg =timepointsH , main= "Embryo Consecutive Step (-) DE")
barplot2(sapply(1:(length(conditH)-1), function(i) {sum(resH[[i]]$log2FoldChange > 1.3, na.rm = T)}) , ylab="# of DE Genes",xlab="Timepoint Index", names.arg =timepointsH , main= "Embryo Consecutive Step (+) DE")




#GENERATE TOTALS AND MATCHES FOR UP
upMatchGenes <-  getDEMatchListMatrix(DEUpHi, DEUpBi)
upMatch <- matrix(sapply(upMatchGenes, FUN = length), nrow = length(DEUpHi),ncol=length(DEUpBi))
totalUpH <-  matrix(sapply(DEUpHi, FUN = length), nrow = length(DEUpHi),ncol=length(DEUpBi))
totalUpB <- matrix(sapply(DEUpBi, FUN = length), nrow = length(DEUpHi),ncol=length(DEUpBi), byrow = TRUE)

#DO THE SAME FOR DOWN
downMatchGenes <-  getDEMatchListMatrix(DEDownHi, DEDownBi)
downMatch <- matrix(sapply(downMatchGenes, FUN = length), nrow = length(DEDownHi),ncol=length(DEDownBi))
totalDownH <-  matrix(sapply(DEDownHi, FUN = length), nrow = length(DEDownHi),ncol=length(DEDownBi))
totalDownB <- matrix(sapply(DEDownBi, FUN = length), nrow = length(DEDownHi),ncol=length(DEDownBi), byrow = TRUE)
#rename columns

colnames (upMatch) <- timepointsB
rownames(upMatch) <- timepointsH
colnames (downMatch) <- timepointsB
rownames(downMatch) <- timepointsH
############DE MATCHING################

  DEUpH <- getDEUpDown(foldChangesH, EBOutResultsH, up=TRUE,plusthreshold =2)
  DEDownH <- getDEUpDown(foldChangesH, EBOutResultsH, up=FALSE, negthreshold = 1/2)
  DEUpB <- getDEUpDown(foldChangesB, EBOutResultsB, up=TRUE,plusthreshold =2)
  DEDownB <- getDEUpDown(foldChangesB, EBOutResultsB,up=FALSE, negthreshold = 1/2)

  DEUpB[[1]]
  foldChangesB[,1]$RealFC[EBOutResultsB[,1]$DEfound]
  
  stageNumH <- c("1","2","3","4","5","6","7","8","9","10","11","12","14","16","19","24","40")
  DEPath <- file.path(path,"DE")
  DEUpH <- sapply(seq(2, length(stageNumH)), function(i){
    read.table( file.path(DEPath, "Up", paste0("Stage_",stageNumH[i] , "_Vs_",stageNumH[1],".Up.txt")), header = TRUE, sep = "\t")$Gene.ID
  })
  DEDownH <- sapply(seq(2, length(stageNumH)), function(i){
    read.table( file.path(DEPath, "Down", paste0("Stage_", stageNumH[i] , "_Vs_",stageNumH[1],".Down.txt")), header = TRUE, sep = "\t")$Gene.ID
  })
  vignette("EBSeq_Vignette")
  
  barplot2(sapply(1:(ncol(ECH)-1), function(i) length(DEUpH[[i]])), ylab="# of DE Genes",xlab="Timepoint Stage", names.arg = timepointsH[-1] , main= "Embryo Vs0 DE")
  
  DEUpH[[14]]
  
  #GENERATE TOTALS AND MATCHES FOR UP
  upMatchGenes <-  getDEMatchListMatrix(DEUpH, DEUpB)
  upMatch <- matrix(sapply(upMatchGenes, FUN = length), nrow = length(DEUpH),ncol=length(DEUpB))
  totalUpH <-  matrix(sapply(DEUpH, FUN = length), nrow = length(DEUpH),ncol=length(DEUpB))
  totalUpB <- matrix(sapply(DEUpB, FUN = length), nrow = length(DEUpH),ncol=length(DEUpB), byrow = TRUE)
  #DO THE SAME FOR DOWN
  downMatchGenes <-  getDEMatchListMatrix(DEDownH, DEDownB)
  downMatch <- matrix(sapply(downMatchGenes, FUN = length), nrow = length(DEDownH),ncol=length(DEDownB))
  totalDownH <-  matrix(sapply(DEDownH, FUN = length), nrow = length(DEDownH),ncol=length(DEDownB))
  totalDownB <- matrix(sapply(DEDownB, FUN = length), nrow = length(DEDownH),ncol=length(DEDownB), byrow = TRUE)
  #rename columns
  colnames (upMatch) <- timepointsB
  rownames(upMatch) <- timepointsH
  colnames (downMatch) <- timepointsB
  rownames(downMatch) <- timepointsH
  colnames (upMatchGenes) <- timepointsB
  rownames(upMatchGenes) <- timepointsH
  colnames (downMatchGenes) <- timepointsB
  rownames(downMatchGenes) <- timepointsH
  #####FOR VS0
  colnames (upMatch) <- timepointsB[-1]
  rownames(upMatch) <- timepointsH[-1]
  colnames (downMatch) <- timepointsB[-1]
  rownames(downMatch) <- timepointsH[-1]
  colnames (upMatchGenes) <- timepointsB[-1]
  rownames(upMatchGenes) <- timepointsH[-1]
  colnames (downMatchGenes) <- timepointsB[-1]
  rownames(downMatchGenes) <- timepointsH[-1]
  
  
  #######CROSSING UP AND DOWN
  xHUBDg <- getDEMatchListMatrix(DEUpH, DEDownB)
  xHDBUg <- getDEMatchListMatrix(DEDownH, DEUpB)
  xHUBD <- matrix(sapply(xHUBDg, FUN = length), nrow = length(DEUpH),ncol=length(DEUpB))
  xHDBU <-  matrix(sapply(xHDBUg, FUN = length), nrow = length(DEUpH),ncol=length(DEUpB))
  colnames (xHUBDg) <- timepointsB
  rownames(xHUBDg) <- timepointsH
  colnames (xHDBUg) <- timepointsB
  rownames(xHDBUg) <- timepointsH
  #for vs0
  colnames (xHUBDg) <- timepointsB[-1]
  rownames(xHUBDg) <- timepointsH[-1]
  colnames (xHDBUg) <- timepointsB[-1]
  rownames(xHDBUg) <- timepointsH[-1]
  
  HG.HUBD <- hypergeoMatch(xHUBD, totalUpB ,totalUpH, numG=nrow(matchH))
  HG.HUBD <- matrix(p.adjust(HG.HUBD, method="BH"), nrow=nrow(HG.HUBD), ncol=ncol(HG.HUBD))
  HG.HDBU <- hypergeoMatch(xHDBU, totalUpB ,totalUpH, numG=nrow(matchH))
  HG.HDBU <- matrix(p.adjust(HG.HDBU, method="BH"), nrow=nrow(HG.HDBU), ncol=ncol(HG.HDBU))
  colnames (HG.HUBD) <- timepointsB[-1]
  rownames(HG.HUBD) <- timepointsH[-1]
  colnames (HG.HDBU) <- timepointsB[-1]
  rownames(HG.HDBU) <- timepointsH[-1]
  pvalPlot(HG.HUBD, "Fisher's DE Overlap \n EmbryoUpBlastemaDown p value")
  pvalPlot(HG.HDBU, "Fisher's DE Overlap \n EmbryoDownBlastemaUp p value")
  heatmapMatch(xHUBD[1:16,1:11], positive =F, cluster = FALSE, main="EBSeq, # Matches \nEmbryoUpBlastemaDown")
  heatmapMatch(xHDBU[1:16,1:11], positive =TRUE, cluster = FALSE, main="EBSeq, # Matches \n EmbryoDownBlastemaUp")

  #Hypergeometric Probability of Overlap significance
  HG.pos <- hypergeoMatch(upMatch, totalUpB ,totalUpH, numG=nrow(matchH))
  HG.pos <- matrix(p.adjust(HG.pos, method="BH"), nrow=nrow(HG.pos), ncol=ncol(HG.pos))
  colnames (HG.pos) <- timepointsB[-1]
  rownames(HG.pos) <- timepointsH[-1]
  heatmapMatch(upMatch[1:16,1:11], positive =TRUE, cluster = FALSE, main="EBSeq, # Matches")
  pvalPlot(HG.pos ,mainT="(+) DE Overlap\n FET P-Value")
  
  
  heatmapMatch(-log( HG.pos,10), positive =TRUE, cluster = FALSE, main="EBSeq, # Matches")
  
  HG.neg <- hypergeoMatch(downMatch, totalDownB ,totalDownH, numG=nrow(matchH))
  HG.neg <- matrix(p.adjust(HG.neg, method="BH"), nrow=nrow(HG.neg), ncol=ncol(HG.neg))
  colnames (HG.neg) <- timepointsB[-1]
  rownames(HG.neg) <- timepointsH[-1]
  pvalPlot(HG.neg, mainT = "(-) DE Overlap\n FET P-Value") 
  heatmapMatch(downMatch[1:16,1:11], positive =F, cluster = FALSE, main="EBSeq, # Matches")
  
  stepnameH <- c("1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12","12-14","14-16","16-19","19-24","24-40")
  stepnameB <- c("0h-3h","3h-6h", "6h-12h","12h-1d","1d-3d","3d-5d","5d-7d","7d-10d","10d-14d", "14d-21d","21d-28d")
  
  
  
  l <- c("EUBU", "EDBD")
  writePath <- "/Users/mschmitz/Desktop/figs/vs0DE/zscores/geneLists"
  xx <- c(list(zscoreDEUp), list(zscoreDEDown))
  xp <- c(list(HG.zUp), list(HG.zDown))
  names(xp) <- l
  names(xx) <- l
  stepnameH <- c("1","2","3","4","5","6","7","8","9","10","11","12","14","16","19","24","40")
  stepnameB <- c("0h","3h", "6h","12h","1d","3d","5d","7d","10d", "14d","21d","28d")
  
  
  writePath <- "/Users/mschmitz/Desktop/figs/vs0DE/vs0DE/geneLists"
  l <- c("EDBU", "EUBD", "EUBU", "EDBD")
  xx <- c(list(xHDBUg), list(xHUBDg), list(upMatchGenes), list(downMatchGenes))
  xp <- c(list(HG.HDBU), list(HG.HUBD), list(HG.pos), list(HG.neg))
  names(xp) <- l
  names(xx) <- l
  stepnameH <- stepnameH[-1]
  stepnameB <- stepnameB[-1]

  for(ll in l){
    x <- xx[ll][[1]]
    pvals <- xp[ll][[1]]
    ind <- which(pvals<.05, arr.ind=TRUE)
    if(nrow(ind)>0){
    for(i in 1:nrow(ind)){
      write.table(as.data.frame(x[ind[i,1],ind[i,2]]), file = paste0(writePath,"/", ll,"/", stepnameH[ind[i,1]],"_", stepnameB[ind[i,2]],".csv"), row.names = F, col.names = F,sep="\t")
    }}
  }
  
  
##################ONCOGENES LIST######################
  
  DEOncoUpH <- sapply(1:length(DEUpH),function(i){
    sum( oncogenes %in% DEUpH[[i]])
  } )
  
  DEOncoDownH <- sapply(1:length(DEDownH),function(i){
    sum( oncogenes %in% DEDownH[[i]])
  } )
  DEOncoUpB <- sapply(1:length(DEUpB),function(i){
    sum( oncogenes %in% DEUpB[[i]])
  } )
  
  DEOncoDownB <- sapply(1:length(DEDownB),function(i){
    sum( oncogenes %in% DEDownB[[i]])
  } )
  
  DEUpH
  barplot2(DEUpH, ylab="# of DE Genes",xlab="Stage", names.arg = timepointsH , main= "Upregulated DE Oncogenes in Embryo")
  barplot2(DEOncoDownH, ylab="# of DE Genes",xlab="Stage", names.arg = timepointsH , main= "Downregulated DE Oncogenes in Embryo")
  barplot2(DEOncoUpB, ylab="# of DE Genes",xlab="Stage", names.arg = timepointsB , main= "Upregulated DE Oncogenes in Blastema")
  barplot2(DEOncoDownB, ylab="# of DE Genes",xlab="Stage", names.arg = timepointsB , main= "Downregulated DE Oncogenes in Blastema")
  
  
  
  
  upEBonco <- getDEMatchRefined(upMatchGenes, oncogenes)
  upEBonco <- upEBonco
  colnames(upEBonco) <- timepointsB
  rownames(upEBonco) <- timepointsH
  downEBonco <- getDEMatchRefined(downMatchGenes, oncogenes)
  downEBonco <- downEBonco
  colnames (downEBonco) <- timepointsB
  rownames(downEBonco) <- timepointsH
  
  heatmapMatch(upMatch[1:16,1:11], positive =TRUE, cluster = FALSE, main="EBSeq Matches")
  heatmapMatch(upEBonco[1:16,1:11], positive =TRUE, cluster = FALSE, main="EBSeq Matches \n Oncogenes")
  heatmapMatch(upEBonco[1:16,1:11]/upMatch[1:16,1:11], positive =TRUE, cluster = FALSE, main="EBSeq Matches \n Oncogenes")

  heatmapMatch(downMatch[1:16,1:11], positive =F, cluster = FALSE, main="EBSeq Matches")
  heatmapMatch(downEBonco[1:16,1:11], positive =F, cluster = FALSE, main="EBSeq Matches \n Oncogenes")
  heatmapMatch(downEBonco[1:16,1:11]/downMatch[1:16,1:11], positive =F, cluster = FALSE, main="EBSeq Matches \n Oncogenes")
  
  #Look at the go terms represented
  
  barplot2(sapply(DEUpH, length))

  index <- unlist(xHUBDg[9,1])
  mypar(5,4)
  for(i in 1:length(index)){
    x <- index[i]
    #plot(ECH[x,1:14], main=paste("Embryo", x),xlab="TimePoint",ylab="EC")
    plot(matchH[x,1:14], main=paste("Embryo", x),xlab="TimePoint",ylab="EC")
    plot(matchB[x,], main= paste("Blastema",x),xlab="TimePoint",ylab="EC")
  }
  index
  terms <- sapply(index, function(i){
    terms <- getGOHuman(index)
  })
  printTopGO(terms$name_1006, 10)
  terms[,2]$name_1006
  
########################Test using Z Scores###################
  
  install.packages("mixOmics")
  library(mixOmics)
  
  MNHB <- MedianNorm(cbind(matchH, matchB))
  NormMatchHB <- GetNormalizedMat(cbind(matchH, matchB),MNHB)

  
  cim(NormMatchHB)

  bc1 <- biclust(cbind(zscoreUpH, zscoreUpB),method=BCBimax(), minr=2, minc=2, number=100)
  bc2 <- biclust(zscoreH,method=BCBimax(), minr=2, minc=2, number=100)
  
  heatmapBC(cbind(zscoreUpH, zscoreUpB),bc1,axes=T)
  bubbleplot(zscoreH, bc1, bc2)
    ?heatmapBC
  
  listA <- NormMatchH
  listB <- NormMatchB
  mypar(1,1)
  upperLim=1.2
  lowerLim=-1.2
  zscoreH <-  (listA-rowMedians(listA))/rowSds(listA)
  zscoreB <-  (listB-rowMedians(listB))/rowSds(listB)
  zscoreUpH <- zscoreH > upperLim
  zscoreUpB <- zscoreB > upperLim
  zscoreDownH <- zscoreH < lowerLim
  zscoreDownB <- zscoreB < lowerLim
  zscoreDEUp <-  zscoreDE(zscoreUpH, zscoreUpB)
  zscoreDEDown <-  zscoreDE(zscoreDownH, zscoreDownB)
  zscoreUpMatch <- matrix(sapply(zscoreDEUp, FUN = length), nrow = nrow(zscoreDEUp),ncol=ncol(zscoreDEUp))
  zscoreDownMatch <- matrix(sapply(zscoreDEDown, FUN = length), nrow = nrow(zscoreDEDown),ncol=ncol(zscoreDEDown))
  colnames (zscoreUpMatch) <- colnames(listB)
  rownames(zscoreUpMatch) <- colnames(listA)
  totalUpBNum <- matrix(  colSums(zscoreUpH, na.rm = T),nrow=ncol(listA),ncol = ncol(listB),byrow = F)
  totalUpHNum <- matrix(colSums(zscoreUpB,na.rm=T),nrow=ncol(listA),ncol = ncol(listB), byrow = T)
  totalDownBNum <- matrix(  colSums(zscoreDownH, na.rm = T),nrow=ncol(listA),ncol = ncol(listB),byrow = F)
  totalDownHNum <- matrix(colSums(zscoreDownB,na.rm=T),nrow=ncol(listA),ncol = ncol(listB), byrow = T)

  HG.zUp <- hypergeoMatch(zscoreUpMatch,totalUpHNum, totalUpBNum ,numG=nrow(listA))
  HG.zUp <- matrix(p.adjust(HG.zUp, method="BH"), nrow=nrow(HG.zUp), ncol=ncol(HG.zUp))
  colnames (HG.zUp) <- colnames(listB)
  rownames(HG.zUp) <- conditH
  pvalPlot(HG.zUp ,"Zscore Up Pvalues")
  heatmapMatch(zscoreUpMatch, positive =TRUE, cluster = FALSE, main="Zscore Matches")
  
  HG.zDown <- hypergeoMatch(zscoreDownMatch,totalDownHNum, totalDownBNum ,numG=nrow(listA))
  HG.zDown <- matrix(p.adjust(HG.zDown, method="BH"), nrow=nrow(HG.zDown), ncol=ncol(HG.zDown))
  colnames (HG.zDown) <- colnames(listB)
  rownames(HG.zDown) <- conditH
  colnames (zscoreDownMatch) <- colnames(listB)
  rownames(zscoreDownMatch) <- conditH
  pvalPlot(HG.zDown ,"Zscore Down Pvalues")
  heatmapMatch(zscoreDownMatch, positive =F, cluster = FALSE, main="Zscore Matches")

  
  
  ind <- which(HG.zUp<.05, arr.ind=TRUE)
  index <- unlist(zscoreDEUp[1,1])
  index
  mypar(5,4)
  for(i in 1:length(index)){
    x <- index[i]
    plot(matchH[x,], main=paste("Embryo", x),xlab="TimePoint",ylab="ExpCounts")
    abline(median(matchH[x,]),0)
    abline(median(matchH[x,])+1.2*sd(matchH[x,]),0, col="red")
    abline(median(matchH[x,])-1.2*sd(matchH[x,]),0, col="red")
    plot(matchB[x,], main= paste("Blastema", x),xlab="TimePoint",ylab="ExpCounts")
    abline(median(matchB[x,]),0)
    abline(median(matchB[x,])+1.2*sd(matchB[x,]),0, col="red")
    abline(median(matchB[x,])-1.2*sd(matchB[x,]),0, col="red")
  }

  o <- getDEMatchRefined(zscoreDEUp, oncogenes)
  rownames(o) <- colnames(matchH)
  colnames(o) <- colnames(matchB)
  heatmapMatch(o, positive =T, cluster = FALSE, main="Zscore Matches \n Oncogenes")
  
  o <- getDEMatchRefined(zscoreDEDown, oncogenes)
  rownames(o) <- colnames(matchH)
  colnames(o) <- colnames(matchB)
  heatmapMatch(o, positive =F, cluster = FALSE, main="Zscore Matches \n Oncogenes")
  
  
  mypar(1,1)
  barplot2(sapply(1:ncol(zscoreUpH), function(i){ sum(na.omit(zscoreUpH[,i]))}), ylab="# of DE Genes",xlab="Timepoint Stage", names.arg = colnames(zscoreUpH) , main= "Embryo Zscore>1.2")
  barplot2(sapply(1:ncol(zscoreUpB), function(i){ sum(na.omit(zscoreUpB[,i]))}), ylab="# of DE Genes",xlab="Timepoint Stage", names.arg = colnames(zscoreUpB) , main= "Blastema Zscore>1.2")
  barplot2(sapply(1:ncol(zscoreDownH), function(i){ sum(na.omit(zscoreDownH[,i]))}), ylab="# of DE Genes",xlab="Timepoint Stage", names.arg = colnames(zscoreDownH) , main= "Embryo Zscore<-1.2")
  barplot2(sapply(1:ncol(zscoreDownB), function(i){ sum(na.omit(zscoreDownB[,i]))}), ylab="# of DE Genes",xlab="Timepoint Stage", names.arg = colnames(zscoreDownB) , main= "Blastema Zscore<-1.2")
  

  corEBz <- cor(na.omit(cbind(zscoreH,zscoreB)),method="spearman")
  corEBs <- cor(NormMatchHB,method="spearman")
  corEBp <- cor(log(NormMatchHB+1,2),method="pearson")

  

  mybreaks <- c(0,seq(.2, .7, length.out = 30) , 1)
  cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(30)
  cols <- c(cols, "black")
  heatmap.2(corEBs, main="Blastema vs Embryo \nSpearman Correlation"
            ,na.rm = TRUE, trace = "none",
            srtRow=0, srtCol=60,
            col = cols, ## using your colors
            dendrogram = "none", Rowv = FALSE, Colv = FALSE,
            breaks = mybreaks)   
  
  heatmap.2(corEBp2, main="Blastema vs Embryo \nPearson Correlation"
            ,na.rm = TRUE, trace = "none",
            srtRow=0, srtCol=60,
            col = cols, ## using your colors
            dendrogram = "none", Rowv = FALSE, Colv = FALSE,
            breaks = mybreaks)    
  
 
  
  
  
  qqplot(cbind(NormMatchH, NormMatchB), NormMatchHB)
  
# MUSCLE index <- upMatchGenes[1,7][[1]]
index <- unlist(zscoreDEUp[9,10])
index

mypar(5,4)
for(i in 1:length(index)){
  A <-log(NormMatchH+1,2) 
  B <- log(matchH +1,2)
  x <- index[i]
  plot(A[x,], main=paste("Norm", x),xlab="TimePoint",ylab="Z Score")
  abline(median(A[x,]),0)
  abline(median(A[x,])+1.2*sd(A[x,]),0, col="red")
  abline(median(A[x,])-1.2*sd(A[x,]),0, col="red")
  plot(B[x,], main= paste("unNorm", x),xlab="TimePoint",ylab="Z Score")
  abline(median(B[x,]),0)
  abline(median(B[x,])+1.2*sd(B[x,]),0, col="red")
  abline(median(B[x,])-1.2*sd(B[x,]),0, col="red")
}

###########GENIE3################
weight.matrixE <- get.weight.matrix(normNonAvgH)
link.listE <- get.link.list(weight.matrixE, report.max=1000)
weight.matrixB <- get.weight.matrix(NormMatchB)
link.listB <- get.link.list(weight.matrixB, report.max=1000)

# weight.matrixEB <- get.weight.matrix(NormMatchHB)
# link.listEB <- get.link.list(weight.matrixEB, report.max=1000)
genieCRNE <- coregnet(weight.matrixE, expressionDATA =NormMatchH, TFlis)
display(crnE)
crnEB <- coregnet(weight.matrixEB, expressionDATA =NormMatchHB)
display(crnEB)

browseVignettes("CoRegNet")


TFs <- read.table(paste0(path,"/my_human_TF_v5_final_11282011.txt"),header=T, sep="\t")
?hLICORN
# regnetEB <- hLICORN( NormMatchHB ,TFlist = TFs$symbol, minCoregSupport=.08,minGeneSupport=0.08)
# regnetE <- hLICORN( NormMatchH,TFlist = TFs$symbol, minCoregSupport=.08,minGeneSupport=0.08)
regnetB <- hLICORN( NormECB,TFlist = TFs$symbol,minCoregSupport=.08,minGeneSupport=0.08)
regnetEna <- hLICORN( normNonAvgH,TFlist = TFs$symbol, minCoregSupport=.08,minGeneSupport=0.08)

# tfaEB <- regulatorInfluence(regnetEB,NormMatchHB)
# tfaE <- regulatorInfluence(regnetE,NormMatchH)
tfaB <- regulatorInfluence(regnetB,NormECB)
tfaEna <- regulatorInfluence(regnetEna,normNonAvgH)


# display(regnetEB, NormMatchHB,TFA=tfaEB)
# display(regnetE, NormMatchH,TFA=tfaE)
display(regnetB, NormECB,TFA=tfaB)
display(regnetEna, normNonAvgH,TFA=tfaEna)

regulatorInfluence(regnetEB,NormMatchHB)
display(regnetEB, NormMatchHB)

genieCRNE <- coregnet(weight.matrixE, expressionDATA =NormECE)

?refine
#########OLD#######################################################################


colnames (totalUpH) <- timepointsB
rownames(totalUpH) <- timepointsH
colnames (totalUpB) <- timepointsB
rownames(totalUpB) <- timepointsH
colnames (totalDownB) <- timepointsB
rownames(totalDownB) <- timepointsH
colnames (totalDownH) <- timepointsB
rownames(totalDownH) <- timepointsH

heatmap.2(totalUpB, xlab = "Blastema Time Point", ylab="Embryo Time Point", main="(+) Differential Expression Blastema"
          ,na.rm = TRUE, trace = "none",colsep=1:11,
          rowsep=1:16,
          col=cols)


heatmap.2(totalDownH, xlab = "Blastema Time Point", ylab="Embryo Time Point", main="(-) Differential Expression Embryo"
          ,na.rm = TRUE, trace = "none",colsep=1:11,
          rowsep=1:16,
          col=cols)

heatmap.2(totalDownB, xlab = "Blastema Time Point", ylab="Embryo Time Point", main="(-) Differential Expression Blastema"
          ,na.rm = TRUE, trace = "none",colsep=1:11,
          rowsep=1:16,
          col=cols)

DEMatchList[8,which(timepointsB=="14d")]
which.max(DEMatchs[,1:11])
DEMatchList[56]
mypar(5,4)
index <- unlist(DEMatchList[8,2])
for(i in 1:length(index)){
  x <- index[i]
  plot(matchH[x,], main=paste("Embryo", x),xlab="TimePoint",ylab="TPM")
  plot(matchB[x,], main= paste("Blastema",x),xlab="TimePoint",ylab="TPM")
}

terms <- sapply(timepointsB, function(i){
  terms <- getGOHuman(DEMatchList[11,which(timepointsB==i)])
})



getGOHuman(DEMatchList[11,which(timepointsB=="14d")])
printTopGO(terms$name_1006, 70)

limbMatch <- getDEMatchRefined(DEMatchList, Limbs)
colnames (limbMatch) <- timepointsB
rownames(limbMatch) <- timepointsH
heatmap.2(limbMatch, xlab = "Blastema Time Point", ylab="Embryo Time Point", main="# of DE\n Limb Gene Matches"
          ,na.rm = TRUE, trace = "none",colsep=1:12,
          rowsep=1:16,
          col=cols)

#############PCA#################

mypar(1,1)
s <- svd(ECB)
plot(s$d^2/sum(s$d^2),ylab="% variance explained",xlab="Principal component", main="PCs & Variance")

#PCA
mypar(3,4)
for(i in 1:12)
  boxplot(split(s$v[,i], factor(colnames(ECB))), main=c("PC ",i))
#What happens if we simply remove the top two PC from the data?
D <- s$d; D[1:2]<-0 ##take out first 2
cleandat <- sweep(s$u,2,D,"*")%*%t(s$v)
plot(cleandat)

#Reconstruct only the 3rd PC
heatmap(s$u%*% s$d %*%t(s$v[,3]), main="3rd PC Heatmap")
heatmap(cleandat)

corHB <- cor(matchH,matchB)
heatmap(corHB, main="Correlations between Blastema and Embryo")

for(i in 1:12){
  ifile <- paste0(i,'_heatmap.pdf')
  pdf(ifile)
  heatmap(cor(matchB[upregs[[i]]$symbol,], matchH[upRegs[[i]]$symbol,]), main=paste("Blastema Point", i, "Genes"))  
  d <- dev.off()
}

#####################TIME WARPING##################
pou5f1Mat <- matrix(0 ,nrow=nrow(matchH), ncol= ncol(matchH))
for(i in 1:nrow(pou5f1Mat)){
  pou5f1Mat[i,] <-  matchH["POU5F1",]
}

ages.H <- c(1,2,3,4,5,6,7,8,9,10,11,12,14,16,19,24,40)
ages.B <- c(0,.125,.25,.5,1,3,5,7,10,14,21,28)
n.age.H <-c(1,2,3,4,5,6,7,8,9,10,11,12,14,16,19,24,40)
n.age.B <- c(1,2,3,4,5,6,7,8,9,10,11,12,14,16,19,24,40)
#n.age.H <-approx(c(ages.H, ages.B), c(matchH[9466,],matchB[9466,]))$x
#n.age.B <- approx(c(ages.H, ages.B), c(matchH[9466,],matchB[9466,]))$x
HvBshift <-  TimeShift::shift(matchH,pou5f1Mat, ages.H, ages.H, n.age.H, n.age.H )
H_Bshift <- dtw(matchH[,6:17],matchB, window.type = "itakura")


corMin <- .5
shiftDistHi <- 2
shiftDistLo <- -2
devRange <- 3:14
rm <- rowMeans(HvBshift)
passFilter <- sapply(1:length(rm), function(i){
  rm[i] < shiftDistHi && rm[i]> shiftDistLo && rm[i]!=0 && cor(matchH[i,devRange],matchB[i,])>corMin
})
sum(passFilter, na.rm = TRUE)
index<- rownames(matchH)[which(unlist(passFilter))]

index
mypar(5,4)
for(i in 1:length(index)){
  x <- index[i]
  plot(matchH[x,], main=paste("Embryo", x),xlab="TimePoint",ylab="TPM")
  plot(matchH["POU5F1",], main= paste("POU5F1"),xlab="TimePoint",ylab="TPM")
}

shiftTerms <- getGOHuman(rownames(matchH)[index])
printTopGO(shiftTerms$name_1006, 42)

#Rank by the shift distance, calculate proportion of de genes that are expressed similarly in the blastema and emb
#Changing the time shift window also finds delayed genes etc... WHAT DOES IT MEAN??!

shiftTerms <- getGOHuman(blastemaGL[[timePoint]][[1]])
printTopGO(shiftTerms$name_1006, 40)

blastemaGL <- sapply(1:length(upRegs), function(i){
  warpWindow(shiftLimHi=3,shiftLimLo=1,corMin =.5, devRangeH = 1:12, devRangeB=1:12 ,geneSet = upRegs[[i]]$symbol)
} )


timePoint <- 12
length(blastemaGL[[timePoint]][[1]])
length(blastemaGL[[timePoint]][[1]])/nrow(upRegs[[timePoint]])
showPlot(timePoint,blastemaGL)
terms <-  getGOTerms(timepoint, blastemaGL)

mypar(5,4)
index <- upRegs[[7]]$symbol
index <- Limbs
for(i in 1:length(index)){
  x <- index[i]
  plot(matchH[x,], main=paste("Embryo", x),xlab="TimePoint",ylab="TPM")
  plot(matchB[x,], main= paste("Blastema",x),xlab="TimePoint",ylab="TPM")
}

matchH["SHH",]
mypar(1,2)
g="SHH"
HvBshift[10431,]
dtw(matchH[g,],matchB[g,], window.type = "itakura")
which(rownames(matchH)==g)
plot(matchH[g,])
plot(matchB[g,])
cor(matchH[g,6:17],matchB[g,])