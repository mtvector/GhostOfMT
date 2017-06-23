source("http://bioconductor.org/biocLite.R")
biocLite("EBSeq")
library(EBSeq)
library(matrixStats)
install.packages("randomForest", repos="http://cran.rstudio.com/")
library(randomForest)
biocLite("networkBMA")
library(networkBMA)

path <-"~/code/data/axo"


blastemaFile <- file.path(path, "tanaka_norrna_rsem_human_gene_expressionForUnixJustJuvTCECs.txt")
ECBAll <- read.table(blastemaFile,header=TRUE)
rownames(ECBAll) <-ECBAll$symbol
ECB<- as.matrix(ECBAll[,9:ncol(ECBAll)])
sizesB <- MedianNorm(ECB)
NormECB <- GetNormalizedMat(ECB, sizesB)

filenameH = file.path(path,"Human_Annotation/ExpectCounts_all_samples.txt")
#Expected Counts All Samples
ECEAll <- read.table(filenameH,header=TRUE)
#Assign rownames as gene names
rownames(ECEAll) <- ECEAll[,1]
ECE <- as.matrix(ECEAll[-1:-2])
sizesE <- MedianNorm(ECE)
NormECE <-  GetNormalizedMat(ECE, sizesE)

factor(colnames(NormECE))
stagesH <- c("Stage_1","Stage_1","Stage_2","Stage_2","Stage_2", "Stage_3","Stage_3","Stage_3","Stage_4","Stage_4","Stage_4","Stage_5","Stage_5","Stage_5","Stage_6",	"Stage_6",	"Stage_6",	"Stage_7",	"Stage_7",	"Stage_7","Stage_8",	"Stage_8",	"Stage_8","Stage_9",	"Stage_9",	"Stage_9",	"Stage_10",	"Stage_10",	"Stage_10",	"Stage_11",	"Stage_11",	"Stage_12","Stage_12","Stage_12","Stage_14","Stage_14","Stage_14","Stage_16","Stage_16","Stage_16","Stage_19","Stage_19","Stage_19","Stage_24","Stage_24","Stage_24","Stage_40","Stage_40", "Stage_40")


#number of reps
n <- 3
equalRepsE <- matrix(nrow=nrow(NormECE),ncol=0)
for(i in 1:length(stagesH)){
  stage <- stagesH[i]
  equalRepsE <- cbind(equalRepsE,NormECE[,i])
  if(table(stagesH)[stage]<n){
    less <- n - table(stagesH)[stage]
    stagesH[i:(length(stagesH)+1)] <- c(stagesH[i], stagesH[i:length(stagesH)])
    equalRepsE <- cbind(equalRepsE, matrix(rowMeans(NormECE[,(i):(i+n-less-1)]),nrow = nrow(equalRepsE),ncol=(n-less-1)))
  }
}
colnames(equalRepsE) <- stagesH

TFs <- read.table(paste0(path,"/my_human_TF_v5_final_11282011.txt"),header=T, sep="\t")

stageNumH <- c("1","2","3","4","5","6","7","8","9","10","11","12","14","16","19","24","40")
DEPath <- file.path(path,"DE")
DEUpE <- sapply(1:(length(stageNumH)-1), function(i){
  read.table( file.path(DEPath, "Up", paste0("Stage_",stageNumH[i+1] , "_Vs_",stageNumH[i],".Up.txt")), header = TRUE, sep = "\t")$Gene.ID
})
DEDownE <- sapply(1:(length(stageNumH)-1), function(i){
  read.table( file.path(DEPath, "Down", paste0("Stage_", stageNumH[i+1] , "_Vs_",stageNumH[i],".Down.txt")), header = TRUE, sep = "\t")$Gene.ID
})

DEGE <- c()
for(i in seq(1,8)){
  DEGE <- c(DEGE, as.character( DEUpE[[i]])) 
  DEGE <- c(DEGE, as.character(DEDownE[[i]]))
}
DEGE1to8 <- unique(DEGE)

DEGE <- c()
for(i in seq(9,16)){
  DEGE <- c(DEGE, as.character( DEUpE[[i]])) 
  DEGE <- c(DEGE, as.character(DEDownE[[i]]))
}
DEGE9to40 <- unique(DEGE)

DEGE <- c()
for(i in seq(1,16)){
  DEGE <- c(DEGE, as.character( DEUpE[[i]])) 
  DEGE <- c(DEGE, as.character(DEDownE[[i]]))
}
DEGE1to40 <- unique(DEGE)

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
# 
# DEGB <- c()
# for(i in seq(1,11)){
#   DEGB <- c(DEGB, as.character( EBOutResultsB[,i]$DEfound)) 
# }
# DEGB <- unique(DEGB)
# write.csv(DEGB, file=file.path(path,"DEGB.csv"))
DEGB <-  read.csv(file.path(path, "DEGB.csv"))[,2]
DEGB <- as.character(DEGB)

btab <-  read.table(file = file.path(path,"rmaB.txt"), header = F, sep = "\t")
btabM <- as.data.frame.matrix(btab[-1,c(-1,-2,-3)])
timepoints <- lapply(btab[1,-(1:3)], as.character)
timepoints <- unlist(timepoints)
btabM <- apply(btabM,c(1,2), as.numeric)
names(timepoints) <- NULL
rownames(btabM) <- btab[-1,1]
head(btab)

n <- 10
tpNames <- c("0", ".5", "1","1.5","2","3","4","5","7","9","10","12","14","16","18","20","22","24","26","28")
equalRepsB <- matrix(nrow=nrow(btabM),ncol=0)
for(i in 1:ncol(btabM)){
  stage <- timepoints[i]
  equalRepsB <- cbind(equalRepsB,btabM[,i])
  if(table(timepoints)[stage]<n){
    less <- n - table(timepoints)[stage]
    timepoints[i:(length(timepoints)+1)] <- c(timepoints[i], timepoints[i:length(timepoints)])
    equalRepsB <- cbind(equalRepsB, matrix(rowMeans(btabM[,(i):(i+n-less-2)]),nrow = nrow(btabM),ncol=less))
  }
}
colnames(equalRepsB) <-  make.names(timepoints, unique = T)


#nbmaBMicroarray <- networkBMA(t(equalRepsB), 20, prior.prob =NULL,verbose=T)
#save.image(file="regnetsBVoss")
#writeEdges( nbmaBMicroarray, threshhold = .7,fileName = file.path(path,"nbmaBVoss.txt"))

#nbmaB <- networkBMA(t(NormECB[DEGB,]), 12, prior.prob = .5,verbose=T)
#writeEdges( nbmaB, threshhold = .7,fileName = file.path(path,"nbmaB.txt"))
#save.image(file="regnetsEB") 



nbmaE1to8 <-  networkBMA(t(log(equalRepsE[ DEGE1to8, c(1:8, 18:25,35:42)]+2,2)),8, priorprob=.3 , verbose=T)
save.image(file=file.path(path,"regnetsEB.Rdata"))
writeEdges( nbmaBMicroarray, threshhold = .7,fileName = file.path(path,"nbmaE1to8.txt"))


nbmaE9to17 <-  networkBMA(t(log(equalRepsE[DEGE9to40,c(9:17,26:34,43:51)]+2,2)),9,priorprob=.3, verbose=T)
writeEdges( nmba9to17, threshhold = .7,fileName = file.path(path,"nbmaE9to17.txt"))
save.image(file=file.path(path,"regnetsEB.Rdata"))


#nbma1to40 <-  networkBMA(t(log(equalRepsE[DEGE1to40,]+2,2)),17, verbose = T)
#writeEdges( nmba9to17, threshhold = .7,fileName = file.path(path,"nbmaE1to40.txt"))
#save.image(file="nbmaE1to40")


