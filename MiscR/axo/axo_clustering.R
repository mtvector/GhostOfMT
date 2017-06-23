library("randomForest")
source("http://bioconductor.org/biocLite.R")
biocLite("Mfuzz")
biocLite("EBSeq")
library(EBSeq)
library(Mfuzz)
install.packages("ape")
library(ape)
biocLite("GOstats")
library(GOstats)
biocLite("topGO")
library(topGO)

?topGO
path <-"~/code/data/axo"

blastemaFile <- file.path(path, "tanaka_norrna_rsem_human_gene_expressionForUnixJustJuvTCECs.txt")
ECBAll <- read.table(blastemaFile,header=TRUE)
rownames(ECBAll) <-ECBAll$symbol
ECB<- as.matrix(ECBAll[,9:ncol(ECBAll)])
sizesB <- MedianNorm(ECB[rowMeans(ECB)>5,])
NormECB <- GetNormalizedMat(ECB[rowMeans(ECB)>5,], sizesB)

filenameH = file.path(path,"Human_Annotation/ExpectCounts_all_samples.txt")
#Expected Counts All Samples
ECEAll <- read.table(filenameH,header=TRUE)
#Assign rownames as gene names
rownames(ECEAll) <- ECEAll[,1]
ECE <- as.matrix(ECEAll[-1:-2])
sizesE <- MedianNorm(ECE[rowMeans(ECE)>5,])
NormECE <-  GetNormalizedMat(ECE[rowMeans(ECE)>5,], sizesE)

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
colnames(equalRepsE) <- make.names(stagesH,unique = T)

btab <-  read.table(file = file.path(path,"rmaB.txt"), header = F, sep = "\t")
btabM <- as.data.frame.matrix(btab[-1,c(-1,-2,-3)])
timepoints <- lapply(btab[1,-(1:3)], as.character)
timepoints <- unlist(timepoints)
btabM <- apply(btabM,c(1,2), as.numeric)
names(timepoints) <- NULL
rownames(btabM) <- btab[-1,1]

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

dis

distB <- dist(equalRepsB)
hcB <-  hclust(distB)
plot(hcB,  hang=-1)
distE <- dist(equalRepsE)
hcE <- hclust(distE)
dendE <- as.dendrogram(hcE)
str(dendE, max = 10000, last.str =  "'")
dim(equalRepsB)
plot(dendE,  hang=-1)
plot(as.phylo(hcE), type = "fan")
?cutree
memb <- cutree(hcE, k = 400)
plot(memb)
cent <- NULL
for(k in 1:10){
  cent <- rbind(cent, colMeans(USArrests[memb == k, , drop = FALSE]))
}
hc1 <- hclust(dist(cent)^2, method = "cen", members = table(memb))
opar <- par(mfrow = c(1, 2))
plot(hc1, labels = FALSE, hang = -1, main = "Re-start from 10 clusters")

interEB <- intersect(rownames(equalRepsE),rownames(NormECB))
distEB <- cbind(NormECE[interEB,], NormECB[interEB,])
save.image(file = file.path(path, "clustersEB.Rdata"))

# mtestB <- partcoef(blastemaV.s)
# mtestE <- partcoef(embryo.s)
# 
# F <- mtestE[[1]];F.n <- mtestE[[2]];F.min <- mtestE[[3]]
# F > 1.01 * F.min

blastemaV.f  <-  filter.std(ExpressionSet(equalRepsB), min.std=0.04)
blastemaV.s <- standardise(blastemaV.f)

?cselection
optimizedc.B <- cselection(blastemaV.s,crange=seq(16,60,4),m=1.25, visu = T)
clB.Voss <- mfuzz(blastemaV.s,c=20, m=1.25)
mfuzz.plot(blastemaV.s,cl=clB.Voss,mfrow=c(4,4),time.labels=colnames(equalRepsB))
save.image(file = file.path(path, "clustersEB.Rdata"))
?mfuzz

embryo.f  <-  filter.std(ExpressionSet(log(NormECB+1,2)), min.std=.3)
embryo.s <- standardise(embryo.f)
optimizedc.E <- cselection(embryo.s,crange=seq(16,60,4),m=1.25, visu=T)
cl.Embryo <- mfuzz(embryo.s,c=40, m=1.25)
mfuzz.plot(embryo.s,cl=cl.Embryo,mfrow=c(4,4),time.labels=colnames(NormECB))
overlap.plot(cl.Embryo,overlap(cl.Embryo))
length(table(cl.Embryo$cluster))
cl.Embryo$cluster
library(org.Hs.eg.db)
goLists <- sapply(1:length(table(cl.Embryo$cluster)),function(i){
params <- new("GOHyperGParams",
              geneIds=ECBAll[names(cl.Embryo$cluster[cl.Embryo$cluster==i]),]$geneid,
              universeGeneIds=ECBAll$geneid,
              ontology="BP",
              annotation = "org.Hs.eg.db",
              pvalueCutoff=.0001,
              conditional=F,
              testDirection="over"
)
hgOver <- hyperGTest(params)
})
topGO::
?hyperGTest
vignette("Mfuzz")
