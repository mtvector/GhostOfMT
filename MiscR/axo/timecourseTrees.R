library("randomForest")
source("http://bioconductor.org/biocLite.R")
source("~/code/axo/genie3.R")
library(matrixStats)
library(rafalib)

path <-"~/code/data/axo"

btab <-  read.table(file = file.path(path,"rmaB.txt"), header = F, sep = "\t")
btabM <- as.data.frame.matrix(btab[-1,c(-1,-2,-3)])
timepoints.o <- lapply(btab[1,-(1:3)], as.character)
timepoints.o <- unlist(timepoints.o)
timepoints <- lapply(btab[1,-(1:3)], as.character)
timepoints <- unlist(timepoints)
btabM <- apply(btabM,c(1,2), as.numeric)
names(timepoints) <- NULL
rownames(btabM) <- btab[-1,1]
head(btabM)

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

rmaNormLogAvg <- sapply(1:nrow(btabM), function(i){
  tapply(btabM[i,], na.omit(timepoints.o), mean, simplify=TRUE)
})

rmaNormLogAvg <- t(rmaNormLogAvg)
rownames(rmaNormLogAvg) <- btab$V2[-1]
B.rmana.match=rmaNormLogAvg[interBB,]
plotMA(B.rmana.match)
plotMA(B.rs.match)

qqplot(log(B.rs.match-rowMeans(B.rs.match)+1,2), B.rmana.match)

interBB <- intersect(rownames(NormECB),btab$V2[-1])
rownames(btabM) <- btab$V2[-1]
B.ma.match=btabM[interBB,]
B.rs.match = NormECB[interBB,]
colnames(B.rs.match) <- c("0h","3h", "6h","12h","1d","3d","5d","7d","10d", "14d","21d", "28d")
zB.ma <- (B.rmana.match-rowMedians(B.rmana.match))/rowSds(B.rmana.match)
zB.rs <- (B.rs.match-rowMedians(B.rs.match))/rowSds(B.rs.match)
colnames(zB.ma) <-  make.names(timepoints.o, unique = T)
blastemaCor <-  cor(cbind(zB.ma,zB.rs),method = "spearman")

 zB.rs

mypar(4,5)
for(i in interBB[1:30]){
  plot(B.ma.match[i,])
  plot(B.rs.match[i,])
}

head(btabM[interBB,1:6])
head(B.ma.match[,1:6])

heatmap.2(blastemaCor, main="Blastema 2"
          ,na.rm = TRUE, trace = "none",
          srtRow=0, srtCol=60,
          col=cols) 

sum(rowSds(equalRepsB)>.06)
num.genes=10000
expr.matrix <- t(btabM[rowSds(btabM)>.06,])
gene.names <-  colnames(expr.matrix)
input.gene.names <- colnames(expr.matrix)
mtry=10000
weight.matrix <- matrix(0.0, nrow=num.genes, ncol=num.genes)
rownames(weight.matrix) <- gene.names
colnames(weight.matrix) <- gene.names

# for (target.gene.idx in seq(from=1, to=num.genes)) {
#   cat(paste("Computing gene ", target.gene.idx, "/", num.genes, "\n", sep=""))
#   flush.console()
#   target.gene.name <- gene.names[target.gene.idx]
#   # remove target gene from input genes
#   these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
#   x <- expr.matrix[,these.input.gene.names]
#   y <- expr.matrix[,target.gene.name]
#   rf <- randomForest(x, y, mtry= round(sqrt(40)), ntree=1000, importance=TRUE)
#   im <- importance(rf)[,"IncNodePurity"]
#   im.names <- names(im)
#   weight.matrix[im.names, target.gene.name] <- im
# }

stageForest <- randomForest(x= expr.matrix,y=timepoints, mtry=round(sqrt(10000)), ntree=5000, importance=T)
stageForest

