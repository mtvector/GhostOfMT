source("http://bioconductor.org/biocLite.R")
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(EACI)
library(matrixStats)
library(dtw)
library(EBSeq)
# library(genefilter)
options("mc.cores"=2)
library(parallel)
library(rafalib)
library(calibrate)
source("~/code/cb/LoadFunctions.R")
#source("~/code/cb/LoadData.R")
#remove.packages("SegReg")
#install.packages("~/code/general/SegRegPar/package/SegReg",repos = NULL,type = "source")
library(SegReg)
library(MASS)
library(ggrepel)
library(zoo)

#Function takes a vector of booleans (value > threshold) and an index (representing the max). 
#Walks out from this index and returns contiguous indices > threshold
stepOut <- function(v,i,skip.num=1){
  inds <- 1:length(v)
  iu <- i+1
  id <- i-1
  u <- 0
  d <- 0
  ret.inds <- c()
  if(v[i]){
    ret.inds <- c(ret.inds,i)
    while(u<=skip.num){
      if(iu %in% inds){
        if(v[iu]){
          ret.inds <- c(ret.inds,iu)
          iu <- iu+1
        }else{
          u <- u+1
        }
      }else{
        u <- u+1
      }
    }
    while(d<=skip.num){
      if(id %in% inds){
        if(v[id]){
          ret.inds <- c(ret.inds,id)
          id <- id-1
        }else{
          d <- d+1
        }
      }else{
        d <- d+1
      }
    }
  }
  ret.inds
}

#Finds the expression area around the max expression for each gene
maxArea <- function(ds,tpds,thresh=1,skip.num=1){
  w.m <- apply(ds,1,which.max)
  #vals <- rowMaxs(ds)
  vals <-  sapply(1:nrow(ds),function(i){
    #Get the sum of values contiguous with max
    sum(ds[i,stepOut(ds[i,]>thresh,which.max(ds[i,]),skip.num = skip.num)])
  })
  vals
}

peakScore <- function(ds,tpds,thresh=seq(.5,5,by = .5)){
  ds <- round.log(ds[,order(tpds)])
  weightedScores <- lapply(1:length(thresh),function(t){
    bscore <- rep(0,nrow(ds))
    escore <- rep(0,nrow(ds))
    utpds <- sort(unique(tpds))
    #count the number of zeroes in the beginning (bscore)
    for(i in 1:(length(utpds)-1)){
      if(sum(tpds== utpds[i])>1){
        binds <-  rowMeans(ds[,tpds == utpds[i]]) < thresh[t]
      }else{
        binds <-  ds[,tpds == utpds[i]] < thresh[t]
      }
      bscore <-  bscore + as.integer(bscore >= i-1 & binds)
    }
    #count the number of zeroes in the end (escore)
    for(i in length(utpds):2){
      if(sum(tpds==utpds[i])>1){
        einds  <-  rowMeans(ds[,tpds == utpds[i]]) < thresh[t]
      }else{
        einds <-  ds[,tpds == utpds[i]] < thresh[t]
      }
      escore <-  escore +  as.integer(escore >= abs(i-length(utpds)) & einds)
    }
    res <- cbind(bescore=(bscore*escore),peakArea=(maxArea(ds,tpds,thresh[t])))
    rownames(res) <- rownames(ds)
    res
  })
  condenseScores <- Reduce(cbind,lapply(weightedScores,function(x){x[,1]*x[,2]}))
  normMat <- GetNormalizedMat(condenseScores , colMaxs(condenseScores))
  #normMat <- GetNormalizedMat(condenseScores , seq(minweight,1,length.out = ncol(condenseScores)))
  mxs <- rowMaxs(normMat)
  names(mxs) <- rownames(ds)
  mxs
}

sort(peakScore(datasets[[1]],tpDatasets[[1]]),T)
plot(datasets[[1]]["LHX9",])


h.inds <- which(grepl("Human",datasetNames))
m.inds <- which(grepl("Mouse",datasetNames))

m.inds <- c(3,5,8,10,14)
h.inds <- c(1,7,11,20)

sort(peakScore(datasets[[1]],tpDatasets[[1]]),T)
plot(datasets[[1]]["LHX9",])

scoresHuman <- lapply(h.inds,function(i){
  print(i)
  peakScore(datasets[[i]],tpDatasets[[i]])
})

scoresMouse <- lapply(m.inds,function(i){
  print(i)
  peakScore(datasets[[i]],tpDatasets[[i]])
})

peakGenesH <-  lapply(scoresHuman,function(x){names(sort(x,T)[1:200])})
peakGenesM <-  lapply(scoresMouse,function(x){names(sort(x,T)[1:200])})

pdf(file = "~/Desktop/HumanPeaks.pdf")
mypar(2,2)
for(g in names(sort(table(Reduce(c,peakGenesH)),T))[1:200]){
  ds <- datasets[h.inds]
  ds2 <- datasets[m.inds]
  tpds <- tpDatasets[h.inds]
  tpds2 <- tpDatasets[m.inds]
  colz <- Reduce(c, lapply(1:length(tpds),function(i){if(g %in% rownames(ds[[i]]))i}))
  plot(Reduce(c, lapply(1:length(tpds),function(i){if(g %in% rownames(ds[[i]]))tpds[[i]]})),Reduce(c, lapply(ds,function(s){if(g %in% rownames(s))s[g,]})),main = g,ylab = "TPM", xlab = "Day",col=colz)
  if(g %in% Reduce(union,lapply(ds2,rownames))) plot(Reduce(c, lapply(1:length(tpds2),function(i){if(g %in% rownames(ds2[[i]]))tpds2[[i]]})),Reduce(c, lapply(ds2,function(s){if(g %in% rownames(s))s[g,]})),main = paste("Mouse",g),ylab = "TPM", xlab = "Day",col=colz)
}
dev.off()

pdf(file = "~/Desktop/MousePeaks.pdf")
mypar(2,2)
for(g in names(sort(table(Reduce(c,peakGenesM)),T))[1:200]){
  ds <- datasets[m.inds]
  ds2 <- datasets[h.inds]
  tpds <- tpDatasets[m.inds]
  tpds2 <- tpDatasets[h.inds]
  colz <- Reduce(c, lapply(1:length(tpds),function(i){if(g %in% rownames(ds[[i]]))i}))
  plot(Reduce(c, lapply(1:length(tpds),function(i){if(g %in% rownames(ds[[i]]))tpds[[i]]})),Reduce(c, lapply(ds,function(s){if(g %in% rownames(s))s[g,]})),main = g,xlim = c(0,45),ylab = "TPM", xlab = "Day",col=colz)
  if(g %in% Reduce(union,lapply(ds2,rownames))) plot(Reduce(c, lapply(1:length(tpds2),function(i){if(g %in% rownames(ds2[[i]]))tpds2[[i]]})),Reduce(c, lapply(ds2,function(s){if(g %in% rownames(s))s[g,]})),main = paste("Human",g),ylab = "TPM", xlab = "Day",col=colz)
}
dev.off()

pdf(file = "~/Desktop/BothPeaks.pdf")
mypar(2,2)
for(g in intersect(names(sort(table(Reduce(c,peakGenesH)),T))[1:400],names(sort(table(Reduce(c,peakGenesM)),T))[1:400])){
  ds <- datasets[m.inds]
  ds2 <- datasets[h.inds]
  tpds <- tpDatasets[m.inds]
  tpds2 <- tpDatasets[h.inds]
  colz <- Reduce(c, lapply(1:length(tpds),function(i){if(g %in% rownames(ds[[i]]))i}))
  plot(Reduce(c, lapply(1:length(tpds),function(i){if(g %in% rownames(ds[[i]]))tpds[[i]]})),Reduce(c, lapply(ds,function(s){if(g %in% rownames(s))s[g,]})),main = g,xlim = c(0,45),ylab = "TPM", xlab = "Day",col=colz)
  if(g %in% Reduce(union,lapply(ds2,rownames))) plot(Reduce(c, lapply(1:length(tpds2),function(i){if(g %in% rownames(ds2[[i]]))tpds2[[i]]})),Reduce(c, lapply(ds2,function(s){if(g %in% rownames(s))s[g,]})),main = paste("Human",g),ylab = "TPM", xlab = "Day",col=colz)
}
dev.off()

#c("DBX1", "NEUROD4", "NEUROG2", "SCUBE3", "ASCL1", "SHH", "WNT7A", "ADAMTS16", "EMILIN2")
pdf(file = "~/Desktop/BothPeaks.pdf")
mypar(2,2)
for(g in intersect(names(sort(table(Reduce(c,peakGenesH)),T)),names(sort(table(Reduce(c,peakGenesM)),T)))){
  ds <- datasets[m.inds]
  ds2 <- datasets[h.inds]
  tpds <- tpDatasets[m.inds]
  tpds2 <- tpDatasets[h.inds]
  colz <- Reduce(c, lapply(1:length(tpds),function(i){if(g %in% rownames(ds[[i]]))i}))
  plot(na.approx(Reduce(merge.zoo,lapply(1:length(tpds),function(i){ aggregate(zoo(ds[[i]][g,],tpds[[i]]),index(zoo(ds[[i]][g,],tpds[[i]])), mean)}))),plot.type = "single",col = 1:10,type = "l",ylab = "TPM", xlab = "Day",main = paste("Mouse",g))
  legend("topright",datasetNames[m.inds] , col = 1:5, lty = 1)
  #plot(Reduce(cbind, lapply(1:length(tpds),function(i){if(g %in% rownames(ds[[i]]))tpds[[i]]})),Reduce(cbind, lapply(ds,function(s){if(g %in% rownames(s))s[g,]})),main = g,ylab = "TPM", xlab = "Day",col=colz,type = 'l')
  plot(na.approx(Reduce(merge.zoo,lapply(1:length(tpds2),function(i){ if(g %in% rownames(ds2[[i]])){aggregate(zoo(ds2[[i]][g,],tpds2[[i]]),index(zoo(ds2[[i]][g,],tpds2[[i]])), mean)}else{zoo(rep(0,length(tpds2[[i]])),tpds2[[i]])}}))),plot.type = "single",col = 1:10,type = "l",ylab = "TPM", xlab = "Day",main = paste("Human",g))
  legend("topright",datasetNames[h.inds] , col = 1:5, lty = 1)
}
dev.off()