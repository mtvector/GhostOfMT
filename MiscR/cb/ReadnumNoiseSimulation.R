library(dtw)
library(SegReg)
options("mc.cores"=2)
source("~/code/cb/LoadFunctions.R")
source("~/code/cb/LoadData.R")
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
library(WGCNA)

str(AllMixSetsEC)
dat <- AllMixSetsEC[[1]]$mm[,1]


tpmreal <- dat/sum(dat)

reads <- sapply(1:length(dat),function(i){
  rep(i, round(dat)[i])
})

genes <- names(dat)

reads <- reads[sapply(reads, length)>0]

nonzerogenes <- names(dat)[sapply(reads, length)>0]
nonzerogeneinds <- 1:length(nonzerogenes)

compReads <- Reduce(c,reads)

n <- length(compReads)
out <- sapply(seq(1000,n,length.out = 50),function(x){
  print(x)
  tally <-  lapply(1:100, function(y){
    s<-table(sample(compReads,x))
    addon <- nonzerogeneinds[!nonzerogeneinds%in% names(s)]
    names(addon) <- nonzerogeneinds[!nonzerogeneinds%in% names(s)]
    c(s,sapply(addon,function(z)0))
  })
  outmat <- Reduce(cbind,tally)
  tpmoutmat <- sweep(outmat,2,colSums(outmat),"/")
  list(rowSds(tpmoutmat), rowMeans(tpmoutmat), outmat)
})

