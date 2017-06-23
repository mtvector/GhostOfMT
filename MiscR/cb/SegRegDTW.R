library(dtw)
library(SegReg)
options("mc.cores"=2)
source("~/code/cb/LoadFunctions.R")
source("~/code/cb/LoadData.R")
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
library(WGCNA)

bm <- blockwiseModules(mat1[1:2000,])

detectEvents <- function(o,tp){
  mclapply(o, function(x){
    #x <- o[[n]]
    EventPoints <-  c(tp[1],x$bp, tp[length(x$fitted)])
    EventPoints <- na.omit(EventPoints)
    if(length(EventPoints)==1){
      re <- as.matrix(c(start= EventPoints[1],end= max(tp) , slp= x$slp[1],slp.pval= x$slp.pval[1],fc = x$slp[1]*(max(tp) -EventPoints[1])))
    }else{
    Event.Table <- sapply(1:(length(EventPoints)-1), function(i){
      re <- c(start= EventPoints[i],end= EventPoints[i+1] , slp= x$slp[i],slp.pval= x$slp.pval[i],fc = x$slp[i]*(EventPoints[i+1] -EventPoints[i]))
    } )
    t(Event.Table)
  }})
}

neural_list <- toupper(sort(as.character(read.csv("~/code/data/cb/markers/FullNeuralGeneList.csv",header = T)[,1])))
neural_related_genes <- make.names(toupper(as.character(read.csv("~/code/data/cb/markers/neuralrelatedgenes.csv",header = T)[,1])))
pluripotent_genes <- toupper(sort(as.character(read.csv("~/code/data/general/marker_lists/ES_Germ/ES.txt",header = T)[,1])))

pc1 <- princomp(mat1)
pc1$scores
pc1$call

#tpm.threshold =3 , gene.list=neural_list ,pvalcut = .1,maxk = 4,min.num.in.seg = 3,cutdiff = .01
set1 = list(datasets[[11]])
tp1 = list(tpDatasets[[11]])
set2 = list(datasets[[12]])
tp2 = list(tpDatasets[[12]])

tpm.threshold =3 
gene.list=neural_list[1:20] 
pvalcut = .1
maxk = 5
min.num.in.seg = 3
cutdiff = .01

all.rn <- as.character(Reduce(intersect, lapply(c(set1,set2),rownames)))
all.rn <- all.rn[all.rn%in%gene.list]
set1 <- lapply(set1, function(n) n[all.rn,])
s1.mx <-  lapply(set1, function(x) rowMaxs(x) )
pass.thresh <-  Reduce( "+" ,lapply(s1.mx, ">",tpm.threshold)) >= length(set1)-1
set1 <- lapply(set1, function(n) n[pass.thresh,])
mat1 <- Reduce(cbind, set1)
mat1 <- round.log(mat1,2)
tp1cat <- Reduce(c , tp1)
tp1.i <- seq(min(tp1cat),max(tp1cat))


set2 <- lapply(set2, function(n) n[all.rn,])
s2.mx <-  lapply(set2, function(x) rowMaxs(x) )
pass.thresh <-  Reduce( "+" ,lapply(s2.mx, ">",tpm.threshold)) >= length(set2)-1
set2 <- lapply(set2, function(n) n[pass.thresh,])
mat2 <- Reduce(cbind, set2)
mat2 <- round.log(mat2,2)
tp2cat <- Reduce(c , tp2)
tp2.i <- seq(min(tp2cat),max(tp2cat))

if(!file.exists("~/code/data/temp/segregtest.rds")){
  sr1 <- segreg(mat1,t.vect = tp1cat,meancut = tpm.threshold,maxk = maxk,min.num.in.seg = min.num.in.seg,pvalcut = pvalcut,cutdiff = cutdiff)
  sr2 <- segreg(mat2,t.vect = tp2cat,meancut = tpm.threshold,maxk = maxk,min.num.in.seg = min.num.in.seg,pvalcut = pvalcut,cutdiff = cutdiff)
}
save(sr1,sr2,file="~/code/data/temp/segregtest.rds")
hist(sapply(sr1,"[[","radj"))



fit.mat.1 <- t(sapply(sr1,function(x) {
  a <- x[["fitted"]]
  r <- approx(tp1cat,a,xout = tp1.i)$y
  names(r) <- tp1.i
  r
  }))

fit.mat.2 <- t(sapply(sr2,function(x) {
  a <- x[["fitted"]]
  r <- approx(tp2cat,a,xout = tp2.i)$y
  names(r) <- tp2.i
  r
}))

gene.list <- intersect(rownames(fit.mat.1),rownames(fit.mat.2))

plotmarker(mat1,tp1cat,fittedres = sr1,listname = c("IFRD1"))
plotmarker(mat2,tp2cat,fittedres = sr2,listname = c("ISL1"))

w.xcors <- sapply(gene.list, function(x){
  a=ccf(fit.mat.1[x,],fit.mat.2[x,],plot=F)$acf
  which.max(a)-(as.integer(length(a)/2)+1)
})
names(w.xcors) <- gene.list

dt.fc <- sapply(gene.list,function(i){
  a=fit.mat.1[i,]
  b=fit.mat.2[i,]
  min((max(a)-min(a)), (max(b)-min(b)))
})

lapply(sr1,function(x)x['radj'])

deEvents1 <- detectEvents(sr1, sort(unique(tp1cat)))
deEvents2 <- detectEvents(sr2, sort(unique(tp2cat)))

gene.names <- unlist(lapply(names(deEvents1), function(n){  rep(n, nrow(deEvents1[[n]]))} ))
de1 <- Reduce(rbind, deEvents1)
de1 <- data.frame(de1,stringsAsFactors = F)
de1 <- cbind(gene.names, de1)
colnames(de1) <- c("gene","start","end", "slp", "slp.pval","fc")

gene.names <- unlist(lapply(names(deEvents2), function(n){ rep(n, nrow(deEvents2[[n]]))} ))
de2 <- Reduce(rbind, deEvents2)
de2 <- data.frame(de2,stringsAsFactors = F)
de2 <- cbind(gene.names, de2)
colnames(de2) <- c("gene","start","end", "slp", "slp.pval","fc")

#write.table(de,file = file.path(output.dir,paste0(fn,".txt")),sep = "\t")




d <- dist(mat1["POU3F2",],mat2["POU3F2",])
d[dim(d)[1],] <- 0
d[,dim(d)[2]] <- 0
dt <- dtw(d,step.pattern = rabinerJuangStepPattern(1,"d"))
dtwPlotTwoWay(dt,mat1["POU3F2",],mat2["POU3F2",])
dt$index1
