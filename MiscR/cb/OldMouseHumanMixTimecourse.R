source("http://bioconductor.org/biocLite.R")
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(plyr)
#biocLite("EACI")
library(EACI)
library(clusterProfiler)
library(matrixStats)
library(zoo)
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
source("~/code/cb/LoadData.R")
biocLite("maSigPro")
library(maSigPro)
library(edgeR)
library(EBSeq)
library(DESeq2)
library(devtools)
library(biomaRt)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(scales)
#devtools::install_github("greenelab/TDM")
library(TDM)

ec <-  read.table("~/code/data/cb/1_ExpectCounts_all_samples.txt",header = T,row.names = 1)[,-1]
indexes <- t(simplify2array(sapply(rownames(ec), function(x){
  strsplit(x,"|",fixed = T)
})))
indexes[,1] <- toupper(indexes[,1])

Hind <- indexes[,2] == "hg19"
Mind <- indexes[,2] == "mm10"
gene.inds <-  lapply(intersect(indexes[Hind,1], indexes[Mind,1]), function(x){
  c( which( indexes[,1]==x & Hind == T), which( indexes[,1]==x & Hind == F) ) 
})
names(gene.inds) <- intersect(indexes[Hind,1], indexes[Mind,1])
HoverM <-  t(sapply(gene.inds, function(gn){
  (ec[gn[[1]],])/((ec[gn[[1]],])+(ec[gn[[2]],]))
}))
rownames(HoverM) <- names(gene.inds)
HoverM <- as.matrix(HoverM)
HoverM[is.na(HoverM)] <- 0
HoverM <- HoverM[,c(1,2,4,5,6,7,3)]

write.pdf({
  x <- HoverM
  nsamp <- dim(x)[2]
  h <- hist(unlist(x[,1]), plot=FALSE,xlim = c(0,1))
  plot(h$mids, h$counts, type="l", col=rainbow(nsamp)[1], main="% of Reads To Human \n Each TP Each Gene",
       xlab="% of reads assigned to Human", ylab="# observations")
  for(i in 2:nsamp){
    h <- hist(unlist(x[,i]), plot=FALSE)
    lines(h$mids, h$counts, col=rainbow(nsamp)[i])
  }
  legend( "topleft",colnames(x),col =rainbow(nsamp)[1:nsamp],lty=1)
},filename = "~/Desktop/Old0to100ControlsHist.pdf")

combinedH <- read.table(file.path(datapath,"CombinedLibraries_H.txt"),row.names = 1,header = T)
combinedM <- read.table(file.path(datapath,"CombinedLibraries_M.txt"),row.names = 1,header = T)
head(combinedM)
indexes <- t(simplify2array(sapply(rownames(combinedM), function(x){
  strsplit(x,"|",fixed = T)
})))
indexes[,1] <- toupper(indexes[,1])
table(indexes[,2])

Hind <- indexes[,2] == "hg19"
Mind <- indexes[,2] == "mm10"
gene.inds <-  lapply(intersect(indexes[Hind,1], indexes[Mind,1]), function(x){
  c( which( indexes[,1]==x & Hind == T), which( indexes[,1]==x & Hind == F) ) 
})
names(gene.inds) <- intersect(indexes[Hind,1], indexes[Mind,1])

HoverM <-  t(sapply(gene.inds, function(gn){
  sum(combinedH[gn[[2]],])/(sum(combinedH[gn[[1]],])+sum(combinedH[gn[[2]],]))
}))
names(HoverM) <- names(gene.inds)
HoverM[is.na(HoverM)] <- 0
HoverM <- HoverM[1,]

HoverM <- HoverM[order(HoverM,decreasing = T)]
hist(HoverM,main = "Percentage of TPMs detected as Mouse",xlab = "mm TPM / (mmTPM + hsTPM)",ylab = "# of Genes")
write.table(HoverM,file = file.path(datapath,"PercentageHtoMTPM.txt"),sep = "\t")
write.table(HoverM[names(HoverM)%in% neural_list],file = file.path(datapath,"PercentageHtoMTPMneural.txt"),sep = "\t")
sortedHData <-  lapply(names(HoverM), function(x) rbind( combinedH[gene.inds[[x]],]))
sortedHDataa <-  Reduce(rbind, sortedHData)
write.table(sortedHDataa,file = file.path(datapath,"HvsCombinedLibraries_sorted.txt"),sep = "\t")
sortedHDataneural <-  lapply(names(HoverM)[names(HoverM)%in% neural_list], function(x) rbind( combinedH[gene.inds[[x]],]))
sortedHDataaneural <-  Reduce(rbind, sortedHDataneural)
write.table(sortedHDataaneural,file = file.path(datapath,"HvsCombinedLibraries_sorted_neural.txt"),sep = "\t")
sortedHDataa[rownames(sortedHDataa)%in% neural_list,]

MoverH <-  t(sapply(gene.inds, function(gn){
  sum(combinedM[gn[[1]],])/(sum(combinedM[gn[[1]],])+sum(combinedM[gn[[2]],]))
}))
names(MoverH) <- names(gene.inds)
MoverH[is.na(MoverH)] <- 0
MoverH <- MoverH[1,]

MoverH <- MoverH[order(MoverH,decreasing = T)]
hist(MoverH,main = "Percentage of Mouse TPMs detected as Human",xlab = "hs TPM / (mmTPM + hsTPM)",ylab = "# of Genes")
write.table(MoverH,file = file.path(datapath,"PercentageMtoHTPM.txt"),sep = "\t")
write.table(MoverH[names(MoverH)%in% neural_list],file = file.path(datapath,"PercentageMtoHTPMneural.txt"),sep = "\t")
sortedMData <-  lapply(names(MoverH), function(x) rbind( combinedM[gene.inds[[x]],]))
sortedMDataa <-  Reduce(rbind, sortedMData)
write.table(sortedMDataa,file = file.path(datapath,"MvsCombinedLibraries_sorted.txt"),sep = "\t")
sortedMDataneural <-  lapply(names(MoverH)[names(MoverH)%in% neural_list], function(x) rbind( combinedM[gene.inds[[x]],]))
sortedMDataaneural <-  Reduce(rbind, sortedMDataneural)
write.table(sortedMDataaneural,file = file.path(datapath,"MvsCombinedLibraries_sorted_neural.txt"),sep = "\t")
sortedMDataa[rownames(sortedMDataa)%in% neural_list,]

#look at mitochondrial stuff

datapath <- "~/code/data/cb"
setwd(datapath)
mitodatafiles <-  file.path(datapath, dir()[ grepl(".tsv", dir())])
mitodatasets <- sapply( mitodatafiles, read.csv, header=T, sep="\t",row.names=1)
head(mitodatasets[[3]])
for(i in 1:length(mitodatasets)){
  rownames(mitodatasets[[i]]) <-  make.names(toupper(rownames(mitodatasets[[i]])),unique = T)
  mitodatasets[[i]] <- log(as.matrix(mitodatasets[[i]][sort(rownames(mitodatasets[[i]])),])+1,2)
  #class(mitodatasets[[i]]) <- "numeric"
}
matplot(t(mitodatasets$`~/code/data/cb/human_TERA_mito.tsv`),type="l")
matplot(t(mitodatasets$`~/code/data/cb/human_TERA_mito_renorm.tsv`),type="l")
matplot(t(mitodatasets$`~/code/data/cb/mouse_TERA_mito.tsv`),type="l")
rownames(mitodatasets$`~/code/data/cb/human_TERA_mito.tsv`)
rownames(mitodatasets$`~/code/data/cb/human_TERA_mito_renorm.tsv`)[which(grepl(pattern = "MT",x = rownames(mitodatasets$`~/code/data/cb/human_TERA_mito_renorm.tsv`)))]

ind=mito.Genes %in% rownames(exp_2)
matplot(t(( exp_2[mito.Genes[ind],]- rowMeans(exp_2[mito.Genes[ind],]))/rowSds(exp_2[mito.Genes[ind],])),type="l")
matplot(t( exp_2[mito.Genes[ind],]),type="l", main = "All Mitochondrial Associated Genes Mix")
ind=mito.Genes %in% rownames(datasets$`~/code/data/cb/human_ALONE.csv`)
matplot(t(( datasets$`~/code/data/cb/human_ALONE.csv`[mito.Genes[ind],]- rowMeans(datasets$`~/code/data/cb/human_ALONE.csv`[mito.Genes[ind],]))/rowSds(datasets$`~/code/data/cb/human_ALONE.csv`[mito.Genes[ind],])),type="l", main = "All Mitochondrial Associated Genes Human Z")
matplot(t( datasets$`~/code/data/cb/human_ALONE.csv`[mito.Genes[ind],]),type="l", main = "All Mitochondrial Associated Genes Human")
ind=mito.Genes %in% rownames(datasets$`~/code/data/cb/mouse_ALONE.csv`)
matplot(t(( datasets$`~/code/data/cb/mouse_ALONE.csv`[mito.Genes[ind],]- rowMeans(datasets$`~/code/data/cb/mouse_ALONE.csv`[mito.Genes[ind],]))/rowSds(datasets$`~/code/data/cb/mouse_ALONE.csv`[mito.Genes[ind],])),type="l", main = "All Mitochondrial Associated Genes Mouse Z")
matplot(t( datasets$`~/code/data/cb/mouse_ALONE.csv`[mito.Genes[ind],]),type="l", main = "All Mitochondrial Associated Genes Mouse")
mito.Genes[ind]

for(i in 1:length(mitodatasets)){
  rownames(mitodatasets[[i]]) <-  make.names(toupper(rownames(mitodatasets[[i]])),unique = T)
  mitodatasets[[i]] <- log(as.matrix(mitodatasets[[i]][sort(rownames(mitodatasets[[i]])),])+1,2)
  #class(mitodatasets[[i]]) <- "numeric"
}


datapath <- "~/code/data/cb"
setwd(file.path(datapath,"mitochondriaGenes"))
mitodatafiles <-  file.path(datapath, "mitochondriaGenes",dir()[ grepl(".txt", dir())])
mitoGenes <- lapply( mitodatafiles, read.csv, header=T, sep="\t",row.names=1)
mitoGenes <- lapply(lapply(mitoGenes,"[[","Approved.Symbol"),as.character)
names(mitoGenes) <- gsub(".txt","" ,gsub("~/code/data/cb/mitochondriaGenes/","",mitodatafiles))
mito.Genes <- unlist(mitoGenes)

