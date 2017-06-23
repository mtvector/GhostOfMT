library(ggplot2)
library(gridExtra)
library(ggrepel)
datapath <- "~/code/data/cb/NBNogData/MixingTest"
setwd(datapath)
hMat <- as.matrix(read.csv2(file.path(datapath,"All_Sub527_hg19_ReplacedWithSub537samples.txt"),sep = "\t",header = T,row.names = 1))
hMat <- hMat[,-c(ncol(hMat))]
hMat <- matrix(as.numeric(hMat) , nrow = nrow(hMat), ncol = ncol(hMat), dimnames =dimnames(hMat))
rownames(hMat) <- toupper(rownames(hMat))
mMat <- read.csv2(file.path(datapath,"All_Sub527_mm10.txt"),sep = "\t",header = T)
mMat <- as.matrix(mMat)
rownames(mMat) <- mMat[,1]
mMat <- mMat[!duplicated(mMat[,1]),]
mMat <- mMat[,-c(1,ncol(mMat))]
mMat <- matrix(as.numeric(mMat) , nrow = nrow(mMat), ncol = ncol(mMat), dimnames =dimnames(mMat))
rownames(mMat) <- toupper(rownames(mMat))

fast <- read.csv2("~/Desktop/hTERAvshContr_Fast.txt",sep = "\t")
slow <- read.csv2("~/Desktop/hTERAvshContr_Slow.txt",sep = "\t")

fastneu <- unlist(lapply(fast,  function(x) x%in% neural_related_genes))
slowneu <- unlist(lapply(slow,  function(x) x%in% neural_related_genes))
sum(fastneu)
mean(fastneu)
sum(slowneu)
mean(slowneu)

neural_selections <- read.csv2(file.path(datapath,"CustomGenes.txt"),sep = "\t",header = T)
neural_selections <- toupper(as.character(neural_selections[,1]))

#sum(neural_selections%in% neural_related_genes)/length(neural_selections%in% neural_related_genes)
#neural_selections[!neural_selections%in% neural_related_genes]

colnames(hMat) <- gsub("_H2BmChIP13","",colnames(hMat))
colnames(hMat) <- gsub("Split_","Split",colnames(hMat))
conditH <- sapply(strsplit(colnames(hMat),"_d"),"[[",1)
tpH <- as.numeric(sapply(strsplit(colnames(hMat),"_d"),"[[",2))

colnames(mMat) <- gsub("_H2BmChIP13","",colnames(mMat))
colnames(mMat) <- gsub("Split_","Split",colnames(mMat))
conditM <- sapply(strsplit(colnames(mMat),"_d"),"[[",1)
tpM <- as.numeric(sapply(strsplit(colnames(mMat),"_d"),"[[",2))

spec <- list("Mouse","Human")
condits <- list(conditM,conditH)
tps <- list(tpM,tpH)
mats <- list(mMat,hMat)
for(j in 1:2){
  for( i in which(!duplicated(condits[[j]]))){
    print(paste(spec[[j]],condits[[j]][i]))
    curMat <- mats[[j]][,which(condits[[j]]==condits[[j]][i])]
    curTP <- tps[[j]][which(condits[[j]]==condits[[j]][i])]
    if(spec[[j]]=="Human"){
      curMat <- cbind(mats[[j]][,which(condits[[j]]== "H9")],curMat)
      curTP <- c(tps[[j]][which(condits[[j]]== "H9")],curTP)
    }
    if(spec[[j]]=="Mouse"){
      curMat <- cbind(mats[[j]][,which(condits[[j]]== "X3")],curMat)
      curTP <- c(tps[[j]][which(condits[[j]]== "X3")],curTP)
    }
    plotlist <- lapply( neural_selections,function(gn){
      if(gn %in% rownames(curMat)){
        gg <-  ggplot(as.data.frame(cbind("tme"=curTP, "value"=curMat[gn,])), aes(x=tme, y=value)) + 
          geom_line() +
          geom_point(size=.8) +
          theme_bw(base_size = 6) +
          theme(legend.position="none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank())+
        #scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
          labs(x = "Day", y = "TPM", 
               title = gn,family="arial",size=6)
        gg
        }
      })
      pdf(paste0("~/Desktop/",spec[[j]],condits[[j]][i],".pdf"), onefile = T)
      print(marrangeGrob(grobs=Filter(Negate(is.null), plotlist) , nrow=4, ncol=4))
      dev.off()
      next
  }
}

condition <- c("X3","H9")
for(j in 1:2){
    curMat <- mats[[j]][,which(grepl(condition[[j]],condits[[j]]) | grepl("Mix" ,condits[[j]]))]
    curTP <- tps[[j]][which(grepl(condition[[j]],condits[[j]]) | grepl("Mix" ,condits[[j]]))]
    plotlist <- lapply( neural_selections,function(gn){
      if(gn %in% rownames(curMat)){
        d <- as.data.frame(cbind("tme"=curTP, "value"=curMat[gn,]))
        d <- cbind(d,"labs"=rownames(d))
        gg <-  ggplot(d, aes(x=tme, y=value)) + 
          geom_point(size=.8) +
          theme_bw(base_size = 6) +
          geom_text_repel(data = d,aes(tme,value,label=labs),segment.size = 0.1,size=.9,box.padding = unit(0.07, "lines"))+
          theme(legend.position="none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank())+
        #scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
          labs(x = "Day", y = "TPM", 
               title = gn,family="arial",size=6)
        gg
        }
      })
      pdf(paste0("~/Desktop/",spec[[j]],condition[[j]],".pdf"), onefile = T)
      print(marrangeGrob(grobs=Filter(Negate(is.null), plotlist) , nrow=2, ncol=2))
      dev.off()
}

###SECOND BIG MIXING EXPERIMENT TESTS

neural_list <- toupper(read.csv2("~/code/data/general/marker_lists/THEneuralList_2017.txt",header = F)[,1])

hMat <- as.matrix(read.csv2("~/code/data/cb/NBNogData/MixingTest/Sub_0563_Developmental_Clock_Human_hg19__d2402ef5dc0c968d/genes.no_mt.ec.tab",sep = "\t",header = T,row.names = 1))
hMat <- hMat[,-c(ncol(hMat))]
hMat <- matrix(as.numeric(hMat) , nrow = nrow(hMat), ncol = ncol(hMat), dimnames =dimnames(hMat))
rownames(hMat) <- toupper(rownames(hMat))
colnames(hMat) <- gsub("X","",colnames(hMat))

mMat <- read.csv2("~/code/data/cb/NBNogData/MixingTest/Sub_0563_Developmental_Clock_Human_mm10__4d7a29b3ee4ac699/genes.no_mt.ec.tab",sep = "\t",header = T)
mMat <- as.matrix(mMat)
rownames(mMat) <- mMat[,1]
colnames(mMat) <- gsub("X","",colnames(mMat))
mMat <- mMat[!duplicated(mMat[,1]),]
mMat <- mMat[,-c(1,ncol(mMat))]
mMat <- matrix(as.numeric(mMat) , nrow = nrow(mMat), ncol = ncol(mMat), dimnames =dimnames(mMat))
rownames(mMat) <- toupper(rownames(mMat))

q <- rn.merge.normalize(hMat,mMat)
hMat <- q[[1]]
mMat <- q[[2]]

tpM <- as.numeric(sapply(strsplit(colnames(mMat),"_d"),"[[",2))
tpH <- as.numeric(sapply(strsplit(colnames(hMat),"_d"),"[[",2))
conditM <- sapply(strsplit(colnames(mMat),"_d"),"[[",1)
conditH <- sapply(strsplit(colnames(hMat),"_d"),"[[",1)

spec <- list("Mouse","Human")
condits <- list(conditM,conditH)
tps <- list(tpM,tpH)
mats <- list(mMat,hMat)

for(j in 1:2){
  curMat <- mats[[j]]
  curTP <- tps[[j]]
  plotlist <- lapply( neural_list,function(gn){
    if(gn %in% rownames(curMat)){
      d <- as.data.frame(cbind("tme"=curTP, "value"=curMat[gn,]))
      d <- cbind(d,"labs"=rownames(d))
      gg <-  ggplot(d, aes(x=tme, y=value)) + 
        geom_point(size=.8) +
        theme_bw(base_size = 6) +
        geom_text_repel(data = d,aes(tme,value,label=labs),segment.size = 0.1,size=.9,box.padding = unit(0.07, "lines"))+
        theme(legend.position="none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank())+
        #scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
        labs(x = "Day", y = "Normalized EC", 
             title = gn,family="arial",size=6)
      gg
    }
  })
  pdf(paste0("~/Desktop/",spec[[j]],".pdf"), onefile = T)
  print(marrangeGrob(grobs=Filter(Negate(is.null), plotlist) , nrow=2, ncol=2))
  dev.off()
}
std.heatmap(cor.compare(hMat,mMat,method="spe"))
