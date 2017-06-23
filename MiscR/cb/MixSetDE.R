#source("~/code/cb/MixSetDE.R")
source("http://bioconductor.org/biocLite.R")
library(ggplot2)
library(Biobase)
library(RColorBrewer)
library(gplots)
library(plyr)
#biocLite("EACI")
#library(EACI)
library(clusterProfiler)
library(matrixStats)
library(zoo)
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
#biocLite("maSigPro")
#library(maSigPro)
#library(edgeR)
library(EBSeq)
library(DESeq2)
library(devtools)
library(biomaRt)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(scales)
library(gridExtra)
source("~/code/cb/LoadFunctions.R")
source("~/code/cb/LoadData.R")
library(parallel)
options("mc.cores"=64)
library(DESeq2)
library(EBSeq)


condit.combos <- expand.grid(c("mm","hs"),unique(conditMixSets[[2]]),stringsAsFactors = F)
condit.combos <- condit.combos[-c(2,7),]

mclapply(c("hs","mm"),function(spec){
spec.condit.combos <- condit.combos[condit.combos[,1]==spec,]
c.spec.cc <- combn(1:nrow(spec.condit.combos),2)
mclapply(1:ncol(c.spec.cc),function(j){
  j <- c.spec.cc[,j]
  x <- unlist(spec.condit.combos[j[1],])
  y <- unlist(spec.condit.combos[j[2],])
  cdir <- paste0("~/code/data/cb/shiftFiles/DifferentExpressions/","DESEQ",paste0(x[1],x[2],y[1],y[2]))
      dir.create(cdir)
      design.df <- data.frame("time"=as.factor(c(tpMixSets[[x[1]]], tpMixSets[[y[1]]])),
                              "condition"=as.factor(c( conditMixSets[[x[1]]],conditMixSets[[y[1]]])))
      design.df$combo <- factor(paste(as.character(design.df$condition),as.character(design.df$time),sep = "."),levels=paste(as.character(design.df$condition),as.character(design.df$time),sep = "."))
      DESDS <-  DESeqDataSetFromMatrix(countData = round(cbind(MixSetsEC[[x[1]]], MixSetsEC[[y[1]]] )),
                                                        design.df,
                                                        design = ~combo) 
      DESOut <- DESeq(DESDS)
      mclapply(levels(design.df$condition),function(c){
        t.vect <- design.df$combo[design.df$condition==c]
        t.vect <- unique(c( ifelse(spec=="mm","0.0","100.0"),levels(t.vect)))
        sapply(1:(length(t.vect)-1),function(ti){
          print(c)
          print(ti)
          DESOutResults <- results(DESOut,contrast = list(unique(paste0("combo", t.vect[1])), 
                                                          unique(paste0("combo",t.vect[ti+1]))))
          DESOutResults <- DESOutResults[order(DESOutResults$log2FoldChange,decreasing = T),]
          write.table(DESOutResults,file = paste0(cdir,"/",c,t.vect[ti],"_" ,t.vect[ti+1] ))
        })
      })
      mclapply(levels(design.df$condition),function(c){
        t.vect <- design.df$combo[design.df$condition==c]
        t.vect <- unique(c( ifelse(spec=="mm","0.0","100.0"),levels(t.vect)))
        sapply(1:(length(t.vect)-1),function(ti){
          print( t.vect[ti])
          print( t.vect[ti+1])
          DESOutResults <- results(DESOut,contrast = list(unique(paste0("combo", t.vect[ti])), 
                                                          unique(paste0("combo",t.vect[ti+1]))))
          DESOutResults <- DESOutResults[order(DESOutResults$log2FoldChange,decreasing = T),]
          write.table(DESOutResults,file = paste0(cdir,"/",c,t.vect[ti],"_" ,t.vect[ti+1] ))
        })
      })

      apply(combn(levels(design.df$condition),2),2, function(c){
        DESOutResults <- results(DESOut,contrast = list(unique(paste0("combo", design.df$combo[design.df$condition==c[1]])), unique(paste0("combo",design.df$combo[design.df$condition==c[2]]))))
        DESOutResults <- DESOutResults[order(DESOutResults$log2FoldChange,decreasing = T),]
        write.table(DESOutResults,file = paste0(cdir,"/",c[1],"_",c[2] ))
      })
      #DESOutResults$padj <- p.adjust(DESOutResults$pvalue ,method = "BH")
})
})  

spec="hs"
spec.condit.combos <- condit.combos[condit.combos[,1]==spec,]
c.spec.cc <- combn(1:nrow(spec.condit.combos),2)
mclapply(1:ncol(c.spec.cc),function(j){
  j <- c.spec.cc[,j]
  x <- unlist(spec.condit.combos[j[1],])
  y <- unlist(spec.condit.combos[j[2],])
  cdir <- paste0("~/code/data/cb/shiftFiles/DifferentExpressions/","DESEQ",paste0(x[1],x[2],y[1],y[2]))
  dir.create(cdir)
  design.df <- data.frame("time"=as.factor(c(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ], tpMixSets[[y[1]]][conditMixSets[[y[1]]]==y[2] ])),
                          "condition"=as.factor(c( conditMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ],conditMixSets[[y[1]]][conditMixSets[[y[1]]]==y[2] ])))
  design.df$combo <- factor(paste(as.character(design.df$condition),as.character(design.df$time),sep = "."),levels=paste(as.character(design.df$condition),as.character(design.df$time),sep = "."))
  DESDS <-  DESeqDataSetFromMatrix(countData = round(cbind(MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ], MixSetsEC[[y[1]]][,conditMixSets[[y[1]]]==y[2] ] )),
                                   design.df,
                                   design = ~combo) 
  DESOut <- DESeq(DESDS)
  apply(combn(levels(design.df$condition),2),2, function(c){
    DESOutResults <- results(DESOut,contrast = list(unique(paste0("combo", design.df$combo[design.df$condition==c[1]])), unique(paste0("combo",design.df$combo[design.df$condition==c[2]]))))
    DESOutResults <- DESOutResults[order(DESOutResults$log2FoldChange,decreasing = T),]
    write.table(DESOutResults,file = paste0(cdir,"/",c[1],"_",c[2] ))
  })
  
  mclapply(levels(design.df$condition),function(c){
    t.vect <- design.df$combo[design.df$condition==c]
    t.vect <- unique(c( ifelse(spec=="mm","0.0","100.0"),levels(t.vect)))
    sapply(1:(length(t.vect)-1),function(ti){
      print( t.vect[ti])
      print( t.vect[ti+1])
      DESOutResults <- results(DESOut,contrast = list(unique(paste0("combo", t.vect[ti])), 
                                                      unique(paste0("combo",t.vect[ti+1]))))
      DESOutResults <- DESOutResults[order(DESOutResults$log2FoldChange,decreasing = T),]
      write.table(DESOutResults,file = paste0(cdir,"/",c,t.vect[ti],"_" ,t.vect[ti+1] ))
    })
  })
  mclapply(levels(design.df$condition),function(c){
    t.vect <- design.df$combo[design.df$condition==c]
    t.vect <- unique(c( ifelse(spec=="mm","0.0","100.0"),levels(t.vect)))
    sapply(1:(length(t.vect)-1),function(ti){
      print(c)
      print(ti)
      DESOutResults <- results(DESOut,contrast = list(unique(paste0("combo", t.vect[1])), 
                                                      unique(paste0("combo",t.vect[ti+1]))))
      DESOutResults <- DESOutResults[order(DESOutResults$log2FoldChange,decreasing = T),]
      write.table(DESOutResults,file = paste0(cdir,"/",c,t.vect[ti],"_" ,t.vect[ti+1] ))
    })
  })
  #DESOutResults$padj <- p.adjust(DESOutResults$pvalue ,method = "BH")
})

spec="hs"
spec.condit.combos <- condit.combos[condit.combos[,1]==spec,]
c.spec.cc <- combn(1:nrow(spec.condit.combos),2)
DE.sets <-  lapply(1:(ncol(c.spec.cc)-1),function(j){
  j <- c.spec.cc[,j]
  print(j)
  x <- unlist(spec.condit.combos[j[1],])
  y <- unlist(spec.condit.combos[j[2],])
  cdir <- paste0("~/code/data/cb/shiftFiles/DifferentExpressions/","DESEQ",paste0(x[1],x[2],y[1],y[2]))
  setwd(cdir)
  design.df <- data.frame("time"=as.factor(c(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ], tpMixSets[[y[1]]][conditMixSets[[y[1]]]==y[2] ])),
                          "condition"=as.factor(c( conditMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ],conditMixSets[[y[1]]][conditMixSets[[y[1]]]==y[2] ])))
  design.df$combo <- factor(paste(as.character(design.df$condition),as.character(design.df$time),sep = "."),levels=paste(as.character(design.df$condition),as.character(design.df$time),sep = "."))
  q <- lapply(levels(design.df$condition),function(c){
    print(c)
    t.vect <- design.df$combo[design.df$condition==c]
    t.vect <- unique(c( ifelse(spec=="mm","0.0","100.0"),levels(t.vect)))
    r <-sapply(1:(length(t.vect)-1),function(ti){
      print(t.vect[1])
      print(t.vect[ti+1])
      print( paste0(cdir,"/",c,t.vect[1],"_" ,t.vect[ti+1] ))
      read.table(file = paste0(cdir,"/",c,t.vect[1],"_" ,t.vect[ti+1] ),header = T,row.names = 1)
    })
    colnames(r) <- sapply(1:(length(t.vect)-1),function(ti){paste0(t.vect[1],"_" ,t.vect[ti+1] )})
    r
  })
  names(q) <- levels(design.df$condition)
  q
})
names(DE.sets) <- c("mm0vs10","mm0vs85")





DE.sets <-  lapply(c("hs","mm"),function(spec){
spec.condit.combos <- condit.combos[condit.combos[,1]==spec,]
c.spec.cc <- combn(1:nrow(spec.condit.combos),2)
lapply(1:(ncol(c.spec.cc)),function(j){
  j <- c.spec.cc[,j]
  print(j)
  x <- unlist(spec.condit.combos[j[1],])
  y <- unlist(spec.condit.combos[j[2],])
  cdir <- paste0("~/code/data/cb/shiftFiles/DifferentExpressions/","DESEQ",paste0(x[1],x[2],y[1],y[2]))
  setwd(cdir)
  design.df <- data.frame("time"=tpMixSets[[x[1]]],
                          "condition"=conditMixSets[[x[1]]],stringsAsFactors = F)
  design.df$combo <- paste(as.character(design.df$condition),as.character(design.df$time),sep = ".")
  l <- lapply(unique(design.df$condition[conditMixSets[[x[1]]]==x[2]| conditMixSets[[y[1]]]==y[2]]),function(c){
    t.vect <- design.df$combo[design.df$condition==c]
    t.vect <- unique(c( ifelse(spec=="mm","0.0","100.0"),t.vect))
    lapply(1:(length(t.vect)-1),function(ti){
        print(design.df$combo)
        print(paste0(cdir,"/",c,t.vect[1],"_" ,t.vect[ti+1] ))
        #CountSet <- MixSetsEC[[spec]][, design.df$combo== t.vect[1] | design.df$combo== t.vect[ti+1] ]
        read.table(file = paste0(cdir,"/",c,t.vect[1],"_" ,t.vect[ti+1] ),header = T,row.names = 1)
    })
  }) 
  names(l) <- unique(design.df$condition[conditMixSets[[x[1]]]==x[2]| conditMixSets[[y[1]]]==y[2]])
  l
})
})

DE.ForGo <-  list("hs10"=DE.sets[[1]][[1]][[1]],"hs100"=DE.sets[[1]][[2]][[2]],"mm0" = DE.sets[[2]][[1]][[1]])
str(DE.ForGo[[1]][[1]])
DE.MixSets.Up <- lapply(DE.ForGo,function(x)sapply(x,function(y) rownames(y)[y$status=="DE" & y$PostFC <.71]))
DE.MixSets.Down <- lapply(DE.ForGo,function(x)sapply(x,function(y) rownames(y)[y$status=="DE" & y$PostFC >1.4]))
DE.ForGo$hs10[[4]][1:100,]
"ASCL1"%in% DE.MixSets.Up$hs10[[7]]
DE.MixSets$hs10[[4]]
DE.MixSets.Up$hs10[[10]]

MixSets[["hs"]]["SOX21",]

str(DE.MixSets)
DE.ForGo[[1]]
xx <- as.list(org.Hs.egALIAS2EG)
names(xx) <- toupper(names(xx))
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="useast.ensembl.org")
GOTermsHS <- getBM(attributes = c('hgnc_symbol', "go_id",
                                  "name_1006"), filters = "entrezgene",
                   values = sapply(xx[rownames(MixSets$hs)],"[",1), mart = mart)

xx <- as.list(org.Mm.egALIAS2EG)
names(xx) <- toupper(names(xx))
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="useast.ensembl.org")
GOTermsMM <- getBM(attributes = c('mgi_symbol', "go_id",
                                  "name_1006"), filters = "entrezgene",
                   values = unlist(sapply(xx[rownames(MixSets$mm)],"[",1)), mart = mart)


GOTermsHS$hgnc_symbol <- toupper(GOTermsHS$hgnc_symbol)
GOTermsMM$mgi_symbol <- toupper(GOTermsMM$mgi_symbol)
colnames(GOTermsMM)[1] <- colnames(GOTermsHS)[1] <- "Symbol"

GOTermsUD <- lapply(1:length(DE.ForGo),function(i){
  print(i)
  if(grepl("mm",names(DE.ForGo)[i])){
    dict <- GOTermsMM
  }else{
    dict <- GOTermsHS
  }
  apply(DE.ForGo[[i]],2,function(r){
    resultsUp <- dict$name_1006[dict$Symbol %in% unlist(r[1])]
    resultsDown <- dict$name_1006[dict$Symbol %in% unlist(r[2])]
    return(list("up"= resultsUp,"down"=resultsDown))
  })
})
save("GOTerms",file = "~/Desktop/GOobject.RData")


###################EBSEQ

GetDEResults <- function (EBPrelim, FDR = 0.05, Method = "robust", FDRMethod = "hard", 
          Threshold_FC = 0.7, Threshold_FCRatio = 0.3, SmallNum = 0.01) 
{
  if (!"PPDE" %in% names(EBPrelim)) 
    stop("The input doesn't seem like an output from EBTest")
  Conditions = EBPrelim$Conditions
  Levels = levels(as.factor(Conditions))
  PPcut = FDR
  GeneMat = EBPrelim$DataNorm
  PP = GetPPMat(EBPrelim)
  if (FDRMethod == "hard") {
    DEfound = rownames(PP)[which(PP[, "PPDE"] >= (1 - PPcut))]
  }
  else {
    SoftThre = crit_fun(PP[, "PPEE"], PPcut)
    DEfound = rownames(PP)[which(PP[, "PPDE"] >= SoftThre)]
  }
  if (Method == "classic") {
    Gene_status = rep("EE", dim(GeneMat)[1])
    names(Gene_status) = rownames(GeneMat)
    Gene_status[DEfound] = "DE"
    NoTest_genes = rownames(GeneMat)[!(rownames(GeneMat) %in% 
                                         rownames(PP))]
    Gene_status[NoTest_genes] = "Filtered: Low Expression"
    PPMatWith0 = EBPrelim$PPMatWith0
    PPMatWith0[NoTest_genes, ] = c(NA, NA)
    return(list(DEfound = DEfound, PPMat = PPMatWith0, Status = Gene_status))
  }
  else {
    PostFoldChange = PostFC(EBPrelim)
    PPFC = PostFoldChange$PostFC
    OldPPFC = PPFC[DEfound]
    OldPPFC[which(OldPPFC > 1)] = 1/OldPPFC[which(OldPPFC > 
                                                    1)]
    FilterFC = names(OldPPFC)[which(OldPPFC > Threshold_FC)]
    NewFC1 = apply(matrix(GeneMat[DEfound, which(Conditions == 
                                                   Levels[[1]])] + SmallNum, nrow = length(DEfound)), 
                   1, median)
    NewFC2 = apply(matrix(GeneMat[DEfound, which(Conditions == 
                                                   Levels[[2]])] + SmallNum, nrow = length(DEfound)), 
                   1, median)
    NewFC = NewFC1/NewFC2
    NewFC[which(NewFC > 1)] = 1/NewFC[which(NewFC > 1)]
    FCRatio = NewFC/OldPPFC
    FCRatio[which(OldPPFC < NewFC)] = 1/FCRatio[which(OldPPFC < 
                                                        NewFC)]
    FilterFCR = names(FCRatio)[which(FCRatio < Threshold_FCRatio)]
    Gene_status = rep("EE", dim(GeneMat)[1])
    names(Gene_status) = rownames(GeneMat)
    Gene_status[DEfound] = "DE"
    NoTest_genes = rownames(GeneMat)[!(rownames(GeneMat) %in% 
                                         rownames(PP))]
    Gene_status[NoTest_genes] = "Filtered: Low Expression"
    Filtered_DEfound = setdiff(DEfound, union(FilterFC, FilterFCR))
    PPMatWith0 = EBPrelim$PPMatWith0
    NAGenes = union(NoTest_genes, union(FilterFC, FilterFCR))
    PPMatWith0[NAGenes, ] = c(NA, NA)
    Gene_status[FilterFC] = "Filtered: Fold Change"
    Gene_status[FilterFCR] = "Filtered: Fold Change Ratio"
    return(list(DEfound = Filtered_DEfound, PPMat = PPMatWith0, 
                Status = Gene_status))
  }
}

runEBTest <- function(dataset,conditions, selected.conditions){
  conditin <- lapply(1:length(conditions),function(q) sapply(conditions[[q]],function(y) y%in%selected.conditions[[q]]))
  ind <- rowSums(data.frame(conditin)) == length(conditions)
  dataset <- dataset[,ind]
  conditions <- lapply(1:length(conditions),function(i)conditions[[i]][ind])
  #condit1 <-apply(data.frame(conditions,stringsAsFactors = F)[ind1,],1,function(x) Reduce(paste0,x) )
  #condit2 <-apply(data.frame(conditions,stringsAsFactors = F)[ind2,],1,function(x) Reduce(paste0,x) )
  condit <- lapply(1:length(selected.conditions),function(i)sapply(unique(selected.conditions[[i]]),function(y) y==conditions[[i]]))
  condit <- factor(apply(condit[[1]],1,which),levels = unique(apply(condit[[1]],1,which))) 
  Sizes= MedianNorm(dataset)
  EBOut <- EBTest(dataset, Conditions = condit,sizeFactors = Sizes, maxround=5)
  EBOut$DataNorm <- dataset
  EBDERes=GetDEResults(EBOut, FDR=0.05,Method = "robust")
  p = PostFC(EBOut)
  EBMat <- data.frame(Reduce( cbind,list(EBDERes$PPMat[names(p$PostFC),],"PostFC"=p$PostFC,p$RealFC)))
  EBMat <- merge(EBMat,EBDERes$Status ,by="row.names",all.x=TRUE)
  rownames(EBMat) <- EBMat[,1] 
  EBMat <- EBMat[,-1] 
  colnames(EBMat) <- c("PPEE","PPDE" ,"PostFC", "RealFC","status")
  EBMat <- EBMat[order(EBMat$PostFC,decreasing = T),]
  EBMat
}

###########
DE.MixSets <- mclapply(c("hs","mm"),function(spec){
  spec.condit.combos <- condit.combos[condit.combos[,1]==spec,]
  c.spec.cc <- combn(1:nrow(spec.condit.combos),2)
  mclapply(1:ncol(c.spec.cc),function(j){
    j <- c.spec.cc[,j]
    x <- unlist(spec.condit.combos[j[1],])
    y <- unlist(spec.condit.combos[j[2],])
    cdir <- paste0("~/code/data/cb/shiftFiles/DifferentExpressions/","EBSEQ",x[1],x[2],y[1],y[2])
    dir.create(cdir)
    print(cdir)
    design.df <- data.frame("time"=tpMixSets[[x[1]]],
                                 "condition"=conditMixSets[[x[1]]],stringsAsFactors = F)
    design.df$combo <- paste(as.character(design.df$condition),as.character(design.df$time),sep = ".")
    #design.df <- data.frame("time"=c(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ], tpMixSets[[y[1]]][conditMixSets[[y[1]]]==y[2] ]),
    #                        "condition"=c( conditMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ],conditMixSets[[y[1]]][conditMixSets[[y[1]]]==y[2] ]),stringsAsFactors = F)
    #design.df$combo <- paste(as.character(design.df$condition),as.character(design.df$time),sep = ".")
    
    apply(combn(unique(design.df$condition[conditMixSets[[x[1]]]==x[2]| conditMixSets[[y[1]]]==y[2]]),2),2, function(c){
      if(!file.exists(paste0(cdir,"/",c[1],"_",c[2] ))){
      print(c)
      DESOutResults <- runEBTest(dataset=MixSetsEC[[spec]], conditions = list(design.df$condition),selected.conditions=  list(c(design.df$condition[design.df$condition==c[1]], design.df$condition[design.df$condition==c[2]])))
      write.table(DESOutResults,file = paste0(cdir,"/",c[1],"_",c[2] ))
      }
    })
    mclapply(unique(design.df$condition[conditMixSets[[x[1]]]==x[2]| conditMixSets[[y[1]]]==y[2]]),function(c){
      print(c)
      t.vect <- design.df$combo[design.df$condition==c]
      t.vect <- unique(c( ifelse(spec=="mm","0.0","100.0"),t.vect))
      sapply(1:(length(t.vect)-1),function(ti){
        if(!file.exists(paste0(cdir,"/",c,t.vect[ti],"_" ,t.vect[ti+1] ))){
        print( t.vect[ti])
        print( t.vect[ti+1])
        DESOutResults <- runEBTest(dataset =  MixSetsEC[[spec]],conditions =  list(design.df$combo),selected.conditions =  list(c(unique( t.vect[ti]),unique(t.vect[ti+1]))))
        write.table(DESOutResults,file = paste0(cdir,"/",c,t.vect[ti],"_" ,t.vect[ti+1] ))
      }
      })
    })
    mclapply(unique(design.df$condition[conditMixSets[[x[1]]]==x[2]| conditMixSets[[y[1]]]==y[2]]),function(c){
      t.vect <- design.df$combo[design.df$condition==c]
      t.vect <- unique(c( ifelse(spec=="mm","0.0","100.0"),t.vect))
      print(t.vect)
      sapply(1:(length(t.vect)-1),function(ti){
        if(!file.exists(paste0(cdir,"/",c,t.vect[1],"_" ,t.vect[ti+1] ))){
          print(design.df$combo)
          print(dim(MixSetsEC[[spec]]))
          #CountSet <- MixSetsEC[[spec]][, design.df$combo== t.vect[1] | design.df$combo== t.vect[ti+1] ]
          DESOutResults <- runEBTest(dataset =   MixSetsEC[[spec]],conditions =  list(design.df$combo),selected.conditions =  list(c(t.vect[1],t.vect[ti+1])))
          write.table(DESOutResults,file = paste0(cdir,"/",c,t.vect[1],"_" ,t.vect[ti+1] ))
        }
      })
    })
    #DESOutResults$padj <- p.adjust(DESOutResults$pvalue ,method = "BH")
  })
})


#inds <- lapply(inds, function(x)x[-length(x)] )
vs0 <- T
DE.MixSets <- mclapply(c("hs","mm"),function(s){
  sapply( 1:length(inds),function(x){
    if(grepl(s,names(inds)[[x]])){
      EBOut <- sapply(1:(length(inds[[x]])-1),function(i){
        cdir <- paste0("~/code/data/cb/shiftFiles/DifferentExpressions/",names(inds)[x])
        if(grepl(s,names(inds)[[x]])){
          if(vs0){
            ind1 <- inds[[x]][1]
            ind2 <-inds[[x]][i+1]
          }else{
            ind1 <- inds[[x]][i]
            ind2 <-inds[[x]][i+1]
          }
          if(!file.exists(paste0(cdir,"/",names(inds)[x],"_" ,ind1,"_",ind2 ))){
            
            dir.create(cdir)
            subs <- MixSetsEC[[s]][,c(ind1,ind2),drop=F]
            tp <- tps[[names(inds)[x]]]
            condit <- factor(c(ind1, ind2))
            Sizes=MedianNorm(subs)
            #         print(i)
            #         print(condit)
            EBOut <- EBTest(subs, Conditions = condit,sizeFactors = Sizes, maxround=5)
            EBDERes=GetDEResults(EBOut, FDR=0.05,Method = "robust")
            p = PostFC(EBOut)
            #         print(head(p))
            #         print(head(EBDERes))
            EBMat <- data.frame(Reduce( cbind,list(EBDERes$PPMat[names(p$PostFC),], p$PostFC,p$RealFC)))
            EBMat <- merge(EBMat,EBDERes$Status ,by="row.names",all.x=TRUE)
            rownames(EBMat) <- EBMat[,1] 
            EBMat <- EBMat[,-1] 
            colnames(EBMat) <- c("PPEE","PPDE" ,"PostFC", "RealFC","status")
            EBMat <- EBMat[order(EBMat$PostFC,decreasing = T),]
            write.table(EBMat,file = paste0(cdir,"/",names(inds)[x],"_" ,ind1,"_",ind2 ))
          }else{
            EBMat <- read.table(file = paste0(cdir,"/",names(inds)[x],"_" ,ind1,"_",ind2 ),header = T,row.names = 1)
          }
          up <- rownames(EBMat)[EBMat$PPEE<.01 & EBMat$PostFC > 4 & !is.na(EBMat$PPEE) ]
          down <- rownames(EBMat)[EBMat$PPEE<.01 & EBMat$PostFC < 1/4 & !is.na(EBMat$PPEE) ]
          return(c(list(up),list(down)))
        }
      })
    }
  })
})



library(DESeq2)
library(EBSeq)
int_rns <-  intersect(rownames(datasets$`~/code/data/cb/TC_2_Mix.csv`), rownames(datasets$`~/code/data/cb/TC_2_H.csv`))
experiments = list("hs100"=MixSetsEC$hs[,77:102],"hs90"= MixSetsEC$hs[,c(77,52:76)] , "hs10"=MixSetsEC$hs[,c(77,27:51)]  , "mm0"=MixSetsEC$mm[,1:26], "mm10"=MixSetsEC$mm[,c(1,27:51)] , "mm90"=MixSetsEC$mm[,c(1,52:76)], "hs33dayMix"=datasets$`~/code/data/cb/TC_2_Mix.csv`[int_rns,], "hs33dayH"=datasets$`~/code/data/cb/TC_2_H.csv`[int_rns,])
tps <- list("hs100"=tpMixSets[77:102],"hs90"= tpMixSets[c(77,52:76)] , "hs10"=tpMixSets[c(77,27:51)]  , "mm0"=tpMixSets[1:26], "mm10"=tpMixSets[c(1,27:51)] , "mm90"=tpMixSets[c(1,52:76)], "hs33dayMix"=tpMix[[2]], "hs33dayH"=tpMix[[2]])
experimentGroups <- factor(c(1,1,1,2,2,2,3,3))
DEGenes <-  sapply(levels(experimentGroups),function(g){
  curGroup = experiments[experimentGroups==g]
  combos <-  combn(names(curGroup), 2)
  sapply(1:ncol(combos), function(c){
    i= combos[1,c]
    j= combos[2,c]
    if(i != j ){
      pathh = paste0("~/code/data/cb/shiftFiles/DifferentExpressions/",make.names(i),"vs.",make.names(j))
      if(!file.exists(file.path(pathh,"DEOutput.txt"))){
        DESDS <-  DESeqDataSetFromMatrix(countData = cbind(round(curGroup[[i]]), round(curGroup[[j]])), data.frame("time"=factor(c(tps[[i]], tps[[j]])), "group" = factor(c(rep(i,ncol(curGroup[[i]])),rep(j,ncol(curGroup[[j]]))))),design = ~time+group) 
        DESOut <- DESeq(DESDS,minReplicatesForReplace = 100)
        DESOutResults <- results(DESOut,contrast = c("group",i,j))
        DESOutResults <- DESOutResults[order(DESOutResults$log2FoldChange,decreasing = T),]
        DESOutResults$padj <- p.adjust(DESOutResults$pvalue ,method = "BH")
        EBMat <- DESOutResults
        #       condit <- factor(c(rep(i,ncol(curGroup[[i]])),rep(j,ncol(curGroup[[j]]))),levels = c(i,j) ) 
        #       Sizes= MedianNorm(cbind(curGroup[[i]], curGroup[[j]]))
        #       EBOut <- EBTest(cbind(curGroup[[i]], curGroup[[j]]), Conditions = condit,sizeFactors = Sizes, maxround=5)
        #       EBDERes=GetDEResults(EBOut, FDR=0.05)
        #       p = PostFC(EBOut)
        #       EBMat <- data.frame(Reduce( cbind,list(EBDERes$PPMat[names(p$PostFC),],"PostFC"=p$PostFC,p$RealFC)))
        #       EBMat <- merge(EBMat,EBDERes$Status ,by="row.names",all.x=TRUE)
        #       rownames(EBMat) <- EBMat[,1] 
        #       EBMat <- EBMat[,-1] 
        #       colnames(EBMat) <- c("PPEE","PPDE" ,"PostFC", "RealFC","status")
        #       EBMat <- EBMat[order(EBMat$PostFC,decreasing = T),]
        hist(DESOutResults$log2FoldChange)
        
        #up <- rownames(EBMat)[EBMat$PPEE<.05 & EBMat$PostFC > 2 & !is.na(EBMat$PPEE) ]
        #down <- rownames(EBMat)[EBMat$PPEE<.05 & EBMat$PostFC < .5 & !is.na(EBMat$PPEE) ]
        dir.create(pathh)
        write.table(EBMat,file = paste0(pathh,"/","DEOutput.txt"),sep = "\t",quote = F )
        
      }else{
        EBMat <- read.table(file = paste0(pathh,"/","DEOutput.txt"),sep = "\t",header = T,row.names = 1)
      }
      up <- rownames(EBMat)[EBMat$padj<.05 & EBMat$log2FoldChange >1.584 & !is.na(EBMat$padj) ]
      down <- rownames(EBMat)[EBMat$padj<.05 & EBMat$log2FoldChange < -1.584 & !is.na(EBMat$padj) ]
      other = setdiff(names(curGroup), c(i,j) )
      
      if(length(up)==0L)up <- rownames(EBMat)[1:5]
      if(length(down)==0L)down <- rownames(EBMat)[(nrow(EBMat)-5):nrow(EBMat) ]
      
      print(up)
      print(down)
      
      for(ggg in c("up","down")){
        write.pdf({
          mypar(3,3)
          gg <-  ifelse(ggg == "up", list(up), list(down))
          for(gn in gg[[1]]){
            if(grepl("hs",i)){
              if(gn %in% rownames(experiments[["mm0"]])){
                x1 <- zoo(experiments[["mm0"]][gn,],tps[["mm0"]])
                x2 <- zoo(experiments[["mm10"]][gn,],tps[["mm10"]])
              }else{
                x1 <- x2 <- zoo(rep(0.0,ncol(experiments[["mm0"]])),tps[["mm0"]])
              }
              nx1 <- "mm0"
              nx2 <- "mm10"
            }else{
              if(gn %in% rownames(experiments[["hs,0%"]])){
                x1 <- zoo(experiments[["hs100"]][gn,],tps[["hs100"]])
                x2 <- zoo(experiments[["hs10"]][gn,],tps[["hs10"]])
              }else{
                x1 <- x2 <- zoo(rep(0.0,ncol(experiments[["hs100"]]),tps[["hs100"]]))
              }
              nx1 <- "hs100"
              nx2 <- "hs10"
            }
            s1 <- zoo(curGroup[[i]][gn,],tps[[i]])  
            s2 <- zoo(curGroup[[j]][gn,],tps[[j]])
            if(length(other)!= 0L){ 
              so <- zoo(curGroup[[other]][gn,],tps[[other]]) 
              a = na.approx(Reduce(merge,list(s1,s2,so,x1,x2)))
              plot.zoo(a,plot.type = "single",col = 1:6, lty = 1:2,main=paste(gn),xlab = "TimePoint",ylab = "Counts" )
              legend("topleft", c(i,j,other,nx1,nx2), col = 1:6, lty = 1:2)
            }else{
              a = na.approx(Reduce(merge,list(s1,s2,x1,x2)))
              plot.zoo(a,plot.type = "single",col = 1:6, lty = 1:2,main=paste(gn),xlab = "TimePoint",ylab = "Counts" )
              legend("topleft", c(i,j,nx1,nx2), col = 1:6, lty = 1:2)
            }
          }},filename = paste0(pathh,"/",ggg,"DE.pdf"))}
      
      
      if(grepl("mm",i)){
        libspec = "mouse"
        xx <- as.list(org.Mm.egALIAS2EG)
        names(xx) <- toupper(names(xx))
      }else{
        libspec = "human"
        xx <- as.list(org.Hs.egALIAS2EG)
        names(xx) <- toupper(names(xx))
      }
      
      u=enrichGO(gene= sapply(xx[up],"[",1), organism= libspec,universe = unlist(sapply(xx[rownames(EBMat)],"[",1)),ont = "BP")
      d=enrichGO(gene= sapply(xx[down],"[",1), organism= libspec,universe = unlist(sapply(xx[rownames(EBMat)],"[",1)),ont = "BP")
      
      write.table(summary(u),file = paste0(pathh,"/","UpGOTerms.txt"),sep = "\t",row.names = F ,quote = F )
      write.table(summary(d),file = paste0(pathh,"/","DownGOTerms.txt"),sep ="\t",row.names = F ,quote = F )
    }
  })
})


