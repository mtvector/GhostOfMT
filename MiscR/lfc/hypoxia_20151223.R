library(EBSeq)
datapath <- file.path("~/code/data/lfc")
EndoMarkers <- as.character(read.table("~/code/data/general/marker_lists/Endoderm/endoderm.txt")[,1])
ESMarkers <- as.character(read.table("~/code/data/general/marker_lists/ES_Germ/ES.txt")[,1])

TPMHypoxiaFile <- file.path(datapath, "Sub347_H1H9DF1911_TPM.txt")
TPMHypoxia.frame <- read.csv(TPMHypoxiaFile, header = T, sep = "\t" )
rownames(TPMHypoxia.frame) <-TPMHypoxia.frame[,1]
TPMHypoxia <- as.matrix(TPMHypoxia.frame[,c(-1, -length(TPMHypoxia.frame))])
TPMHypoxia <- TPMHypoxia[,  -which(colnames(TPMHypoxia)=="H1_LowO2_rep2" ) ]
sizes.Hypoxia <- MedianNorm(TPMHypoxia)
TPMHypoxiaNorm <- log(GetNormalizedMat(TPMHypoxia, sizes.Hypoxia)+1,2)


#load dataset
ECHypoxiaFile <- file.path(datapath, "Sub347_H1H9DF1911_EC.txt")
ECHypoxia.frame <- read.csv(ECHypoxiaFile, header = T, sep = "\t" )
rownames(ECHypoxia.frame) <-ECHypoxia.frame[,1]
ECHypoxia <- as.matrix(ECHypoxia.frame[,-1])
#Remove outlier
ECHypoxia <- ECHypoxia[,  -which(colnames(ECHypoxia)=="H1_LowO2_rep2" ) ]
sizes.Hypoxia <- MedianNorm(ECHypoxia)

#EndoMarkers <- c("CXCR4", "SOX17", "FOXA1","FOXA2", "KRT19", "CER1", "KIT", "HNF1B", "CDX2", "VWF", "GATA4", "CD14", "CD36") #"CD31","GD3", "VEGF"

line <- sapply(colnames(ECHypoxia), function(x){
  strsplit(x, split = "_")[[1]][1]
})

cond <- sapply(colnames(ECHypoxia), function(x){
  strsplit(x, split = "_")[[1]][2]
})

day <- sapply(colnames(ECHypoxia), function(x){
  s <-  strsplit(x, split = "_")[[1]][3]
  if(grepl("rep", s)){
    s <- "d0"
  }
  return(s)
})

cond.day <-  paste(cond,day , sep = "_")

#O2.Contrast.DE <-  EBMultiTest(ECHypoxia, Conditions = as.factor(cond.day) , sizeFactors =sizes.Hypoxia , Print = T,maxround = 7)

each.line.contrasts.same.times <-  sapply(unique(line), function(l){
  lapply(unique(day),function(d){
    col.ind <-  line == l & day == d
    EH <-  ECHypoxia[,col.ind]
    EBMultiTest(EH, Conditions= as.factor(cond.day[col.ind]),sizeFactors=MedianNorm(EH),Print =T, maxround = 7)
  })
})
p <- as.matrix(each.line.contrasts.same.times)
rownames(p) <- unique(day)

DE.Results.same.times <-  apply( p, 1:2, function(x){
  GetMultiPP( unlist(x,recursive = F) )
} )

FC.Results.same.times <- apply( p, 1:2, function(x){
  GetMultiFC(unlist(x,recursive = F))
}  )

for(i in 1:ncol(FC.Results.same.times)){
  for(j in 1:nrow(FC.Results.same.times)){
    EM <- EndoMarkers[ EndoMarkers %in% rownames(FC.Results.same.times[[j,i]]$Log2PostFCMat)]
    ESM <- ESMarkers[ ESMarkers %in% rownames(FC.Results.same.times[[j,i]]$Log2PostFCMat)]
    myplot( expressionHeatmap(FC.Results.same.times[[j,i]]$Log2PostFCMat[EM,]) ,file.path(datapath, paste0(colnames(FC.Results.same.times)[i], rownames(FC.Results.same.times)[j],"EndoMarkers", ".pdf")))
    myplot( expressionHeatmap(FC.Results.same.times[[j,i]]$Log2PostFCMat[ESM,]), file.path(datapath, paste0(colnames(FC.Results.same.times)[i], rownames(FC.Results.same.times)[j],"ESMarkers", ".pdf")))
  }
}

each.line.contrasts.d0.d3 <- sapply(unique(line), function(l){
  lapply(unique(cond),function(c){
    col.ind <-  line == l & cond == c
    EH <-  ECHypoxia[,col.ind]
    EBTest(EH, Conditions= as.factor(cond.day[col.ind]),sizeFactors=MedianNorm(EH), maxround = 7)
  })
})
q <-  as.matrix(each.line.contrasts.d0.d3)
rownames(q) <- unique(cond)

DE.Results.d0.d3 <-  apply( q, 1:2, function(x){
  GetDEResults( unlist(x,recursive = F) , Method = "classic", SmallNum = .01)
} )

FC.Results.d0.d3 <- apply( q, 1:2, function(x){
 PostFC( unlist(x,recursive = F))
}  )

apply(FC.Results, 1:2 , function(x){
  x[[1]]$PostFC[EndoMarkers]
  #median(x[[1]]$PostFC)
})

apply(FC.Results, 1:2 , function(x){
  x[[1]]$PostFC["SOX17"]
})

DE.Attack <- function(ECMat, names, line=NULL, time=NULL, condition=NULL){
  
  each.line.contrasts.diff.conditions <-  sapply(unique(line), function(l){
    lapply(unique(time),function(t){
      col.ind <-  line == l & time == t
      EH <-  ECMat[,col.ind]
      EBMultiTest(EH, Conditions= as.factor(names[col.ind]),sizeFactors=MedianNorm(EH),Print =T, maxround = 7)
    })
  })
  
  each.line.contrasts.diff.times <- sapply(unique(line), function(l){
    lapply(unique(cond),function(c){
      col.ind <-  line == l & cond == c
      EH <-  ECMat[,col.ind]
      EBTest(EH, Conditions= as.factor(names[col.ind]),sizeFactors=MedianNorm(EH), maxround = 7)
    })
  })
  
  
  p <- as.matrix(each.line.contrasts.same.times)
  rownames(p) <- unique(time)
  
  DE.Results.same.times <-  apply( p, 1:2, function(x){
    GetMultiPP( unlist(x,recursive = F) )
  } )
  
  FC.Results.same.times <- apply( p, 1:2, function(x){
    GetMultiFC(unlist(x,recursive = F))
  }  )
  
  
  each.line.contrasts.d0.d3 <- sapply(unique(line), function(l){
    lapply(unique(cond),function(c){
      col.ind <-  line == l & cond == c
      EH <-  ECHypoxia[,col.ind]
      EBTest(EH, Conditions= as.factor(names[col.ind]),sizeFactors=MedianNorm(EH), maxround = 7)
    })
  })
  q <-  as.matrix(each.line.contrasts.d0.d3)
  rownames(q) <- unique(cond)
  
  DE.Results.d0.d3 <-  apply( q, 1:2, function(x){
    GetDEResults( unlist(x,recursive = F) , Method = "classic", SmallNum = .01)
  } )
  
  FC.Results.d0.d3 <- apply( q, 1:2, function(x){
    PostFC( unlist(x,recursive = F))
  }  )
  
  
}

