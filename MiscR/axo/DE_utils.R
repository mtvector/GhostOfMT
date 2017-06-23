#source("http://bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("AnnotationDbi")
#biocLite("IRanges")
#biocLite("GenomeInfoDb")
#biocLite("genefilter")
#install.packages("devtools")
#install.packages("rafalib")
#install.packages("RColorBrewer")
library(devtools)
library(rafalib)
library(RColorBrewer)
library(genefilter)
library(Biobase)
#biocLite("EBSeq")
#library(foreign)
#install(file.path(path ,"TimeShift"))
#library(TimeShift)
library(EBSeq)
#install.packages("rgl")
#install.packages("shapes")
#library(dplyr)
#library(shapes)
#biocLite("biomaRt")
library(biomaRt)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("gplots")
library(gplots)
#install.packages("dtw")
#library("dtw")
#install.packages("matrixStats")
library(matrixStats)
#biocLite("DESeq2")
#library(DESeq2)
#install.packages("randomForest")
library(randomForest)
#install.packages("biclust")
#library(biclust)
#biocLite("CoRegNet")
#library(CoRegNet)
library(contrast)
library(limma)
library(RColorBrewer)
library(RCurl)

summary <- function(x) {
  funs <- c(mean, median, sd, mad, IQR)
  lapply(funs, function(f) f(x, na.rm = TRUE))
}

?summary

datapath <- "~/code/data/axo"

#Get the GO terms associated with differentially expressed genes
getGOHuman <- function(geneList){
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  results <- getBM(attributes = c( 'hgnc_symbol', 'go_id', 'name_1006'), filters = 'hgnc_symbol',
                   values = geneList, mart = mart)
  return(results)
}

#Takes values and a set of genes and some params
#prints top 42 GO Terms and gene

printTopGO <- function(terms, num){
  .a <- environment()
  GOFactor <- factor(terms)
  tt <- data.frame(table(GOFactor))
  tt <- tt[order(tt$Freq, decreasing = TRUE),]
  tt <- tt[1:num,]
  tt$GOFactor <- as.character(tt$GOFactor)
  tt$GOFactor <- factor(tt$GOFactor, levels = unique(tt$GOFactor))
  mypar(1,1)
  return(qplot(tt$GOFactor,tt$Freq ,stat="identity", geom = "bar",environment=.a) + coord_flip())
}

fractionCharToNumeric <- function(string){
  string <- as.character(string)
  a <-  unlist(strsplit(string, "/"))
  as.numeric(a[1])/as.numeric(a[2])
}

showPlot <- function(timePoint, GL){
  mypar(10,4)
  for(i in 1:length(GL[[timePoint]][[1]])){
    gH <- as.matrix(GL[[timePoint]][[2]][[i]])
    gB <- as.matrix(GL[[timePoint]][[3]][[i]])
    plot(gH[1,], main=paste("Embryo", rownames(gH)[1]),xlab="TimePoint" ,ylab="TPM")
    plot(gB[1,], main=paste("Blastema", rownames(gB)[1]),xlab="TimePoint" ,ylab="TPM") 
  }
}

addPrefix <- function(prefix,list){
  lapply(list, function(x){
    paste0(prefix,x)
  })
}

contrastString <- function(nameList1,nameList2){
  string <- "("
  lapply(nameList1, function(x){
    string <<- paste0(string,x,"+")
  })
  string <-  substr(string ,1, nchar(string)-1)
  string <-  paste0(string,")/",as.character(length(nameList1))," - ")
  string2 <- "("
  lapply(nameList2, function(x){
    string2 <<- paste0(string2,x,"+")
  })
  string2 <-  substr(string2 ,1, nchar(string2)-1)
  string2 <- paste0(string2,")/",as.character(length(nameList2)))
  paste0(string,string2)
}

Amby002gene <- function(probe){
  idx <- chipAnnotations$Probe.Set.ID == probe
  if(chipAnnotations[idx,]$Gene!="" | chipAnnotations[idx,]$Gene!="Err:512"){
    return(chipAnnotations[idx,]$Gene)
  }
  else{
   return(chipAnnotations[idx,]$Source.Seq) 
  }
}

AgilentAmbyGene <- function(probe){
  tab =  fData(GSE36451)
  idx <- tab$ID == probe
  g <- as.character(tab[idx,]$GENE_SYMBOL)
  if( g!="" | g !="N/A" | g !="no hit"){
    return(g)
  }
  else{
    return(tab[idx,]$ID ) 
  }
}


expressionHeatmap <- function(expression, maint="Cors"){
  cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(20)
  heatmap.2(expression, main=maint
            ,na.rm = TRUE, trace = "none",
            srtRow=0, srtCol=60, cexCol = .45, cexRow = .45,
            dendrogram = "column", Rowv = FALSE, #Colv = FALSE,
            col = cols)
} 

hist.normalized <- function(x){
  nsamp <- dim(x)[2]
  h <- hist(x[,1], plot=FALSE)
  plot(h$mids, h$density, type="l", col=rainbow(nsamp)[1], main="",
       xlab="Log normalized expression value", ylab="Proportion of molecules")
  for(i in 2:nsamp){
    h <- hist(x[,i], plot=FALSE)
    lines(h$mids, h$density, col=rainbow(nsamp)[i])
  }
}

getSingleGo <- function(goreport){
  decent=F
  i<<-1
  term <- ""
  badterms <- c( "macromolecule biosynthetic process","cellular macromolecular biosynthetic process","biological_process","biological process", "single-organism process", "cellular process", "cellularprocess","homeostatic process", "RNA processing","regulation of RNA metabolic process" ,"regulation of RNA biosynthetic process", "regulation of metabolic process","regulation of cellular biosynthetic process","regulation of cellular macromolecule biosynthetic process")
  while(!decent){
    term <- as.data.frame(summary(goreport))[i, "Description"]
    if(term %in% badterms){
      i <<- i+1
    }else{
      decent <<- T
      break
    }
  }
  return(term)
}

getDEMatchListMatrix <- function(listA, listB){
  sMat <- matrix(list(),nrow=length(listA), ncol=length(listB))
  for(i in 1:length(listA)){
    DEGenesA <-  listA[[i]]
    for(j in 1:length(listB)){
      for(k in 1:length(listB[[j]])){
        print(i)
        print(j)
        print(k)
        DEGenesB <- listB[[j]]
        if(DEGenesB[k]%in%DEGenesA){
          sMat[i,j][[1]] <- c(sMat[i,j][[1]],as.character(DEGenesB[k][1]))
        }
      }
    }
  } 
  return(sMat)
}



getDEMatchListMatrix <- function(listA, listB){
  sMat <- matrix(list(),nrow=length(listA), ncol=length(listB))
  for(i in 1:length(listA)){
    for(j in 1:length(listB)){
      sMat[i,j] <- list(intersect(listA[[i]], listB[[j]]))
    }
  } 
  return(sMat)
}

getDEMatchRefined <- function(listMatrix, testlist=NULL){
  if(is.null(testlist)){
    return(matrix(sapply(listMatrix, FUN = length), nrow = nrow(listMatrix),ncol=ncol(listMatrix)))
  }
  sMat <- matrix(0,nrow=nrow(listMatrix), ncol=ncol(listMatrix))
  for(i in 1:nrow(listMatrix)){
    for(j in 1:ncol(listMatrix)){
      listMatrixGenes <-  listMatrix[i,j][[1]]
      for(k in 1:length(testlist)){
          if(testlist[k]%in%listMatrixGenes){
            sMat[i,j] <- sMat[i,j] + 1
          }
      }
    }
  } 
  return(sMat)
}

##Returns matrix of [timepoint in EB, timepoint in blastema it is DE] = number of genes from EB that are expressed at blastema timepoint
zscoreDE <- function(mat1, mat2){
  sMat <- matrix(list(),nrow=ncol(mat1), ncol=ncol(mat2))
  for(i in 1:ncol(mat1)){
    for(j in 1:ncol(mat2)){
      for(k in 1:nrow(mat1)){
        if(!is.na(mat1[k,i])&& !is.na(mat2[k,j])){
          if(mat1[k,i]==TRUE &&mat2[k,j]==TRUE)  
            sMat[i,j][[1]] <- c(sMat[i,j][[1]], rownames(mat1)[k])
          }
        }
      }
    }
  return(sMat)
}

#A is H, Blastema is B
heatmapMatch <- function(Match ,main, positive=TRUE, cluster=TRUE, xLab = "Blastema Time Point", yLab="Embryo Stage",...){  
  cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)
    if(!cluster){
     heatmap.2(Match, xlab = xLab, ylab=yLab, main=paste(main)
               ,na.rm = TRUE, trace = "none",colsep=1:ncol(Match),
               rowsep=1:nrow(Match),srtRow=0, srtCol=60,
                dendrogram = "none", Rowv = FALSE, Colv = FALSE,
                col=cols,...)     
    }else{
      heatmap.2(Match, xlab = xLab, ylab=yLab, main=paste("(+) DE\n",main)
                ,na.rm = TRUE, trace = "none",colsep=1:ncol(Match),
                rowsep=1:nrow(Match),srtRow=0, srtCol=60,
                col=cols,...)     
    }
}

fisherTestMatch <- function(intersectMat, aMat, bMat, numG = 5000){
  matrix(sapply(1:length(intersectMat), function(i){
    A <- aMat[[i]]
    B <- bMat[[i]]
    X <- intersectMat[[i]]
    fisher.test(matrix(X, B-X ,  A-X , numG-A-B-X))$p.value
  }), nrow = nrow(aMat),ncol=ncol(bMat))
}

hypergeoMatch <- function(intersectMat, aMat, bMat, numG=10000){
  aMat <- matrix(sapply(aMat, FUN = length), nrow = length(aMat),ncol=length(bMat))
  bMat <- matrix(sapply(bMat, FUN = length), nrow = length(aMat),ncol=length(bMat))
  matrix(sapply(1:length(intersectMat), function(i){
    A <- aMat[[i]]
    B <- bMat[[i]]
    AuB <- sum(aMat[,1])+sum(bMat[1,])
    phyper(intersectMat[[i]]-1, A, numG-A, B, lower.tail=F)
  }), nrow = nrow(aMat),ncol=ncol(bMat))
}


getZScoreMatches <- function(matchH, matchB,upperLim=1.5, lowerLim = -1.5, up=TRUE){
  zscoreH <-  (matchH-rowMeans(matchH))/rowSds(matchH)
  zscoreB <-  (matchB-rowMeans(matchB))/rowSds(matchB)
  zscoreUpH <- zscoreH > upperLim
  zscoreUpB <- zscoreB > upperLim
  zscoreDownH <- zscoreH < lowerLim
  zscoreDownB <- zscoreB < lowerLim
  zscoreDEUp <-  zscoreDE(zscoreUpB, zscoreUpH)
  zscoreDEDown <-  zscoreDE(zscoreDownB, zscoreDownH)
  
  zscoreUpMatch <- matrix(sapply(zscoreDEUp, FUN = length), nrow = nrow(zscoreDEUp),ncol=ncol(zscoreDEUp))
  zscoreDownMatch <- matrix(sapply(zscoreDEDown, FUN = length), nrow = nrow(zscoreDEDown),ncol=ncol(zscoreDEDown))
  if(up){
    return(zscoreUpMatch)
  }else{
    return(zscoreDownMatch)
  }
}

pvalPlot <- function(pMatrix , mainT,...) {
    pCol <- c( "red", "orange", "yellow", "green","blue", "gray60", "gray15")
  # Defining breaks for the color scale
  myBreaks <- c(0, 0.0000005, 0.00005, .0005, .005, 0.05, 0.1, 1)
  heatmap.2(pMatrix, scale="none",... ,
            main = mainT ,col = pCol, ## using your colors
            breaks = myBreaks, ## using your breaks
            na.rm = T, trace = "none",colsep=1:ncol(pMatrix),
            rowsep=1:nrow(pMatrix),srtRow=0, srtCol=60,
            dendrogram = "none", Rowv = F, Colv = F,
            cexRow=0.8, cexCol=0.8, key=F, keysize=1.5)
  legend(x="topright", fill = pCol, cex=0.6,
         legend = c("0 to 5e-6", "5e-6 to 5e-5", "5e-5 to 5e-4", "5e-4 to 5e-3", "5e-3 to 5e-2", "5e-2 to 0.1", "0.1 to 1"))
}






warpWindow <- function(shiftLimHi=.1,shiftLimLo = -.1 ,devRangeH = 6:17, devRangeB=1:12, corMin=.50,geneSet=NULL){
  if(!is.null(geneSet)){
    gB <- matchB[geneSet,]
    gH <- matchH[geneSet,]
  }else{
    gB <- matchB
    gH <- matchH
  }
  rm <- rowMeans(HvBshift)
  passFilter <- sapply(1:length(geneSet), function(i){
    rm[i] < shiftLimHi && rm[i]> shiftLimLo && rm[i]!=0 && cor(gH[i,devRangeH],gB[i,devRangeB]) >corMin
  })
  index<- which(unlist(passFilter))
  if(sum(passFilter, na.rm = TRUE)==0){
    return(NULL)
  }
  print(sum(passFilter, na.rm = TRUE))  
  Bplots <- list()
  Hplots <- list()
  for(i in 1:length(index)){
    x <- index[i]
    .e <- environment()
    plotH <- as.data.frame(t(gH[x,]))
    plotB <- as.data.frame(t(gB[x,]))
    colnames(plotH)<- colnames(gH)
    colnames(plotB)<- colnames(gB)
    rownames(plotH)<- rownames(gH)[x] 
    rownames(plotB)<- rownames(gB)[x] 
    Hplots[[length(Hplots)+1]] <- plotH
    Bplots[[length(Bplots)+1]] <- plotB
  }
  #Need to execute functions in GO analysis section first
  return(list(rownames(gH)[index], Hplots,Bplots))
}




##Goterms is a list of terms
getRepresentativeGo <- function(goterms){
  library(tm)
  stop_words <- stopwords("SMART")
  
  #goterms <-  sumGO[, c("Description")]
  goterms <-  paste(unlist(goterms), collapse =" " )
  
  reviews <- goterms
  
  reviews <- gsub(" of ", "-of-", reviews)
  reviews <- gsub("cell ", "cell-", reviews)
  reviews <- gsub(" cell ", "-cell", reviews)
  reviews <- gsub(" II", "-II", reviews)
  reviews <- gsub(" 3", "-3", reviews)
  reviews <- gsub(" system", "-system", reviews)
  #reviews <- gsub("RNA ", "RNA-", reviews)
  #reviews <- gsub("DNA ", "DNA-", reviews)
  reviews <- gsub(" regulation", "-regulation", reviews)
  reviews <- gsub(" response", "-response", reviews)
  reviews <- gsub(" process", "-process", reviews)
  reviews <- gsub(" type", "-type", reviews)
  reviews <- gsub(" activity", "-activity", reviews)
  reviews <- gsub(" beta", "beta-", reviews)
  reviews <- gsub(" pathway", "-pathway", reviews)
  reviews <- gsub(" transport", "-transport", reviews)
  reviews <- gsub(" complex", "-complex", reviews)
  reviews <- gsub(" protein", "-protein", reviews)
  reviews <- gsub(" complex", "-complex", reviews)
  reviews <- gsub("gene ", "gene-", reviews)
  reviews <- gsub(" phase", "-phase", reviews)
  reviews <- gsub(" binding", "-binding", reviews)
  reviews <- gsub("cellular ", "cellular-", reviews)
  reviews <- gsub(" signaling", "-signaling", reviews)
  reviews <- gsub("histone ", "histone-", reviews)
  reviews <- gsub(" region", "-region", reviews)
  reviews <- gsub(" involved", "-involved", reviews)
  reviews <- gsub(" part", "-part", reviews)
  reviews <- gsub(" acid", "-acid", reviews)
  reviews <- gsub(" factor", "-factor", reviews)
  reviews <- gsub(" compound", "-compound", reviews)
  reviews <- gsub(" to ", "-to-", reviews)
  reviews <- gsub(" and ", "-and-", reviews)
  reviews <- gsub("process", "", reviews)
  reviews <- gsub("biological", "", reviews)
  reviews <- gsub( "organism", "", reviews)
  reviews <- gsub("process", "", reviews)
  reviews <- gsub("process", "", reviews)
  reviews <- gsub("'", "", reviews)  # remove apostrophes
  reviews <- gsub("[[:punct:]]", " ", reviews)  # replace punctuation with space
  reviews <- gsub("[[:cntrl:]]", " ", reviews)  # replace control characters with space
  reviews <- gsub("^[[:space:]]+", "", reviews) # remove whitespace at beginning of documents
  reviews <- gsub("[[:space:]]+$", "", reviews) # remove whitespace at end of documents
  reviews <- tolower(reviews)
  
  doc.list <- strsplit(reviews, "[[:space:]]+")
  
  doc.list[grep("^ion",doc.list)] <- "*ion*"
  doc.list[grep("-ion-",doc.list)] <- "*ion*"
  doc.list[grep("-ion",doc.list)] <- "*ion*"
  doc.list[grep("^ribosom",doc.list)] <- "ribosom*"
  doc.list[grep("^epithel",doc.list)] <- "epithel*"
  doc.list[grep("^neur",doc.list)] <- "neur*"
  doc.list[grep("^dna",doc.list)] <- "dna*"
  doc.list[grep("^skeletal",doc.list)] <- "skeletal*"
  doc.list[grep("^structur",doc.list)] <- "structur*"
  doc.list[grep("^transcription",doc.list)] <- "transcription*"
  doc.list[grep("^transmembrane-transport",doc.list)] <- "transmembrane-transport*"
  
  term.table <- table(unlist(doc.list))
  term.table <- sort(term.table, decreasing = TRUE)
  
  vocab <- names(term.table)
  
  get.terms <- function(x) {
    index <- match(x, vocab)
    index <- index[!is.na(index)]
    rbind(as.integer(index - 1), as.integer(rep(1, length(index))))
  }
  documents <- lapply(doc.list, get.terms)
  
  
  
  
  D <- length(documents)  # number of documents (2,000)
  W <- length(vocab)  # number of terms in the vocab (14,568)
  doc.length <- sapply(documents, function(x) sum(x[2, ]))  # number of tokens per document [312, 288, 170, 436, 291, ...]
  N <- sum(doc.length)  # total number of tokens in the data (546,827)
  term.frequency <- as.integer(term.table)
  
  K <- 3
  G <- 5000
  alpha <- 0.02
  eta <- 0.02
  
  # Fit the model:
  library(lda)
  set.seed(357)
  t1 <- Sys.time()
  fit <- lda.collapsed.gibbs.sampler(documents = documents, K = K, vocab = vocab, 
                                     num.iterations = G, alpha = alpha, 
                                     eta = eta, initial = NULL, burnin = 0,
                                     compute.log.likelihood = TRUE)
  t2 <- Sys.time()
  t2 - t1
  fit$document_sums
  str(fit)
  fit$topics
  fit$topic_sums
  theta <- t(apply(fit$document_sums + alpha, 2, function(x) x/sum(x)))
  phi <- t(apply(t(fit$topics) + eta, 2, function(x) x/sum(x)))
  
  MovieReviews <- list(phi = phi,
                       theta = theta,
                       doc.length = doc.length,
                       vocab = vocab,
                       term.frequency = term.frequency)
  library(LDAvis)
  
  # create the JSON object to feed the visualization:
  json <- createJSON(phi = MovieReviews$phi, 
                     theta = MovieReviews$theta, 
                     doc.length = MovieReviews$doc.length, 
                     vocab = MovieReviews$vocab, 
                     term.frequency = MovieReviews$term.frequency)
  install.packages("servr")
  library(servr)  
  serVis(json, out.dir = 'vis', open.browser = T)
}


SymbolToEntrez <- function(list, handle=NULL){
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
  results <- getBM(attributes = c( 'entrezgene'), filters = 'hgnc_symbol', curl = handle,
                   values = list, mart = mart)
  return(results$entrezgene)
}

if(!file.exists(file.path(datapath,"AllProbeDict.txt"))){
  AllProbeDict <- as.data.frame(rbind(cbind(as.character(chipAnnotations$Probe.Set.ID), as.character(chipAnnotations$Gene)),  cbind(as.character(fData(GSE36451)$ID), as.character(fData(GSE36451)$GENE_SYMBOL)) ), stringsAsFactors=F)
  AllProbeDict <-rbind(AllProbeDict, cbind(union(rownames(TPMEAll),rownames(TPMBAll)),union(rownames(TPMEAll),rownames(TPMBAll))))
  colnames(AllProbeDict) <- c("Probe", "Gene")
  for(i in 1:nrow(AllProbeDict)){
    if(AllProbeDict[i,]$Gene == ""| AllProbeDict[i,]$Gene =="no hit" | AllProbeDict[i,]$Gene =="N/A"| AllProbeDict[i,]$Gene =="Err:512"){
      AllProbeDict[i,]$Gene = AllProbeDict[i,]$Probe
    }
  }
  
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
  symbol.entrez <- getBM(attributes = c( 'hgnc_symbol','entrezgene'), filters = 'hgnc_symbol', curl = handle,
                         values = AllProbeDict$Gene, mart = mart)
  AllProbeDict <-  cbind(AllProbeDict, as.character(symbol.entrez[match(AllProbeDict$Gene, symbol.entrez$hgnc_symbol),]$entrezgene))
  colnames(AllProbeDict) <- c("Probe", "Gene","Entrez")
  rownames(AllProbeDict) <- AllProbeDict$Probe
  unmapped <- as.data.frame(cbind( paste0(rownames(NormTPME.unmapped.contigs)),paste0(rownames(NormTPME.unmapped.contigs)),rep(NA, length(NormTPME.unmapped.contigs[,1]))))
  colnames(unmapped) <- c("Probe", "Gene","Entrez")
  rownames(unmapped) <- unmapped$Probe
  AllProbeDict <-  rbind(AllProbeDict,unmapped)
  rownames(AllProbeDict) <- AllProbeDict$Probe
  
  write.table(AllProbeDict ,file.path(datapath,"AllProbeDict.txt"), row.names = T,col.names = T, quote = F)
}else{
  AllProbeDict <- read.table(file.path(datapath,"AllProbeDict.txt"), header = T )
}

