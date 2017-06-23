library(EACI)
library(clusterProfiler)
source("http://bioconductor.org/biocLite.R")
source("~/code/general/wgcna_functions.R")
#biocLite("WGCNA")
library(WGCNA)
options("mc.cores"=2)
enableWGCNAThreads(2)
library(doParallel)
library(mclust)
registerDoParallel()
options(stringsAsFactors = FALSE)


datapath <- "~/code/data/cb/shiftFiles/"


euclidianDendroModules <- function(set, dpsplit=T,simplecut=F,k = 20){
  Q <- list()
  d <- dist(set)
  tree <- hclust(d)
  maxht <-  max(tree$height)/1.1
  if(simplecut){
    modules <- cutree(tree ,k=20)
  }else{
    modules <- cutreeDynamicTree(tree,maxTreeHeight = maxht,deepSplit = dpsplit)
  }
  Q$colors <-labels2colors( modules,zeroIsGrey = T)
  names(Q$colors) <- rownames(set)
  Q$dendrograms <- tree
  Q$MEs <-sapply( unique(Q$colors) , function(i){
    ind <- names(Q$colors)[Q$colors == i]
    getME(set[ind,])
  } )
  return(Q)
}


goterms <- lapply(unique(labels),function(c){
  ind <- c == labels
  ind <- as.numeric(ind)
  names(ind) <- names(labels)
  eacitest(ind,"org.Hs.eg","SYMBOL",sets = "GO",minsetsize = 30)$setscores
})

dtw.clustering <-  function(output.dir,gene.list=neural_list,nameofclust=nameofclust,sc=F,k=10){
  load(file.path(output.dir, "fits.RData"))
  gene.list <- gene.list[(gene.list%in%names(fits[[1]]))]
  exprsmat <- t(sapply(fits[[1]], function(q) q$mu$y))
  shifts <- read.table(file.path(output.dir,"shift.txt"),row.names = 1,header = T)
  colnames(shifts) <- gsub("V","",colnames(shifts))
  clust <- euclidianDendroModules(shifts,dpsplit = F,simplecut = sc,k = k)
  clust$shifts <- shifts
  clust$exprs <- exprsmat
  clust$params <- list( maxcutheight="max/1.5" , corFnc = "euclidian", deepSplit = T)
  clust$name <- nameofclust
  clust$clusters <- lapply( unique(clust$colors) , function(c){
    ind = clust$colors == c
    rownames(clust$shifts)[ind]
  })
  names(clust$clusters) <- unique(clust$colors)
  plotDendroAndColors(clust$dendrograms ,clust$colors, dendroLabels = F,main="Clusters and colors",autoColorHeight = T)
  #allgenes <- rownames(exprsmat)
  modulesGO <- dtwModGO(clust,all.genes)
  clust$clusterGO <- modulesGO 
  return(clust)
  }

ods <- c( "~/code/data/cb/shiftFiles/MousevsHInVitro2cvsplinerj7","~/code/data/cb/shiftFiles/MixvsHInVitro2cvsplinerj7", "~/code/data/cb/shiftFiles/MousevsMixInVitro2cvsplinerj7")
#DTW.Mods <- lapply(ods, dtw.clustering)
c1 <- dtw.clustering(output.dir = ods[1],nameofclust = "MousevsHInVitro")
c2 <- dtw.clustering(output.dir = ods[2],nameofclust = "MixvsHInVitro")
c3 <- dtw.clustering(output.dir = ods[3],nameofclust = "MousevsMixInVitro")
cs1 <- dtw.clustering(output.dir = ods[1],nameofclust = "MousevsHInVitro",sc = T,k=25)
cs2 <- dtw.clustering(output.dir = ods[2],nameofclust = "MousevsHInVitro",sc = T,k=25)
cs3 <- dtw.clustering(output.dir = ods[3],nameofclust = "MousevsHInVitro",sc = T,k=25)
save(cs1,cs2,cs3,c1,c2,c3,file ="~/code/clusters.RData")

dtwModGO <- function(clust, allgenes){
   matMembership <- sapply(clust$clusters,function(x){
    sapply(allgenes,function(n){
      if(n%in%x){
        1.0
      }else{
        0.0
      }
    })
  })
  return(apply(matMembership,2,function(q){eacitest(q,"org.Hs.eg","SYMBOL",sets = "GO",minsetsize = 30)$setscores}))
}

makePPT.shift(c1,"MousevsH")
makePPT.shift(c2,"MixvsH")
makePPT.shift(c3,"MousevsMix")
makePPT.shift(cs1,"MousevsHBlunt")
makePPT.shift(cs2,"MixvsHBlunt")
makePPT.shift(cs3,"MousevsMixBlunt")
load("~/code/data/cb/clusters.RData")
makePPT.shift <-  function(clustas,name){
  #Load all the datasets
  library(ReporteRs)
  mydoc = pptx()
  mydoc = addSlide( mydoc, "Title and Content" )
  mydoc = addTitle(mydoc, name)
#   paramstring = pot(paste(clustas$params, collapse = " ") )
#   mydoc = addParagraph( mydoc, value = paramstring)
  mydoc = addPlot(mydoc, function(){plotDendroAndColors(clustas$dendrograms ,clustas$colors, dendroLabels = F,addGuide = T, guideAll = T,main="Clusters and colors",autoColorHeight = F)})
  apath <- file.path(datapath,paste0(name,".pptx"))
  mydoc = addSlide( mydoc, "Title and Content" )
  mydoc = addTitle(mydoc, "Avg Shift at each TP")
  mydoc = addPlot(mydoc, function(){
    matplot(sapply(clustas$clusters,function(i){ colMeans(clustas$shifts[i,])}), col = names(clustas$clusters),type = "l",ylab = "Avg Shift",xlab = "TP",xaxt="n") 
    legend(x = "topright",legend = names(clustas$clusters),fill=names(clustas$clusters),col=names(clustas$clusters),cex = .3)
    axis(1, at = 1:ncol(clustas[["shifts"]]), labels = gsub("X","Day ", colnames(clustas[["shifts"]])), cex.axis = 0.7 )
  },vector.graphic=F)
  mydoc = addSlide( mydoc, "Title and Content" )
  mydoc = addTitle(mydoc, "# Neural Genes in Each Cluster")
  tab <-  t(as.data.frame(colSums(sapply(clustas$clusters, function(x)neural_list %in% x ))))
  rownames(tab) <- "#"
  mydoc = addFlexTable(mydoc,flextable = FlexTable(tab))
  for(c in names(clustas[["clusters"]])){
    print(c)
    mydoc = addSlide(mydoc,"Title and Content")
    mydoc = addTitle(mydoc, paste("Shifts for",c))
    mydoc = addPlot(mydoc, function() {matplot( t(clustas$shifts[clustas[["clusters"]][[c]],]),type = "l",ylab = "Shifts",xlab = "TP",xaxt="n") 
      axis(1, at = 1:ncol(clustas[["shifts"]]), labels = gsub("X","Day ", colnames(clustas[["exprs"]])), cex.axis = 0.7 )
    }, vector.graphic = F)
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle(mydoc, paste("Gene Names",c))
    remove("nm")
    nm =  pot(paste(clustas[["clusters"]][[c]],collapse = " "))
    mydoc = addParagraph( mydoc, value = nm)
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle(mydoc, paste("GO Enrichment",c))
    mydoc = addFlexTable(mydoc,flextable = FlexTable(clustas$clusterGO[[c]][1:50,c("Term","set.size","pval")]))
    mydoc = addSlide( mydoc, "Two Content" )
    mydoc = addTitle( mydoc, c)
    mydoc = addPlot(mydoc, function() expressionHeatmap(cor(clustas[["exprs"]][clustas[["clusters"]][[c]],]), paste0(c ," Pearson\n Shift Cors"),clust = T), vector.graphic = F)
    mydoc = addPlot(mydoc, function() expressionHeatmap(cor(t(clustas[["exprs"]][clustas[["clusters"]][[c]],])), paste0(c ," Pearson\n Gene Expression Cors"),clust = T), vector.graphic = F)  
    
  }
  writeDoc( mydoc, apath )
}


ods <- c( "~/code/data/cb/shiftFiles/MousevsHInVitro2cvsplinerj7","~/code/data/cb/shiftFiles/MixvsHInVitro2cvsplinerj7", "~/code/data/cb/shiftFiles/MousevsMixInVitro2cvsplinerj7")
shiftSets <-  lapply(ods,function(x){
  shifts <- read.table(file.path(x,"shift.txt"),row.names = 1,header = T)
})
comS <- Reduce(intersect, sapply(shiftSets,rownames))
comT <- Reduce(intersect, sapply(shiftSets,colnames))
nl <- neural_list[neural_list%in% comS]
hist(unlist(shiftSets[[1]][comS,comT] - shiftSets[[2]][comS,comT] +shiftSets[[3]][comS,comT]))

a=as.matrix(shiftSets[[1]][comS,comT])
b=as.matrix(shiftSets[[2]][comS,comT])
c=as.matrix(shiftSets[[3]][comS,comT])

heatmap.2(a[nl,], F, F, trace = "none",col = cols)
heatmap.2(b[nl,], F, F, trace = "none",col=cols)
heatmap.2(c[nl,], F, F, trace = "none",col=cols)

combo <- Reduce(cbind, lapply(1:ncol(a), function(i){
  ttt=Reduce(cbind,list(a[,i],b[,i],c[,i]))
}))
colnames(combo) <- rep(colnames(a), each=3)

heatmap.2(combo[nl,], F, F, trace = "none",col = cols)

