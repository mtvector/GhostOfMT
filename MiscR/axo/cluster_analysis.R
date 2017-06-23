source("http://bioconductor.org/biocLite.R")
#biocLite("WGCNA")
source("~/code/axo/LoadDatasets.R")
source("~/code/axo/DE_utils.R")
source("~/code/axo/DE_MA.R")
datapath <- "~/code/data/axo"
library(Mfuzz)
library(WGCNA)
enableWGCNAThreads(32)
#biocLite("clusterProfiler")
library(clusterProfiler)
#biocLite("RDAVIDWebService")
#library(RDAVIDWebService)
#biocLite("DAVIDQuery")
#library(DAVIDQuery)
library(org.Hs.eg.db)
library(biomaRt)
library(gridExtra)
options("mc.cores"=32)
library(parallel)
library(doParallel)
registerDoParallel()
#library(EACI)
#library(allez)

#ALL IS THE LIST OF EXPERIMENT EXPRESSION SETS
each.cors <- function(all, genelist){
  c.l <-  lapply(all, function(i){
    nms <-  intersect( AllProbeDict[genelist, ]$Gene, AllProbeDict[rownames(i), ]$Gene) 
    ind <- AllProbeDict[rownames(i), ]$Gene %in% nms
    if(sum(ind)>1){
    M <- as.matrix(i[ind,])
    }else{
      M <- matrix(NA,0,0) 
    }
    #Get the average, absolute value of cor for each set
    if(dim(M)[1]>0){
      Q <- abs(cor(t(M)))
      diag(Q) <- NA
      return(list(mean(Q, na.rm=T), sum(ind) ))
    }else{list(0,0)}
  })
  return(c.l)
}

myplot <- function(plot, filename){
  pdf(filename)
  plot
  dev.off()
}

mygrid <- function(plot, filename){
  theme <- ttheme_default(
    core = list(fg_params=list(cex = .6)),
    colhead = list(fg_params=list(cex = 0.3)),
    rowhead = list(fg_params=list(cex = 0.3)))
  pdf(filename , width = 50)
  grid.newpage()
  g <- plot
  grid.draw(g)
  dev.off()
}

DE.Rep <- function(gnz, DElist){
  DE.Pres <-  sapply(names(DElist), function(l){
    gnz <- unique(AllProbeDict[gnz,]$Gene)
    up= gnz %in% DE.Universe[[l]]$upGenes
    down= gnz %in% DE.Universe[[l]]$downGenes
      list(sum(up),sum(down))
    })
}

getCorMat <- function(clustdirs){
  lapply(clustdirs, function(dirrr){
    setwd(dirrr)
    blokfile <- list.files("./", pattern = "block")
    wgcna.outputs <-  lapply(blokfile ,readRDS)
    names(wgcna.outputs) <- blokfile
    
    all.avg.cors <- lapply(wgcna.outputs, function(x){
      sapply(x$clusters, function(y){
        each.cors(experiments, y)
      })
    })
  })
}

getModMatches <- function(clustdirs){
  lapply(clustdirs, function(dirrr){
    setwd(dirrr)
    blokfile <- list.files("./", pattern = "block")
    wgcna.outputs <-  lapply(blokfile ,readRDS)
    names(wgcna.outputs) <- blokfile
    
    all.avg.cors <- lapply(wgcna.outputs, function(x){
      sapply(x$clusters, function(y){
        each.match(wgcna.outputs, y)
      })
    })
  })
}


makeDAVIDReports <-  function(clustdirs){
  for(dirrr in clustdirs){
  setwd(dirrr)
  blokfile <- list.files("./", pattern = "block")
  wgcna.outputs <-  lapply(blokfile ,readRDS)
  names(wgcna.outputs) <- blokfile
  for(x in 1:length(wgcna.outputs)){
    if(is.null(wgcna.outputs[x][[1]]$clusterDAVID)){
    wgcna.outputs[x][[1]]$clusterDAVID <- lapply( unique(wgcna.outputs[x][[1]]$colors) , function(c){
      try({
        dq <-  DAVIDQuery(AllProbeDict[ wgcna.outputs[x][[1]]$clusters[c][[1]], ]$Entrez, tool="chartReport", type="ENTREZ_GENE_ID",annot ="GOTERM_BP_ALL",details = T)$DAVIDQueryResult
        return(dq)
        }) 
    })
    names(wgcna.outputs[x][[1]]$clusterDAVID) <- unique(wgcna.outputs[x][[1]]$colors)
    }
  }
}

  lapply(names(wgcna.outputs), function(x){
    saveRDS(wgcna.outputs[[x]], file.path(dirrr, x))
  })

reportpath <- file.path( dirrr ,"Reports")
  dir.create(reportpath)
  lapply(names(wgcna.outputs), function(exp){
    print(exp)
      sapply(names(wgcna.outputs[[exp]][["clusters"]]), function(clust){
        print(clust)
        try({
        print(wgcna.outputs[[exp]][["clusterDAVID"]][[clust]])
        mygrid(tableGrob(wgcna.outputs[[exp]][["clusterDAVID"]][[clust]][1:35, c("Category", "Term", "Count","%", "PValue", "Benjamini", "FDR")], theme = theme), file.path(reportpath,exp,clust,"DAVIDanalysis.pdf"))
          })
          })
        })
    }

# mean(rowMads(NormTPME[grepl("HOX", rownames(NormTPME)),]) > median(rowMads(NormTPME)))
# 
# matplot(t(TPME[grepl("HOX", rownames(TPME)),]),type="l")
# 
# hist(rowSds(NormTPME), seq(0,7,.05))
# min(rowMads(TPME))
# median(rowSds(NormTPME))/3



# sapply(
# sapply(wgcna.outputs$blockEmbryoSeq$clusters, function(w){
#   sapply(w, function(x){
#     grepl("HOX",AllProbeDict[x,]$Gene )
#   })
# } )
# , sum)
# 
# 
# length(wgcna.outputs$blockEmbryoSeq$clusters$turquoise)
# 
# ec <- eacitest(score=v,lib="org.Hs.eg",idtype="SYMBOL",locallist=NULL)
# 
# 
# 
# v <- as.vector(cor( t(wgcna.outputs$blockEmbryoSeq$exprs[wgcna.outputs$blockEmbryoSeq$clusters$turquoise,]), unlist(wgcna.outputs$blockEmbryoSeq$MEs$turquoise$blockEmbryoSeq) ))
# names(v) <- rownames(cor( t(wgcna.outputs$blockEmbryoSeq$exprs[wgcna.outputs$blockEmbryoSeq$clusters$turquoise,]), unlist(wgcna.outputs$blockEmbryoSeq$MEs$turquoise$blockEmbryoSeq) ))


makeReports <-  function(clustdirs){
    for(dirrr in clustdirs){
    setwd(dirrr)
    blokfile <- list.files("./", pattern = "block")
    wgcna.outputs <-  lapply(blokfile ,readRDS)
    names(wgcna.outputs) <- blokfile
        for(x in 1:length(wgcna.outputs)){
          if(is.null(wgcna.outputs[x][[1]]$clusterGO)){
            print(x)
            wgcna.outputs[x][[1]]$clusterGO <- lapply( unique(wgcna.outputs[x][[1]]$colors), FUN= function(c){
            enrichGO( AllProbeDict[ wgcna.outputs[x][[1]]$clusters[c][[1]], ]$Entrez, ont="BP", universe = unique(AllProbeDict$Entrez))
          })
          names(wgcna.outputs[x][[1]]$clusterGO) <- unique(wgcna.outputs[x][[1]]$colors)
          }
        }
        
        lapply(names(wgcna.outputs), function(x){
          saveRDS(wgcna.outputs[[x]], file.path(dirrr, x))
        })
        
        new1 = T
        if(is.null(wgcna.outputs[[1]][["avg.cors"]])){
        new1 =T
        print("all.avg.cor")
         cors.sums<- lapply(wgcna.outputs, function(x){
          sapply(x$clusters, function(y){
            each.cors(experiments, y)
          })
        })
        
        #Turn list of matrix of lists into 2 lists of matrices
        all.avg.cors <- 
          lapply(cors.sums, function(i){
            matrix(
              mapply(as.matrix(i),FUN=function(j){
                unlist(j)[1]
              }), nrow = nrow(i), ncol = ncol(i), dimnames = dimnames(i))
          })
        all.matches <- 
          lapply(cors.sums, function(i){
            matrix(
              mapply(as.matrix(i),FUN=function(j){
                unlist(j)[2]
              }), nrow = nrow(i), ncol = ncol(i), dimnames = dimnames(i))
          })
        }
    sapply(1:length(wgcna.outputs), function(i){
      if(new1){wgcna.outputs[[i]]$avg.cors <<- all.avg.cors[[i]]}
      if(new1) {wgcna.outputs[[i]]$each.matches <<- all.matches[[i]]}
    })

    lapply(names(wgcna.outputs), function(x){
      saveRDS(wgcna.outputs[[x]], file.path(dirrr, x))
    })    
    
    reportpath <- file.path( dirrr ,"Reports")
    dir.create(reportpath)
    mclapply(names(wgcna.outputs), function(exp){
      print(exp)
      dir.create(file.path(reportpath,exp))
      #try({
      theme <- ttheme_default(
        core = list(fg_params=list(cex = .6)),
        colhead = list(fg_params=list(cex = 0.3)),
        rowhead = list(fg_params=list(cex = 0.3)))
      #myplot(matplot(wgcna.outputs[[exp]][["MEs"]], type=c("l"), col =gsub("ME","",  names(wgcna.outputs[[exp]][["MEs"]]) )),file.path(reportpath,exp,"MEs.pdf") )
      ###
      mygrid(tableGrob(round(wgcna.outputs[[exp]][["avg.cors"]], digits = 3), theme = theme), file.path(reportpath,exp,paste0("crosscortable.pdf") ))

      mygrid(tableGrob(round(wgcna.outputs[[exp]][["each.matches"]], digits = 3), theme = theme),file.path(reportpath,exp,paste0("cormatchtable.pdf") ) )

      ###
      sapply(names(wgcna.outputs[[exp]][["clusters"]]), function(clust){
        print(clust)
        dir.create(file.path(reportpath,exp,clust))
        #myplot(matplot(wgcna.outputs[[exp]][["MEs"]][[paste0("ME",clust)]], type=c("l")),file.path(reportpath,exp,clust,paste0( clust," ME.pdf") ) )
        ###
        myplot({barplot(DE.Rep(wgcna.outputs[[exp]][["clusters"]][[clust]],DE.Universe),
                                         main =paste("DE Representation \n",exp,clust),  ylab = "Representation count" ,col = c("red","blue"),cex.names = .5,las=3 )  
                                 legend("topleft", c("DE Up", "DE Down"  ) , pch = 1, horiz = F , text.width = 1,cex=.60,
                                        col = c("red","blue") )}, file.path(reportpath,exp,clust,"DEReport.pdf"))

        myplot(expressionHeatmap(cor(t(wgcna.outputs[[exp]][["exprs"]][wgcna.outputs[[exp]][["clusters"]][[clust]],])),paste0(clust ," Pearson\n Gene Cors")),file.path(reportpath,exp,clust,"cor.pdf") )
        ###
        mygrid(tableGrob(as.data.frame(summary(wgcna.outputs[[exp]][["clusterGO"]][[clust]]))[1:30, c("ID", "Description","GeneRatio", "BgRatio", "p.adjust")], theme = theme), file.path(reportpath,exp,clust,"GOanalysis.pdf"))
      }) #})
      graphics.off()
  })
}}


deepsplit <- c(3)
corfn <- c("cor","bicor")
sign <- c("unsigned","signed hybrid")
clustpaths <- vector(mode = "list")
for(dpsplt in deepsplit ){
  for(crfn in corfn){
    for(sgn in sign){
      clustpath <- file.path(datapath,"hClustering", paste0("MEAN",crfn ,gsub(" ", "",sgn), dpsplt))
      clustpaths <- c(clustpaths,clustpath)
      }
  }
}



# About .61 is the cor for probes blasting to the same gene
# by.gene.cors = sapply(unique(AllProbeDict$Gene), function(g){
#   p = AllProbeDict[AllProbeDict$Gene==g,]$Probe
#   ind = rownames(allAmby002Data) %in% p
#   if(sum(ind)>1){mean(cor(t(allAmby002Data[ind,])),na.rm = T)}
#   else{NA}
# })
# mean(by.gene.cors,na.rm = T)

# plotDendroAndColors(wgcna.outputs$blockEmbryoSeq$dendrograms[[2]], colors = labels2colors(wgcna.outputs$blockEmbryoSeq$colors), dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Gene dendrogram and module colors")

# names(DElist)[g %in% DElist[[i]]]
# 
# mean(cor(t(wgcna.outputs$blockEmbryoSeq$exprs[wgcna.outputs$blockEmbryoSeq$clusters$red,])))
# all.avg.cors$blockEmbryoSeq
# sapply(all.avg.cors, mean ,na.rm=T)
# colMeans( all.avg.cors$blockallAmby002, na.rm = T)
# mean(cor(t(wgcna.outputs$blockEmbryoSeq$exprs)))
# 
# sets <-  lapply(wgcna.outputs, function(x){
#   o = list(data = as.data.frame(t(x$exprs)))
#   names(o$data) <- rownames(x$exprs)
#   rownames(o$data) <- colnames(x$exprs)
#   o
# })
# 
# colors <- lapply(wgcna.outputs,function(x){
#   x$colors
# })
# 
# sets <- fixDataStructure(sets)
# fixDataStructure(sets)
# checkSets(sets)
# modulePreservation(sets, colors)
# expressionHeatmap(cor(t(wgcna.outputs$blockallAmby002$exprs[ wgcna.outputs$blockallAmby002$clusters$darkred,])))
# 
# enrichMap(wgcna.outputs$blockEmbryoSeq$clusterGO$red,15)
# summary(wgcna.outputs$blockEmbryoSeq$clusterGO$red)
# matplot(wgcna.outputs$blockEmbryoSeq$MEs, type=c("l"), col =gsub("ME","",  names(wgcna.outputs$blockEmbryoSeq$MEs) ))
# legend("right", gsub("ME","",names(wgcna.outputs$blockEmbryoSeq$MEs)) , pch = 1, horiz = F , text.width = 1,cex=.40,
#         col = gsub("ME","",names(wgcna.outputs$blockEmbryoSeq$MEs)) )
# 
# 
# 
# matplot(wgcna.outputs$blockallAmby002$MEs, type=c("l"), col =gsub("ME","",  names(wgcna.outputs$blockallAmby002$MEs) ))
# 
# ?matplot()
# 
# 
# 
# s= NormTPME
# E = s[rowMads(s)>median(rowMads(s)),]
# mergedColors = labels2colors(x$colors)
# # Plot the dendrogram and the module colors underneath
# table(x$colors)
# plotDendroAndColors(x$dendrograms[[1]], mergedColors[x$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# 
# ind = x$colors =="midnightblue"
# matplot(x$MEs[,1], type=c("l"), col = 1:20)
# matplot(rowMeans(t(NormTPME[ind,])), type=c("l"), col = 1:20)
# matplot(t(NormTPME[ind,]), type=c("l"), col = 1:20)
# expressionHeatmap(cor(t(NormTPME[ind,])))
# 
# 
# 
# 
# 
# names(wgcna.Entrez) <- unique(x$colors)
# 
# 
# std.E <- standardise(ExpressionSet(E))
# #optim.E <- cselection(std.E,crange=seq(32,44,4),m=1.25, visu = T)
# cl.Embryo <- mfuzz(std.E,c=32, m=1.25)
# mfuzz.plot(std.E,cl=cl.Embryo,mfrow=c(4,4),time.labels=colnames(ECEtp))
# 
# cl.Embryo.GO <-  sapply(1:32, function(i){
#   enrichGO( SymbolToEntrez( names(cl.Embryo$cluster)[cl.Embryo$cluster==i]),ont = "BP",readable=T)
# })
# 
# cl.Embryo.entrez <-  sapply(1:32, function(i){
#  SymbolToEntrez( names(cl.Embryo$cluster)[cl.Embryo$cluster==i])
# })
# 
# 
# 
# 
# cl.Embryo.GO <-  sapply(1:32, function(i){
#   enrichGO( as.character(unlist(cl.Embryo.entrez[[i]])),ont = "BP",readable=T)
# })
# enrichMap(cl.Embryo.GO[[3]], n=1)
# summary(cl.Embryo.GO[[8]])[1:10,]
# ?enrichment()
# wgcna.Embryo.entrez <- 
#   
# 
# cl.Embryo.DAVID <-  sapply(1:32, function(i){
#   DAVIDQuery( SymbolToEntrez(names(cl.Embryo$cluster)[cl.Embryo$cluster==i]),tool="chartReport", type="ENTREZ_GENE_ID",annot ="GOTERM_BP_ALL",details = T)$DAVIDQueryResult
# })
# 
# wgcna.Embryo.GO <-  sapply(unique(x$colors), function(i){
#   ind = x$colors == i
#   enrichGO( SymbolToEntrez(names(x$colors)[ind]) ,ont = "BP",readable=T)
# })
# 
# wgcna.Embryo.DAVID <-  sapply(unique(x$colors), function(i){
#   ind = x$colors == i
#   DAVIDQuery(SymbolToEntrez(names(x$colors)[ind]),tool="chartReport", type="ENTREZ_GENE_ID",annot ="GOTERM_BP_ALL",details = T)$DAVIDQueryResult  
# })
# 
# enrichMap(y, n=10)
# ?enrichMap
# y=enrichGO( SymbolToEntrez( names(cl.Embryo$cluster)[cl.Embryo$cluster==2]),ont = "BP",readable=T)
# head(summary(y))
# z=enrichKEGG(SymbolToEntrez( names(cl.Embryo$cluster)[cl.Embryo$cluster==15]))
# head(summary(z))
# 
# a
# DAVIDQuery( as.character(names(cl.Embryo$cluster)[cl.Embryo$cluster==15]) ,tool="chartReport",URLlengthLimit = 40000, type="GENE_SYMBOL",annot ="GOTERM_BP_ALL")$DAVIDQueryResult
# 
# 
# ?"DAVIDQuery"
# NormTPME[ind,]
# SymbolToEntrez(c("HPS5","TRIO"))
# Amby002gene(rownames(allAmby002Data))
# rownames(NormTPME)
# c(Amby002gene(rownames(allAmby002Data)), rownames(NormTPME), rownames(NormECB) )
# 
# 
# 
# symbols = union( union(as.character(Amby002gene(rownames(allAmby002Data))), AgilentAmbyGene(rownames(exprsGSE36451))) , union(rownames(NormTPME), rownames(NormECB)) )
# names(symbols) <- symbols
# SymbolEntrez <-  sapply( symbols, USE.NAMES = T , SymbolToEntrez)
# 
# 
# ?lapply
# ?"DAVIDQuery"
# ECEtpmat <-  matrix(ECEtp,dim(x$MEs)[1],dim(x$MEs)[2])
# Data <- list(y=x$MEs, x=ECEtpmat)
# ?lasso
# data(chin07)
# dim(chin07$cn)
# chin07$cn
# ?lasso
# ?cv
# ?rowWeightedMads()
# ?rowMads
# hist(rowWeightedMads(allAmby002Data),breaks= seq(0,9,.2))
# 
# clusts$blockGSE67118$MEsOK
# 
# newGSE67118 <-  recutBlockwiseTrees(t(exprsGSE67118), clusts$blockGSE67118$goodSamples, clusts$blockGSE67118$goodGenes, clusts$blockGSE67118$blocks, file.path(datapath,"hClustering","TOMGSE67118"), clusts$blockGSE37198$dendrograms , minModuleSize = 50 )
# 
# 
# 
# dissTOM = readRDS(file.path(datapath, "hClustering", "TOMGSE67118"))
# geneTree = readRDS(file.path(datapath, "hClustering", "TOMclustGSE67118"))
# 
# minModuleSize = 60
# # Module identification using dynamic tree cut:
# dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
#                             deepSplit = 2, pamRespectsDendro = FALSE,
#                             minClusterSize = minModuleSize)
# table(dynamicMods)
# 
# dynamicColors = labels2colors(dynamicMods)
# table(dynamicColors)
# # Plot the dendrogram and colors underneath
# sizeGrWindow(8,6)
# plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Gene dendrogram and module colors")
# 
# 
# 
# 
# 
# 
# 
# 
# #Consensus
# consMEs = net$multiMEs;
# moduleLabels = net$colors;
# # Convert the numeric labels to color labels
# moduleColors = labels2colors(moduleLabels)
# consTree = net$dendrograms[[1]];
# #A quick way to take a look at the results is to plot the gene dendrogram and the corresponding module colors:
#   sizeGrWindow(8,6);
# plotDendroAndColors(consTree, moduleColors,
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Consensus gene dendrogram and module colors")
# 
# ?plotClusterTreeSamples
# ?gzfile()
# ?TOMsimilarity()
# 
# # ?mfrow
# # plot(hcd, main="Main")
# # plot(cut(hcd, h=75)$upper, 
# #      main="Upper tree of cut at h=75")
# # plot(cut(hcd, h=75)$lower[[2]], 
# #      main="Second branch of lower tree with cut at h=75")
# # plot(cutree(dendE,h=100), xlab="", sub="", main = "Gene clustering on Euclidian dissimilarity",
# #      labels = FALSE, hang = 0.04)
# # 
# # str(dendE, max = 10000, last.str =  "'")
# # 
# # plot(dendE,  hang=-1)
# # plot(as.phylo(hcE), type = "fan")
# # ?cutree
# # memb <- cutree(hcE, k = 400)
