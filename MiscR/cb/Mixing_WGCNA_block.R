source("http://bioconductor.org/biocLite.R")
#biocLite("WGCNA")
library(WGCNA)
options("mc.cores"=64)
enableWGCNAThreads(64)
library(doParallel)
registerDoParallel()
options(stringsAsFactors = FALSE)
source("http://bioconductor.org/biocLite.R")
library(ggplot2)
library(Biobase)
library(RColorBrewer)
library(gplots)
library(plyr)
library(EACI)
library(clusterProfiler)
library(matrixStats)
library(zoo)
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
#biocLite("maSigPro")
library(maSigPro)
library(edgeR)
library(EBSeq)
library(DESeq2)
library(devtools)
library(biomaRt)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(scales)
library(gridExtra)
#devtools::install_github("greenelab/TDM")
library(TDM)
library(reshape2)
source("~/code/cb/LoadFunctions.R")
source("~/code/cb/LoadData.R")



#' Lag Penalized Correlation
#'
#' This functions calculates a distance matrix of peptides i and j using penalized correlation
#'
#' @param data a matrix of time series data of peptides
#' @param C a numeric value of constant
#' @param maxlag a numeric value of the maximum lags to be included
#' @param timepoints a vector of time points of the data
#' @return a list of two matrix a matrix of lag selected and a
#' matrix of distance of peptide i and j
#' @author Thevaa Chandereng
#' @details
#' This function computes the distance between peptides using penalized correlation
#' @seealso \code{dist}
#' @export
#'

lag.penalized.corr <- function(data, C = NULL, max.lag = 2, timepoints = NULL){
  #check if all the inputs are right
  stopifnot(is.matrix(as.matrix(data)), max.lag < length(timepoints), is.numeric(max.lag),
            is.vector(timepoints), length(timepoints) > 1, dim(data)[2] == length(timepoints))
  #creating an empty matrix for dist and lag
  dist <- array(NA, c(dim(data)[1], dim(data)[1]))
  lag <- array(NA, c(dim(data)[1], dim(data)[1]))
  #making sure the dataset is a matrix
  data <- as.matrix(data)
  #all the different lags
  lags <- -max.lag:max.lag
  #lloping to store the dist data for (i > j)
  for(i in 1:dim(data)[1]){
    for(j in 1:dim(data)[1]){
      if(i > j){
        corr <- rep(NA, length(lags))
        for(k in max.lag:1){
          corr[max.lag - k + 1] <- exp(-C * (timepoints[k + 1]- timepoints[1])^2) *
            cor(data[i, 1:(length(timepoints) - k)], 
                data[j, (k + 1):length(timepoints)])
        }
        corr[max.lag + 1] <- cor(data[i, ], data[j, ])
        for(m in 1:max.lag){
          corr[m + max.lag + 1] <- exp(-C * (timepoints[m + 1]- timepoints[1])^2) *
            cor(data[j, 1:(length(timepoints) - m)], 
                data[i, (m + 1):length(timepoints)])
        }
        #picking the max distance and related lag
        dist[i, j] <- corr[which.max(abs(corr))]
        lag[i, j] <- lags[which.max(abs(corr))]
      }
    }
  }
  #making sure the distance is in the class of as.dist
  dist <- as.dist(dist)
  lag <- as.dist(lag)
  #returning the result
  return(list(lag = lag, distance = dist))
}


corDist <- lag.penalized.corr(datasets[[1]][neural_list[neural_list%in%rownames(datasets[[1]])],],max.lag = 6,C = 10,timepoints = tpDatasets[[1]])
corDist

condit.combos <- expand.grid(c("mm","hs"),unique(conditMixSets[[2]]),stringsAsFactors = F)
condit.combos <- condit.combos[-c(2,7),]
condit.combos <- condit.combos[c(1,6,2,3,4,5),]
condit.combos <- condit.combos[c(1,3,6),]
x <- unlist(condit.combos[6,])
set <- MixSets[[x[1]]][1:8000,conditMixSets[[x[1]]]==x[2]]
set <- apply(set[1:8000,],1, aggregate,list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), mean )
set <- t(Reduce(cbind,lapply(set,"[",2)))
rownames(set) <- rownames(MixSets[[x[1]]][1:8000,conditMixSets[[x[1]]]==x[2]])
set <- log(set+1,2)
set <- set[rowSds(set)>1,]
tp <- unique(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2]])
corMat <- lag.penalized.corr(set,C = 2,max.lag = 2,timepoints = tp)
corMat <- (as.matrix(corMat$distance)+1)/2
rownames(corMat) <- colnames(corMat) <- rownames(set)
power <- pickSoftThreshold.fromSimilarity(corMat)$powerEstimate
adj <-  adjacency.fromSimilarity(corMat,power = power)
TOM <- TOMsimilarity(adj,TOMType = "unsigned")
dissTOM <- as.dist(1-TOM)
geneTree = hclust(dissTOM, method = "average")
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = as.matrix(dissTOM),
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

clusters <- lapply( unique(dynamicColors) , function(c){
  ind = dynamicColors == c
  rownames(set)[ind]
})
names(clusters) <- unique(dynamicColors)

matplot(t(set[clusters$brown,]),type = "l")



apply(condit.combos,1,function(x){
  SegRegEvents(set1 = list(log(GetNormalizedMat( MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ], MedianNorm(MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]))+1,2)), tp1 = list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]),output.dir = file.path(datapath,"segregevents"), fn=paste0(x[1],x[2],"Events")) 
})
condit.combos



dpsplt <- 3
corfn <- c("bicor")
sign <- c("unsigned", "signed hybrid")
clustpaths <- vector(mode = "list")
foreach(crfn = corfn)%dopar%{
foreach(sgn = sign)%dopar%{
clustpath <- file.path(datapath,"hClustering", paste0("MEAN",crfn ,gsub(" ", "",sgn), dpsplt,"SD"))
clustpaths <<- list(clustpaths,clustpath)
dir.create(clustpath)

powerValues <- lapply(experiments, function(x){
  # Choose a set of soft-thresholding powers
  powers = c(c(3:10), seq(from = 12, to=24, by=2))
  idx <- rowSds(x)>median(rowSds(x))
  x <- x[idx,]
  # Call the network topology analysis function
  sft = pickSoftThreshold(t(x), powerVector = powers,RsquaredCut=.85,removeFirst = T ,corFnc = crfn, networkType = sgn, verbose = 5)
})

powerValues <- sapply(powerValues, function(i){
  if(i$powerEstimate >2 | is.na(i$powerEstimate)){
      i$powerEstimate <-  i$fitIndices$Power[which.max(i$fitIndices$SFT.R.sq)] 
  }
  i$powerEstimate
})

foreach(es = 1:length(experiments), .errorhandling='pass', .export='errval')%dopar%{
if(!is.null(powerValues[[es]]) & !is.na(powerValues[[es]])){
  softPower = powerValues[[es]]
  x <- experiments[[es]]
  idx <- rowSds(x)>median(rowSds(x))
  x <- x[idx,]
  clust <- blockwiseModules(t(x), maxBlockSize = 10000, power=softPower, networkType = sgn, corFnc = crfn, maxPOutliers = 0.1 , minKMEtoStay = .3, TOMDenom = "mean", deepSplit = dpsplt,verbose = 5)
  clust$exprs <- x
  clust$clusters <- lapply( unique(clust$colors) , function(c){
    ind = clust$colors == c
    rownames(clust$exprs)[ind]
  })
  names(clust$clusters) <- unique(clust$colors)
  saveRDS( clust, file =file.path(clustpath, paste0("block",names(experiments)[es])))
}}
}}

deepsplit <- c(3)
corfn <- c( "bicor")
sign <- c("unsigned","signed hybrid")
clustpaths <- vector(mode = "list")
for(dpsplt in deepsplit ){
  for(crfn in corfn){
    for(sgn in sign){
      clustpath <- file.path(datapath,"hClustering", paste0("MEAN",crfn ,gsub(" ", "",sgn), dpsplt, "SD"))
      clustpaths <- c(clustpaths,clustpath)
    }
  }
}

# makeReports(clustpaths)
# source("~/code/axo/clustering/cluster_meta.R")

intToBits(74)

nSets = 3
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("GSE35255","GSE67118","GSE37198")
shortLabels = c("Terrestrial", "Blastema","Nerve")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(experiments[[1]])))
names(multiExpr[[1]]$data) = rownames(experiments[[1]])
rownames(multiExpr[[1]]$data) = colnames(exprsGSE35255)
multiExpr[[2]] = list(data = as.data.frame(t(experiments[[2]])))
names(multiExpr[[2]]$data) = rownames(experiments[[2]])
rownames(multiExpr[[2]]$data) = colnames(exprsGSE67118)
multiExpr[[3]] = list(data = as.data.frame(t(experiments[[3]])))
names(multiExpr[[3]]$data) = rownames(experiments[[3]])
rownames(multiExpr[[3]]$data) = colnames(exprsGSE37198)

net$exprs <- multiExpr
net$clusters <- lapply( unique(net$colors) , function(c){
  ind = net$colors == c
  rownames(net$exprs)[ind]
})
names(net$clusters) <- unique(net$colors)

net = blockwiseConsensusModules(
  multiExpr, power = 6, minModuleSize = 30, deepSplit = 3,
  pamRespectsDendro = F, checkMissingData = F,
  networkType = "signed", corFnc = "bicor", maxPOutliers = 0.1 , minKMEtoStay = .2, TOMDenom = "mean",
  mergeCutHeight = 0.25, numericLabels = TRUE,
  maxBlockSize = 22000,
  saveTOMs = TRUE, verbose = 5)

saveRDS( net, file =file.path(datapath,"hClustering",paste0("consensusAmby002")))

nSets = 5
collapseExperiments <- list()
for(i in 1:nSets){
  collapseExperiments[[i]] <- collapseRows(experiments[[i]],AllProbeDict[rownames(experiments[[i]]),]$Gene, rownames(experiments[[i]]))  
}

for(i in 1:nSets){
  collapseExperiments[[i]] <- collapseExperiments[[i]]$datETcollapsed 
}

for(i in 1:nSets){
  rownames(collapseExperiments[[i]]) <- AllProbeDict[rownames(collapseExperiments[[i]]),]$Gene
}

RN <- list()
for(i in 1:nSets){
  RN[[i]] <- rownames(collapseExperiments[[i]])
}

RNcommon <- Reduce(intersect, RN)

for(i in 1:nSets){
  collapseExperiments[[i]] <- collapseExperiments[[i]][RNcommon,]
}

# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("GSE35255","GSE67118","GSE37198","EmbryoSeq","BlastemaSeq")
shortLabels = c("Terrestrial", "Blastema","Nerve","EmbryoSeq","BlastemaSeq")
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(collapseExperiments[[1]])))
names(multiExpr[[1]]$data) = rownames(collapseExperiments[[1]])
rownames(multiExpr[[1]]$data) = colnames(exprsGSE35255)
multiExpr[[2]] = list(data = as.data.frame(t(collapseExperiments[[2]])))
names(multiExpr[[2]]$data) = rownames(collapseExperiments[[2]])
rownames(multiExpr[[2]]$data) = colnames(exprsGSE67118)
multiExpr[[3]] = list(data = as.data.frame(t(collapseExperiments[[3]])))
names(multiExpr[[3]]$data) = rownames(collapseExperiments[[3]])
rownames(multiExpr[[3]]$data) = colnames(exprsGSE37198)
multiExpr[[4]] = list(data = as.data.frame(t(collapseExperiments[[4]])))
names(multiExpr[[4]]$data) = rownames(collapseExperiments[[4]])
rownames(multiExpr[[4]]$data) = colnames(collapseExperiments[[4]])
multiExpr[[5]] = list(data = as.data.frame(t(collapseExperiments[[5]])))
names(multiExpr[[5]]$data) = rownames(collapseExperiments[[5]])
rownames(multiExpr[[5]]$data) = colnames(collapseExperiments[[5]])

net = blockwiseConsensusModules(
  multiExpr, power = 6, minModuleSize = 30, deepSplit = 3,
  pamRespectsDendro = F, checkMissingData = F,
  networkType = "signed", corFnc = "bicor", maxPOutliers = 0.1 , minKMEtoStay = .2, TOMDenom = "mean",
  mergeCutHeight = 0.25, numericLabels = TRUE,
  maxBlockSize = 22000,
  saveTOMs = TRUE, verbose = 5)

saveRDS( net, file =file.path(datapath,"hClustering",paste0("consensusAmby002plusSeq")))

net$exprs <- multiExpr
   net$clusters <- lapply( unique(net$colors) , function(c){
     ind = net$colors == c
     rownames(net$exprs)[ind]
   })
   names(net$clusters) <- unique(net$colors)

saveRDS( net, file =file.path(datapath,"hClustering",paste0("consensusAmby002plusSeq")))

