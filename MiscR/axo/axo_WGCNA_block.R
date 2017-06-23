source("http://bioconductor.org/biocLite.R")
#biocLite("WGCNA")
library(WGCNA)
options("mc.cores"=64)
enableWGCNAThreads(64)
library(doParallel)
registerDoParallel()
options(stringsAsFactors = FALSE)
source("~/code/axo/LoadDatasets.R")
source("~/code/axo/clustering/cluster_analysis.R")

# dpsplt <- 3
# corfn <- c("bicor")
# sign <- c("unsigned", "signed hybrid")
# clustpaths <- vector(mode = "list")
# foreach(crfn = corfn)%dopar%{
# foreach(sgn = sign)%dopar%{
# clustpath <- file.path(datapath,"hClustering", paste0("MEAN",crfn ,gsub(" ", "",sgn), dpsplt,"SD"))
# clustpaths <<- list(clustpaths,clustpath)
# dir.create(clustpath)
# 
# powerValues <- lapply(experiments, function(x){
#   # Choose a set of soft-thresholding powers
#   powers = c(c(3:10), seq(from = 12, to=24, by=2))
#   idx <- rowSds(x)>median(rowSds(x))
#   x <- x[idx,]
#   # Call the network topology analysis function
#   sft = pickSoftThreshold(t(x), powerVector = powers,RsquaredCut=.85,removeFirst = T ,corFnc = crfn, networkType = sgn, verbose = 5)
# })
# 
# powerValues <- sapply(powerValues, function(i){
#   if(i$powerEstimate >2 | is.na(i$powerEstimate)){
#       i$powerEstimate <-  i$fitIndices$Power[which.max(i$fitIndices$SFT.R.sq)] 
#   }
#   i$powerEstimate
# })
# 
# foreach(es = 1:length(experiments), .errorhandling='pass', .export='errval')%dopar%{
# if(!is.null(powerValues[[es]]) & !is.na(powerValues[[es]])){
#   softPower = powerValues[[es]]
#   x <- experiments[[es]]
#   idx <- rowSds(x)>median(rowSds(x))
#   x <- x[idx,]
#   clust <- blockwiseModules(t(x), maxBlockSize = 10000, power=softPower, networkType = sgn, corFnc = crfn, maxPOutliers = 0.1 , minKMEtoStay = .3, TOMDenom = "mean", deepSplit = dpsplt,verbose = 5)
#   clust$exprs <- x
#   clust$clusters <- lapply( unique(clust$colors) , function(c){
#     ind = clust$colors == c
#     rownames(clust$exprs)[ind]
#   })
#   names(clust$clusters) <- unique(clust$colors)
#   saveRDS( clust, file =file.path(clustpath, paste0("block",names(experiments)[es])))
# }}
# }}

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

