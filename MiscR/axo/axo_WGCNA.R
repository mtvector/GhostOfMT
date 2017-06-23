source("http://bioconductor.org/biocLite.R")
#biocLite("WGCNA")
library(WGCNA)
enableWGCNAThreads(32)
source("~/code/axo/LoadDatasets.R")

powerValues <- lapply(experiments, function(x){
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(t(x), powerVector = powers, verbose = 5)
})
# ?pickSoftThreshold
# sizeGrWindow(9, 5)
# par(mfrow = c(1,1))
# cex1 = .9
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red")
# 
# cutreeDynamic()
# 
# minModuleSize = 30;
# # Module identification using dynamic tree cut:
# dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
#                             deepSplit = 2, pamRespectsDendro = FALSE,
#                             minClusterSize = minModuleSize);
# table(dynamicMods)
# 

for(es in 1:length(experiments)){
  softPower = powerValues[[es]]$powerEstimate
  adjacency = adjacency(t(experiments[[es]]), power = softPower)
  TOM = TOMsimilarity(adjacency)
  dissTOM = 1-TOM
  geneTree = hclust(as.dist(dissTOM), method = "average")
  saveRDS( geneTree, file =file.path(datapath,"hClustering",paste0("TOMclust",names(experiments)[es])))
  saveRDS( TOM, file =file.path(datapath,"hClustering",paste0("TOM",names(experiments)[es])))
}
for(es in 1:length(experiments)){
  dist <- dist(experiments[[es]])
  hc <-  hclust(dist)
  saveRDS( dist, file =file.path(datapath,"hClustering",paste0("dist",names(experiments)[es])))
  saveRDS( hc, file =file.path(datapath,"hClustering",paste0("hc",names(experiments)[es])))
}
