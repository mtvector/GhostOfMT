source("http://bioconductor.org/biocLite.R")
#biocLite("WGCNA")
source("~/code/axo/LoadDatasets.R")
source("~/code/axo/DE_utils.R")
datapath <- "~/code/data/axo"
library(Mfuzz)
library(WGCNA)
enableWGCNAThreads(32)
#biocLite("clusterProfiler")
library(clusterProfiler)
#biocLite("RDAVIDWebService")
library(RDAVIDWebService)
#biocLite("DAVIDQuery")
library(DAVIDQuery)
library(org.Hs.eg.db)
library(biomaRt)
library(doParallel)
registerDoParallel()
library(parallel)
options("mc.cores"=32)
clustpath <- file.path(datapath, "hClustering","nofilter_blk")

setwd(clustpath)
blokfile <- list.files("./", pattern = "block")
wgcna.outputs <-  lapply(blokfile ,readRDS)
names(wgcna.outputs) <- blokfile

# Convert labels to colors for plotting

names(experiments) <-  lapply(names(experiments), function(x){
  paste0("block",x)
})

for(x in names(experiments)){
  matt <- experiments[x][[1]]
  idx <- rowMads(matt)>median(rowMads(matt))
  wgcna.outputs[x][[1]]$exprs <- matt[idx,]
}

names(experiments) <-  lapply(names(experiments), function(x){
  gsub("block","",x)
})

for(x in 1:length(wgcna.outputs)){
  wgcna.outputs[x][[1]]$clusters <- lapply( unique(wgcna.outputs[x][[1]]$colors) , function(c){
    ind = wgcna.outputs[x][[1]]$colors == c
    rownames(wgcna.outputs[x][[1]]$exprs)[ind]
  })
  names(wgcna.outputs[x][[1]]$clusters) <- unique(wgcna.outputs[x][[1]]$colors)
}

for(x in 1:length(wgcna.outputs)){
  wgcna.outputs[x][[1]]$clusterGO <- lapply( unique(wgcna.outputs[x][[1]]$colors), FUN= function(c){
    print(c)
    enrichGO( AllProbeDict[ wgcna.outputs[x][[1]]$clusters[c][[1]], ]$Entrez, ont="BP", universe = unique(AllProbeDict$Entrez))
  })
  names(wgcna.outputs[x][[1]]$clusterGO) <- unique(wgcna.outputs[x][[1]]$colors)
}


foreach(x = 1:length(wgcna.outputs))%dopar%{
    wgcna.outputs[x][[1]]$clusterDAVID <- lapply( unique(wgcna.outputs[x][[1]]$colors) , function(c){
      tryCatch({
      DAVIDQuery( AllProbeDict[ wgcna.outputs[x][[1]]$clusters[c][[1]], ]$Entrez,tool="chartReport", type="ENTREZ_GENE_ID",annot ="GOTERM_BP_ALL",details = T)$DAVIDQueryResult
      }, error=function(cond) {
        message("object 'secondStageResult' not found")
        return(NA)
      })
    })
    names(wgcna.outputs[x][[1]]$clusterDAVID) <- unique(wgcna.outputs[x][[1]]$colors)
}

lapply(names(wgcna.outputs), function(x){
  saveRDS(wgcna.outputs[[x]], file.path(clustpath, paste0(x)))
})
