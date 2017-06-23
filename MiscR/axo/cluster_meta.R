# source("http://bioconductor.org/biocLite.R")
# #biocLite("WGCNA")
# #source("~/code/axo/LoadDatasets.R")
source("~/code/axo/DE_utils.R")
# #source("~/code/axo/DE_MA.R")
# datapath <- "~/code/data/axo"
# library(Mfuzz)
# library(WGCNA)
# enableWGCNAThreads(32)
# #biocLite("clusterProfiler")
# library(clusterProfiler)
# #biocLite("RDAVIDWebService")
# #library(RDAVIDWebService)
#biocLite("DAVIDQuery")
# library(DAVIDQuery)
# library(org.Hs.eg.db)
# library(biomaRt)
# library(gridExtra)
options("mc.cores"=64)
# library(parallel)
# library(doParallel)
# registerDoParallel()

# deepsplit <- c(3)
# corfn <- c("cor","bicor")
# sign <- c("unsigned","signed hybrid")
# clustpaths <- vector(mode = "list")
# for(dpsplt in deepsplit ){
#   for(crfn in corfn){
#     for(sgn in sign){
#       clustpath <- file.path(datapath,"hClustering", paste0("MEAN",crfn ,gsub(" ", "",sgn), dpsplt))
#       clustpaths <- c(clustpaths,clustpath)
#     }
#   }
# }

addFullExpression <- function(clustdirs){
  for(dirrr in clustdirs){
    setwd(dirrr)
    blokfile <- list.files("./", pattern = "block")
    wgcna.outputs <-  lapply(blokfile ,readRDS)
    names(wgcna.outputs) <- blokfile
    lapply(names(wgcna.outputs), function(set){
      print(set)
      wgcna.outputs[[set]][["full.exprs"]] <<- experiments[[gsub("block","",set)]]
    })
    lapply(names(wgcna.outputs), function(x){
      saveRDS(wgcna.outputs[[x]], file.path(dirrr, x))
    }) 
  }
}

addFullExpression(clustpaths)

addDE <-  function(clustdirs){
  for(dirrr in clustdirs){
    setwd(dirrr)
    blokfile <- list.files("./", pattern = "block")
    wgcna.outputs <-  lapply(blokfile ,readRDS)
    names(wgcna.outputs) <- blokfile
    lapply(names(wgcna.outputs), function(set){
        print(set)
        lapply(names(wgcna.outputs[[set]][["clusters"]]), function(clust){
          #if(is.null(wgcna.outputs[[set]][["DE.Profile"]])){
          wgcna.outputs[[set]][["DE.Profile"]][[clust]] <<-  DE.Rep(wgcna.outputs[[set]][["clusters"]][[clust]],DE.Universe)
          #str(wgcna.outputs[[set]][["DE.Profile"]][[clust]])
          #}
        })
    })
  lapply(names(wgcna.outputs), function(x){
    saveRDS(wgcna.outputs[[x]], file.path(dirrr, x))
  }) 
  }
}
addDE(clustpaths)

getME <- function(exprset, genelist){
  nms <-  intersect( AllProbeDict[genelist, ]$Gene, AllProbeDict[rownames(exprset), ]$Gene) 
  ind <- AllProbeDict[rownames(exprset), ]$Gene %in% nms
  M <- as.matrix(exprset[ind,])
  if(sum(ind)>1){
  s <- svd((M- rowMeans(M))/rowSds((M)))
  r <- as.vector(as.vector(s$u[1,])%*% diag(c(s$d[1],rep(0,length(s$d)-1))) %*% t(s$v))
  #print(cor( colMeans((M- rowMeans(M))/rowSds((M))),r ))
  #print((s$d[1:4])^2/sum((s$d)^2))
  if(cor( colMeans((M- rowMeans(M))/rowSds((M))),r ) < 0 ){
    r <- -r
  }
  names(r) <- colnames(exprset)
  }else{r <- list(rep(0.0,ncol(M)))}
  return(r)
}

addME <- function(clustdirs){
  for(dirrr in clustdirs){
    setwd(dirrr)
    blokfile <- list.files("./", pattern = "block")
    wgcna.outputs <-  lapply(blokfile ,readRDS)
    names(wgcna.outputs) <- blokfile
    lapply(names(wgcna.outputs), function(exp){
      wgcna.outputs[[exp]][["MEs"]] <<- NULL
    })
    lapply(names(wgcna.outputs), function(exp){
      print(exp)
      sapply(names(wgcna.outputs[[exp]][["clusters"]]), function(clust){
        print(clust)
        sapply(names(wgcna.outputs), function(set){
          print(set)
          wgcna.outputs[[exp]][["MEs"]][[clust]][[set]] <<-  getME(wgcna.outputs[[set]][["full.exprs"]],wgcna.outputs[[exp]][["clusters"]][[clust]])
        })
      })
    })
    lapply(names(wgcna.outputs), function(x){
      saveRDS(wgcna.outputs[[x]], file.path(dirrr, x))
    }) 
  }
}
addME(clustpaths[[1]])

findModules <-  function(clustdirs, colMeanCor, numGO, set = NULL, DEcond=NULL, updown=NULL, DElim=NULL){
  dir.create(file.path(datapath,"checkout"))
  eg <- lapply( clustdirs, function(dirrr){
    setwd(dirrr)
    print(dirrr)
    blokfile <- list.files("./", pattern = "block")
    wgcna.outputs <-  lapply(blokfile ,readRDS)
    names(wgcna.outputs) <- blokfile
    if(is.null(set)){
      mgenes <-  sapply(names(wgcna.outputs), function(set){
        try(findSpecificModules(dirrr,wgcna.outputs, colMeanCor, numGO, set, DEcond, updown, DElim))
      })
    }else{
      mgenes <-  try(findSpecificModules(dirrr,wgcna.outputs, colMeanCor, numGO, set, DEcond, updown, DElim))
    }
  })
  if(is.null(set)){
    ll <<- list()
    for(x in eg){
      lapply(names(x), function(n){
        ll[[n]] <<- c(ll[[n]],x[[n]])
      } )
    }
    jj <<-  list()
    for(x in ll){
      lapply(names(x), function(n){
        jj[[n]] <<- cbind(jj[[n]],x[[n]])
      } )
    }
    return(jj)
  }else{
    return(unlist(eg, recursive = F))
  }
}

findSpecificModules<-  function(dirrr,wgcna.outputs, colMeanCor, numGO = 0 , set=NULL, DEcond=NULL, updown=NULL, DElim=NULL){
  #     moduleGenes <<-  lapply(names(wgcna.outputs),function(x){
  #       matrix(nrow = ncol(wgcna.outputs[[x]][["exprs"]]), ncol = 0)
  #     })
  moduleGenes <<- list()
  #names(moduleGenes) <<- names(wgcna.outputs)
  for(c in names(wgcna.outputs[[set]][["clusters"]])){
    sumGO <- summary(wgcna.outputs[[set]][["clusterGO"]][[c]])
    cmc <- mean(wgcna.outputs[[set]][["avg.cors"]][,c], na.rm=T )
    deq = T
    if(!is.null(DEcond)){
      deq <- wgcna.outputs[[set]][["DE.Profile"]][[c]][updown,DEcond] > DElim
    }
    if( length(sumGO) == 9 & !is.na(cmc) & deq  ){
      if(sum(sumGO[,"p.adjust"]<.05) >= numGO & cmc >= colMeanCor){
        print(paste(set,c))
        #bestgo <- as.data.frame(summary(wgcna.outputs[[set]][["clusterGO"]][[c]]))[1, "Description"]
        bestgo <- getSingleGo(wgcna.outputs[[set]][["clusterGO"]][[c]])
        apath <- file.path(datapath,"checkout",gsub(pattern = " ","",bestgo) ,paste0(strsplit(dirrr, "/")[[1]][6],"_",gsub("block","", set),"_",c ))
        dir.create(file.path(datapath,"checkout",gsub(pattern = " ","",bestgo)))
        try({
          system(paste("cp -r", file.path(dirrr, "Reports", set ,c),  apath))
          system(paste("cp ", file.path(dirrr, "Reports", set,"*.pdf"), apath))
          #dscolnames <-  c( colnames(TPME), colnames(TPMB), allAmby002OriginalSamplenames, knappOriginalSamplenames)
          nms <-  AllProbeDict[wgcna.outputs[[set]][["clusters"]][[c]], ]$Gene
          o <-  wgcna.outputs[[set]][["exprs"]][ wgcna.outputs[[set]][["clusters"]][[c]],]
          o <- cbind( rownames(wgcna.outputs[[set]][["exprs"]][ wgcna.outputs[[set]][["clusters"]][[c]],]) ,AllProbeDict[wgcna.outputs[[set]][["clusters"]][[c]], ]$Gene,o)
          colnames(o)[1:2] <-c("probe", "symbol")
          write.table(o,file.path(apath, paste0(c,"_","ExpressionSet.txt")),sep="\t",quote = F,row.names = F,col.names =T)
          for(datset in names(wgcna.outputs)){
            myplot(plot =barplot(wgcna.outputs[[set]][["MEs"]][[c]][[datset]] , main=paste(c,datset), cex.names = .25,las=3 ) ,filename = file.path(apath,paste0(datset, "ExpressionOfME.pdf")))
            moduleGenes[[datset]] <<- cbind(moduleGenes[[datset]],  wgcna.outputs[[set]][["MEs"]][[c]][[datset]] )
            colnames(moduleGenes[[datset]])[ncol(moduleGenes[[datset]])] <<- paste0(c,"\n", bestgo)
          }
        })
      }
    }
    
  }
  return(moduleGenes)
}

eigenGeneReports <- function(eg){
  sapply(names(eg),function(x){
    dir.create(file.path(datapath, "checkout",x))
    myplot(
      {matplot(eg[[x]] , type = "l", col =1:ncol(eg[[x]]), ylab = "PC1", cex=.20, xlab = "TP", main=paste("Module Eigengenes\n",x),xaxt = "n")
        legend("left", colnames(eg[[x]]), col = sapply(strsplit(colnames(eg[[x]]),split = "\n"), function(x)x[1]), pch = 1, horiz = F , text.width = 1,cex=.50)
        axis(1, 1:nrow(eg[[x]]), rownames(eg[[x]]),cex=.20, las=3)}
      ,filename = file.path(datapath, "checkout",x, "ME_chart.pdf"))
    
    cc <- cor(eg[[x]])
    cc[is.na(cc)] <- 0
    myplot(expressionHeatmap(cc,maint = "ME Pearson Cor"), filename = file.path(datapath, "checkout",x,"ME_Cor.pdf"))
    myplot(plot(hclust(dist(t(eg[[x]])), method="ward.D2"), cex = .6, xlab="Module", main=paste0("Distance Between Eigengenes\nof Modules in ", x) ),filename = file.path(datapath, "checkout",x,"ME Distances.pdf") )
  })
}

# eg <- findModules("~/code/data/axo/hClustering/MEANbicorsignedhybrid3", .40 , 15)#, set = "blockEmbryoSeq")#,which(names(DE.Universe) == "BlastemavsNot"),1,5)
# eigenGeneReports(eg)
wgcna.outputs$blockEmbryoSeq$MEs$black

makePPT <-  function(dirrr, colMeanCor=0.0, numGO = 0, set= NULL , maxGO.BG = 1.0 , DEcond=NULL, updown=NULL, DElim=NULL){
  #Load all the datasets
  library(ReporteRs)
  mydoc <<- pptx()
  mydoc <<- addSlide( mydoc, "Title and Content" )
  mydoc <<- addTitle(mydoc, "Params")
  setwd(dirrr)
  blokfile <- list.files("./", pattern = "block")
  wgcna.outputs <-  lapply(blokfile ,readRDS)
  names(wgcna.outputs) <- blokfile
  expnames <-  names(wgcna.outputs)
  #check if doing all
  if(!is.null(set)){
    expnames <- set
  }
  mydoc <<- addParagraph( mydoc, value = paste(dirrr,paste0("CorThreshold =",colMeanCor), paste0("GOThreshold =",numGO)), "Modules from sets:" ,expnames )
  apath <- file.path(datapath,"checkout",paste0(strsplit(dirrr, "/")[[1]][6],set,".pptx"))
  
  for(set in expnames){
    mydoc <<- addSlide( mydoc, "Title and Content" )
    mydoc <<- addTitle(mydoc, "Module Correlations Across All Sets")
    mydoc <<- addFlexTable(mydoc,flextable = FlexTable( as.data.frame(wgcna.outputs[[set]][["avg.cors"]]),add.rownames = T))
  for(c in names(wgcna.outputs[[set]][["clusters"]])[which(names(wgcna.outputs[[set]][["clusters"]])!="grey")]){
    print(c)
    sumGO <- summary(wgcna.outputs[[set]][["clusterGO"]][[c]])
    SGok <- F
    #If the genes are all unknown, there will be no GO summary

    if(length(sumGO)==9){ 
      sumGO <- sumGO[ sapply(sumGO[,"BgRatio"], fractionCharToNumeric) < maxGO.BG , ]
      SGok <-  sum(sumGO[,"p.adjust"]<.05, na.rm = T) >= numGO }
    else{
      sumGO <- NA
      SGok <- T}
    
    cmc <- mean(wgcna.outputs[[set]][["avg.cors"]][,c], na.rm=T )
    
    deq = T
    if(!is.null(DEcond)){
      deq <- wgcna.outputs[[set]][["DE.Profile"]][[c]][updown,DEcond] >= DElim
    }
    
    if(  deq  ){ #length(sumGO) == 9 & !is.na(cmc) &
      if( SGok & cmc >= colMeanCor){
        #bestgo <- as.data.frame(summary(wgcna.outputs[[set]][["clusterGO"]][[c]]))[1, "Description"]
          if(length(sumGO)!=9){
            bestgo <- NA
          }else{
            bestgo <- getSingleGo(wgcna.outputs[[set]][["clusterGO"]][[c]])
          }
          nms <-  AllProbeDict[wgcna.outputs[[set]][["clusters"]][[c]], ]$Gene
          print("Names")
          mydoc <<- addSlide( mydoc, "Title and Content" )
          mydoc <<- addTitle(mydoc, paste("Gene Names",set,c,bestgo))
          mydoc <<- addParagraph( mydoc, value = paste(nms, collapse=', ') )
          mydoc <<- addSlide( mydoc, "Title and Content" )
          print("GO")
          mydoc <<- addTitle(mydoc, paste("GO Enrichment",set,c,bestgo))
          if(length(sumGO)==9){mydoc <<- addFlexTable(mydoc,flextable = FlexTable(sumGO[1:50, c("ID", "Description","GeneRatio", "BgRatio", "p.adjust")]))}
          mydoc <<- addSlide( mydoc, "Two Content" )
          print("GRAPHS")
          mydoc <<- addTitle( mydoc, paste(set,c,bestgo))
          mydoc <<- addPlot(mydoc, function() expressionHeatmap(cor(t(wgcna.outputs[[set]][["exprs"]][wgcna.outputs[[set]][["clusters"]][[c]],])), paste0(c ," Pearson\n Gene Cors")), vector.graphic = F)
          mydoc <<- addPlot(mydoc, function(){barplot(DE.Rep(wgcna.outputs[[set]][["clusters"]][[c]],DE.Universe),  main =paste("DE Representation \n",set,c),  ylab = "Representation count" ,col = c("red","blue"),cex.names = .5,las=3 )  
             legend("topleft", c("DE Up", "DE Down"  ) , pch = 1, horiz = F , text.width = 1,cex=.60,
                    col = c("red","blue") )
          }, vector.graphic = F
          )
            mydoc <<- addSlide(mydoc,"Title and Content")
            mydoc <<- addTitle(mydoc, paste("ME for",c,bestgo))
            mydoc <<- addPlot(mydoc, function() {
              mypar(3,3)
              for(datset in names(wgcna.outputs)){
                barplot(as.numeric(unlist(wgcna.outputs[[set]][["MEs"]][[c]][[datset]])),names.arg =names(wgcna.outputs[[set]][["MEs"]][[c]][[datset]]) , main=paste(c,datset), cex.names = .25,las=3 )
                }
          }, vector.graphic = F)
            mydoc <<- addSlide(mydoc,"Title and Content")
            mydoc <<- addTitle(mydoc, paste("All Genes for",c,bestgo))
            mydoc <<- addPlot(mydoc, function() {
              mypar(3,3)
              for(datset in names(wgcna.outputs)){
                mm <-  wgcna.outputs[[datset]][["full.exprs"]]
                nms <-  intersect( AllProbeDict[wgcna.outputs[[set]][["clusters"]][[c]], ]$Gene, AllProbeDict[rownames(mm), ]$Gene) 
                if(length(nms)>1){
                ind <- AllProbeDict[rownames(mm), ]$Gene %in% nms
                M <- as.matrix(mm[ind,])
                matplot(t( (M-rowMeans(M))/rowSds(M) ),type = "l", main=paste(c,datset),ylab = "Row Z Score"  )
                }
              }
            }, vector.graphic = F)
      }
   }
  }
  writeDoc( mydoc, apath )
  }
}

makePPT(clustpaths[[1]], .0, 0, "blockEmbryoSeq",maxGO.BG = .2)
makePPT(clustpaths[[2]], .0, 0, "blockEmbryoSeq",maxGO.BG = .2)

makePPT(clustpaths[[1]], .0, 0, "blockGSE67118",maxGO.BG = .2)
makePPT(clustpaths[[2]], .0, 0, "blockGSE67118",maxGO.BG = .2)

grepl("stage_16",colnames(NormTPME) )

q <- (rowMeans(NormTPME[,grepl("Stage_16",colnames(NormTPME))])+1) - (rowMeans(NormTPME[,grepl("Stage_14",colnames(NormTPME))])+1)
matplot(t(NormTPME[names(q[order(q, decreasing = T)][1:75]), ]), type="l")
expressionHeatmap(cor(NormTPME[names(q[order(q, decreasing = T)]),35:45 ]))
expressionHeatmap(cor(NormTPME.unmapped.contigs[,35:45 ]))
expressionHeatmap(cor(NormTPME.unmapped.contigs[rownames(screen.unmapped.contigs)[(rownames(screen.unmapped.contigs) %in% rownames(NormTPME.unmapped.contigs))],]), maint= "Log2 TPM of Unmapped Contigs By Stage")

NormTPME.unmapped.contigs
expressionHeatmap(cor(screen.u[names(q[order(q, decreasing = T)]),35:45 ]))
expressionHeatmap(t(NormTPME[names(q[order(q, decreasing = T)][1:500]), ]))
expressionHeatmap(t(screen.unmapped.contigs[,-1]))
expressionHeatmap(t(NormTPME.unmapped.contigs[rownames(screen.unmapped.contigs)[(rownames(screen.unmapped.contigs) %in% rownames(NormTPME.unmapped.contigs))],]), maint= "Log2 TPM of Unmapped Contigs By Stage")

rownames(screen.unmapped.contigs)[(rownames(screen.unmapped.contigs) %in% rownames(NormTPME.unmapped.contigs))]
expressionHeatmap(t(NormTPME.unmapped.contigs[rownames(screen.unmapped.contigs)[(rownames(screen.unmapped.contigs) %in% rownames(NormTPME.unmapped.contigs))],]), maint= "Log2 TPM\n of Unmapped Contigs\n By Stage")
expressionHeatmap(t(log(screen.unmapped.contigs[,-1]+1,2)), maint = "Log2 TPM\n of Unmapped Contigs\n By Stage Averaged")

q[order(q, decreasing = T)][1:100]
MysteryBunch <- names(q[order(q, decreasing = T)][1:150])
NormTPME[MysteryBunch,]


screen.unmapped.contigs = file.path(datapath,"Human_Annotation/1_Screened_TPM_50_and_FoldChange_10_Genes.txt")
screen.unmapped.contigs <- read.table(screen.unmapped.contigs,header=TRUE)
rownames(screen.unmapped.contigs) <- paste0("contig",screen.unmapped.contigs[,1])
head(screen.unmapped.contigs)
gsub("contig","",MysteryBunch)%in%rownames(screen.unmapped.contigs)
gsub("contig","",MysteryBunch)
rownames(screen.unmapped.contigs) <- screen.unmapped.contigs$Gene.ID
screen.unmapped.contigs[gsub("contig","",MysteryBunch),]
screen.unmapped.contigs["603309",]
NormTPME["contig603309",]
NormTPME[MysteryBunch,]

sapply(wgcna.outputs$blockEmbryoSeq$clusters, function(x) sum(x %in% MysteryBunch))
x <- NormTPME
idx <- rowSds(x)> median(rowSds(x))
x <- x[idx,]
mean(idx)
sum(MysteryBunch %in% rownames(x))

?grepl
dirrr
dirrr = clustpaths[[2]]
dirrr = "~/code/data/axo/hClustering/corsignedhybrid3"
mydoc <<- pptx()
setwd(dirrr)
blokfile <- list.files("./", pattern = "block")
wgcna.outputs <-  lapply(blokfile ,readRDS)
names(wgcna.outputs) <- blokfile
apath <- file.path(datapath,"checkout",paste0(strsplit(dirrr, "/")[[1]][6],".pptx"))
set = names(wgcna.outputs)[[3]]
c = names(wgcna.outputs$blockallAmby002$clusters)[[1]]
sumGO <- summary(wgcna.outputs[[set]][["clusterGO"]][[c]])
bestgo <- getSingleGo(wgcna.outputs[[set]][["clusterGO"]][[c]])
nms <-  AllProbeDict[wgcna.outputs[[set]][["clusters"]][[c]], ]$Gene
mydoc <<- addSlide( mydoc, "Title and Content" )
mydoc <<- addTitle(mydoc, paste("Gene Names",set,c,bestgo))
mydoc <<- addParagraph( mydoc, value = paste(nms, collapse=', ') )
mydoc <<- addSlide( mydoc, "Title and Content" )
mydoc <<- addFlexTable(mydoc,flextable = FlexTable(sumGO[1:30, c("ID", "Description","GeneRatio", "BgRatio", "p.adjust")]))
mydoc <<- addSlide( mydoc, "Two Content" )
mydoc <<- addTitle( mydoc, paste(set,c,bestgo))
mydoc <<- addPlot(mydoc, function() expressionHeatmap(cor(t(wgcna.outputs[[set]][["exprs"]][wgcna.outputs[[set]][["clusters"]][[c]],])), paste0(c ," Pearson\n Gene Cors")), vector.graphic = F)
writeDoc( mydoc, apath )

library("timeline")
data(ww2)
ww2
ww2.events
# person n start end
# 1         sam 6     0   6
# 2        greg 5     6  11
# 3     teacher 4    11  15
# 4         sam 4    15  19
# 5        greg 5    19  24
# 6       sally 5    24  29
# 7        greg 4    29  33
# 8         sam 3    33  36
# 9       sally 5    36  41
# 10 researcher 6    41  47
# 11       greg 6    47  53

Experiment = "blockGSE67118"
tps = as.integer(GSE67118tp)

startstop<<-data.frame("name"=character(),goName = character(), "start" =integer(), "end" = integer(),stringsAsFactors=FALSE)
sapply(names(wgcna.outputs[[Experiment]][["MEs"]]), function(men){
#sapply("blue", function(men){
  me = wgcna.outputs[[Experiment]][["MEs"]][[men]][[Experiment]]
  start = 0
  end = 0
  sumGO = summary(wgcna.outputs[[Experiment]][["clusterGO"]][[men]])
  if(length(sumGO)==9){
  sumGO <- sumGO[ sapply(sumGO[,"BgRatio"], fractionCharToNumeric) < .15 , ]
    goName <- paste( sumGO[1,"Description"],sumGO[3,"Description"],sep = "\n")
    }
  else{goName="XXX"}
  if(grepl("XXX",goName))goName= men
  if(grepl("NA\nNA",goName))goName= men
  for(i in 1:length(me)){
    if(me[i]>0.3 & start == 0){
      #start = as.integer(strsplit( names(wgcna.outputs[[Experiment]][["MEs"]][[men]][[Experiment]])[i], "_")[[1]][2])
      start = tps[i]
    }
    if(me[i]<=0.3 & end == 0 & start !=0){
      #end = as.integer(strsplit( names(wgcna.outputs[[Experiment]][["MEs"]][[men]][[Experiment]])[i], "_")[[1]][2])
      end = tps[i]
      startstop[nrow(startstop)+1,] <<- c("name"=men,"goName"=goName, "start" =start, "end" = end)
      start = 0
      end = 0
    }
    if(i == length(me) & start != 0){
      #end = as.integer(strsplit( names(wgcna.outputs[[Experiment]][["MEs"]][[men]][[Experiment]])[i], "_")[[1]][2])
      end = tps[i]
      startstop[nrow(startstop)+1,] <<- c("name"=men,"goName"=goName, "start" =start, "end" = end)
      start = 0
      end = 0
    }
  }
})
startstop$start <- as.integer(startstop$start)
startstop$end <- as.integer(startstop$end)
startstop
#library(ggplot2)
ggplot(startstop, aes(colour=goName)) + 
  geom_segment(aes(x=start, xend=end, y=goName, yend=goName), size=3) + 
  ggtitle("Upregulation of Modules\n in Axolotl Limb Regeneration")+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #labs(x=names(wgcna.outputs$blockEmbryoSeq$MEs$blue$blockEmbryoSeq))+
  xlab("Time (days)")
  
