datapath <- "~/code/data/cb/NBNogData"
setwd(datapath)
datafiles <-  file.path(datapath, dir()[ grepl(".csv", dir())])
datasets <- sapply( datafiles, read.csv, header=T, sep=",",stringsAsFactors=F)

HumanTFs <- toupper(read.table(paste0("~/code/data/general/marker_lists/my_human_TF_v5_final_11282011.txt"),header=T, sep="\t")$symbol)
MouseTFs <-  toupper(read.table(paste0("~/code/data/general/marker_lists/RIKEN_Mouse_TFdb.txt"),header=T, sep="\t")$Symbol)
oncogenes <- toupper(read.table( file.path("~/code/data/general/marker_lists/MSKCCOncogenes.txt"), header = T, sep = "\t")$Gene.Symbol)
neural_related_genes <- make.names(toupper(as.character(read.csv("~/code/data/cb/markers/neuralrelatedgenes.csv",header = T)[,1])))

for(i in 1:length(datasets)){
  #rownames(datasets[[i]]) <-  make.names(toupper(datasets[[i]][,1]),unique = T)
  datasets[[i]] <- datasets[[i]][!duplicated(toupper(datasets[[i]][,1])),]
  rownames(datasets[[i]]) <-  toupper(datasets[[i]][!duplicated(datasets[[i]][,1]),1])
  colnames(datasets[[i]]) <-  gsub("X","",colnames(datasets[[i]]))
  datasets[[i]] <- as.matrix(datasets[[i]][sort(rownames(datasets[[i]])),-1])
  #class(datasets[[i]]) <- "numeric"
}
dfnames <- gsub(".csv", "", dir()[ grepl(".csv", dir())])
for(i in 1:length(datasets)){
  colnames(datasets[[i]]) <-  paste(dfnames[i] ,colnames(datasets[[i]]))
}
datasetNames <- c("HumanTControl","HumanTera","MouseTControl","MouseTera","17DayMouse","17DayHumanAndMix","33DayHuman","33DayMouse","33DayMix","LongMouse","NBNogHuman100","NBNogHuman10","NBNogHuman90","NBNogMouse0","NBNogMouse10", "NBNogMouse90")

mixDir <- file.path(datapath, "MixSets")
setwd(mixDir)
AllMixSetsEC <- lapply(dir(), function(d){
  setwd(file.path(mixDir,d,  dir(file.path(mixDir,d))[which(grepl("results.hg19_and_mm10.mm10",dir(file.path(mixDir,d))))]))
  mixdatafiles <- list()
  mixdatafiles[["mm"]] <-  file.path(getwd(), dir()[ grepl("genes.no_mt.ec.tab", dir())])
  setwd(file.path(mixDir, d, dir(file.path(mixDir,d))[which(grepl("results.hg19_and_mm10.hg19",dir(file.path(mixDir,d))))]))
  mixdatafiles[["hs"]] <-  file.path(getwd(), dir()[ grepl("genes.no_mt.ec.tab", dir())])
  MS <-  lapply( mixdatafiles, read.csv2, header=T, sep="\t",row.names=1)
  for(i in 1:length(MS)){
    #rownames(MS[[i]]) <- make.names(toupper(rownames(MS[[i]])),unique = T)
    MS[[i]] <- MS[[i]][!duplicated(toupper(rownames(MS[[i]]))),]
    rownames(MS[[i]]) <- toupper(rownames(MS[[i]]))
    colnames(MS[[i]]) <-  gsub("X","",colnames(MS[[i]]))
    MS[[i]] <- as.matrix(MS[[i]][,-ncol(MS[[i]])])
    #MS[[i]] <- matrix(as.numeric(MS[[i]]) , nrow = nrow(MS[[i]]), ncol = ncol(MS[[i]]), dimnames =dimnames(MS[[i]]))
    MS[[i]] <- matrix(as.numeric(MS[[i]]) , nrow = nrow(MS[[i]]), ncol = ncol(MS[[i]]), dimnames =dimnames(MS[[i]]))
  }
  MS
})

MixSetsEC <- lapply(AllMixSetsEC, function(y){
  lapply(y, function(x){
    rownames(x) <- gsub("\\|HG19", "", rownames(x))
    rownames(x) <- gsub("\\|MM10", "", rownames(x))
    x
  })
} )

MixSetsEC <- lapply(1:length(MixSetsEC),function(x) Reduce(rn.merge,lapply(MixSetsEC,"[[",x)))
names(MixSetsEC) <- c("mm","hs")
MixSetsEC <- lapply(MixSetsEC, function(y){
  colnames(y) <-  gsub("_rep\\d+", "", colnames(y))
  y[,!grepl(pattern = "control",x = colnames(y))]
})

names(MixSetsEC) <- c("mm","hs")


mixDir <- file.path(datapath, "MixSets")
setwd(mixDir)
AllMixSets <- lapply(dir(), function(d){
  setwd(file.path(mixDir, d, dir(file.path(mixDir,d))[which(grepl("results.hg19_and_mm10.mm10",dir(file.path(mixDir,d))))]))
  mixdatafiles <- list()
  mixdatafiles[["mm"]] <-  file.path(getwd(), dir()[ grepl("genes.no_mt.tpm.rescale.tab", dir())])
  setwd(file.path(mixDir, d, dir(file.path(mixDir,d))[which(grepl("results.hg19_and_mm10.hg19",dir(file.path(mixDir,d))))]))
  mixdatafiles[["hs"]] <-  file.path(getwd(), dir()[ grepl("genes.no_mt.tpm.rescale.tab", dir())])
  MS <-  lapply( mixdatafiles, read.csv2, header=T, sep="\t",row.names=1)
  for(i in 1:length(MS)){
    #rownames(MS[[i]]) <- make.names(toupper(rownames(MS[[i]])),unique = T)
    MS[[i]] <- MS[[i]][!duplicated(toupper(rownames(MS[[i]]))),]
    rownames(MS[[i]]) <- toupper(rownames(MS[[i]]))
    colnames(MS[[i]]) <-  gsub("X","",colnames(MS[[i]]))
   MS[[i]] <- as.matrix(MS[[i]][,-ncol(MS[[i]])])
   MS[[i]] <- matrix(as.numeric(MS[[i]]) , nrow = nrow(MS[[i]]), ncol = ncol(MS[[i]]), dimnames =dimnames(MS[[i]]))
  }
MS
})

MixSets <- lapply(AllMixSets, function(y){
  lapply(y, function(x){
    rownames(x) <- gsub(".HG19", "", rownames(x))
    rownames(x) <- gsub(".MM10", "", rownames(x))
    x
  })
} )
MixSets <- lapply(1:length(MixSets),function(x) Reduce(rn.merge,lapply(MixSets,"[[",x)))

MixSets <- lapply(MixSets, function(y){
  colnames(y) <-  gsub("_rep\\d+", "", colnames(y))
  y[,!grepl(pattern = "control",x = colnames(y))]
})

MixSets <-lapply(1:length(MixSets),function(i) MixSets[[i]][,lapply(lapply(MixSetsEC,colSums),">",5e5)[[i]]])
MixSetsEC <-lapply(1:length(MixSetsEC),function(i) MixSetsEC[[i]][,lapply(lapply(MixSetsEC,colSums),">",5e5)[[i]]])
tpMixSets <- lapply(MixSets, function(x){
  as.numeric(sapply(colnames(x), function(y) strsplit(y,"d")[[1]][2]))
})
names(tpMixSets) <- c("mm","hs")

conditMixSets <- lapply(MixSets, function(x){
  sapply(colnames(x), function(y) strsplit(y,"_")[[1]][1])
})
names(conditMixSets) <- c("mm","hs")


orderAll <- lapply(1:length(conditMixSets),function(z) order(as.numeric(conditMixSets[[z]]),tpMixSets[[z]]))
conditMixSets <- lapply(1:length(conditMixSets),function(z)conditMixSets[[z]][orderAll[[z]]])
tpMixSets <- lapply(1:length(tpMixSets),function(z)tpMixSets[[z]][orderAll[[z]]])
MixSets <- lapply(1:length(MixSets),function(z)MixSets[[z]][,orderAll[[z]]])
MixSetsEC <- lapply(1:length(MixSetsEC),function(z)MixSetsEC[[z]][,orderAll[[z]]])
names(MixSets) <-names(tpMixSets) <-names(conditMixSets) <- c("mm","hs")#<-names(MixSetsEC)
 
#conditMixSets[["mm"]][1:3] <- NA
#conditMixSets[["hs"]][conditMixSets[["hs"]]=="0"& tpMixSets[["hs"]] ==0 ] <- NA

condit.combos <- expand.grid(c("mm","hs"),unique(conditMixSets[[2]]),stringsAsFactors = F)
condit.combos <- condit.combos[-c(2,7),]
condit.combos <- condit.combos[c(6,3,5,1,2,4),]

datasets[11:16] <- lapply(1:6, function(i){ MixSets[[condit.combos[i,1]]][,conditMixSets[[condit.combos[i,1]]] == condit.combos[i,2]] })
MixSetsECMat <-  lapply(1:6, function(i){ 
  m <- MixSetsEC[[condit.combos[i,1]]][,conditMixSets[[condit.combos[i,1]]] == condit.combos[i,2]] 
  #GetNormalizedMat(m,MedianNorm(m))
  })

datasets <- lapply(datasets,function(dd)dd[!is.na(rowSums(dd,na.rm = F)),])

all.genes <- Reduce(union, sapply(MixSets,rownames))
#neural_related_genes[!neural_related_genes%in% rownames(MixSets$hs)]

# datapath <- "~/code/data/cb"
# setwd(datapath)
# mitodatafiles <-  file.path(datapath, dir()[ grepl(".tsv", dir())])
# mitodatasets <- sapply( mitodatafiles, read.csv, header=T, sep="\t",row.names=1)
# head(mitodatasets[[3]])
# for(i in 1:length(mitodatasets)){
#   rownames(mitodatasets[[i]]) <-  make.names(toupper(rownames(mitodatasets[[i]])),unique = T)
#   mitodatasets[[i]] <- log(as.matrix(mitodatasets[[i]][sort(rownames(mitodatasets[[i]])),])+1,2)
#   #class(mitodatasets[[i]]) <- "numeric"
# }
# matplot(t(mitodatasets$`~/code/data/cb/human_TERA_mito.tsv`),type="l",main="Human MT Genes in Teratoma")
# matplot(t(mitodatasets$`~/code/data/cb/human_TERA_mito_renorm.tsv`[rownames(mitodatasets$`~/code/data/cb/mouse_TERA_mito.tsv`),]),type="l")
# matplot(t(mitodatasets$`~/code/data/cb/mouse_TERA_mito.tsv`),type="l",main="Mouse MT Genes in Teratoma")
# rownames(mitodatasets$`~/code/data/cb/human_TERA_mito.tsv`)
#neural_list <- sort(as.character(read.csv("~/code/data/cb/markers/marker_union.csv",header = T)[,2]))

samples_1 <- colnames(datasets$`~/code/data/cb/NBNogData/TC_1.csv`)
samples_1_M <- colnames(datasets$`~/code/data/cb/NBNogData/TC_1_M.csv`)
tp_1_M <- as.numeric(gsub('.1',"",sapply(samples_1_M, function(x) strsplit(x,"d")[[1]][2])))
tp_1 <- as.numeric(sapply(samples_1, function(x) strsplit(x,"d")[[1]][2]) )
tp_2 <- colnames(datasets$`~/code/data/cb/NBNogData/TC_2_Mix.csv`)
tp_2 <- as.numeric(gsub("TC_2_Mix ","", tp_2))

neural_list <- toupper(sort(as.character(read.csv("~/code/data/cb/markers/FullNeuralGeneList.csv",header = T)[,1])))
#neural_list <- toupper(sort(as.character(read.table("~/code/data/cb/markers/FavNeuralGeneList3.txt",header = F,sep = "\t")[,1])))
cc_list <- sort(as.character(read.csv2("~/code/data/cb/markers/cell_cycle_list.txt",header = T,sep = "\t",row.names = 1)[,2]))

tp.m.long <- as.numeric( sapply(strsplit(colnames(datasets$`~/code/data/cb/NBNogData/TC_L_150415.csv`)," "),"[",2) )

tp.h.control <- as.numeric(sapply(colnames(datasets$`~/code/data/cb/NBNogData/human_ALONE.csv`), function(x) strsplit(x,"d")[[1]][2]))
tp.m.control <- as.numeric(sapply(colnames(datasets$`~/code/data/cb/NBNogData/mouse_ALONE.csv`), function(x) strsplit(x,"d")[[1]][2]))
tp.m.control[is.na(tp.m.control)] <- 0

tp.h.tera <- as.numeric(sapply(colnames(datasets$`~/code/data/cb/NBNogData/human_TERA.csv`), function(x) strsplit(x," ")[[1]][2]))
tp.m.tera <- as.numeric(sapply(colnames(datasets$`~/code/data/cb/NBNogData/mouse_TERA.csv`), function(x) strsplit(x," ")[[1]][2]))
tp.h.tera[is.na(tp.h.tera)] <- 0
tp.m.tera[is.na(tp.m.tera)] <- 0

ngm <- function(mat){
  ng <- neural_list[neural_list%in%rownames(mat)]
  return(mat[ng,])
}
collapseMats <- function(mat1, mat2){
  g <- intersect(rownames(mat1), rownames(mat2))
  return(list(mat1[g,], mat2[g,]))
}

tpDatasets <- list(tp.h.control,tp.h.tera,tp.m.control, tp.m.tera,tp_1_M,tp_1, tp_2,tp_2,tp_2,tp.m.long)
tpDatasets <- c(tpDatasets, lapply(1:6, function(i){ tpMixSets[[condit.combos[i,1]]][conditMixSets[[condit.combos[i,1]]] == condit.combos[i,2]] }))

tpDatasets[11:16] <- lapply(1:6, function(i){ if(!0%in% tpMixSets[[condit.combos[i,1]]][conditMixSets[[condit.combos[i,1]]] == condit.combos[i,2]]){
  c(tpMixSets[[condit.combos[i,1]]][tpMixSets[[condit.combos[i,1]]]==0],tpDatasets[[10+i]])
}else{
  tpDatasets[[10+i]]
}})
datasets[11:16] <- lapply(1:6, function(i){ if(!0%in% tpMixSets[[condit.combos[i,1]]][conditMixSets[[condit.combos[i,1]]] == condit.combos[i,2]]){
  cbind(MixSets[[condit.combos[i,1]]][,tpMixSets[[condit.combos[i,1]]]==0],datasets[[10+i]])
}else{
  datasets[[10+i]]
}})


datapathTest <- "~/code/data/cb/NBNogData/MixingTest"
setwd(datapathTest)
hMat <- as.matrix(read.csv2(file.path(datapathTest,"All_Sub527_hg19_ReplacedWithSub537samples.txt"),sep = "\t",header = T,row.names = 1))
hMat <- hMat[,-c(ncol(hMat))]
hMat <- matrix(as.numeric(hMat) , nrow = nrow(hMat), ncol = ncol(hMat), dimnames =dimnames(hMat))
rownames(hMat) <- toupper(rownames(hMat))
mMat <- read.csv2(file.path(datapathTest,"All_Sub527_mm10.txt"),sep = "\t",header = T)
mMat <- as.matrix(mMat)
rownames(mMat) <- mMat[,1]
mMat <- mMat[!duplicated(mMat[,1]),]
mMat <- mMat[,-c(1,ncol(mMat))]
mMat <- matrix(as.numeric(mMat) , nrow = nrow(mMat), ncol = ncol(mMat), dimnames =dimnames(mMat))
rownames(mMat) <- toupper(rownames(mMat))
colnames(hMat) <- gsub("_H2BmChIP13","",colnames(hMat))
colnames(hMat) <- gsub("Split_","Split",colnames(hMat))
conditH <- sapply(strsplit(colnames(hMat),"_d"),"[[",1)
tpHTest <- as.numeric(sapply(strsplit(colnames(hMat),"_d"),"[[",2))

colnames(mMat) <- gsub("_H2BmChIP13","",colnames(mMat))
colnames(mMat) <- gsub("Split_","Split",colnames(mMat))
conditM <- sapply(strsplit(colnames(mMat),"_d"),"[[",1)
tpMTest <- as.numeric(sapply(strsplit(colnames(mMat),"_d"),"[[",2))
tpHTest <- tpHTest[grepl('H9',conditH)]
hMat <- hMat[,grepl('H9',conditH)]
tpMTest <- tpMTest[grepl('X3',conditM)]
mMat <- mMat[,grepl('X3',conditM)]
conditM <- conditM[grepl('X3',conditM)]
conditH <- conditH[grepl('H9',conditH)]
spec <- list("Mouse","Human")
condits <- list(conditM,conditH)
tps <- list(tpMTest,tpHTest)
mats <- list(mMat,hMat)

datasets <- c(datasets,mats)
tpDatasets <- c(tpDatasets,tps)
datasetNames <- c(datasetNames,"MouseProtocolTest","HumanProtocolTest")
# for(j in 1:2){
#   for( i in which(!duplicated(condits[[j]]))){
#     print(paste0(spec[[j]],condits[[j]][i]))
#     curMat <- mats[[j]][,which(condits[[j]]==condits[[j]][i])]
#     curTP <- tps[[j]][which(condits[[j]]==condits[[j]][i])]
#     if(spec[[j]]=="Human"){
#       curMat <- cbind(mats[[j]][,which(condits[[j]]== "H9")],curMat)
#       curTP <- c(tps[[j]][which(condits[[j]]== "H9")],curTP)
#     }
#     if(spec[[j]]=="Mouse"){
#       curMat <- cbind(mats[[j]][,which(condits[[j]]== "X3")],curMat)
#       curTP <- c(tps[[j]][which(condits[[j]]== "X3")],curTP)
#     }
#     datasets[[length(datasets)+1]]  <- curMat
#     tpDatasets[[length(tpDatasets)+1]] <-  unlist(curTP)
#     datasetNames[[length(datasetNames)+1]] <- paste0(spec[[j]],condits[[j]][i])
#     }
# }


datapathIncTemp <- "~/code/data/cb/NBNogData/IncTempComparison/"
setwd(datapathIncTemp)

data <- read.csv2(paste0(datapathIncTemp,"genes.no_mt.tpm.rescale.tab"), header=T, sep="\t",row.names=1)
data <- data[,-ncol(data)]
data <- as.matrix(data)
storage.mode(data) <- "numeric"
colnames(data) <- gsub("X","",colnames(data))

conditions <- sapply(colnames(data),function(x) unlist(strsplit(x,'_')))
conditions[3,] <- as.integer(gsub('d','',conditions[3,]))
datasetsIT <- lapply(unique(conditions[2,]), function(x){
  data[,conditions[2,]==x]
} )
names(datasetsIT) <- unique(conditions[2,])
tpsIT <- lapply(unique(conditions[2,]), function(x){
  q=conditions[3,]
  q[conditions[2,]==x]
} )
names(tpsIT) <- unique(conditions[2,])
datasets <- c(datasets,datasetsIT)
tpDatasets <- c(tpDatasets,tpsIT)
datasetNames <- c(datasetNames, paste0("Human",names(tpsIT)))

data <- read.csv2("~/code/data/cb/NBNogData/GSE56796.results.hg19/genes.no_mt.tpm.renorm.tab", header=T, sep="\t",row.names=1)
data <- data[,-ncol(data)]
data <- as.matrix(data)
storage.mode(data) <- "numeric"

colnames(data) <- gsub("X[0-9]+_","",colnames(data))
colnames(data) <- gsub("d","",colnames(data))
colnames(data) <- gsub("Day.","",colnames(data))
colnames(data) <- gsub(".[A-Z+]","",colnames(data))
datasets <- c(datasets,list(data[,1:24],data[,25:ncol(data)]))
tpDatasets <- c(tpDatasets, list(colnames(data)[1:24]), list(colnames(data)[25:ncol(data)]))
datasetNames <- c(datasetNames, "HumanCORTECONHighDensity","HumanCORTECONLowDensity")

datapathRobot <- "~/code/data/cb/NBNogData/robotHiResolution/"
setwd(datapathRobot)

data <- read.csv2(paste0(datapathRobot, "/Sub_0146_Developmental_Clock_Human_hg19__890bf46a388bf6b4/","genes.no_mt.tpm.rescale.tab"), header=T, sep="\t",row.names=1)
data <- data[,-ncol(data)]
data <- as.matrix(data)
storage.mode(data) <- "numeric"
colnames(data) <- gsub("X","",colnames(data))
data <- data[,colMeans(cor(data,method = "spe"))>.93]

conditions <- sapply(colnames(data),function(x) unlist(strsplit(x,'_')))
tpsRob <- as.numeric(gsub('min|\\.1','',conditions[2,]))/(60*24)

datasets <- c(datasets,list(data))
tpDatasets <- c(tpDatasets,list(tpsRob))
datasetNames <- c(datasetNames, "NBNogHiResRobotHuman")


data <- read.csv2(paste0(datapathRobot, "/Sub_0267_Developmental_Clock_Human_hg19__da2c173ac11fc21c/","genes.no_mt.tpm.rescale.tab"), header=T, sep="\t",row.names=1)
data <- data[,-ncol(data)]
data <- as.matrix(data)
storage.mode(data) <- "numeric"
colnames(data) <- gsub("X","",colnames(data))
data <- data[,colMeans(cor(data,method = "spe"))>.93]
data <- data[,-c(3,4)]


conditions <- sapply(colnames(data),function(x) unlist(strsplit(x,'_')))
tpsRobE8 <- as.integer(gsub('mins|min','',conditions[2,]))/(60*24)

datasets <- c(datasets,list(data))
tpDatasets <- c(tpDatasets,list(tpsRobE8))
datasetNames <- c(datasetNames, "E8HiResRobotHuman")


data <- read.csv2(paste0(datapathRobot, "/Sub_0545_830c1_mm10_00128a5650c3660a/","genes.no_mt.tpm.rescale.tab"), header=T, sep="\t",row.names=1)
data <- data[,-ncol(data)]
data <- as.matrix(data)
storage.mode(data) <- "numeric"
rownames(data) <- toupper(rownames(data))
colnames(data) <- gsub("X","",colnames(data))
data <- data[,colMeans(cor(data,method = "pear"))>.92]

conditions <- sapply(colnames(data),function(x) unlist(strsplit(x,'_')))
tpsMRob <- as.integer(gsub('mins|min','',conditions[2,]))/(60*24)

datasets <- c(datasets,list(data))
tpDatasets <- c(tpDatasets,list(tpsMRob))
datasetNames <- c(datasetNames, "NBNogHiResRobotMouse")


datapathMega <- "~/code/data/lfc/MegaTimeSeries/"
setwd(datapathMega)

data <- read.csv2(paste0(datapathMega, "HUMAN_normEC_96h.csv"), header=T, sep=",", row.names=1)
data <- as.matrix(data)
storage.mode(data) <- "numeric"
rownames(data) <- toupper(rownames(data))
colnames(data) <- gsub("X","",colnames(data))
#data <- data[,colMeans(cor(data,method = "spe"))>.93]
tpsHMega <- as.integer(gsub('h|\\_p','',colnames(data)))/24
datasets <- c(datasets,list(data))
tpDatasets <- c(tpDatasets,list(tpsHMega))
datasetNames <- c(datasetNames, "HumanMegaEndoderm")

data <- read.csv2(paste0(datapathMega, "MOUSE_normEC.csv"), header=T, sep=",", row.names=1)
data <- as.matrix(data)
storage.mode(data) <- "numeric"
colnames(data) <- gsub("X","",colnames(data))
rownames(data) <- toupper(rownames(data))
#data <- data[,colMeans(cor(data,method = "spe"))>.93]
tpsMMega <- as.integer(gsub('R1_|\\_p','',colnames(data)))/24
rownames(data) <- toupper(rownames(data))
datasets <- c(datasets,list(data))
tpDatasets <- c(tpDatasets,list(tpsMMega))
datasetNames <- c(datasetNames, "MouseMegaEndoderm")


tpDatasets <- lapply(tpDatasets,as.numeric)
names(tpDatasets) <- datasetNames
names(datasets) <- datasetNames


carnegie.equivalents <- read.table(file.path(datapath,"carnegie_equivalents.txt"), sep = "\t",header = T)

if(T){
BrainSpan <- read.csv2(file = "~/code/data/cb/NBNogData/brainspan_devcourse/expression_matrix.csv",sep = ",",header = F,row.names = 1)
BrainSpanCols <-  read.csv2(file = "~/code/data/cb/NBNogData/brainspan_devcourse/columns_metadata.csv",sep = ",",header = T,row.names = 1)
BrainSpanRows <-  read.csv2(file = "~/code/data/cb/NBNogData/brainspan_devcourse/rows_metadata.csv",sep = ",",header = T,row.names = 1)
BrainSpan <- as.matrix(BrainSpan)
BrainSpan <- matrix(as.numeric(BrainSpan) , nrow = nrow(BrainSpan), ncol = ncol(BrainSpan), dimnames =dimnames(BrainSpan))
#rownames(BrainSpan) <-  make.names(BrainSpanRows$gene_symbol,unique = T)
colnames(BrainSpan) <- paste(BrainSpanCols$structure_acronym,BrainSpanCols$age)
BrainSpanFPKM <- BrainSpan
BrainSpan <- (sweep(BrainSpan, 2, colSums(BrainSpan), "/") * 10^6) 
rownames(BrainSpan) <- rownames(BrainSpanFPKM) <-toupper(rownames(BrainSpan))
BrainSpanCols$structure_name <- as.character(BrainSpanCols$structure_name)

MHOrthos <-  read.table(file.path(datapath,"MouseHumanOrthologs.txt"),header = F,sep = "\t")[,1:4]
colnames(MHOrthos) <- c("HSGENE","HSENTREZ","MMENTREZ","MMGENE")
MHOrthos$MMGENE <- toupper(MHOrthos$MMGENE)
MHOrthos$HSGENE <- toupper(MHOrthos$HSGENE)

MtoH.Orthologs <- function(l , MtoH=F, entrez=F){
  if(entrez){
    x1=MHOrthos$HSENTREZ[match(x=l , MHOrthos$HSENTREZ)]
    x2=MHOrthos$MMENTREZ[match(x=l , MHOrthos$HSENTREZ)]
    n1="hs"
    n2="mm"
  }else{
    x1=MHOrthos$HSGENE[match(x=l , MHOrthos$HSGENE)]
    x2=MHOrthos$MMGENE[match(x=l , MHOrthos$HSGENE)]
    n1="hs"
    n2="mm"
  }
  if(MtoH){
    x1=MHOrthos$MMGENE[match(x=l , MHOrthos$MMGENE)]
    x2=MHOrthos$HSGENE[match(x=l , MHOrthos$MMGENE)]
    n2="hs"
    n1="mm"
  }
  print(length(x1))
  print(length(x2))
  r <- na.omit(data.frame(n1 = x1 , n2 = x2,stringsAsFactors = F))
  colnames(r) <- c(n1,n2)
  return(r)
}
}
