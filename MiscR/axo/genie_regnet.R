source("http://bioconductor.org/biocLite.R")
library(EBSeq)
library(matrixStats)
library(limma)
#install.packages("randomForest", repos="http://cran.rstudio.com/")
library(randomForest)
library(devtools)
datapath <-"~/code/data/axo"


options("cores"=16)
source("~/code/axo/genie3.R")
source("~/code/axo/LoadDatasets.R")
set.seed(42)

blastemaRNAseqNetGenie <-  get.weight.matrix( NormTPMB, input.idx = intersect(fullRegulatorList, rownames(NormTPMB)),K="all")
saveRDS(blastemaRNAseqNetGenie ,file = file.path(datapath,"genie_outputs", "regOnlyConsideredNorm","blastemaRNAseqNetGenie.rds"))
remove(blastemaRNAseqNetGenie)
#get.link.list(blastemaRNAseqNetGenie,100)
# get.link.list(readRDS(file.path(datapath,"genie_outputs", "regOnlyConsidered","blastemaRNAseqNetGenie.rds")),100)

embryoRNAseqNetGenie <-  get.weight.matrix( NormTPME, input.idx = intersect(fullRegulatorList, rownames(NormTPME)),K="all") 
saveRDS(embryoRNAseqNetGenie ,file = file.path(datapath,"genie_outputs", "regOnlyConsideredNorm","embryoRNAseqNetGenie.rds"))
remove(embryoRNAseqNetGenie)

# woundData <- normalizeQuantiles(cbind(exprsGSE35255[,GSE35255samplenames=="Aquatic"],exprsGSE37198[,GSE37198samplenames=="Flankwound"]))
# woundAmby002NetGenie <-  get.weight.matrix(woundData,K="all", input.idx = regulatorProbes)
# saveRDS(woundAmby002NetGenie , file = file.path(datapath,"genie_outputs", "regOnlyConsideredNorm","woundAmby002NetGenie.rds"))
# remove(woundAmby002NetGenie)
# 
# EBMatchRNAseqNetGenie <-  get.weight.matrix( NormEBMatch, input.idx = intersect(fullRegulatorList, rownames(NormEBMatch)),K="all") 
# saveRDS(EBMatchRNAseqNetGenie ,file = file.path(datapath,"genie_outputs", "regOnlyConsidered","NormEBMatch.rds"))
# remove(EBMatchRNAseqNetGenie)
# 
# blastemaData <- normalizeQuantiles(exprsGSE67118[,GSE67118samplenames!="T0"])
# blastemaAmby002NetGenie <-  get.weight.matrix(blastemaData, input.idx = regulatorProbes,K="all")
# saveRDS(blastemaAmby002NetGenie , file = file.path(datapath,"genie_outputs", "regOnlyConsideredNorm","blastemaAmby002NetGenie.rds"))
# remove(blastemaAmby002NetGenie)
# 
# allAmby002Data <- normalizeQuantiles(cbind(exprsGSE67118,exprsGSE35255,exprsGSE37198))
# allAmby002NetGenie <-  get.weight.matrix(allAmby002Data, input.idx = regulatorProbes,K="all")
# saveRDS(allAmby002NetGenie ,file = file.path(datapath,"genie_outputs", "regOnlyConsideredNorm","allAmby002NetGenie.rds"))
# remove(allAmby002NetGenie)
# 
# blastemawoundGSE36451 <- get.weight.matrix(exprsGSE36451, input.idx = KnappRegulatorProbes,K="all")
# saveRDS(blastemawoundGSE36451 ,file = file.path(datapath,"genie_outputs", "regOnlyConsideredNorm","blastemawoundGSE36451.rds"))


source("~/code/axo/genie3.R")
setwd(file.path("~/code/data/axo/genie_outputs/regOnlyConsideredNorm"))
genie <- list.files("./", pattern = "")
networks <-  lapply(genie ,readRDS)
names(networks) <- gsub(".rds","",genie)
for(i in 1:length(networks)){
  tab <- get.link.list(networks[[i]], 5000)
  write.table(tab ,row.names = F, col.names = T, quote = F , file = file.path(datapath,"genie_outputs", "regOnlyConsideredNorm",paste0(names(networks)[i],".txt"))) 
}
