#MTSchmitz Meta Analysis of Axolotl Blastema
#loads/writes Blastema dataset files
#GSE35255 , GSE37198 , GSE67118, GSE36451,GSE34394 + Embryo and Blastema RNA-Seq
source("http://bioconductor.org/biocLite.R")
#biocLite("EBSeq")
#biocLite("GEOquery")
#biocLite("matrixStats")
#biocLite("limma")
library(EBSeq)
library(GEOquery)
library(matrixStats)
library(limma)

#### SET THIS TO THE PATH OF metaAxo
datapath <- "~/code/data/axo"
#Set working directory for download
setwd(datapath)

#GSE35255 , GSE37198 , GSE67118, GSE36451,GSE34394
#Blastema timecourse Voss et al
GSE67118 <- getGEO("GSE67118",destdir=datapath, GSEMatrix =TRUE)
if (length(GSE67118) > 1) idx <- grep("GPL15153", attr(GSE67118, "names")) else idx <- 1
GSE67118 <- GSE67118[[idx]]
#Wound Healing, normal and metamorphosized axolotl
GSE35255 <- getGEO("GSE35255",destdir=datapath, GSEMatrix =TRUE)
if (length(GSE35255) > 1) idx <- grep("GPL15153", attr(GSE35255, "names")) else idx <- 1
GSE35255 <- GSE35255[[idx]]
#Aquatic axolotl full thickness epithelial flank wounds, innervated limbs, and denervated limbs collected over seven days
GSE37198 <- getGEO("GSE37198",destdir=datapath, GSEMatrix =TRUE)
if (length(GSE37198) > 1) idx <- grep("GPL15153", attr(GSE37198, "names")) else idx <- 1
GSE37198 <- GSE37198[[idx]]

# ##Load Blastema RNA Seq Stewart Expected counts !REMOVES LOWEST EXPRESSORS, log2
# blastemaFile <- file.path(datapath, "tanaka_norrna_rsem_human_gene_expressionForUnixJustJuvTCECs.txt")
# ECBAll <- read.table(blastemaFile,header=TRUE)
# rownames(ECBAll) <-ECBAll$symbol
# ECB<- as.matrix(ECBAll[,9:ncol(ECBAll)])
# colnames(ECB) <- timepointsB <- c("0h","3h", "6h","2h","1d","3d","5d","7d","10d", "14d","21d", "28d")
# sizesB <- MedianNorm(ECB[rowMeans(ECB)>5,])
# NormECB <- log(GetNormalizedMat(ECB[rowMeans(ECB)>5,], sizesB)+1,2)
# 
# #Load Axoltl Embryo RNA-Seq Timecourse
# filenameH = file.path(datapath,"Human_Annotation/ExpectCounts_all_samples.txt")
# #Expected Counts All Samples
# ECEAll <- read.table(filenameH,header=TRUE)
# #Assign rownames as gene names
# rownames(ECEAll) <- ECEAll[,1]
# ECE <- as.matrix(ECEAll[-1:-2])
# sizesE <- MedianNorm(ECE[rowMeans(ECE)>5,])
# NormECE <-  log(GetNormalizedMat(ECE[rowMeans(ECE)>5,], sizesE)+1,2)
# 
# ECEtp <-  lapply(colnames(NormECE),function(x){
#   s = unlist(strsplit(x, "_"))
#   as.numeric( s[2])
# } )

##Load Blastema RNA Seq Stewart Expected counts !REMOVES LOWEST EXPRESSORS, log2
blastemaFile <- file.path(datapath, "Stewart_Gene_Expression_Across_TimecourseTPMs_0.txt")
TPMBAll <- read.csv(blastemaFile, header=TRUE , sep = "\t")
rownames(TPMBAll) <-TPMBAll$symbol
TPMB<- as.matrix(TPMBAll[,10:ncol(TPMBAll)])
colnames(TPMB) <- timepointsB <- c("0h","3h", "6h","2h","1d","3d","5d","7d","10d", "14d","21d", "28d")
sizesB <- MedianNorm(TPMB[rowMeans(TPMB)>1,])
NormTPMB <- log(GetNormalizedMat(TPMB[rowMeans(TPMB)>1,], sizesB)+1,2)


#Load Axoltl Embryo RNA-Seq Timecourse
filenameH = file.path(datapath,"Human_Annotation/TPM_all_samples.txt")
filenameE.unmapped.contigs = file.path(datapath,"Human_Annotation/1_TPM_all_samples.txt")
TPME.unmapped.contigs <- read.table(filenameE.unmapped.contigs,header=TRUE)
rownames(TPME.unmapped.contigs) <- paste0("contig", TPME.unmapped.contigs[,1])
TPME.unmapped.contigs <- as.matrix(TPME.unmapped.contigs[,c(-1,-2)])
#Expected Counts All Samples
TPMEAll <- read.table(filenameH,header=TRUE)
#Assign rownames as gene names
rownames(TPMEAll) <- TPMEAll[,1]
TPME <- as.matrix(TPMEAll[,c(-1,-2)])
TPME <- rbind(TPME, TPME.unmapped.contigs)
sizesE <- MedianNorm(TPME[rowMeans(TPME)>5,])
NormTPME <-  log(GetNormalizedMat(TPME[rowMeans(TPME)>5,], sizesE)+1,2)

filenameE.unmapped.contigs = file.path(datapath,"Human_Annotation/1_TPM_all_samples.txt")
TPME.unmapped.contigs <- read.table(filenameE.unmapped.contigs,header=TRUE)
rownames(TPME.unmapped.contigs) <- paste0("contig", TPME.unmapped.contigs[,1])
TPME.unmapped.contigs <- as.matrix(TPME.unmapped.contigs[,c(-1,-2)])
sizesE.unmapped.contigs <- MedianNorm(TPME.unmapped.contigs)
NormTPME.unmapped.contigs <-  log(GetNormalizedMat(TPME.unmapped.contigs, sizesE.unmapped.contigs)+1,2)

TPMEtp <-  lapply(colnames(NormTPME),function(x){
  s = unlist(strsplit(x, "_"))
  as.numeric( s[2])
} )



#Load axolotl embryo + blastema together and normalize !REMOVES LOW EXPRESSORS
TPMEmatch <- TPME[intersect(rownames(TPMB),rownames(TPME)),]
TPMBmatch <- TPMB[intersect(rownames(TPMB),rownames(TPME)),]
idx <- rowMads(TPMEmatch)> median(rowMads(TPMEmatch)) & rowMads(TPMBmatch)> median(rowMads(TPMBmatch))
TPMEmatch <- TPMEmatch[idx,]
TPMBmatch <- TPMBmatch[idx,]
unionnnorm <- cbind(TPMEmatch,TPMBmatch)
mnunion <-  MedianNorm(unionnnorm)
NormEBMatch <- log(GetNormalizedMat(unionnnorm,mnunion)+1,2)
#hist.normalized(NormEBMatch)
#expressionHeatmap(cor(NormEBMatch, method = "spearman"))


#Tripartite expression Knapp et all Quantile Normalized, Knapp et al.
GSE36451 <- getGEO("GSE36451",destdir=datapath, GSEMatrix =TRUE)
if (length(GSE36451) > 1) idx <- grep("GPL15342", attr(GSE36451, "names")) else idx <- 1
GSE36451 <- GSE36451[[idx]]

##Identify probes that correspond to putative regulators (blast to oncogenes, transcription factors, phosphatases and protein kinases)
chipAnnotations <- read.csv(file.path(datapath,"Amby002-Annotations-Final.tsv"), sep = "\t",header = T)
TFs <- read.table(paste0(datapath,"/my_human_TF_v5_final_11282011.txt"),header=T, sep="\t")
oncogenes <- read.table( file.path(datapath,"/MSKCCOncogenes.txt"), header = TRUE, sep = "\t")
regulators <- read.table(paste0(datapath,"/regulators.csv"),header=T, sep="\t")
fullRegulatorList <-union(union(TFs$symbol,oncogenes$Gene.Symbol),regulators$Gene)
fullRegulatorList <- fullRegulatorList[fullRegulatorList !=""]
Amby002Regulators <-  intersect(fullRegulatorList,chipAnnotations$Gene)
KnappRegulators <-  intersect(fullRegulatorList,fData(GSE36451)$GENE_SYMBOL)
KnappRegulators <- KnappRegulators[!is.null(KnappRegulators)]
index <- chipAnnotations$Gene %in% Amby002Regulators
regulatorProbes <-  chipAnnotations[index,]$Probe.Set.ID
annot <- fData(GSE36451)
index <-  annot$GENE_SYMBOL %in% KnappRegulators
KnappRegulatorProbes <- rownames(annot)[index]

#Get expression sets, exclude standardizing probes
exprsGSE37198 <- exprs(GSE37198)
exprsGSE67118 <- exprs(GSE67118)
exprsGSE35255 <- exprs(GSE35255)
index <-rownames(exprsGSE37198) %in% chipAnnotations$Probe.Set.ID
exprsGSE37198 <- exprsGSE37198[index,]
exprsGSE67118 <- exprsGSE67118[index,]
exprsGSE35255 <- exprsGSE35255[index,]
#tripartite expression set
exprsGSE36451 <- exprs(GSE36451)

allAmby002OriginalSamplenames <- list(pData(GSE67118)$source_name_ch1, pData(GSE35255)$source_name_ch1,pData(GSE37198)$source_name_ch1)
knappOriginalSamplenames <- pData(GSE36451)$source_name_ch1

#wound, deinervated, inervated
GSE37198samplenames <- unlist(lapply(pData(GSE37198)$source_name_ch1, function(x){
  wordVector <- strsplit(as.character(x),"\ ")
  paste0(wordVector[[1]][1],wordVector[[1]][2])
}))

GSE37198tp <- unlist(lapply(pData(GSE37198)$source_name_ch1, function(x){
  wordVector <- strsplit(as.character(x),"\ ")
  time <-  paste0(wordVector[[1]][4],wordVector[[1]][5])
  gsub("time|day|days","",time)
}))

#Blastema time course
GSE67118tp <- unlist(lapply(pData(GSE67118)$source_name_ch1, function(x){
  wordVector <- strsplit(as.character(x),"\ ")
  time <- paste(wordVector[[1]][13])
  gsub("T","", time)
}))

GSE67118samplenames <- unlist(lapply(pData(GSE67118)$source_name_ch1, function(x){
  wordVector <- strsplit(as.character(x),"\ ")
  paste0(wordVector[[1]][5],wordVector[[1]][7])
}))

#Wound healing terrestrial aquatic
GSE35255samplenames <- unlist(lapply(pData(GSE35255)$source_name_ch1, function(x){
  wordVector <- strsplit(as.character(x),"\ ")
  paste(wordVector[[1]][1])
}))

GSE35255tp <- unlist(lapply(pData(GSE35255)$source_name_ch1, function(x){
  wordVector <- strsplit(as.character(x),"\ ")
  time <- paste0(wordVector[[1]][3],wordVector[[1]][4])
  gsub("time|day","",time)
}))

#Tripartite expression
GSE36451samplenames <- unlist(lapply(pData(GSE36451)$source_name_ch1, function(x){
  wordVector <- strsplit(as.character(x),"\ ")
  paste(wordVector[[1]][1])
}))

GSE36451tp<- unlist(lapply(pData(GSE36451)$source_name_ch1, function(x){
  wordVector <- strsplit(as.character(x),"\ ")
  paste(wordVector[[1]][3])
}))

colnames(exprsGSE37198) <- make.unique(GSE37198samplenames)
colnames(exprsGSE67118) <- make.unique(GSE67118samplenames)
colnames(exprsGSE35255) <- make.unique(GSE35255samplenames)
colnames(exprsGSE36451) <- make.unique(GSE36451samplenames)

#Put together all 3 datasets from Voss & gang on Amby002 Platform
allAmby002Data <-  normalizeBetweenArrays(cbind(exprsGSE67118,exprsGSE35255,exprsGSE37198))
allAmby002Samplenames <- c(colnames(exprsGSE67118),colnames(exprsGSE35255),colnames(exprsGSE37198))

experiments <- list(exprsGSE35255,exprsGSE67118,exprsGSE37198,NormTPME,NormTPMB, allAmby002Data,exprsGSE36451)
names(experiments) <- c("GSE35255","GSE67118","GSE37198","EmbryoSeq","BlastemaSeq", "allAmby002","GSE36451")

write.table(cbind(make.names(unlist(allAmby002Samplenames), unique = T),as.character(unlist(allAmby002OriginalSamplenames))) , "~/SampleKey.txt",quote = F,row.names = F,col.names = F)


# experiments <- list(NormTPME)
# names(experiments) <- c("EmbryoSeq")


#write transposed datasets to files, for CLI tools
# for(es in 1:length(experiments)){
#  write.table( experiments[[es]], file =file.path(datapath ,"ExpressionSets", paste0(names(experiments)[es],".txt")), quote = F ,col.names=T, row.names = T,sep="\t")
# }

#normWound <- cbind(exprsGSE35255[,GSE35255samplenames=="Aquatic"],exprsGSE37198[,GSE37198samplenames=="Flankwound"]) - rowMeans(cbind(exprsGSE35255[,GSE35255samplenames=="Aquatic"],exprsGSE37198[,GSE37198samplenames=="Flankwound"]))
#normBlastema <- exprsGSE67118[,GSE67118samplenames!="T0"] - rowMeans(exprsGSE67118[,GSE67118samplenames!="T0"])
#normAllAmby <- cbind(exprsGSE67118,exprsGSE35255,exprsGSE37198) - rowMeans(cbind(exprsGSE67118,exprsGSE35255,exprsGSE37198))
#write.table(regulatorProbes, file =file.path(datapath,"ExpressionSets", "RegulatorProbes.txt"),sep="\n",quote=F,col.names=F,row.names = F)
#write.table(intersect(fullRegulatorList,rownames(NormTPMB) ), file =file.path(datapath,"ExpressionSets", "BlastemaSeqRegulators.txt"),sep="\n",col.names=F,row.names = F)
# write.table( t(cbind(exprsGSE35255[,GSE35255samplenames=="Aquatic"],exprsGSE37198[,GSE37198samplenames=="Flankwound"])), file =file.path(datapath,"ExpressionSets","wound.txt"),quote = F ,col.names=T,row.names = F,sep="\t")
# write.table( t(exprsGSE67118[,GSE67118samplenames!="T0"]), file =file.path(datapath,"ExpressionSets","blastemaAmby002.txt"),quote = F ,col.names=T,row.names = F,sep="\t")
# write.table( t(cbind(exprsGSE67118,exprsGSE35255,exprsGSE37198)), file =file.path(datapath,"ExpressionSets","AllAmby002.txt"),quote = F ,col.names=T,row.names = F,sep="\t")
#write.table( t(normalizeBetweenArrays(cbind(exprsGSE35255[,GSE35255samplenames=="Aquatic"],exprsGSE37198[,GSE37198samplenames=="Flankwound"]))), file =file.path(datapath,"ExpressionSets","normwound.txt"),quote = F ,col.names=T,row.names = F,sep="\t")
#write.table( t(normalizeBetweenArrays(exprsGSE67118[,GSE67118samplenames!="T0"])), file =file.path(datapath,"ExpressionSets","normblastemaAmby002.txt"),quote = F ,col.names=T,row.names = F,sep="\t")
#write.table( t(normalizeBetweenArrays(cbind(exprsGSE67118,exprsGSE35255,exprsGSE37198))), file =file.path(datapath,"ExpressionSets","normAllAmby002.txt"),quote = F ,col.names=T,row.names = F,sep="\t")
# write.table( t(log(NormTPMB+1)), file =file.path(datapath,"ExpressionSets","blastemaRNASeq.txt"),quote = F ,col.names=T,row.names = F,sep="\t")
# write.table( t(log(NormTPME+1)), file =file.path(datapath,"ExpressionSets","embryoRNASeq.txt"),quote = F ,col.names=T,row.names = F,sep="\t")
