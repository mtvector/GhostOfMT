#install the core bioconductor packages, if not already installed
source("http://bioconductor.org/biocLite.R")

# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")
biocLite("affy")
biocLite("gcrma")
biocLite("IRanges")
biocLite("makecdfenv")
biocLite("hugene10stv1cdf")
biocLite("hugene10stv1probe")
biocLite("hugene10stprobeset.db")
biocLite("hugene10sttranscriptcluster.db")
biocLite("limma")

#Load the necessary libraries
library(makecdfenv)
library(GEOquery)
library(affy)
library(gcrma)
library(limma)
library(devtools)
library(matrixStats)
library(gplots)
library(ggplot2)
library(RColorBrewer)
#Set working directory for download
setwd("~/code/data/axo")
datapath <- "~/code/data/axo"
#GSE35255 , GSE37198 , GSE67118, GSE36451,GSE34394
#Blastema timecourse Voss et al
GSE67118 <- getGEO("GSE67118",destdir=datapath, GSEMatrix =TRUE)
if (length(GSE67118) > 1) idx <- grep("GPL15153", attr(GSE67118, "names")) else idx <- 1
GSE67118 <- GSE67118[[idx]]
#
GSE35255 <- getGEO("GSE35255",destdir=datapath, GSEMatrix =TRUE)
if (length(GSE35255) > 1) idx <- grep("GPL15153", attr(GSE35255, "names")) else idx <- 1
GSE35255 <- GSE35255[[idx]]
# Quantile Normalized Expression data from aquatic axolotl full thickness epithelial flank wounds, innervated limbs, and denervated limbs collected over seven days
GSE37198 <- getGEO("GSE37198",destdir=datapath, GSEMatrix =TRUE)
if (length(GSE37198) > 1) idx <- grep("GPL15153", attr(GSE37198, "names")) else idx <- 1
GSE37198 <- GSE37198[[idx]]

#Tripartite expression Knapp et all Quantile Normalized
GSE36451 <- getGEO("GSE36451",destdir=datapath, GSEMatrix =TRUE)
if (length(GSE36451) > 1) idx <- grep("GPL15342", attr(GSE36451, "names")) else idx <- 1
GSE36451 <- GSE36451[[idx]]
#Stewart et all blastema and mouse digit
GSE34394 <- getGEO("GSE34394",destdir=datapath, GSEMatrix =TRUE)


<<<<<<< HEAD
EBSeq::RankNorm()
nq <- cor(cbind(ECE[intersect(rownames(ECE),rownames(ECB)),],ECB[intersect(rownames(ECE),rownames(ECB)),]))
heatmap.2(nq)

qqplot(nq[,1:49], nq[,50:61])
which.max(nq)
rownames(nq) <- rownames(ECB[intersect(rownames(ECE),rownames(ECB)),])
which.max(rowMeans(nq))
abline(0,1)

featureData(GSE37198)
fData(GSE37198)
pData(GSE37198)

plotMA( exprs(GSE36451))
plotMA( exprs(GSE37198))
plotMA( exprs(GSE67118))
plotMA( exprs(GSE35255))
qqplot(exprs(GSE37198),exprs(GSE35255))

head(exprs(GSE36451))
qqplot(normalizeQuantiles(exprs(GSE36451)),exprs(GSE37198))
