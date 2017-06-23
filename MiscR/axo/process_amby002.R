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
datapath= "~/code/data/axo"
setwd(datapath)
dataset="GSE67118"
process_amby002 <-  function(dataset){
  setwd("~/code/data/axo")
  #Download the CEL file package for this dataset (by GSE - Geo series id)
  getGEOSuppFiles(dataset,makeDirectory = F)
  dir.create(dataset)
  untar( paste0(dataset,"_RAW.tar"), exdir=file.path(".",dataset))
  path <- file.path("~/code/data/axo",dataset)
  #Unpack the CEL files
  setwd(path)
  cels = list.files("./", pattern = "CEL")
  sapply(paste0("./", cels), gunzip)
  cdfiles = list.files("./", pattern = "cdf")
  gunzip(cdfiles[[1]],gsub("[.]gz$", "", cdfiles[[1]]), remove=T)
  make.cdf.package(gsub("[.]gz$", "", cdfiles[[1]]),packagename = "ambymexcdf", species = "Ambystoma_mexicanum", compress = F, unlink = T)
  install("ambymexcdf")
  library(ambymexcdf)
  tabfile= list.files(".",pattern = "probe_tab")
  gunzip(tabfile)
  probe_tab <- read.table(gsub("[.]gz$", "", tabfile), header = TRUE, sep = "\t")
  
  raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname = "ambymex") #From bioconductor
  #perform RMA normalization (I would normally use GCRMA but it did not work with this chip)
  data.rma.norm=rma(raw.data)
  rmadata <- rma(raw.data)
  gcrmadata <- gcrma(raw.data)
  #Get the important stuff out of the data - the expression estimates for each array
  rma=exprs(data.rma.norm)
  #Format values to 5 decimal places
  rma=format(rma, digits=5)
  #Map probe sets to gene symbols or other annotations
  #To see all available mappings for this platform
  
  
  #Extract probe ids, entrez symbols, and entrez ids
  probes=row.names(rma)
  rmaM <- matrix(as.numeric(rma),nrow=nrow(rma),ncol = ncol(rma))
  rownames(rmaM) <- rownames(rma)
  colnames(rmaM) <- colnames(rma)
  rmaNorm <- normalizeQuantiles(rmaM)
  ?
  GSE67118 <- getGEO('GSE67118', destdir=".")
  GSE67118 <-  GSE67118$GSE67118_series_matrix.txt.gz
  
  #GET THE TIMEPOINTS FROM LONG STRING EXPERIMENT NAMES
  timepoints <- unlist(strsplit(as.character(GSE67118$description),split = "zeugopod\ "))[seq(2,2*length(timepoints),2)]
  timepoints <- unlist(strsplit(timepoints,split = "\ post"))[seq(1,2*length(timepoints),2)]
  timepointsFactor <- factor(timepoints, levels = unique(timepoints))
 
  logrmaNorm <-  backgroundCorrect.matrix(log(rmaNorm,2))
  chipAnnotations <- read.csv(file.path("","Amby002-Annotations-Final.csv"), sep = "\t",header = T)
  rmaTotal <-  cbind(as.character(chipAnnotations$Probe.Set.ID.1),as.character(chipAnnotations$Gene), as.character(chipAnnotations[,4]), logrmaNorm[chipAnnotations$Probe.Set.ID.1,])
  rmaTotal <- rbind(c("Probe.Set.ID", "Gene", "Annotation",as.character(timepoints)), rmaTotal)
  head(rmaTotal)
  
  path="~/code/data/axo"
  #Write RMA-normalized, mapped data to file
  write.table(rmaTotal, file = file.path(path,"rmaB.txt"),quote = T, row.names = F,col.names = F, sep="\t")
  btab <-  read.table(file = file.path(path,"rmaB.txt"), header = F, sep = "\t")
  
  }
