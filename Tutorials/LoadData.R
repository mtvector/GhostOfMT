getConditionsFromColnames <- function(data,split="_",ind=1 ){
  splitList <- strsplit(x = colnames(data),split = split)
  splitMatrix <- sapply(splitList, "[", 1:max(sapply(splitList,length)) )
  return(splitMatrix[ind,])
}

rn.merge <- function(x,y,fill=0){
  rn <- intersect(rownames(x),rownames(y))
  zerosx <- setdiff(rownames(x),rownames(y))
  zerosy <- setdiff(rownames(y),rownames(x))
  out <- cbind(x[rn,],y[rn,])
  if(length(zerosx)!=1  & length(zerosy)!=1){
    zx <- matrix(fill, nrow=length(zerosx), ncol =ncol(y), dimnames = list(zerosx,NULL))
    zy <- matrix(fill, nrow=length(zerosy), ncol =ncol(x), dimnames = list(zerosy,NULL))
    zx <- cbind(x[zerosx,],zx)
    zy <- cbind(zy,y[zerosy,])
  }else if(length(zerosx)==1){
    zx <- rep(fill, ncol(y))
    zy <- matrix(fill, nrow=length(zerosy), ncol =ncol(x), dimnames = list(zerosy,NULL))
    zx <- c(x[zerosx,],zx)
    zy <- cbind(zy,y[zerosy,])
  }else if(length(zerosy)==1){
    zx <- matrix(fill, nrow=length(zerosx), ncol =ncol(y), dimnames = list(zerosx,NULL))
    zy <- rep(fill, ncol(x))
    print(zx)
    print(zy)
    zx <- cbind(x[zerosx,],zx)
    zy <- c(zy,y[zerosy,])
  }
  out <- rbind(out,rbind(zx,zy))
  return(out)
  
}

#Load in RNAseq to test for developmental genes only now
pathToFileH <- "~/code/GhostOfMT/Data/GSE90053_H1_EC.txt"
#NOTE THE PARAMETERS
countsH <- read.csv2(file=pathToFileH,header = T,row.names = 1,sep = "\t")
#countsH <- countsH[,-ncol(countsH)] #uncomment this line if the last column is gene name description text
#It removes the last column of the data from the data set
countsH <- as.matrix(countsH)
storage.mode(countsH) <- "numeric"
normCountsH <- round.log(GetNormalizedMat(countsH, MedianNorm(countsH)))
splitListH <- strsplit(x = colnames(normCountsH),split = "_d")
splitMatrixH <- sapply(splitListH, "[", 1:max(sapply(splitListH,length)) )
tpsH <- as.numeric(splitMatrixH[2,])
colnames(normCountsH) <- tpsH

pathToFileM <- "~/code/GhostOfMT/Data/GSE90053_mEpi_EC.txt"
#NOTE THE PARAMETERS
countsM <- read.csv2(file=pathToFileM,header = T,row.names = 1,sep = "\t")
#countsM <- countsM[,-ncol(countsM)] #uncomment this line if the last column is gene name description text
#It removes the last column of the data from the data set
countsM <- as.matrix(countsM)
storage.mode(countsM) <- "numeric"
#Make the genenames for mouse all uppercase.
#It will be very annoying to compare mouse and human if they are left with different conventional capitalizations
rownames(countsM) <- toupper(rownames(countsM))
normCountsM <- round.log(GetNormalizedMat(countsM, MedianNorm(countsM)))
splitListM <- strsplit(x = colnames(normCountsM),split = "_d")
splitMatrixM <- sapply(splitListM, "[", 1:max(sapply(splitListM,length)) )
tpsM <- as.numeric(splitMatrixM[2,])
colnames(normCountsM) <- tpsM


pathToFile <- "~/code/GhostOfMT/Data/GSE90053_H1_EC.txt"
#NOTE THE PARAMETERS
hTPMs <- read.csv2(file=pathToFile,header = T,row.names = 1,sep = "\t")
#hTPMs <- hTPMs[,-ncol(hTPMs)] #uncomment this line if the last column is gene name description text (or modify it if it is a different column)
#It removes the last column of the data from the data set
hTPMs <- as.matrix(hTPMs)
storage.mode(hTPMs) <- "numeric"
rownames(hTPMs) <- toupper(rownames(hTPMs))

pathToFileM <- "~/code/GhostOfMT/Data/GSE90053_mEpi_EC.txt"
#NOTE THE PARAMETERS
mTPMs <- read.csv2(file=pathToFileM,header = T,row.names = 1,sep = "\t")
#mTPMS <- mTPMs[,-ncol(mTPMs)] #uncomment this line if the last column is gene name description text (or modify it if it is a different column)
#It removes the last column of the data from the data set
mTPMs <- as.matrix(mTPMs)
storage.mode(mTPMs) <- "numeric"
rownames(mTPMs) <- toupper(rownames(mTPMs))


pathToFile <- "~/code/GhostOfMT/Data/GSE90053_H1_TPM.txt"
#NOTE THE PARAMETERS
hTPMs <- read.csv2(file=pathToFile,header = T,row.names = 1,sep = "\t")
#hTPMs <- hTPMs[,-ncol(hTPMs)] #uncomment this line if the last column is gene name description text (or modify it if it is a different column)
#It removes the last column of the data from the data set
hTPMs <- as.matrix(hTPMs)
storage.mode(hTPMs) <- "numeric"
rownames(hTPMs) <- toupper(rownames(hTPMs))

pathToFileM <- "~/code/GhostOfMT/Data/GSE90053_mEpi_TPM.txt"
#NOTE THE PARAMETERS
mTPMs <- read.csv2(file=pathToFileM,header = T,row.names = 1,sep = "\t")
#mTPMS <- mTPMs[,-ncol(mTPMs)] #uncomment this line if the last column is gene name description text (or modify it if it is a different column)
#It removes the last column of the data from the data set
mTPMs <- as.matrix(mTPMs)
storage.mode(mTPMs) <- "numeric"
rownames(mTPMs) <- toupper(rownames(mTPMs))

splitList <- strsplit(x = colnames(hTPMs),split = "_d")
splitMatrix <- sapply(splitList, "[", 1:max(sapply(splitList,length)) )
tpsH <- as.numeric(splitMatrix[2,])
conditionsH <- splitMatrix[1,]
splitList <- strsplit(x = colnames(mTPMs),split = "_d")
splitMatrix <- sapply(splitList, "[", 1:max(sapply(splitList,length)) )
tpsM <- as.numeric(splitMatrix[2,])
conditionsM <- gsub("_Sorted","",splitMatrix[1,])

tpsH <- as.numeric(getConditionsFromColnames(hTPMs,"_d",2))
conditionsH <- getConditionsFromColnames(hTPMs,"_d",1)
tpsM <- as.numeric(getConditionsFromColnames(mTPMs,"_d",2))
conditionsM <- gsub("_Sorted","",getConditionsFromColnames(mTPMs,"_d",1))

pathToList <- "~/code/GhostOfMT/Data/GeneLists/THEneuralList_2017.txt"
#NOTE THE PARAMETERS
rawList <- read.csv2(file=pathToList,header = F,sep = "\t",stringsAsFactors = F)
neuralList <- unique(toupper(rawList[,1]))

MHOrthos <-  read.table("~/code/GhostOfMT/Data/GeneLists/MouseHumanOrthologs.txt",header = F,sep = "\t")[,1:4]
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
