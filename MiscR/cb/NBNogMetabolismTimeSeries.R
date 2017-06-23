library(reshape2)

seahorsePath <- "~/code/data/cb/SeahorseData/NBNogData/"
setwd(seahorsePath)
paths <- dir()
condits <- gsub(pattern = "MTS_MitoStress_|\\.txt","",paths)
species <- sapply(strsplit(condits,"_"),"[[",1)
tp <- sapply(strsplit(condits,"_"),"[[",2)
tp <- as.integer(gsub("d","",tp))


tableList <-  lapply(paths,read.csv2,sep="\t",stringsAsFactors = F)
names(tableList) <- condits

tables <- lapply(names(tableList),function(n){
  print(n)
  metabolismTab <- tableList[[n]]
  rowInds <- sort(unlist(lapply(2:5, function(x)seq(x,to = nrow(metabolismTab),by = 8))))
  if(grepl("mEpi3_",n) ){
    ocrInds <- 3:5
    ecarInds <- 11:13
    pprInds <- 19:21
    rns <-  c(NA,3,6,8,1,4,7,9,2,5,NA,10,11,NA,16,19,12,14,17,20,13,15,18,NA)[1:12]
  }else if(grepl("mEpi3rep2",n) ){
    ocrInds <- 3:8
    ecarInds <- 11:16
    pprInds <- 19:24
    rns <-  c(NA,3,6,8,1,4,7,9,2,5,NA,10,11,NA,16,19,12,14,17,20,13,15,18,NA)
  } else{
    ocrInds <- 3:8
    ecarInds <- 11:16
    pprInds <- 19:24
    rns <-  c(NA,3,6,8,1,4,7,9,2,5,NA,10,11,NA,16,19,12,14,17,20,13,15,18,NA)
  }
  
  mt <- metabolismTab[rowInds,]
  ocr <- lapply(0:11,function(i)mt[i*4+c(1,2,3,4),ocrInds])
  ocr <- sapply(ocr,unlist)
  storage.mode(ocr) <- "numeric"
  rownames(ocr) <- rns
  
  ecar <- lapply(0:11,function(i)mt[i*4+c(1,2,3,4),ecarInds])
  ecar <- sapply(ecar,unlist)
  storage.mode(ecar) <- "numeric"
  rownames(ecar) <- rns
  
  
  ppr <- lapply(0:11,function(i)mt[i*4+c(1,2,3,4),pprInds])
  ppr <- sapply(ppr,unlist)
  storage.mode(ppr) <- "numeric"
  rownames(ppr) <- rns
  noutliers <- rowMeans(ocr[,1:3])>rowMeans(ocr[,4:6]) & rowMeans(ocr[,7:9])>rowMeans(ocr[,4:6])&rowMeans(ocr[,7:9])>rowMeans(ocr[,10:12])&rowMeans(ocr[,1:3])>rowMeans(ocr[,10:12])&rowMeans(ocr[,4:6])>rowMeans(ocr[,10:12])
  #noutliers <- colMeans(cor(t(ocr),method = "spe"),na.rm = T)>.7 & rowSds(ocr)>0
  list("ocr"=ocr[noutliers,], "ecar"=ecar[noutliers,],"ppr"=ppr[noutliers,])
})
names(tables) <- condits
seahorsePath <- "~/code/data/cb/SeahorseData/NBNogDNAData/"
setwd(seahorsePath)
dnadata1 <- read.csv2(paste0(seahorsePath,"Plate1.csv"),sep = ",",header = F,stringsAsFactors = F)
dnadata1[,1] <- as.character(dnadata1[,1])
dnadata1[,2] <- as.numeric(dnadata1[,2])
dnadata2 <- read.csv2(paste0(seahorsePath,"Plate2.csv"),sep = ",",header = F,stringsAsFactors = F)
dnadata2[,1] <- as.character(dnadata2[,1])
dnadata2[,2] <- as.numeric(dnadata2[,2])
dnadata3 <- read.csv2(paste0(seahorsePath,"Plate3.csv"),sep = ",",header = F,stringsAsFactors = F)
dnadata3[,1] <- as.character(dnadata3[,1])
dnadata3[,2] <- as.numeric(dnadata3[,2])

dnadata <- rbind(dnadata1,dnadata2,dnadata3)
dnadata[,1] <- gsub("mEpi","mEpi3", dnadata[,1])
dnadata <- dnadata[dnadata[,2]>2& dnadata[,2]<15,]

ocrAll <- lapply(tables,"[[","ocr")
names(ocrAll) <- condits

ocrAllMat <- Reduce(rbind,ocrAll)
rownames(ocrAllMat) <- Reduce(c,sapply( condits, function(x)paste(rep(x,nrow(ocrAll[[x]])),rownames(ocrAll[[x]]),sep="_")))
chemAdded <- factor(rep(c("baseline","oligomycin","FCCP","rotenone/antimycin"),each=3),levels = unique(c("baseline","oligomycin","FCCP","rotenone/antimycin")))
Respirtypes <-  list(basal=c(1,4), ATP=c(1,2),ProtonLeak=c(2,4),MaxRespiration=c(3,4),NonMito=c(4,NA),RawBasline=c(1,NA))

write.table(dnadata,file = "~/Desktop/DNA_Table.txt",sep = "\t",quote = F,row.names = F,col.names = F)

#ocrAllMat <- ocrAllMat[dnadata[dnadata[,1]%in%rownames(ocrAllMat),1],]/dnadata[dnadata[,1]%in%rownames(ocrAllMat),2]

colnames(ocrAllMat) <- as.vector(chemAdded)
rownames(ocrAllMat) <- gsub("d","",rownames(ocrAllMat))

###WITHDNA 
normMat <-  sapply(names(Respirtypes),function(r){
  dat <- ocrAllMat
  conditions <- Respirtypes[[r]]
  ag <-  t(apply(dat,1,function(q)aggregate(q,by=list(chemAdded),mean)[,2]))
  if(is.na(conditions[2])){
    ag[,conditions[1]]
  }else{
    abs(ag[,conditions[1]]-ag[,conditions[2]])
  }
})
rownames(normMat) <- gsub("d","",rownames(normMat))
normMat
#normMat <- (normMat-rowMeans(normMat))/rowSds(normMat)
normMat <- normMat[dnadata[dnadata[,1]%in%rownames(normMat),1],]/dnadata[dnadata[,1]%in%rownames(normMat),2]
#meltNormMat <- melt(normMat)
ocrAllMatA <- ocrAllMat[dnadata[dnadata[,1]%in%rownames(ocrAllMat),1],]/dnadata[dnadata[,1]%in%rownames(ocrAllMat),2]
meltNormMat <- melt(ocrAllMatA)
meltNormMat$Var1<- as.character(meltNormMat$Var1)
meltNormMat <- cbind(meltNormMat, species= gsub("rep2","",sapply(strsplit(x = meltNormMat$Var1,split ="_" ),"[",1 )))
meltNormMat <- cbind(meltNormMat, tp= as.numeric(gsub("d","",sapply(strsplit(x = meltNormMat$Var1,split ="_" ),"[",2 ))))
#lmList(value~species|tp,meltNormMat)
#summary(lmList(value~species|tp,meltNormMat[meltNormMat$Var2=="basal",]))
tip="baseline"
normMat.sub <- meltNormMat[meltNormMat$Var2==tip,]
splitforstrip <- split(normMat.sub$value,paste(normMat.sub$species,normMat.sub$tp ))
splitforstrip <- splitforstrip[c(1,6,7,8,3,4,5,11,9,12,13,10,14,15)]
stripchart(splitforstrip, vertical=TRUE, pch=1, method="jitter", las=2,ylab = "(mpH/min)/(ngDNA/uL)",main=paste(tip,"Acidification"))
stripchart(splitforstrip, vertical=TRUE, pch=1, method="jitter", las=2,ylab = "(ppMolO2/min)/(ngDNA/uL)",main=paste(tip,"OCR"))

plot(dnadata[dnadata[,1]%in%rownames(ocrAllMat),2],ocrAllMat[dnadata[dnadata[,1]%in%rownames(ocrAllMat),1],1],ylab = "Unnormalized ECAR",xlab = "ng/uL DNA")


l <- lm(value~species,normMat.sub[normMat.sub$tp=="3",])
summary(l)




### WITHOUT DNA
ocrAllMatMelt <- melt(ocrAllMat)
ocrAllMatMelt$Var1<- as.character(ocrAllMatMelt$Var1)
ocrAllMatMelt <- cbind(ocrAllMatMelt, species= sapply(strsplit(x = ocrAllMatMelt$Var1,split ="_" ),"[",1 ))
ocrAllMatMelt <- cbind(ocrAllMatMelt, tp= as.numeric(gsub("d","",sapply(strsplit(x = ocrAllMatMelt$Var1,split ="_" ),"[",2 ))))

tip="baseline"
ocrAllMatMelt.sub <- ocrAllMatMelt[ocrAllMatMelt$Var2==tip,]
stripchart(split(ocrAllMatMelt.sub$value,paste(ocrAllMatMelt.sub$species,ocrAllMatMelt.sub$tp )), vertical=TRUE, pch=1, method="jitter", las=2,ylab = "(mpH/min)/(ngDNA/uL)",main=paste(tip,"OCR"))
l <- lm(value~species,ocrAllMatMelt.sub[ocrAllMatMelt.sub$tp=="3",])
summary(l)


RespirTables <-  lapply(tables, function(x){
  dat <- x$ocr
  colnames(dat) <- chemAdded
  dat <- (dat-rowMeans(dat))/rowSds(dat)
  sapply(names(Respirtypes),function(r){
    conditions <- Respirtypes[[r]]
    ag <-  t(apply(dat,1,function(q)aggregate(q,by=list(chemAdded),mean)[,2]))
    if(is.na(conditions[2])){
      ag[,conditions[1]]
    }else{
      abs(ag[,conditions[1]]-ag[,conditions[2]])
    }
  })
})
names(RespirTables) <- condits
RespirMat <- Reduce(rbind,RespirTables)
rownames(RespirMat) <- Reduce(c,sapply( condits, function(x)paste(rep(x,nrow(RespirTables[[x]])),rownames(RespirTables[[x]]),sep="_")))
colnames(RespirMat) <- names(Respirtypes)
RespirMatMelt <- melt(RespirMat)
RespirMatMelt$Var1<- as.character(RespirMatMelt$Var1)
RespirMatMelt <- cbind(RespirMatMelt, species= sapply(strsplit(x = RespirMatMelt$Var1,split ="_" ),"[",1 ))
RespirMatMelt <- cbind(RespirMatMelt, tp= as.numeric(gsub("d","",sapply(strsplit(x = RespirMatMelt$Var1,split ="_" ),"[",2 ))))

tip="basal"
RespirMatMelt.sub <- RespirMatMelt[RespirMatMelt$Var2==tip,]
stripchart(split(RespirMatMelt.sub$value,paste(RespirMatMelt.sub$species,RespirMatMelt.sub$tp )), vertical=TRUE, pch=1, method="jitter", las=2,ylab = "(mpH/min)/(ngDNA/uL)",main=paste(tip,"OCR"))
l <- lm(value~species,RespirMatMelt.sub[RespirMatMelt.sub$tp=="3",])
summary(l)



pdf(file = "~/Desktop/metabolismTimeGraphsNormalized.pdf",width = 10,height = 10)
plts <- lapply(unique(meltNormMat$Var2),function(i){
  subPlot <- meltNormMat[meltNormMat[,"Var2"]==i,]
  plt <- ggplot(subPlot)+
    geom_point(aes(x=tp, y=value,color=species),show.legend = T) +
    #geom_text(aes(x=tp, y=value,label=Var1))+
    #geom_line(aes(x=tp, y=value,color=species),show.legend = T) +
    #geom_errorbar(aes(x=tp,ymin=value-se, ymax=value+se), width=.03) +
    labs(title = paste0(i," Metabolic Trends"),x="day", y = "Zscore Value") +
    theme_classic(base_size = 16)
  print(plt)
} )
dev.off()




zscoreTables <-  lapply(tables, function(x){
  dat <- x$ocr
  colnames(dat) <- chemAdded
  dat <- (dat-rowMeans(dat))/rowSds(dat)
  sapply(names(Respirtypes),function(r){
    conditions <- Respirtypes[[r]]
    ag <-  t(apply(dat,1,function(q)aggregate(q,by=list(chemAdded),mean)[,2]))
    if(is.na(conditions[2])){
      ag[,conditions[1]]
    }else{
      abs(ag[,conditions[1]]-ag[,conditions[2]])
    }
  })
})

std <- function(x) sd(x)/sqrt(length(x))

forPlotting <- t(sapply(zscoreTables,function(x) apply(x,2,mean)))
forPlotting <-cbind(forPlotting,"tp"=tp)
storage.mode(forPlotting) <- "numeric"
forPlotting <-as.data.frame(forPlotting)
forPlotting <-cbind(forPlotting,"species"=species)
forPlotting <-cbind(forPlotting,"condits"=condits)

forPlottingSE <- t(sapply(zscoreTables,function(x) apply(x,2,std)))
forPlottingSE <-cbind(forPlottingSE,"tp"=tp)
storage.mode(forPlottingSE) <- "numeric"
forPlottingSE <-as.data.frame(forPlottingSE)
forPlottingSE <-cbind(forPlottingSE,"species"=species)
forPlottingSE <-cbind(forPlottingSE,"condits"=condits)

pdf(file = "~/Desktop/metabolismTimeGraphs.pdf",width = 10,height = 10)
plts <- lapply(1:length(Respirtypes),function(i){
  subPlot <- melt(forPlotting,id.vars = c("tp","species","condits"))
  subPlot <- subPlot[subPlot[,"variable"]==names(Respirtypes)[i],]
  subPlotSE <- melt(forPlottingSE,id.vars = c("tp","species","condits"))
  subPlotSE <- subPlotSE[subPlotSE[,"variable"]==names(Respirtypes)[i],]
  subPlot <- cbind(subPlot,se=subPlotSE$value)
  plt <- ggplot(subPlot)+
    geom_point(aes(x=tp, y=value,color=species),show.legend = T) +
    geom_text(aes(x=tp, y=value,label=condits))+
    geom_line(aes(x=tp, y=value,color=species),show.legend = T) +
    geom_errorbar(aes(x=tp,ymin=value-se, ymax=value+se), width=.03) +
    labs(title = paste0(names(Respirtypes)[i]," Metabolic Trends"),x="day", y = "Zscore Value") +
    theme_classic(base_size = 16)
  print(plt)
} )
dev.off()

pdf(file = "~/Desktop/metabolismRoughGraphs.pdf")
mypar(2,2)
for(i in 1:14){
  matplot(t( tables[[i]]$ocr), ylab=" (pmol O2/min)",type = "l",main=paste(condits[i],"OCR"))
  matplot(t( (tables[[i]]$ocr-rowMeans(tables[[i]]$ocr))/(rowSds(tables[[i]]$ocr))  ),type = "o",main=paste(condits[i],"OCR ZScore"))
  matplot(t( tables[[i]]$ecar), ylab="(mpH/min)",type = "l",main=paste(condits[i],"ECAR"))
  matplot(t( (tables[[i]]$ecar-rowMeans(tables[[i]]$ecar))/(rowSds(tables[[i]]$ecar))  ),type = "o",main=paste(condits[i],"ECAR ZScore"))
}
dev.off()



