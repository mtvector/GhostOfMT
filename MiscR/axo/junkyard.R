
#Load Data File for Zebrafish Alignments
filename = "/Users/mschmitz/Documents/Code/ph525/axo/Zebrafish_Annotation/AVE_TPM_each_stage.txt"
#TPM for each stage
TPMZAll <- read.table(filename,header=TRUE)
#Assign rownames as gene names
rownames(TPMZAll) <- TPMZAll[,1]
TPMZFrame <- TPMZAll[-1:-2]
TPMZ <- as.matrix(TPMZFrame)


EBOutTPMZ <- sapply(1:(ncol(TPMZ)-1), function(i){
  condit <- factor(colnames(TPMZAll[(i):(i+1)]))  
  subTPMZ <- TPMZ[,(i):(i+1)]
  Sizes=MedianNorm(subTPMZ)
  EBOut <- EBTest(subTPMZ, Conditions = condit, sizeFactors = Sizes, maxround=5)
  return(EBOut)
})

EBOutResultsZ <- sapply(1:(ncol(TPMZ)-1), function(i){
  EBDERes=GetDEResults(EBOutTPMZ[,i], FDR=0.05)
})




##Returns matrix of [timepoint in EB, timepoint in blastema it is DE] = number of genes from EB that are expressed at blastema timepoint
getDEStageShift <- function(EB, testlist){
  sMat <- matrix(0,nrow=ncol(EB), ncol=length(testlist))
  for(i in 1:ncol(EB)){
    EBGenes <-  EB[,i]$DEfound
    for(j in 1:length(testlist)){
      for(k in 1:length(testlist[[j]])){
        if(testlist[[j]][k]%in%EBGenes){
          sMat[i,j] <- sMat[i,j] + 1
        }
      }
    }
  } 
  return(sMat)
}

getDETotalsMatrix <- function(EB, testlist){
  sMat <- matrix(0,nrow=ncol(EB), ncol=length(testlist))
  for(i in 1:ncol(EB)){
    EBGenes <-  EB[,i]$DEfound
    for(j in 1:length(testlist)){
      for(k in 1:length(testlist[[j]])){
        sMat[i,j] <- sMat[i,j] + 1
      }
    }
  } 
  return(sMat)
}




DEPropH <- function(EB, testlist, matchs){
  sMat <- matrix(0,nrow=ncol(EB), ncol=length(testlist))
  for(i in 1:ncol(EB)){
    EBGenes <-  length(EB[,i]$DEfound)
    for(j in 1:ncol(matchs)){
      sMat[i,j] <- matchs[i,j]/EBGenes
      if(EBGenes==0){
        sMat[i,j] <- 0
      }
    }
  }
  return(sMat)
}


?GetDEResults



DEMatchs <- getDEMatchMatrix(EBOutResultsH, getListGenes(upRegs))
DETotals <- getDETotalsMatrix(EBOutResultsH, getListGenes(upRegs))
DEMatchList <- getDEMatchListMatrix(EBOutResultsH, getListGenes(upRegs))
DEProp <- DEMatchs/DETotals
DEProp_H <- DEPropH(EBOutResultsH, getListGenes(upRegs), DEMatchs)

cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)
colnames (DEProp) <- timepointsB
rownames(DEProp) <- timepointsH
heatmap.2(DEProp, xlab = "Blastema Time Point", ylab="Embryo Time Point", main="Differential Expression Matches\n (Proportion of total blastemal DE)"
          ,na.rm = TRUE, trace = "none",colsep=1:12,
          rowsep=1:16,
          col=cols)

cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)
colnames (DEProp_H) <- timepointsB
rownames(DEProp_H) <- timepointsH
heatmap.2(DEProp_H[,1:11], xlab = "Blastema Time Point", ylab="Embryo Time Point", main="Differential Expression Matches\n (Proportion of total embryo DE)"
          ,na.rm = TRUE, trace = "none",colsep=1:12,
          rowsep=1:16,
          col=cols)

colnames (DEMatchs) <- timepointsB
rownames(DEMatchs) <- timepointsH
heatmap.2(DEMatchs[,1:11], xlab = "Blastema Time Point", ylab="Embryo Time Point", main="# of DE Gene Matches"
          ,na.rm = TRUE, trace = "none",colsep=1:12,
          rowsep=1:16,
          col=cols)



##Read in files of upregulated genes from the paper
setwd("/Users/mschmitz/Documents/Code/ph525/axo/juvTC_UP_ALL_FDR_LT_0.05")
files <-list.files()
filelist = list.files(pattern = ".*.txt")
upRegs <-lapply(filelist, read.table, sep="\t", header=TRUE)

setwd("/Users/mschmitz/Documents/Code/ph525/axo/juvTC_DOWN_ALL_FDR_LT_0.05")
files <-list.files()
filelist = list.files(pattern = ".*.txt")
downRegs <- lapply(filelist, read.table, sep="\t", header=TRUE)

setwd("/Users/mschmitz")


##########PLOTTING################FAILED PLOTTING#########################



multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#   plotH <- ggplot(gH[1,],aes(x=factor(colnames(gH), levels =colnames(gH)),y=gH[1,]), environment = .e) + 
#     geom_point(size = 4,stat = "identity",y=gH[1,])+ theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#     xlab("TimePoint") +
#     ylab("TPM") +
#     ggtitle(paste("Embryo", rownames(gH)))
#   
#   plotB <- ggplot(gB[1,],aes(x=factor(colnames(gB), levels =colnames(gB)),y=gB[1,]), environment = .e) + 
#     geom_point(size = 4,stat = "identity",y=gB[1,])+ theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#     xlab("TimePoint") +
#     ylab("TPM") +
#     ggtitle(paste("Blastema", rownames(gB)))
#   return(multiplot(plotH, plotB, cols=2))

################## CURVE FITTING##################################
library(splines)
library(stats)
mypar(1,2)
fit <- loess(formula(TPMZ["aven",]~seq(1,ncol(TPMZ))), degree = 2)
plot(TPMZ["aven",], main="AVEN Expression")
fit <- loess(formula(TPMZ["aven",]~seq(1,ncol(TPMZ))), degree = 2)
fitted=predict(fit,newdata=data.frame(seq(1,ncol(TPMZ))))
lines(fitted, col="red")
fit2= smooth.spline(TPMZ["aven",], df=7)
fitted2=predict(fit2)
lines(fitted2, col="dodgerblue")
plot(TPMBMatrix["AVEN",])
fit3= smooth.spline(TPMBMatrix["AVEN",], df=7)
fitted3=predict(fit3)
lines(fitted3, col="orange")

indexZH <- 7:17
indexB <- 2:12

splinesZ <-  sapply(1:nrow(TPMZ), function(i){
  fit2= smooth.spline(TPMZ[i,indexZH],df=7)
  fitted2=predict(fit2)
})

splinesH <-  sapply(1:nrow(TPMH), function(i){
  fit2= smooth.spline(TPMH[i,indexZH],df=7)
  fitted2=predict(fit2)
})

splinesB <-  sapply(1:nrow(TPMBMatrix), function(i){
  fit2= smooth.spline(TPMBMatrix[i,indexB], df=7)
  fitted2=predict(fit2)
})


IntegrateCoords <- function(spln){
  Q <- c()
  Q <- sapply( 1:length(spln$x), function(i){
    Q <- c(Q, c(spln$x[i], spln$y[i]))
  })
  return(as.matrix(t(Q)))
}
#IntegrateCoords(splinesZ[,1])


splineZCoords <- sapply( 1:ncol(splinesZ),function(i){
  IntegrateCoords(splinesZ[,i])
}, simplify = 'array' )
splineZCoords

A <- digit3.dat[,,1]
B <- digit3.dat[,,5]
A <- splineZCoords[,,17]
B <- splineBCoords[,,21]
procOPA(A,B, scale = TRUE, reflect = FALSE)
plot(splineZCoords[7:17,,17])
plot(splineBCoords[2:12,,21])
plot(TPMB[21,])

test <- sapply(1:nrow(TPMB), function(i){
  A <- splineHCoords[,,i]
  B <- splineBCoords[,, which(rownames(TPMB)==TPMHFrame[i,1])]
  return(c(dim(A)[1] !=dim(B)[1],dim(A)[2] !=dim(B)[2]))
})
test[2,]
sum(test[1,])

A <- splineHCoords[,,4]
B <- splineBCoords[,, which(rownames(TPMB)==TPMHFrame[4,1])]


OPAHout <- sapply(1:nrow(TPMB), function(i){
  A <- splineHCoords[,,i]
  B <- splineBCoords[,, which(rownames(TPMB)==TPMHFrame[i,1])]
  if(length(B)>0){
    return(procOPA(A,B, scale = FALSE, reflect = FALSE))
  }else{
    NULL
  }
})

PQPQ <- 23
plot(splineHCoords[,,PQPQ])
plot(splineBCoords[,, which(rownames(TPMB)==TPMHFrame[PQPQ,1])])
splineBCoords[,, which(rownames(TPMB)==TPMHFrame[PQPQ,1])]

passFilter <- sapply(1:length(OPAHout), function(i){
  if(!is.null(OPAHout[[i]])){
    riemdist(OPAHout[[i]]$Ahat,OPAHout[[i]]$Bhat) > 1.5
  }else{
    FALSE
  }
})
#riemdist(OPAHout[[2]]$Ahat,OPAHout[[2]]$Bhat)
sum(unlist(passFilter))
ww <- which(unlist(passFilter))
mypar(6,2)
for(i in 1:6){
  x <- ww[i]
  plot(splineHCoords[,,x], main= TPMHFrame[x,1])
  plot(splineBCoords[,, which(rownames(TPMB)==TPMHFrame[x,1])], main= TPMHFrame[x,1])
}




#############JUNKYARD#############
install.packages("devtools")
library(devtools)

numericGenes <- as.numeric(as.matrix(TPMZ[-1:-2]))
biocLite("RamiGO")
library(RamiGO)
list2score(EBOutResultsZ[,1]$DEfound, lib="GO.db")
ls("package:GO.db")
install.packages("dtw")
library("dtw")

#Sizes=MedianNorm(TPMZ[,1:3])
#subTPMZ=TPMZ[,1:3]
#head(TPMZ)
#dim(TPMZ)
#EBOut <- EBMultiTest(subTPMZ, Conditions = factor(colnames(subTPMZ)),sizeFactors = Sizes, maxround=5)

?"TimeShift"
help(shift)
data("PrimateBrain")
typeof(ages.h)
typeof(ages.c)
typeof(n.age.h)
length(n.age.c)
head(human)
head(chimp)
head(ages.h)
head(ages.c)
head(n.age.h)
head(n.age.h)
sh = shift(human, chimp, ages.h, ages.c, n.age.h, n.age.c)
summary(sh)
?cbind



regressBlastema <- sapply(1:nrow(TPMBMatrix), function(i){
  fit <- loess(formula(TPMBMatrix[i,]~seq(1,ncol(TPMBMatrix))), degree = 2)
  fitted=predict(fit,newdata=data.frame(seq(1,ncol(TPMBMatrix))))
})

regressZ <- sapply(1:nrow(TPMZ), function(i){
  fit <- loess(formula(TPMZ[i,]~seq(1,ncol(TPMZ))), degree = 2)
  fitted=predict(fit,newdata=data.frame(seq(1,ncol(TPMZ))))
})