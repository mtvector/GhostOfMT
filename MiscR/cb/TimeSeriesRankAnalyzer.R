library(zoo)
library(EBSeq)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(matrixStats)
source("~/code/cb/LoadFunctions.R")
source("~/code/cb/LoadData.R")

setdiff(rownames(datasets[[11]]),rownames(datasets[[14]]))
setdiff(rownames(datasets[[14]]),rownames(datasets[[11]]))[999:2000]



getMeanExprs <- function(x){
  coef(x)[1,]
}

write.pdf <- function(plot, filename){
  pdf(filename,width = 20, height = 20)
  plot
  dev.off()
}

cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)

#Takes two model fits and calculates their bhattacharya distance, as well as a sum of variance
dist.fit <- function(f1,f2){
  cf1 = coef(f1)
  cf2= coef(f2)
  cv1=1.96*sqrt(diag(vcov(f1)))
  cv2=1.96*sqrt(diag(vcov(f2)))
  #anova(f1,f2)
  x=ncol(cf1)
  y=ncol(cf2)
  d <- matrix(nrow = x, ncol = y)
  v <- matrix(nrow = x, ncol = y)
  for(i in 1:x){
    for(j in 1:y){
      if(cf1[1,i]>cf2[1,j]){
        d[i,j] <-dist(cf1[1,i]+cv1[i], cf2[1,j]-cv2[j])
      }else{
        d[i,j] <- dist(cf1[1,i]-cv1[i], cf2[1,j]+cv2[j])
      }
      v[i,j] <- cv1[i]+cv2[j]
    }
  }
  return(list(d,v))
}

sd.centered <-function(x,mu) sqrt(sum((x -mu)^2) / (length(x)))

smooth.model <- function(m,tp.i=NULL){
  if(is.null(tp.i)){
    tp.i <- seq(min(m$tme), max(m$tme))
  }
  a=anova(lm(formula=y~as.character(tme), data=m))
  ss=smooth.spline(m$tme, m$y,keep.data = F)
  mu=predict(ss,tp.i)
  sd = sapply(sort(unique(m$tme)),function(x){
    sd.centered(as.numeric(m[m$tme==x,]$y), mu$y[which(mu$x==x)])
  })
  names(sd) <- sort(unique(m$tme))[sort(unique(m$tme))%in% tp.i]
  sd = sd[!is.na(names(sd))]
  ind = tp.i%in%names(sd)
  smooth.sd <- rep(mean(sd),length(tp.i))
  names(smooth.sd) <- tp.i
  names(mu$y) <- mu$x
  smooth.sd[names(sd)] <- sd
  list("mu"=mu , "anova"=a,"data"= m,"sd"= sd,"smooth.sd"=smooth.sd,"spline"=ss)
}

range01 <- function(x){ (x-min(x))/(max(x)-min(x))}

getF.P <- function(m){
  sapply(m,function(x){
    x$anova$`Pr(>F)`
  })
}


plot.model <- function(f,show.conifdence.bands =T,...){
  plot(c(f$mu$x,f$data$tme ) ,c(f$mu$y,f$data$y),...)
  if(show.conifdence.bands){
    upper.band <- spline(x = as.numeric(f$mu$x), 
                         y = f$mu$y + 1.96 * f$smooth.sd, 
                         method = "natural", n = length(f$mu$x)/2)
    lower.band <- spline(x = as.numeric(f$mu$x), 
                         y = f$mu$y - 1.96 * f$smooth.sd, 
                         method = "natural", n =  length(f$mu$x)/2)
    
    col.meanCurve.rgb <- col2rgb("red")
    polygon(x = c(upper.band$x, rev(lower.band$x)), y = c(upper.band$y, 
                                                          rev(lower.band$y)), col = rgb(col.meanCurve.rgb[1], 
                                                                                        col.meanCurve.rgb[2], col.meanCurve.rgb[3], alpha = 125, 
                                                                                        maxColorValue = 255), border = NA)
  }
}


sdThresh <- 1
tp <- tpDatasets[inds]
tp.i <- Reduce(intersect, tp)
set <- datasets[inds]
genes <- lapply(set,rownames)
set <- lapply(set,round.log)
allg <- Reduce(intersect, genes)
set <- lapply(set, md, allg)
allg <-allg[Reduce("|",lapply(lapply(set, rowSds), ">",sdThresh))]
set <- lapply(set, md, allg)

frame <- lapply(1:length(set), function(x){
  q <-  lapply(1:nrow(set[[x]]), function(k){
      f<- t(sapply(1:ncol(set[[x]]),function(j){
        c("y"=set[[x]][k,j], "tme"=tp[[x]][j])
      }))
      as.data.frame(f)
  })
  names(q) <- allg
  q
})

fits <- lapply(frame, function(x) lapply(x, smooth.model,tp.i))
nGenes <- lapply(fits,nrow)
geneCors <- lapply(allg,function(g){
  as.dist(cor(sapply(fits, function(x){
    x[[g]]$mu$y
  })))
})
names(geneCors) <- allg

geneDists <- lapply(allg,function(g){
  dist(t(sapply(fits, function(x){
    x[[g]]$mu$y
  })))
})
names(geneDists) <- allg

sort(unlist(lapply(geneCors, which.min)))

mypar(1,3)
for(i in 1:3){
  plot.model(fits[[i]][["HOXB4"]])
}

geneDists[["HOXB4"]]

TOMs <- WGCNA::blockwiseIndividualTOMs(as.vector(lapply(set,function(x) list(data=t(x))),mode = "list"),checkMissingData = T,checkPower = T,saveTOMs = F)
set[[1]][which(is.na(TOMs$blocks)),]
dim(TOMs$TOMSimilarities[[1]])
WGCNA::


c <- cor.compare(a,b,method="spea")
c <- range01(c)
c[c<median(c)] <- 0
std.heatmap(c)

x.cm <- sapply(1:ncol(a), function(i){
  (1/sum(c[i,]))*sum(sapply(1:ncol(b),function(j) c[i,j]*tpb[j]))
})

y.cm <- sapply(1:ncol(b), function(i){
  (1/sum(c[,i]))*sum(sapply(1:ncol(a),function(j) c[j,i]*tpa[j]))
})
plot(tpa,x.cm)
plot(tpb,y.cm)


ar <- apply(a,2,rank,ties.method="min")
plot(tpa,a["TERT",])
plot(tpa,a["CHRNA5",])
TRRUST.interactions[TRRUST.interactions$to=="TERT",]
matplot(t(rbind(a[TRRUST.interactions[TRRUST.interactions$to=="TERT",1],], a["TERT",] )),type="l")

cbind(cor(t(a[TRRUST.interactions[TRRUST.interactions$to=="TERT",1],]), a["TERT",] ),TRRUST.interactions$type)



hist(rowSds(ar))
hist(rowSds(ar[neural_related_genes[neural_related_genes%in%rownames(ar)],]))

c <- cor(ar)
c <- cor(a,method = "spe")

#Maybe this idea is ill-guided...



