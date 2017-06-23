source("http://bioconductor.org/biocLite.R")
library(ggplot2)
library(RColorBrewer)
library(gplots)
#biocLite("EACI")
library(EACI)
library(matrixStats)
library(dtw)
library("amap")
# library(genefilter)
library(randomForest)
library(C50)
#library(RWeka)
library(e1071)
options("mc.cores"=32)
library(parallel)
library(doParallel)
library(rafalib)
#library(EMCluster)
registerDoParallel()
library(calibrate)
source("~/code/cb/LoadData.R")
datapath ="~/code/data/cb"
library(ebdbNet)
library(G1DBN)

getMeanExprs <- function(x){
  coef(x)[1,]
}

write.pdf <- function(plot, filename){
  pdf(filename,width = 20, height = 20)
  plot
  dev.off()
}

cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)

bhattacharyya.dist <- 
function (mu1, mu2, Sigma1, Sigma2) 
{
  aggregatesigma <- (Sigma1 + Sigma2)/2
  d1 <- mahalanobis(mu1, mu2, aggregatesigma)/8
  d2 <- log(det(as.matrix(aggregatesigma))/sqrt(det(as.matrix(Sigma1)) * 
                                                  det(as.matrix(Sigma2))))/2
  out <- d1 + d2
  out
}

polyarea <- function (x, y) 
{
  if (length(x) == 0 && length(y) == 0) 
    return(0)
  if (!(is.numeric(x) || is.complex(x)) || !(is.numeric(y) || 
                                             is.complex(y))) 
    stop("Arguments 'x' and 'y' must be real or complex.")
  if (is.null(dim(x))) 
    x <- matrix(x, length(x), 1)
  if (is.null(dim(y))) 
    y <- matrix(y, length(y), 1)
  if (any(dim(x) != dim(y))) 
    stop("Matrices 'x' and 'y' must be of same size.")
  n <- nrow(x)
  m <- ncol(x)
  z <- numeric(m)
  for (i in 1:m) {
    xi <- x[, i]
    yi <- y[, i]
    p1 <- sum(xi[1:(n - 1)] * yi[2:n]) + xi[n] * yi[1]
    p2 <- sum(xi[2:n] * yi[1:(n - 1)]) + xi[1] * yi[n]
    z[i] <- 0.5 * (p1 - p2)
  }
  return(z)
}


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

smooth.model <- function(m,tp.i){
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

mc.pval <- function(dist,null.dist,upper=F,two.tail=F){
  if(!upper & !two.tail){
    return(sapply(dist,function(x){
      (sum(x>=null.dist)+1)/(length(null.dist)+1)
    }))
  }
  if(upper & !two.tail){
    return(sapply(dist,function(x){
      (sum(x<=null.dist)+1)/(length(null.dist)+1)
    }))
  }
  if(two.tail){
    m.d = mean(dist)
    n.m.d = mean(null.dist)
    return(sapply(dist,function(x){
      (sum( abs(x-m.d) >= abs(null.dist-n.m.d))+1)/(length(null.dist)+1)
    }))
  }
}

mc.pval2 <- function(dist,null.dist,upper=F){
  if(!upper){
    return(sapply(1:length(dist),function(x){
      (sum(dist[x]>=unlist(null.dist[[x]]))+1)/(length(null.dist[[x]])+1)
    }))
  }
  if(upper){
    return(sapply(1:length(dist),function(x){
      (sum(dist[x]<=unlist(null.dist[[x]]))+1)/(length(null.dist[[x]])+1)
    }))
  }
}


t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- list(m1-m2, se, t, 2*pt(-abs(t),df,lower.tail = F))    
  names(dat) <- c("Dmeans", "Std Error", "t", "pvalue")
  return(dat) 
}

AUDTWC <- function(d,e){
  extrasy= seq(max(e$index2),max(d$index1))
  extrasx= rep(max(e$index1),length(extrasy))
  dfun = approxfun(d$index1, d$index2)
  efun = approxfun(c(e$index2,extrasy), c(e$index1,extrasx))
  list(integrate(efun, min(d$index1),max(d$index1)),integrate(dfun,min(d$index1),max(d$index1)))
}

warpArea2 <- function (d) {
  ii <- approx(x = d$index1, y = d$index2, 1:d$N)
  dg <- seq(from = 1, to = d$M, len = d$N)
  ad <- ii$y - dg
  sum(ad)
}

bhatta.dist.mat <- function(my.list, my.list.2) {
  n <- length(my.list$mu$x)
  m <- length(my.list.2$mu$x)
  mat <- matrix(0, ncol = m, nrow = n)
  
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      mat[i,j] <- bhattacharyya.dist(my.list$mu$y[i], my.list.2$mu$y[j],my.list$smooth.sd[i], my.list.2$smooth.sd[j])
    }}
  return(mat)
}

AUDTWC.dist <- function(my.list) {
  n <- length(my.list)
  mat <- matrix(0, ncol = n, nrow = n)
  colnames(mat) <- rownames(mat) <- names(my.list)
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      a = AUDTWC(my.list[[i]],my.list[[j]])
      mat[i,j] <- a[[2]]$value-a[[1]]$value
    }}
  return(as.dist(mat))
}

dtwPlotTwoWay <- function (d, xts = NULL, yts = NULL, offset = 0, ts.type = "l", 
          pch = 21, match.indices = NULL, match.col = "gray70", match.lty = 3, 
          xlab = "Index", ylab = "Query value", ...) 
{
  if (is.null(xts) || is.null(yts)) {
    xts <- d$query
    yts <- d$reference
  }
  if (is.null(xts) || is.null(yts)) 
    stop("Original timeseries are required")
  ytso <- yts + offset
  maxlen <- max(length(xts), length(ytso))
  length(xts) <- maxlen
  length(ytso) <- maxlen
  if (offset != 0) {
    par(mar = c(5, 4, 4, 4) + 0.1)
  }
  matplot(cbind(xts, ytso), type = ts.type, pch = pch, xlab = xlab, 
          ylab = ylab, axes = FALSE, ...)
  box()
  axis(1)
  axis(2, at = pretty(xts))
  if (offset != 0) {
    rightTicks <- pretty(yts)
    axis(4, at = rightTicks + offset, labels = rightTicks)
  }
  if (is.null(match.indices)) {
    ml <- length(d$index1)
    idx <- 1:ml
  }
  else if (length(match.indices) == 1) {
    idx <- seq(from = 1, to = length(d$index1), length.out = match.indices)
  }
  else {
    idx <- match.indices
  }
  segments(d$index1[idx], xts[d$index1[idx]], d$index2[idx], 
           ytso[d$index2[idx]], col = match.col, lty = match.lty)
}

point.line.dist <- function(p1,p2,x){
  ((p2[2]-p1[2])*x[1] - (p2[1]-p1[1])*x[2] + p2[1]*p1[2] - p2[2]*p1[1] )/sqrt((p2[2]-p1[2])^2 + (p2[1]-p1[1])^2)
}



mc.readtest <- sapply(seq(0,1,by=.01), function(pH){
  print(pH)
  G=3200
  Gm=2400
  Gh=2900
  N= 80000000
  O= 1000000
  pM = 1-pH
  
  trueM <- sample(sample(1:G,Gm) , N,replace=T,prob= sample(seq(0,1,by = .01), Gm,replace=T))
  trueH <- sample(sample(1:G,Gh) , N,replace=T, prob = sample(seq(0,1,by = .01), Gh,replace=T))
  obs <-  data.frame(species= sample(c("M","H"),O,T, c(pM,pH)))
  obs[obs[,"species"]=="H","gene"] <-  sample(trueH,sum(obs[,"species"]=="H"))
  obs[obs[,"species"]=="M","gene"] <-  sample(trueM,sum(obs[,"species"]=="M"))
  
  ratioMtrue <- table(trueM)/length(trueM)
  ratioM <- table(obs[obs[,"species"]=="M","gene"])/sum(obs[,"species"]=="M")
  sectM <- names(ratioMtrue)%in% names(ratioM)
  ratioHtrue <- table(trueH)/length(trueH)
  ratioH <- table(obs[obs[,"species"]=="H","gene"])/sum(obs[,"species"]=="H")
  sectH <- names(ratioHtrue)%in% names(ratioH)
  perc.ErrorM = abs(c(ratioMtrue[sectM]- ratioM, ratioMtrue[!sectM]))/c(ratioMtrue[sectM], ratioMtrue[!sectM])
  perc.ErrorH = abs(c(ratioHtrue[sectH]- ratioH, ratioHtrue[!sectH]))/c(ratioHtrue[sectH], ratioHtrue[!sectH])
  c("H.dropped"= sum(!sectH)/length(sectH),"M.dropped"= sum(!sectM)/length(sectM), "perc.ErrorH"=mean(perc.ErrorH),"perc.ErrorM"=mean(perc.ErrorM) )
})

MixSets$mm[1,103:137]
MixSets$hs[1,103:137]
mc.readtest10000 = mc.readtest
mc.readtest100000

plot(mc.readtest[1,])
plot(mc.readtest[2,])
plot(mc.readtest[3,])
plot(mc.readtest[4,])

