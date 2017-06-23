source("http://bioconductor.org/biocLite.R")
library(Biobase)
library(ggplot2)
library(RColorBrewer)
library(gplots)
#biocLite("EACI")
#library(EACI)
library(matrixStats)
library(dtw)
library("amap")
# library(genefilter)
#library(randomForest)
#library(C50)
#library(RWeka)
#library(e1071)
options("mc.cores"=32)
library(parallel)
library(rafalib)
#library(EMCluster)
library(calibrate)
library(clusterProfiler)
library(biomaRt)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

source("~/code/cb/LoadFunctions.R")
source("~/code/cb/LoadData.R")
datapath ="~/code/data/cb"

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


########
#set1 = datasets[c(1,2)]
#set1 =  list(exp_1[,17:24],exp_2)
#set2 = datasets[c(3,4)]
#set1=test.setM
#set2=test.setH
#tp1 = tpDatasets[c(1,2)]
#tp2 =tpDatasets[c(3,4)]
#output.dir = "~/code/data/cb/shiftFiles/MixvsHInVitrorj7zequalized"
#gene.list = neural_related_genes
# null.length = 1
# cor.Threshold = 0
#fc.threshold = 2
#tpm.threshold = 2
#threshold.num =3
#name1="Mouse"
#name2="Human"
# #tp2=tp.h.control
# tp1.i =seq(0,max(tp_2))
# tp2.i =seq(0,max(tp_2))
#tp2.i=seq(0,max(tp.h.control))
# open.start = F
# stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T)
# z = F
######
####DT output from warp is the direction and magnitude of shift for set2 at each tp
#Larger aad.Thresh means more genes get in
#Set1 should be the longer of the two sequences. Called Experimental, query by dtw. "Forward shift" is set1 shift positively

run.DTW.Genes <-  function(set1, set2, output.dir , tp1, tp2,name1,name2,n.null=200 , gene.list=neural_list ,open.start=F, power=1, threshold.num=1,fc.threshold= 2, equalize = F, tpm.threshold=2,ppt=F, derivative = 0,bhatta=F , stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T), z =F){
  pairs <-  c(2,1,4,3,6,5)
  cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
  all.rn <- as.character(Reduce(intersect, c(lapply(set1,rownames),lapply(set2, rownames))))
  set1 <- lapply(set1, function(n) n[all.rn,])
  set2 <- lapply(set2, function(n) n[all.rn,])
  s1.mx <-  lapply(set1, function(x) rowMaxs(x) )
  s2.mx <-  lapply(set2, function(x) rowMaxs(x) )
  pass.thresh <-  Reduce( "+" ,lapply(c(s1.mx,s2.mx), ">",tpm.threshold)) >= sum(length(set1),length(set2))
  set1 <- lapply(set1, function(n) n[pass.thresh,])
  set2 <- lapply(set2, function(n) n[pass.thresh,])
  
  tp <- list(tp1,tp2)
  set <- list(set1,set2)
  frame <- lapply(1:length(set), function(x){
    q <-  lapply(1:nrow(set[[x]][[1]]), function(k){
      as.data.frame(Reduce(rbind,lapply(1:length(set[[x]]),function(i){
        f<- t(sapply(1:ncol(set[[x]][[i]]),function(j){
          c("y"=set[[x]][[i]][k,j], "tme"=tp[[x]][[i]][j])#,"ind"=i)
        }))
        #if(z) f[,"y"] <- f[,"y"]/max(f[,"y"])
        f
      })))
    })
    names(q) <- rownames(set[[1]][[1]])
    q
  })
#   all.exprs <- lapply(1:length(frame),function(s){
#     unlist(sapply(frame[[s]],"[","y"))
#   })
  tp.i=NULL
  if(equalize){
    output.dir <- paste0(output.dir,"equalized")
    tp.i = seq(max(min(frame[[1]][[1]]$tme),min(frame[[2]][[1]]$tme)),min(max(frame[[1]][[1]]$tme), max(frame[[2]][[1]]$tme)))
  }

  dir.create(output.dir)
  dir.create(paste0(output.dir,"/diagnostics"))
  if(file.exists(file.path(output.dir,"fits.RData"))){
    load(file.path(output.dir, "fits.RData"))
  }else{
    fits1 <- mclapply(frame[[1]], smooth.model, tp.i)
    fits2 <- mclapply(frame[[2]], smooth.model,tp.i)
    fits <- list(fits1, fits2)
    save(fits,name1,name2,file=file.path(output.dir,"fits.RData"))
  }
  gene.list <- gene.list[(gene.list%in%names(fits[[1]]))]

  rancors <- sapply(1:length(fits[[1]]), function(x){
    ind <- intersect(fits[[1]][[x]]$mu$x,fits[[pairs[1]]][[x]]$mu$x)
    a <- fits[[1]][[x]]$mu$y[ind]
    b <- fits[[pairs[1]]][[x]]$mu$y[ind]
    cor(a,b  ,method="spearman")
  })
  names(rancors) <- names(fits[[1]])
  
  f.values <- lapply(1:length(fits),function(s){
    ft <- sapply(1:length(fits[[s]]), function(x){
      fits[[s]][[x]]$anova$`Pr(>F)`[1]
    })
    names(ft) <- names(fits[[1]])
    ft
  })

  xcors <- sapply(1:length(fits[[1]]), function(x){
    max(ccf(fits[[1]][[x]]$mu$y,fits[[pairs[1]]][[x]]$mu$y,plot=F)$acf)
  })
  names(xcors) <- names(fits[[1]])

  w.xcors <- sapply(1:length(fits[[1]]), function(x){
    a=ccf(fits[[1]][[x]]$mu$y,fits[[pairs[1]]][[x]]$mu$y,plot=F)$acf
    which.max(a)-(as.integer(length(a)/2)+1)
  })
  names(w.xcors) <- names(fits[[1]])
  
  dt.fc <- lapply(1:length(fits[[1]]),function(i){
    gn <- names(fits[[1]])[i]
    a=fits[[1]][[gn]]$mu$y
    b=fits[[2]][[gn]]$mu$y
    min((max(a)-min(a)), (max(b)-min(b)))
  })
  f.ind <-  f.values[[1]] < sqrt(.25) & f.values[[2]] < sqrt(.25)
  f.ind[is.na(f.ind)] <- T
  print("pass F test:")
  print(sum(f.values[[1]] < sqrt(.05)))
  print(sum(f.values[[2]] < sqrt(.05)))
  print(sum(dt.fc >fc.threshold))
  fitind <- dt.fc >fc.threshold &  f.ind
  fits[[1]] <- fits[[1]][fitind]
  fits[[2]] <- fits[[2]][fitind]
  print("# null")
  print(n.null)
  print(format(Sys.time(),"%D %H:%M:%S"))
  if(file.exists(file.path(output.dir,"dtw.RData"))){
    load(file.path(output.dir,"dtw.RData"))
  }else{
    print("DTWing")
    dt.dist <- list()
    dt.likelihood <- list()
    dt.aligned.dist <- list()
    sumErrors <- list()
    dt.full <- list()
    dt.objects <- list()
    null.full <- list()
    null <- list()
    null.area <- list()
    for(s in 1:length(fits)){
      print(s)
      dt.dist[[s]] <- vector()
      dt.objects[[s]] <- list()
      null[[s]] <- list()
      dt.full[[s]] <- matrix(nrow = length(fits[[s]]), ncol=length(fits[[s]][[1]]$mu$x),dimnames = list(names(fits[[s]]),fits[[s]][[1]]$mu$x))
      null.full[[s]] <- list()
      null.area[[s]] <- list()
      dtt <-  mclapply(1:length(fits[[s]]),function(x){
        #cat("\r",x)
        print(x)
        gn <- names(fits[[s]])[x]
        f1<-fits[[s]][[x]]
        f2<-fits[[pairs[s]]][[x]]
        cm = 1
        if( !bhatta & !z){
          cm <- dist(f1$mu$y,f2$mu$y)^power
        }
        if(z) {
          cm <- dist(f1$mu$y/max(f1$mu$y) ,f2$mu$y/max(f2$mu$y))^power
        }
        if(derivative != 0){
          cm <- cm + dist(predict(f1$spline,x = f1$mu$x,deriv = derivative)$y,predict(f2$spline,x = f1$mu$x,deriv = derivative)$y  )^power
        }
        if(bhatta){
          cm <- bhatta.dist.mat(f1,f2)^power
        }

        dtq1 <-  dtw(cm, open.end = T, open.begin = open.start, step.pattern = stepPattern,keep.internals = F)
        distx <- dtq1$normalizedDistance
        ans <- warp(dtq1,index.reference = T)

        rsind1 <-  sapply(1:n.null, function(x) sample(1:nrow(cm),nrow(cm),replace = F))
        rsind2 <-  sapply(1:n.null, function(x) sample(1:ncol(cm),ncol(cm),replace = F))
        an <- list()
        aarea <- list()
        n <- sapply(1:ncol(rsind1),function(x){
          dtn <-  dtw(cm[rsind1[,x],rsind2[,x]], open.end = T , open.begin = open.start, step.pattern = stepPattern,keep.internals = F)
          if(x%% as.integer(n.null/10) == 0){
            an[[length(an)+1]] <<- warp(dtn,index.reference = T)
            aarea[[length(aarea)+1]] <<- warpArea2(dtn)
          }
          return(dtn$normalizedDistance)
        })
        
        return(list(distx,ans,n,t(simplify2array(an)),dtq1,aarea))
      })
      for(gn in 1:length(fits[[s]])){
        dt.objects[[s]][[gn]] <- dtt[[gn]][[5]]
        dt.full[[s]][gn,] <- dtt[[gn]][[2]]
        dt.dist[[s]][[gn]] <- dtt[[gn]][[1]]
        null[[s]][[gn]] <- dtt[[gn]][[3]]
        null.full[[s]][[gn]] <-  dtt[[gn]][[4]]
        null.area[[s]][[gn]] <- dtt[[gn]][[6]]
      }
      null.full[[s]] <- Reduce(rbind,null.full[[s]])
      null.area[[s]] <- unlist(null.area[[s]])
      rownames(dt.full[[s]]) <- names(dt.objects[[s]])<- names(dt.dist[[s]]) <-  names(null[[s]]) <- names(null.area[[s]]) <-  names(fits[[s]])
      colnames(dt.full[[s]])<- fits[[s]][[1]]$mu$x
      format(Sys.time(),"%D %H:%M:%S")
    }
    print(format(Sys.time(),"%D %H:%M:%S"))
    save(dt.objects,dt.full,dt.dist,null,null.full,null.area,file = file.path(output.dir,"dtw.RData") )
  }

  dt.warp.area <- lapply(1:length(dt.objects),function(s){
    sapply(dt.objects[[s]],function(x){
      warpArea2(x)
    })
  })

  p.dist <- mc.pval2(dt.dist[[1]],null[[1]],F)
  p.dist2 <- mc.pval2(dt.dist[[2]],null[[2]],F)
  names(p.dist) <- names(p.dist2) <- names(dt.dist[[1]])
  
  post.rancors <- sapply(1:length(fits[[1]]), function(x){
    a <- fits[[1]][[x]]$mu$y[colnames(dt.full[[1]])]
    b <- fits[[pairs[1]]][[x]]$mu$y[dt.full[[1]][x,]]
    cor(a,b  ,method="spearman")
  })
  names(post.rancors) <- names(fits[[1]])
  
  post.rancors2 <- sapply(1:length(fits[[1]]), function(x){
    a <- fits[[1]][[x]]$mu$y[dt.full[[2]][x,]]
    b <- fits[[pairs[1]]][[x]]$mu$y[colnames(dt.full[[2]])]
    cor(a,b  ,method="spearman")
  })
  names(post.rancors2) <- names(fits[[1]])
  post.rancors[is.na(post.rancors)] <- 0
  post.rancors2[is.na(post.rancors2)] <- 0
  
  write.pdf({
    mypar(3,4)
    n =names(p.dist)
    for(gn in n[1:96] ){
      plot.model(fits[[1]][[gn]],main=paste(gn,post.rancors[gn],"\n",post.rancors2[gn]))
      plot.model(fits[[2]][[gn]],main=paste(p.dist[gn],"\n",p.dist2[gn]))
      dtwPlotTwoWay(dt.objects[[1]][[gn]], fits[[1]][[gn]]$mu$y,fits[[2]][[gn]]$mu$y,main=paste(dt.dist[[1]][gn],"\n", dt.dist[[2]][gn]),match.col = "black")
      dtwPlotTwoWay(dt.objects[[2]][[gn]], fits[[2]][[gn]]$mu$y, fits[[1]][[gn]]$mu$y, main =paste(warpArea2(dt.objects[[1]][[gn]]), "\n",warpArea2(dt.objects[[2]][[gn]])),match.col = "black")
    }
    for(gn in n[(length(n)-96):length(n)] ){
      plot.model(fits[[1]][[gn]],main=paste(gn,post.rancors[gn],"\n",post.rancors2[gn]))
      plot.model(fits[[2]][[gn]],main=paste(p.dist[gn],"\n",p.dist2[gn]))
      dtwPlotTwoWay(dt.objects[[1]][[gn]], fits[[1]][[gn]]$mu$y,fits[[2]][[gn]]$mu$y,main=paste(dt.dist[[1]][gn],"\n", dt.dist[[2]][gn]),match.col = "black")
      dtwPlotTwoWay(dt.objects[[2]][[gn]], fits[[2]][[gn]]$mu$y, fits[[1]][[gn]]$mu$y,main=paste(warpArea2(dt.objects[[1]][[gn]]), "\n",warpArea2(dt.objects[[2]][[gn]])),match.col = "black")
    }
  },filename = file.path(output.dir,"diagnostics", "check96.pdf")
  )

  write.pdf({
    mypar(3,4)
    n =names(p.dist)[order(p.dist,decreasing = T)]
    for(gn in n[1:36] ){
      plot.model(fits[[1]][[gn]],main=paste(gn,post.rancors[gn],"\n",post.rancors2[gn]))
      plot.model(fits[[2]][[gn]],main=paste(p.dist[gn],"\n",p.dist2[gn]))
    }
    for(gn in n[(length(n)-36):length(n)] ){
      plot.model(fits[[1]][[gn]],main=paste(gn,dt.dist[[1]][gn],"\n",dt.dist[[2]][gn]))
      plot.model(fits[[2]][[gn]],main=paste(p.dist[gn],"\n",p.dist2[gn]))
    }
  },filename = file.path(output.dir,"diagnostics", "ordereddistpvals.pdf")
  )
  
  
  write.pdf({
    mypar(3,4)
    n =names(post.rancors)[order(post.rancors,decreasing = T)]
    for(gn in n[1:48] ){
      plot.model(fits[[1]][[gn]],main=paste(gn,post.rancors[gn],"\n",post.rancors2[gn]))
      plot.model(fits[[2]][[gn]],main=paste(p.dist[gn],"\n",p.dist2[gn]))
    }
    for(gn in n[(length(n)-48):length(n)] ){
      plot.model(fits[[1]][[gn]],main=paste(gn,post.rancors[gn],"\n",post.rancors2[gn]))
      plot.model(fits[[2]][[gn]],main=paste(p.dist[gn],"\n",p.dist2[gn]))
    }
  },filename = file.path(output.dir,"diagnostics", "orderedrankcorpost.pdf")
  )
  
  optimum <- sapply(1:length(p.dist), function(i){
    if(p.dist[i] == p.dist2[i]){
      dt.dist[[1]][[i]] <= dt.dist[[2]][[i]] 
    }else{
      p.dist[i] <= p.dist2[i]
    }
  })
  warpAreaPvals1f <- p.adjust(mc.pval( dt.warp.area[[1]] ,null.area[[1]],upper = T),method = "BH")
  warpAreaPvals2f <-  p.adjust(mc.pval( dt.warp.area[[2]] ,null.area[[2]],upper = F),method = "BH")
  warpAreaPvals1b <-  p.adjust(mc.pval( dt.warp.area[[1]] ,null.area[[1]],upper = F),method = "BH")
  warpAreaPvals2b <-  p.adjust(mc.pval( dt.warp.area[[2]] ,null.area[[2]],upper = T),method = "BH")
  
  names(optimum) <-  names(p.dist)
  opti <-  t(sapply(1:length(optimum),function(x){
    if(optimum[x]){
      c(mean(f.values[[1]][x],f.values[[2]][x],na.rm = T),p.dist[x],post.rancors[x],rancors[x],p.dist2[x],dt.dist[[1]][x],dt.dist[[2]][x],xcors[x],dt.fc[x], abs(w.xcors[x]+rowMeans(dt.full[[1]])[x]),warpAreaPvals1f[x],warpAreaPvals1b[x],T)
    }else{
      c(mean(f.values[[1]][x],f.values[[2]][x],na.rm = T),p.dist2[x],post.rancors2[x],rancors[x],p.dist[x],dt.dist[[1]][x],dt.dist[[2]][x],xcors[x],dt.fc[x],  abs(rowMeans(dt.full[[2]])[x]-w.xcors[x]),warpAreaPvals2f[x],warpAreaPvals2b[x],F)
    }
  }))
  rownames(opti) <- names(dt.dist[[1]])
  opti <- matrix(unlist(opti),nrow = nrow(opti),ncol=ncol(opti),dimnames = dimnames(opti))

  colnames(opti) <- c("f.pval","dist.pval","post.s.cor","pre.s.cor","other.dist.p","dist1","dist2","pre.cor.x","dt.fc","diff","areaPValF", "areaPValB","fwd")
  
  comparable <-   opti[,3] >= .7 & opti[,2] <= .05
  print("Estimating appropriate size of null distribution")
  #Estimate the null distribution
  number.null.obs <-  mclapply(1:length(null.full), function(s){
    sq <- as.integer(seq(from = 2, to = nrow(null.full[[s]]),length.out = 1000 ))
    is <-  sapply(sq,function(i){
      d <<- null.full[[s]][sample(1:nrow(null.full[[s]]),i),]
      bhattacharyya.dist(mean(null.full[[s]]), mean(d),sd(null.full[[s]]),sd(d))
    })
    sq[min(which(is <= 1E-5 & diff(is<=1E-5)==0 & diff(is<=1E-5,lag = 2)==0))]
  })

  gene.list <- gene.list[gene.list %in% names(fits[[1]])]
  
  write.pdf({
    mypar(3,4)
    n =gene.list
    for(gn in n ){
      plot.model(fits[[1]][[gn]],main=paste(gn,post.rancors[gn],"\n",post.rancors2[gn]))
      plot.model(fits[[2]][[gn]],main=paste(p.dist[gn],"\n",p.dist2[gn]))
      dtwPlotTwoWay(dt.objects[[1]][[gn]], fits[[1]][[gn]]$mu$y,fits[[2]][[gn]]$mu$y,main=gn,match.col = "black")
      dtwPlotTwoWay(dt.objects[[2]][[gn]], fits[[2]][[gn]]$mu$y, fits[[1]][[gn]]$mu$y, main =gn,match.col = "black")
    }
  },filename = file.path(output.dir,"diagnostics", "checkNEU.pdf")
  )
  
  dt = list()
  dt[[1]] <- dt.full[[1]][comparable,]
  dt[[2]] <- dt.full[[2]][comparable,]
  opti <- opti[comparable,]
  dt[[3]] <- null.full[[1]][sample(1:nrow(null.full[[1]]),number.null.obs[[1]],replace = F),]
  dt[[4]] <- null.full[[2]][sample(1:nrow(null.full[[2]]),number.null.obs[[2]],replace = F),]
  
  avg.shift <-  colMeans(dt[[3]])
  rev.avg.shift <- colMeans(dt[[4]])
  dt.adj <- sweep(dt[[1]],2, avg.shift)
  nd.adj <- sweep(dt[[3]],2, avg.shift)
  nd <- t(sapply(1:nrow(dt[[3]]),function(gn){
      sapply(1:length(dt[[3]][gn,]), function(x){
        as.numeric(dt[[3]][gn,x])-as.numeric(fits[[1]][[1]]$mu$x[x])
      })
  }))
  nd.rev <- t(sapply(1:nrow(dt[[4]]),function(gn){
    sapply(1:length(dt[[4]][gn,]), function(x){
      as.numeric(dt[[4]][gn,x])-as.numeric(fits[[2]][[1]]$mu$x[x])
    })
  }))
    
  rev.dt.adj <- sweep(dt[[2]],2, rev.avg.shift)
  rev.nd.adj <- sweep(dt[[4]],2, rev.avg.shift)
  rm1 <- rowMeans(dt.adj)
  rm2 <- rowMeans(rev.dt.adj)
  rm3 <- rowMeans(nd.adj)
  rm4 <- rowMeans(rev.nd.adj)
  
  print("W Testing")
  if(file.exists(file.path(output.dir,"metrics.txt"))){
    opti <- read.table(file.path(output.dir,"metrics.txt"),header=T,row.names=1,sep = "\t")
    dt.combined <- read.table(file.path(output.dir,"shift.txt"),header=T,row.names = 1,sep="\t")
  }else{
  p.sig <- unlist(mclapply(1:nrow(opti),function(i){
    #cat("\r",i)
    print(i)
    gn <- rownames(opti)[i]
    if(opti[i,ncol(opti)]){
      wilcox.test(dt.full[[1]][gn,],dt[[3]],alternative = "g")$p.value
    }else{
      wilcox.test(dt.full[[2]][gn,],dt[[4]],alternative = "l")$p.value
    }
  }))
  
  mindex <- intersect(colnames(dt.full[[1]]), colnames(dt.full[[2]]))
  dt.combined <-lapply(1:nrow(opti),function(i){
    gn <- rownames(opti)[i]
    if(opti[i,ncol(opti)]){
      w = warp(dt.objects[[1]][[gn]],T)
      sapply(1:length(w), function(x){
        as.numeric(fits[[1]][[gn]]$mu$x[w[x]])-as.numeric(fits[[1]][[gn]]$mu$x[x])
      })
    }else{
      w = warp(dt.objects[[2]][[gn]],T)
      fits[[2]][[gn]]$mu$x
      -sapply(1:length(w), function(x){
        as.numeric(fits[[2]][[gn]]$mu$x[w[x]])-as.numeric(fits[[2]][[gn]]$mu$x[x])
      })
    }
  })
  
  dt.combined <- t(simplify2array(dt.combined))
  colnames(dt.combined) <- mindex
  rownames(dt.combined) <- rownames(opti)
  opti <- as.data.frame(opti)
  opti[,"p"] <-  p.sig
  opti[,"p.fwd.adj"] <- p.adjust(p.sig,"BH")
  opti[,"p.bck.adj"] <- p.adjust(1-p.sig,"BH")
  write.table(opti,file.path(output.dir,"metrics.txt"),sep="\t")
  write.table(dt.combined,file = file.path(output.dir,"shift.txt"),sep="\t")
  }
  forward <- rownames(opti[(opti[,"p.fwd.adj"]<.05),])
  backward <-rownames(opti[(opti[,"p.bck.adj"]<.05),])
  forwardWA <- rownames(opti[(opti[,"areaPValF"]<.05),])
  backwardWA <-rownames(opti[(opti[,"areaPValB"]<.05),])
  dt.combined <- as.matrix(dt.combined)

  write.pdf(  barplot(c("Comparable"=sum(comparable),"NonComparable"=sum(!comparable)),main = "Distance Pvalue <.05",
                      xlab = "Class", ylab = "# of Genes"), file.path(output.dir,"shift_comparables.pdf"))

  write.pdf({
    mypar(3,4)
    for(gn in sample(names(fits[[1]])[comparable],48,replace = T)){
      plot.model( fits[[1]][[gn]],main=paste(gn,paste(name1, "Log2 TPM")))
      plot.model( fits[[2]][[gn]],main=paste(gn,paste(name2,"Log2 TPM")))
    }
  }  , file.path(output.dir,"diagnostics","comparables.pdf"))
  
  write.pdf({
    mypar(3,4)
    for(gn in sample(names(fits[[1]])[!comparable],48,replace = T)){
      plot.model( fits[[1]][[gn]],main=paste(gn,paste(name1, "Log2 TPM")))
      plot.model( fits[[2]][[gn]],main=paste(gn,paste(name2, "Log2 TPM")))
    }
  }  , file.path(output.dir,"diagnostics","uncomparables.pdf"))
  
  exprsmat1 <- t(sapply(fits[[1]],function(x) x$mu$y))
  exprsmat2 <- t(sapply(fits[[2]],function(x) x$mu$y))
  exprsmat1 <-exprsmat1[comparable,]
  exprsmat2 <-exprsmat2[comparable,]
  c=cor(exprsmat1,exprsmat2)
  write.pdf({
  heatmap.2(c,Rowv = F,Colv = F,trace="none",col = cols ,main=paste("Spearman Correlation\n",name1,"(rows)",name2, "(columns)"))
  },filename = file.path(output.dir,"diagnostics","TPCorrelations.pdf"))
  write.pdf({
    plot(apply(c, 2,which.max),main = "Timepoint Equivalence\n (Spearman Correlation)", type = "l",xlab = paste(name2,"Days"),ylab= paste(name1,"Days"),axes = F)
    axis(side=1,at=1:length(colnames(c)) ,labels= colnames(c), tck=F)
    axis(side=2,at=1:length(rownames(c)) ,labels= rownames(c), tck=F)
  },filename = file.path(output.dir,"diagnostics","TPEquivalence.pdf"))
  

  #######
  write.csv(set1,file = file.path(output.dir,"diagnostics","set1exprs.csv"),row.names = T)
  write.csv(set2,file = file.path(output.dir,"diagnostics","set2exprs.csv"), row.names = T)
  
  write.csv(dt.full[[1]],file = file.path(output.dir,"diagnostics","fwd_shift.csv"))
  write.csv(dt.full[[2]],file = file.path(output.dir,"diagnostics","rev_shift.csv"))

  print(paste("Mean Null shift",mean(rowMeans(nd))))
  print(paste("Mean shift experimental:", mean(rowMeans(dt[[1]]))))
  
  write.pdf(  barplot(colMeans(dt.combined),main = paste("Experiment Mean Shift at each Timepoint, All Genes",name1,"Forward = Positive" ),
                      xlab = "TP (day)", ylab = "Shift (# TPs)"), file.path(output.dir,"exp_Shift_tps.pdf"))
  
  write.pdf( hist(dt.combined,main = paste("Average Experiment Shift (across TPs for each gene), All Genes",name1,"Forward = Positive" ),
                    xlab = "shift (# TPs)", ylab = "Number of Observations"), file.path(output.dir,"exp_Shift_Hist.pdf"))  
  
  write.pdf( hist(c(dt.warp.area[[1]], -dt.warp.area[[2]]) ,main = paste("Area Between Equivalence Curve and Diagonal" ),
                  xlab = "Area (days^2)", ylab = "Number of Observations"), file.path(output.dir,"diagnostics","exp_area_Hist.pdf"))  
  
  write.pdf( hist(c(null.area[[1]], -null.area[[2]]) ,main = paste("Area Between Equivalence Curve and Diagonal" ),
                  xlab = "Area (days^2)", ylab = "Number of Observations"), file.path(output.dir,"null_area_Hist.pdf"))  
  
  
  write.pdf( plot(colMeans(rbind(dt[[3]] ,dt[[4]])),main = paste("Average Null TP Equivalence"),
                  xlab = name1, ylab = name2), file.path(output.dir,"Null_equivalence.pdf"))
  write.pdf( plot(colMeans(dt[[1]] ),main = paste("Average TP Equivalence"),
                  xlab = name1, ylab = name2), file.path(output.dir,"diagnostics","Exp1_equivalence.pdf"))
  write.pdf( plot(colMeans(dt[[2]] ),main = paste("Average TP Equivalence"),
                  xlab = name2, ylab = name1), file.path(output.dir,"diagnostics","Exp2_equivalence.pdf"))
  
  write.pdf( plot(colMeans(dt.combined ),main = paste("Average TP Equivalence"),
                  xlab = name1, ylab = name2), file.path(output.dir,"diagnostics","combined_equivalence.pdf"))
  
  
  write.pdf(  barplot(colMeans(nd),main = paste("Null Mean Shift at each Timepoint, All Genes",name1,"Forward = Positive" ),
                      xlab = "TP (day)", ylab = "Shift (# TPs)"), file.path(output.dir,"Null_Shift_tps.pdf"))
  
  write.pdf(   hist(nd,main = paste("Average Null  Shift (across TPs for each gene), All Genes",name1,"Forward = Positive" ),
                    xlab = "shift (# TPs)", ylab = "Number of Observations"), file.path(output.dir,"Null_Shift_Hist.pdf"))
  
  #####Run T-Tests on overall data
  t <-  t.test(rowMeans(dt[[3]]),rowMeans(dt[[1]]))
  #Significance for each gene
  print("W Test")
 
  write.pdf(barplot( c("Forward Shifted"=length(forward),"Backward Shifted"=length(backward)) ,las=1,xpd = F ,main = paste("# Of of Experimental Genes Shifted Fwd vs Back \nMann-Whitney-Wilcoxon Test (5% FDR):",name1),
                     ylab = "# Of Genes"), file.path(output.dir,"exp_fwd_vs_back.pdf"))
  
  write.pdf(barplot( c("Forward Shifted"=length(forwardWA),"Backward Shifted"=length(backwardWA)) ,las=1,xpd = F ,main = paste("# Of of Experimental Genes Shifted Fwd vs Back \nBy Warp Area (5% FDR):",name1),
                     ylab = "# Of Genes"), file.path(output.dir,"exp_fwd_vs_back_warpArea.pdf"))
  
  write.pdf({
    mypar(3,4)
    for(gn in sample(forward,48,replace = T)){
      plot.model( fits[[1]][[gn]],main=paste(gn,paste(name1, "Log2 TPM")))
      plot.model( fits[[2]][[gn]],main=paste(gn,paste(name2,"Log2 TPM")))
    }
  },file=file.path(output.dir,"diagnostics","forward_shifted_viz.pdf"))
  
  write.pdf({
    mypar(3,4)
    for(gn in sample(backward,48,replace = T)){
      plot.model( fits[[1]][[gn]],main=paste(gn,paste(name1, "Log2 TPM")))
      plot.model( fits[[2]][[gn]],main=paste(gn,paste(name2,"Log2 TPM")))
    }
  },file=file.path(output.dir,"diagnostics","back_shifted_viz.pdf"))
  
  if(ppt){
    fwdGO <-read.table(file =file.path(output.dir,"ForwardEnrichment.tsv") ,sep = "\t",header = T)  
    bckGO <- read.table(file =file.path(output.dir,"BackwardEnrichment.tsv") ,sep = "\t",header = T)  
  }else{
    if(grepl("mm",name1) & grepl("mm",name2)){
      libspec = "mouse"
      xx <- as.list(org.Mm.egALIAS2EG)
      names(xx) <- toupper(names(xx))
    }else{
      libspec = "human"
      xx <- as.list(org.Hs.egALIAS2EG)
      names(xx) <- toupper(names(xx))
    }
    u=enrichGO(gene= sapply(xx[forward],"[",1), organism= libspec,universe = unlist(sapply(xx[all.genes],"[",1)),ont = "BP")
    print(u)
    d=enrichGO(gene= sapply(xx[backward],"[",1), organism= libspec,universe = unlist(sapply(xx[all.genes],"[",1)),ont = "BP")
    print(d)
    if(!is.null(u))write.table(summary(u),file = paste0(output.dir,"/","ForwardEnrichment.tsv"),sep = "\t",row.names = F ,quote = F )
    if(!is.null(d))write.table(summary(d),file = paste0(output.dir,"/","BackwardEnrichment.tsv"),sep ="\t",row.names = F ,quote = F )
    
  }
  pos.shift.pvals <-  sapply(1:ncol(dt[[1]]) , function(i){
    print(paste("Pval tp",i))
    sapply(1:nrow(dt[[1]]),function(j){
      (sum( dt[[1]][j,i] <= dt[[3]][,i])+1)/(length(dt[[3]][,i])+1)
    })
  } )

  neg.shift.pvals <-  sapply(1:ncol(dt[[2]]) , function(i){
    print(paste("Pval tp",i))
    sapply(1:nrow(dt[[2]]),function(j){
      (sum( dt[[2]][j,i] <= dt[[4]][,i])+1)/(length(dt[[4]][,i])+1)
    })
  } )

  dimnames(pos.shift.pvals) <- dimnames(dt[[1]])
  dimnames(neg.shift.pvals) <- dimnames(dt[[2]])
  
  gene.list <- gene.list[gene.list%in%rownames(dt.combined)]

  write.csv(pos.shift.pvals,file = file.path(output.dir,"diagnostics",paste0(name1,"forward_shift_pvals.csv")))
  write.csv(neg.shift.pvals,file = file.path(output.dir,"diagnostics",paste0(name1,"back_shift_pvals.csv")))
  
  write.pdf( heatmap.2( dt.combined[gene.list,] ,col = cols,
                        trace = "none",Rowv = T,Colv = F,cexRow = .35,cexCol = .45,
                        main = paste("Shift (# of TPs)",name1,"forward = positive")), file.path(output.dir,"Shifts.pdf"))

  write.pdf( heatmap.2( pos.shift.pvals[gene.list,],col = rev(cols),
                        trace = "none",Rowv = F,Colv = F,cexRow = .35,cexCol = .45,main = "Forward Shift Pvals"), file.path(output.dir,"Selected_Forward_Pvals.pdf"))
  
  write.pdf( heatmap.2( neg.shift.pvals[gene.list,],col = rev(cols),
                        trace = "none",Rowv = F,Colv = F,cexRow = .35,cexCol = .45,main = "Backward Shift Pvals"), file.path(output.dir,"Selected_Backward_Pvals.pdf"))
  
  write.pdf(heatmap.2( cbind(exprsmat1[gene.list,], exprsmat2[gene.list,]),col = cols,
                       trace = "none",Rowv = F,Colv = F,cexRow = .35,cexCol = .45,main = "Expression of Selected Genes (log2 TPM)"), file.path(output.dir,"Selected_Expression.pdf"))
  
  write.pdf(  barplot(colMeans(pos.shift.pvals < .05),main = "% of Experimental Genes Shifted Ahead \n at each TP\n(5% FDR)",
                      xlab = "TP (day)", ylab = "% of Genes"), file.path(output.dir,"exp_forward_shift_significant.pdf"))
  
  write.pdf(  barplot(colMeans(neg.shift.pvals < .05 ),main = "% of Experimental Genes Shifted Back \n at each TP\n(5% FDR)",
                      xlab = "TP (day)", ylab = "% of Genes"), file.path(output.dir,"exp_back_shift_significant.pdf"))
  
  t.tst <-  data.frame(t$estimate[1], t$estimate[2], t$p.value, row.names = "Value")
  colnames(t.tst) <- c("Control Avg Shift", "Experimental Avg Shift", "P.value")
  write.csv(t.tst,file = file.path(output.dir,"overall_shift_T_Test.csv"))
  if(ppt){
    library(ReporteRs)
    mydoc = pptx()
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle(mydoc, "DTW Analysis")
    mydoc = addParagraph( mydoc, value = paste(output.dir,"Log FC Threshold",fc.threshold , "Step Pattern", paste(stepPattern,collapse = ",") ))
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle( mydoc, "Overall T Test")
    mydoc = addParagraph(mydoc, paste("Mean Null:",t$estimate[1], "Mean",name1,"-",name2 ,":",t$estimate[2], "\nP Value",  t$p.value))
    mydoc = addSlide( mydoc, "Two Content" )
    mydoc = addTitle( mydoc, "Shifts Observed")
    mydoc = addPlot(mydoc, function() hist(dt[[1]],main = "Average Experiment Shift (across TPs for each gene), All Genes",
                                           xlab = "shift (# TPs)", ylab = "Number of Observations"), vector.graphic = F)
    mydoc = addPlot(mydoc, function()  hist(nd,main = "Average Null  Shift (across TPs for each gene), All Genes",
                                            xlab = "shift (# TPs)", ylab = "Number of Observations"), vector.graphic = F)
    
    
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle( mydoc, "Classification")
    mydoc = addPlot(mydoc, function() barplot(c("Comparable"=sum(comparable),"NonComparable"=sum(!comparable)),main = "Classifications By Random Forest",
                                              xlab = "Class", ylab = "# of Genes"), vector.graphic = F)
    
    mydoc = addSlide( mydoc, "Two Content" )
    mydoc = addTitle( mydoc, "Sanity Check\n Comparable vs NonComparable")
    mydoc = addPlot(mydoc, function() {
      mypar(3,4)
      for(gn in sample(names(fits[[1]])[comparable],12,replace = T)){
        plot.model( fits[[1]][[gn]],main=paste(gn,paste(name1, "Log2 TPM")))
        plot.model( fits[[2]][[gn]],main=paste(gn,paste(name2,"Log2 TPM")))
      }
    }, vector.graphic = F)
    mydoc = addPlot(mydoc, function() {
      mypar(3,4)
      for(gn in sample(names(fits[[1]])[!comparable],12,replace = T)){
        plot.model( fits[[1]][[gn]],main=paste(gn,paste(name1, "Log2 TPM")))
        plot.model( fits[[2]][[gn]],main=paste(gn,paste(name2,"Log2 TPM")))
      }
    }, vector.graphic = F)
    
    mydoc = addSlide( mydoc, "Two Content" )
    mydoc = addTitle( mydoc, "Shifts at each TP")
    mydoc = addPlot(mydoc, function() barplot(colMeans(dt[[1]]),main = "Experiment Mean Shift at each Timepoint, All Genes",
                                              xlab = "TP (day)", ylab = "Shift (# TPs)"), vector.graphic = F)
    mydoc = addPlot(mydoc, function() barplot(colMeans(nd),main = "Null Mean Shift at each Timepoint, All Genes",
                                              xlab = "TP (day)", ylab = "Shift (# TPs)"), vector.graphic = F)
    
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle( mydoc, "Mann-Whitney-Wilcox Test")
    mydoc = addPlot(mydoc, function() barplot( c("Forward Shifted"=length(forward),"Backward Shifted"=length(backward)) ,las=1,xpd = F ,main = paste("# Of of Experimental Genes Shifted Fwd vs Back \nMann-Whitney-Wilcoxon Test (5% FDR):",name1),
                                               ylab = "# Of Genes"), vector.graphic = F)
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle(mydoc, "GO Enrichment Fwd Shifted")
    mydoc = addFlexTable(mydoc,flextable = FlexTable(fwdGO[1:50,c("Term","set.size","pval")]))
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle(mydoc, "GO Enrichment Back Shifted")
    mydoc = addFlexTable(mydoc,flextable = FlexTable(bckGO[1:50,c("Term","set.size","pval")]))
    mydoc = addSlide( mydoc, "Two Content" )
    mydoc = addTitle( mydoc, "Shifts & Expression")
    mydoc = addPlot(mydoc, function()  heatmap.2( dt.combined[gene.list,] ,col = cols,
                                                  trace = "none",Rowv = T,Colv = F,cexRow = .35,cexCol = .45,
                                                  main = "Shift (# of TPs)"), vector.graphic = F)
    mydoc = addPlot(mydoc, function()  heatmap.2( cbind(exprsmat1[gene.list,], exprsmat2[gene.list,]),col = cols,
                                                  trace = "none",Rowv = F,Colv = F,cexRow = .35,cexCol = .45,main = paste("Expression of Selected Genes (log2 TPM)\n",name1,"left",name2,"right")), vector.graphic = F)
    
    mydoc = addSlide( mydoc, "Two Content" )
    mydoc = addTitle( mydoc, "P Values")
    mydoc = addPlot(mydoc, function()  heatmap.2( pos.shift.pvals[gene.list,],col = rev(cols),
                                                  trace = "none",Rowv = F,Colv = F,cexRow = .35,cexCol = .45,main = "Forward Shift Pvals"), vector.graphic = F)
    mydoc = addPlot(mydoc, function()  heatmap.2( neg.shift.pvals[gene.list,],col = rev(cols),
                                                  trace = "none",Rowv = F,Colv = F,cexRow = .35,cexCol = .45,main = "Backward Shift Pvals"), vector.graphic = F)
    writeDoc( mydoc, file.path(output.dir,"summary.pptx") )
  }
}

condit.combos <- expand.grid(c("mm","hs"),unique(conditMixSets[[2]]),stringsAsFactors = F)
condit.combos <- condit.combos[-c(2,7),]

apply(combn(1:nrow(condit.combos),2),2,function(c){
    x <- log(MixSets[[condit.combos[c[1],1]]][,conditMixSets[[condit.combos[c[1],1]]]==condit.combos[c[1],2]]+1,2)
    y <- log(MixSets[[condit.combos[c[2],1]]][,conditMixSets[[condit.combos[c[2],1]]]==condit.combos[c[2],2]]+1,2)
    if(ncol(x)> 3 & ncol(y)>3){
      run.DTW.Genes(set1 = list(x) ,set2 = list(y), output.dir = paste0("~/code/data/cb/shiftFiles/", condit.combos[c[1],1],condit.combos[c[1],2],"vs",condit.combos[c[2],1],condit.combos[c[2],2],"rj1cZ"),
                    tp1 =list(tpMixSets[[condit.combos[c[1],1]]][conditMixSets[[condit.combos[c[1],1]]]==condit.combos[c[1],2]]),
                    tp2=list(tpMixSets[[condit.combos[c[2],1]]][conditMixSets[[condit.combos[c[2],1]]]==condit.combos[c[2],2]]),
                    name1 = paste0(condit.combos[c[1],1],condit.combos[c[1],2]),name2 = paste0(condit.combos[c[2],1],condit.combos[c[2],2]),power=2,derivative = 1,tpm.threshold = 3,fc.threshold = 2,
                    n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = T,
                    equalize=T,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="c",smoothed=T))
    }
})

run.DTW.Genes(set1 =  list(datasets[[2]]),set2 = list(datasets[[4]]), output.dir = "~/code/data/cb/shiftFiles/HTeravsMTerarj1Z",tp1 =list(tpDatasets[[2]]), tp2=list(tpDatasets[[4]]),name1 = datasetNames[[2]],name2 = datasetNames[[4]],n.null=100,gene.list = neural_list, open.start = F,ppt = F,power = 2 ,derivative = 1,equalize=T,stepPattern =rabinerJuangStepPattern(type=1,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 =  list(datasets[[1]]),set2 = list(datasets[[3]]), output.dir = "~/code/data/cb/shiftFiles/HTeraControlvsMTeraControlrj1Z",tp1 =list(tpDatasets[[1]]), tp2=list(tpDatasets[[3]]),name1 = datasetNames[[2]],name2 = datasetNames[[4]],n.null=100,gene.list = neural_list, open.start = F,ppt = F,power = 2 ,derivative = 1,equalize=T,stepPattern =rabinerJuangStepPattern(type=1,slope.weighting="c",smoothed=T))


run.DTW.Genes(set1 = setMix[-1] ,set2 = setH[2], output.dir = "~/code/data/cb/shiftFiles/Mix33vsH33InVitrorj1cubehybrid",tp1 =tpMix[-1], tp2=tpH[2],name1 = "Mix33day",name2 = "Human",n.null=100,gene.list = neural_list, open.start = F,ppt = F,derivative = 1,power = 3 ,equalize=T,stepPattern =rabinerJuangStepPattern(type=1,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/Mix33vsHInVitrorj7cubehybrid",tp1 =tpMix[-1], tp2=tpH,name1 = "Mix33day",name2 = "Human",n.null=100,gene.list = neural_list, open.start = F,ppt = F,power = 3 ,derivative = 1,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setH[2], output.dir = "~/code/data/cb/shiftFiles/Mix33vsH33InVitrorj7cubehybrid",tp1 =tpMix[-1], tp2=tpH[2],name1 = "Mix33day",name2 = "Human",n.null=100,gene.list = neural_list, open.start = F,ppt = F,power = 3 ,derivative = 1,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))

run.DTW.Genes(set1 = list(MixSets$mm[,1:26]) ,set2 = list(MixSets$hs[,77:102]), output.dir = "~/code/data/cb/shiftFiles/NewMousevsNewHumanrj7zcubehybrid",tp1 =list(tpMixSets[1:26]), tp2=list(tpMixSets[77:102]),name1 = "NewMouse",name2 = "NewHuman",power=3,derivative = 1,tpm.threshold = 2,fc.threshold = 2,n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = list(MixSets$hs[,27:51]) ,set2 = list(MixSets$hs[,77:102]), output.dir = "~/code/data/cb/shiftFiles/NewMix10vsNewHumanrj7zcubehybrid",tp1 =list(tpMixSets[27:51]), tp2=list(tpMixSets[77:102]),name1 = "NewMix10",name2 = "NewHuman",power=3,derivative = 1,tpm.threshold = 2,fc.threshold = 2,n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = list(MixSets$hs[,52:76]) ,set2 = list(MixSets$hs[,77:102]), output.dir = "~/code/data/cb/shiftFiles/NewMix85vsNewHumanrj7zcubehybrid",tp1 =list(tpMixSets[52:76]), tp2=list(tpMixSets[77:102]),name1 = "NewMix85",name2 = "NewHuman",power=3,derivative = 1,tpm.threshold = 2,fc.threshold = 2 ,n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = list(MixSets$mm[,27:51]) ,set2 = list(MixSets$mm[,1:26]), output.dir = "~/code/data/cb/shiftFiles/NewMixMouse10vsNewMouserj7zcubehybrid",tp1 =list(tpMixSets[27:51]), tp2=list(tpMixSets[1:26]),name1 = "NewMixMouse10",name2 = "NewMouse",derivative = 1,power=3,tpm.threshold = 2,fc.threshold = 2,n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))



experiments = list("hs,100%"=MixSets$hs[,77:102],"hs,90%"= MixSets$hs[,c(77,52:76)] , "hs,10%"=MixSets$hs[,c(77,27:51)]  , "mm,0%"=MixSets$mm[,1:26], "mm,10%"=MixSets$mm[,c(1,27:51)] , "mm,90%"=MixSets$mm[,c(1,52:76)], "hs33dayMix"=datasets$`~/code/data/cb/TC_2_Mix.csv`[int_rns,], "hs33dayH"=datasets$`~/code/data/cb/TC_2_H.csv`[int_rns,])
experimentGroups <- factor(c(1,1,1,2,2,2,3,3))

int_rns <-  intersect(rownames(datasets$`~/code/data/cb/TC_2_Mix.csv`), rownames(datasets$`~/code/data/cb/TC_2_H.csv`))
DifferentGenes <-  sapply(levels(experimentGroups),function(g){
  curGroup = experiments[experimentGroups==g]
  sapply(names(curGroup), function(i){
    sapply( names(curGroup), function(j){
      IntegralOrder = NULL
      if(i != j ){
        contrast <- (rowSums(curGroup[[i]])-rowSums(curGroup[[j]]))
        contrastOrder <- -contrast[order(contrast)]
        major <- (rowMeans(curGroup[[i]]) <= 2 | rowMeans(curGroup[[j]]) <=2 ) & (rowMaxs(curGroup[[i]]) >= 4 | rowMaxs(curGroup[[j]]) >= 4 )
        contrastMajor <- contrast[major]
        IntegralOrder <- -contrastMajor[order(contrastMajor)]
        other = setdiff(names(curGroup), c(i,j) )
        dir.create(paste0("~/code/data/cb/shiftFiles/DifferentExpressions/",make.names(j),"minus",make.names(i)))
        write.pdf({
        mypar(3,3)
        for(gn in names(IntegralOrder[1:27])){
          if(grepl("hs",i)){
            if(gn %in% rownames(experiments[["mm,0%"]])){
            x1 <- experiments[["mm,0%"]][gn,]
            x2 <- experiments[["mm,10%"]][gn,]
            }else{
              x1 <- x2 <- rep(0.0,ncol(curGroup[[i]]))
            }
            nx1 <- "mm,0%"
            nx2 <- "mm,10%"
          }else{
            if(gn %in% rownames(experiments[["hs,0%"]])){
            x1 <- experiments[["hs,100%"]][gn,]
            x2 <- experiments[["hs,10%"]][gn,]
            }else{
              x1 <- x2 <- rep(0.0,ncol(curGroup[[i]]))
            }
            nx1 <- "hs,100%"
            nx2 <- "hs,10%"
          }
          s1 <- zoo(curGroup[[i]][gn,])  
          s2 <- zoo(curGroup[[j]][gn,])
          if(length(other)!= 0L){ so <- zoo(curGroup[[other]][gn,]) 
          plot.zoo(merge(s1,s2,so,x1,x2),plot.type = "single",col = 1:6, lty = 1:3,main=paste(gn, IntegralOrder[gn]),xlab = "TimePoint",ylab = "log2 TPM" )
          legend("topleft", c(i,j,other,nx1,nx2), col = 1:6, lty = 1:3)
          }else{
            plot.zoo(merge(s1,s2,x1,x2),plot.type = "single",col = 1:6, lty = 1:3,main=paste(gn, IntegralOrder[gn]),xlab = "TimePoint",ylab = "log2 TPM" )
            legend("topleft", c(i,j,nx1,nx2), col = 1:6, lty = 1:3)
          }
        }},filename = paste0("~/code/data/cb/shiftFiles/DifferentExpressions/",make.names(j),"minus",make.names(i),"/","FCIntegral.pdf"))
      
#       if(grepl("mm",i)){
#         libspec = "org.Mm.eg"
#       }else{
#         libspec = "org.Hs.eg"
#       }
      #libspec = "org.Hs.eg"
      #GOTerms <- eacitest(range01(contrastOrder),lib =libspec ,idtype = "SYMBOL",sets="GO")
      #write.table(GOTerms$setscores[1:50,],file = paste0("~/code/data/cb/shiftFiles/DifferentExpressions/",make.names(j),"minus",make.names(i),"/","GOTerms.txt") )
      }
      IntegralOrder
    })
  })
})


write.pdf({
  mypar(3,3)
  for(s in 1:length(setH)){
    for(gn in names(percent10IntegralOrder[1:100])){
      h10 <- zoo(MixSets$hs[gn,27:51],tpMixSets[27:51])  
      h85 <- zoo(MixSets$hs[gn,52:76],tpMixSets[52:76])
      h100 <- zoo(MixSets$hs[gn,77:102],tpMixSets[77:102]) 
  }
    plot.zoo(merge(h10,h85,h100),plot.type = "single",col = 1:4, lty = 1,main=paste(gn,percent10IntegralOrder[gn]),xlab = "days",ylab = "log2 TPM" )
    legend("topleft", c("hs,10%","hs,85%","hs,100%"), col = 1:3, lty = 1)
  }}
  ,filename = "~/Desktop/MostUpreg10.pdf")

write.pdf({
  mypar(3,3)
  for(gn in names(percent10IntegralOrder[(length(percent10IntegralOrder)-100):length(percent10IntegralOrder)])){
    h10 <- zoo(MixSets$hs[gn,27:51],tpMixSets[27:51])  
    h85 <- zoo(MixSets$hs[gn,52:76],tpMixSets[52:76])
    h100 <- zoo(MixSets$hs[gn,77:102],tpMixSets[77:102]) 
    plot.zoo(merge(h10,h85,h100),plot.type = "single",col = 1:4, lty = 1,main=paste(gn,percent10IntegralOrder[gn]),xlab = "days",ylab = "log2 TPM" )
    legend("topleft", c("hs,10%","hs,85%","hs,100%"), col = 1:3, lty = 1)
  }}
  ,filename = "~/Desktop/MostDownreg10.pdf")


SetNames <- c("17day100","33day100","TeraControl100","New100")
write.pdf({
  mypar(3,3)
  for(gn in neural_list){
    zooz <-  sapply(1:length(setH), function(i){
      if(gn %in% rownames(setH[[i]])){
        if(max(setH[[i]][gn,])>0){
          zoo( (setH[[i]][gn,])/max(setH[[i]][gn,])  ,tpH[[i]]) 
        }else{
          zoo( (setH[[i]][gn,])  ,tpH[[i]]) 
        }
      }else{
        zoo(rep( 0,ncol(setH[[i]])),tpH[[i]])
      }
    })
    plot.zoo(na.approx(Reduce(merge,zooz)),plot.type = "single",col = 1:6, lty = 1:2,main=paste(gn),xlab = "days",ylab = "% of Max" )
    legend("topleft", SetNames, col = 1:6, lty = 1:2)
  }}
  ,filename = "~/Desktop/HumanNeuralGenesControls.pdf")

sapply(setM, colnames)
SetNames <- c("17day0","33day0","TeraControl0","LongMouse0","New0")
write.pdf({
  mypar(3,3)
  for(gn in neural_list){
    zooz <-  sapply(1:length(setM), function(i){
      if(gn %in% rownames(setM[[i]])){
        if(max(setM[[i]][gn,])>0){
          zoo( setM[[i]][gn,]/max(setM[[i]][gn,]) ,tpM[[i]]) 
        }else{
          zoo(setM[[i]][gn,],tpM[[i]]) 
        }
        
      }else{
        zoo(rep( 0,ncol(setM[[i]])),tpM[[i]])
      }
    })
    plot.zoo(na.approx(Reduce(merge,zooz)),plot.type = "single",col = 1:6, lty = 1:2,main=paste(gn),xlab = "days",ylab = "% of Max" )
    legend("topleft", SetNames, col = 1:6, lty = 1:2)
  }}
  ,filename = "~/Desktop/MouseNeuralGenesControls.pdf")




run.DTW.Genes(set1 = setMTera ,set2 = setHTera, output.dir = "~/code/data/cb/shiftFiles/MouseTeravsHTerarj7zcube",tp1 =tpMTera, tp2=tpHTera,name1 = "MouseTeratoma",name2 = "HumanTeratoma",tpm.threshold = 0,fc.threshold = 0,n.null = 100,gene.list = neural_list, open.start = F,ppt = F,power=3 ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))

run.DTW.Genes(set1 = setMTera ,set2 = setM, output.dir = "~/code/data/cb/shiftFiles/MouseTeravsMouseInVitrorj7z",tp1 =tpMTera, tp2=tpM,name1 = "MouseTeratoma",name2 = "Mouse",tpm.threshold = 0,fc.threshold = 0,n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setM ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MousevsHInVitrorj7z",tp1 =tpM, tp2=tpH,name1 = "Mouse",name2 = "Human",n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="d",smoothed=T))
run.DTW.Genes(set1 = setMix ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixvsHInVitrorj7z",tp1 =tpMix, tp2=tpH,name1 = "Mix",name2 = "Human",n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setM ,set2 = setMix, output.dir = "~/code/data/cb/shiftFiles/MousevsMixInVitrorj7z",tp1 =tpM, tp2=tpMix,name1 = "Mouse",name2 = "Mix",n.null = 100,gene.list = neural_list, open.start = F,ppt = F,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-2] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixShortvsHInVitrorj7z",tp1 =tpMix[-2], tp2=tpH,name1 = "Mix17day",name2 = "Human",n.null = 100,gene.list = neural_list, open.start = F,ppt = F,z = T ,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixLongvsHInVitrorj7z",tp1 =tpMix[-1], tp2=tpH,name1 = "Mix33day",name2 = "Human",n.null=100,gene.list = neural_list, open.start = F,ppt = F ,z = T,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setMix[-2], output.dir = "~/code/data/cb/shiftFiles/MixLongvsMixShortInVitrorj7z",tp1 =tpMix[-1], tp2=tpMix[-2],name1 = "Mix33day",name2 = "Mix17day",n.null=100,gene.list = neural_list, open.start = F,ppt = F ,z = T,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setM[c(-1,-2)] ,set2 = setM[c(-3,-4)], output.dir = "~/code/data/cb/shiftFiles/M34vsM12InVitrorj7z",tp1 =tpM[c(-1,-2)], tp2=tpM[c(-3,-4)],name1 = "Mouse34",name2 = "Mouse12",n.null=100,gene.list = neural_list, open.start = F,z = T,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))



run.DTW.Genes(set1 = setH ,set2 = list(MixSets$hs[,77:102]), output.dir = "~/code/data/cb/shiftFiles/OldHumanvsNewHumanrj7z",tp1 =tpH, tp2=list(tpMixSets[77:102]),name1 = "OldHuman",name2 = "NewHuman",tpm.threshold = 2,fc.threshold = 2,n.null = 100,gene.list = neural_list, open.start = T,ppt = F ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = list(MixSets$mm[,1:26]) ,set2 = setM, output.dir = "~/code/data/cb/shiftFiles/NewMousevsOldMouserj7z",tp1 =list(tpMixSets[1:26]), tp2=tpM,name1 = "NewMouse",name2 = "OldMouse",tpm.threshold = 2,fc.threshold = 2,n.null = 100,gene.list = neural_list, open.start = T,ppt = F ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))

run.DTW.Genes(set1 = setH ,set2 = list(MixSets$hs[,77:102]), output.dir = "~/code/data/cb/shiftFiles/OldHumanvsNewHumanrj1z",tp1 =tpH, tp2=list(tpMixSets[77:102]),name1 = "OldHuman",name2 = "NewHuman",tpm.threshold = 2,fc.threshold = 2,n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = list(MixSets$mm[,1:26]) ,set2 = setM, output.dir = "~/code/data/cb/shiftFiles/NewMousevsOldMouserj1z",tp1 =list(tpMixSets[1:26]), tp2=tpM,name1 = "NewMouse",name2 = "OldMouse",tpm.threshold = 2,fc.threshold = 2,n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = T,equalize=T,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="c",smoothed=T))



a+q + asdfklj

#run.DTW.Genes(set1 = test.setM ,set2 = test.setH, output.dir = "~/code/data/cb/shiftFiles/MousevsHInvitro",tp1 =tpM, tp2=tpH,gene.list = neural_list, open.start = F,ppt = T ,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F),z = F)

run.DTW.Genes(set1 = setM ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MousevsHInVitrorj7",tp1 =tpM, tp2=tpH,name1 = "Mouse",name2 = "Human",n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixvsHInVitrorj7",tp1 =tpMix, tp2=tpH,name1 = "Mix",name2 = "Human",n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setM ,set2 = setMix, output.dir = "~/code/data/cb/shiftFiles/MousevsMixInVitrorj7",tp1 =tpM, tp2=tpMix,name1 = "Mouse",name2 = "Mix",n.null = 100,gene.list = neural_list, open.start = F,ppt = F,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-2] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixShortvsHInVitrorj7",tp1 =tpMix[-2], tp2=tpH,name1 = "Mix17day",name2 = "Human",n.null = 100,gene.list = neural_list, open.start = F,ppt = F,z = F ,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixLongvsHInVitrorj7",tp1 =tpMix[-1], tp2=tpH,name1 = "Mix33day",name2 = "Human",n.null=100,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setMix[-2], output.dir = "~/code/data/cb/shiftFiles/MixLongvsMixShortInVitrorj7",tp1 =tpMix[-1], tp2=tpMix[-2],name1 = "Mix33day",name2 = "Mix17day",n.null=100,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setM[c(-1,-2)] ,set2 = setM[c(-3,-4)], output.dir = "~/code/data/cb/shiftFiles/M34vsM12InVitrorj7",tp1 =tpM[c(-1,-2)], tp2=tpM[c(-3,-4)],name1 = "Mouse34",name2 = "Mouse12",n.null=100,gene.list = neural_list, open.start = F,z = F,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))

run.DTW.Genes(set1 = setM ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MousevsHInVitrorj7bhatta",tp1 =tpM, tp2=tpH,name1 = "Mouse",name2 = "Human",n.null = 100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixvsHInVitrorj7bhatta",tp1 =tpMix, tp2=tpH,name1 = "Mix",name2 = "Human",n.null = 100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setM ,set2 = setMix, output.dir = "~/code/data/cb/shiftFiles/MousevsMixInVitrorj7bhatta",tp1 =tpM, tp2=tpMix,name1 = "Mouse",name2 = "Mix",n.null = 100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-2] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixShortvsHInVitrorj7bhatta",tp1 =tpMix[-2], tp2=tpH,name1 = "Mix17day",name2 = "Human",n.null = 100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F ,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixLongvsHInVitrorj7bhatta",tp1 =tpMix[-1], tp2=tpH,name1 = "Mix33day",name2 = "Human",n.null=100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setMix[-2], output.dir = "~/code/data/cb/shiftFiles/MixLongvsMixShortInVitrorj7bhatta",tp1 =tpMix[-1], tp2=tpMix[-2],name1 = "Mix33day",name2 = "Mix17day",n.null=100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setM[c(-1,-2)] ,set2 = setM[c(-3,-4)], output.dir = "~/code/data/cb/shiftFiles/M34vsM12InVitrorj7bhatta",tp1 =tpM[c(-1,-2)], tp2=tpM[c(-3,-4)],name1 = "Mouse34",name2 = "Mouse12",n.null=100,bhatta = T,gene.list = neural_list, open.start = F,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))

run.DTW.Genes(set1 = setM ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MousevsHInVitrorj7deriv",tp1 =tpM, tp2=tpH,name1 = "Mouse",name2 = "Human",n.null = 100, derivative = 1,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixvsHInVitrorj7deriv",tp1 =tpMix, tp2=tpH,name1 = "Mix",name2 = "Human",n.null = 100,derivative = 1,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setM ,set2 = setMix, output.dir = "~/code/data/cb/shiftFiles/MousevsMixInVitrorj7deriv",tp1 =tpM, tp2=tpMix,name1 = "Mouse",name2 = "Mix",n.null = 100,derivative = 1,gene.list = neural_list, open.start = F,ppt = F,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-2] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixShortvsHInVitrorj7deriv",tp1 =tpMix[-2], tp2=tpH,name1 = "Mix17day",name2 = "Human",n.null = 100,derivative = 1,gene.list = neural_list, open.start = F,z = T,ppt = F ,equalize=T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixLongvsHInVitrorj7deriv",tp1 =tpMix[-1], tp2=tpH,name1 = "Mix33day",name2 = "Human",n.null=100,derivative = 1,gene.list = neural_list, open.start = F,z = T,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setMix[-2], output.dir = "~/code/data/cb/shiftFiles/MixLongvsMixShortInVitrorj7deriv",tp1 =tpMix[-1], tp2=tpMix[-2],name1 = "Mix33day",name2 = "Mix17day",n.null=100,derivative = 1,gene.list = neural_list,z = T, open.start = F,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))
run.DTW.Genes(set1 = setM[c(-1,-2)] ,set2 = setM[c(-3,-4)], output.dir = "~/code/data/cb/shiftFiles/M34vsM12InVitrorj7deriv",tp1 =tpM[c(-1,-2)], tp2=tpM[c(-3,-4)],name1 = "Mouse34",name2 = "Mouse12",n.null=100,derivative = 1,gene.list = neural_list,z = T, open.start = F,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T))


run.DTW.Genes(set1 = setM ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MousevsHInVitrorj1",tp1 =tpM, tp2=tpH,name1 = "Mouse",name2 = "Human",n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))
run.DTW.Genes(set1 = setMix ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixvsHInVitrorj1",tp1 =tpMix, tp2=tpH,name1 = "Mix",name2 = "Human",n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))
run.DTW.Genes(set1 = setM ,set2 = setMix, output.dir = "~/code/data/cb/shiftFiles/MousevsMixInVitrorj1",tp1 =tpM, tp2=tpMix,name1 = "Mouse",name2 = "Mix",n.null = 100,gene.list = neural_list, open.start = F,ppt = F,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))
run.DTW.Genes(set1 = setMix[-2] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixShortvsHInVitrorj1",tp1 =tpMix[-2], tp2=tpH,name1 = "Mix17day",name2 = "Human",n.null = 100,gene.list = neural_list, open.start = F,ppt = F ,equalize=T,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixLongvsHInVitrorj1",tp1 =tpMix[-1], tp2=tpH,name1 = "Mix33day",name2 = "Human",n.null=100,gene.list = neural_list, open.start = F,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setMix[-2], output.dir = "~/code/data/cb/shiftFiles/MixLongvsMixShortInVitrorj1",tp1 =tpMix[-1], tp2=tpMix[-2],name1 = "Mix33day",name2 = "Mix17day",n.null=100,gene.list = neural_list, open.start = F,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))


run.DTW.Genes(set1 = setM ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MousevsHInVitrorj1bhatta",tp1 =tpM, tp2=tpH,name1 = "Mouse",name2 = "Human",n.null = 100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))
run.DTW.Genes(set1 = setMix ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixvsHInVitrorj1bhatta",tp1 =tpMix, tp2=tpH,name1 = "Mix",name2 = "Human",n.null = 100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F ,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))
run.DTW.Genes(set1 = setM ,set2 = setMix, output.dir = "~/code/data/cb/shiftFiles/MousevsMixInVitrorj1bhatta",tp1 =tpM, tp2=tpMix,name1 = "Mouse",name2 = "Mix",n.null = 100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F,z = F,equalize=T,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))
run.DTW.Genes(set1 = setMix[-2] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixShortvsHInVitrorj1bhatta",tp1 =tpMix[-2], tp2=tpH,name1 = "Mix17day",name2 = "Human",n.null = 100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F ,equalize=T,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixLongvsHInVitrorj1bhatta",tp1 =tpMix[-1], tp2=tpH,name1 = "Mix33day",name2 = "Human",n.null=100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))
run.DTW.Genes(set1 = setMix[-1] ,set2 = setMix[-2], output.dir = "~/code/data/cb/shiftFiles/MixLongvsMixShortInVitrorj1bhatta",tp1 =tpMix[-1], tp2=tpMix[-2],name1 = "Mix33day",name2 = "Mix17day",n.null=100,bhatta=T,gene.list = neural_list, open.start = F,ppt = F ,equalize=T,stepPattern =rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))

derivs <- lapply(1:length(fits),function(s){
  t(sapply(fits[[s]],function(f){
    predict(f$spline,x = f$mu$x,deriv = 1)$y
    
  }))
})
?diff()
a=derivs[[1]][order(rowMeans(derivs[[1]]),decreasing = T),]
heatmap.2(a[1:1000,],T,F,trace="none")

bhattacharyya.dist(1,3,.1,.1)

dnorm(.8 , 0, sqrt(.05))
(1/(.005*sqrt(2*pi)))*exp(-((.8^2) /(2*.005)))

hist(sapply(dt.objects[[1]],warpArea2))
dt.full[[1]]["TNIK",]
dtwPlot(dtw(fits[[1]][["SOX11"]]$mu$y,fits[[2]][["SOX11"]]$mu$y,step.pattern =  rabinerJuangStepPattern(1,"c")))#, fits[[1]][["SOX11"]]$mu$y,fits[[2]][["SOX11"]]$mu$y)
dtwPlot(dtw(diff(fits[[1]][["SOX11"]]$mu$y),diff(fits[[2]][["SOX11"]]$mu$y),step.pattern =   rabinerJuangStepPattern(1,"c")))#, diff(fits[[1]][["SOX11"]]$mu$y),diff(fits[[2]][["SOX11"]]$mu$y))
warpArea2(dt.objects[[1]][["TNIK"]])


cor(c(as.dist(outer(setMix[-1][[1]]["PAX6",], setMix[-1][[1]]["PAX6",], `-`) )),
c(as.dist(outer (setMix[-1][[1]]["NEUROG2",], setMix[-1][[1]]["PAX6",], `-`) )))
as.dist(outer (setMix[-1][[1]]["PAX6",], setMix[-1][[1]]["PAX6",], `-`) )
cor(setMix[-1][[1]]["PAX6",],setMix[-1][[1]]["NEUROG2",])

lapply(1:nrow(setMix[-1][[1]]),function(x){
  c(internal.dist(setMix[-1][[1]][x,]))
})

i=3
j=5


plot.model(fits[[1]][["ZIC1"]])
internal.dist <- function(x){
  as.dist(outer (x, x, `-`))
}

b=as.dist(outer (set1[[1]][1:7000,1],set1[[1]][1:7000,1],`-` ))
b2=as.dist(outer (set1[[1]][1:7000,7],set1[[1]][1:7000,7],`-` ))
cor(b,b2)
cor(set1[[1]][1:7000,1],set1[[1]][1:7000,7])
access.diag.mat <- function(i,j,mat){
  n=nrow(mat)
  if(i<j){
    mat[n*(i-1)-(i*(i+1)/2)+j]
  }else{
    print("wa")
   -(mat[n*(j-1)-(j*(j+1)/2)+i ])
  }
}

distdex<-function(i,j,n) #given row, column, and n, return index
  n*(i-1) - i*(i-1)/2 + j-i

rowcol<-function(ix,n) { #given index, return row and column
  nr=ceiling(n-(1+sqrt(1+4*(n^2-n-2*ix)))/2)
  nc=n-(2*n-nr+1)*nr/2+ix+nr
  cbind(nr,nc)
}



setMix[-1][[1]]["PAX6",]
pr_DB$get_entries()
setMash1 = list(exp_1_M, exp_2, datasets$`~/code/data/cb/human_ALONE.csv`)
setMash2 = list( exp_1[,17:24],exp_2_H,datasets$`~/code/data/cb/mouse_ALONE.csv`,exp_1[,10:16])
tpMash1 =list(tp_1_M,tp_2,tp.h.control)
tpMash2 =list(tp_1[17:24],tp_2,tp.m.control,tp_1[10:16])
run.DTW.Genes(set1 = setH[-3] ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/HvsHInVitrorj1",tp1 =tpH[-3], tp2=tpH,name1 = "HumanShort",name2 = "Human",gene.list = neural_list, open.start = F,ppt = F,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="b",smoothed=F))
run.DTW.Genes(set1 = setMash1 ,set2 = setMash2, output.dir = "~/code/data/cb/shiftFiles/Mashuprj1",tp1 =tpMash1, tp2=tpMash2,name1 = "Mash1",name2 = "Mash2",gene.list = neural_list, open.start = F,ppt = F,stepPattern = rabinerJuangStepPattern(type=1,slope.weighting="b",smoothed=F))


a=(1-range01(rowMeans(t(sapply(fits[[1]], function(x) x[["coefficients"]][1,1:6]))))^.3)
b=(1-range01(rowMeans(t(sapply(fits[[2]], function(x) x[["coefficients"]][1,1:10]))))^.3)
c=(1-range01(sumErrors[[3]]))
d=(range01(rowMeans(dt[[1]][,1:16])))
e=(1-range01(dt.dist[[3]]))
measure <- a*b*c*d*e
measure <- measure[ind]
om <-  measure[order(measure,decreasing = T)]
om <- om
write.pdf({
  mypar(3,4)
  for(i in 1:24){
    gn=names(om)[i]
    print(paste(gn, sapply(list(a,b,c,d,e), function(x) x[gn])))
    plot.model( fits[[1]][[gn]],main=paste(gn,"Mouse Log2 TPM"))
    plot.model( fits[[2]][[gn]],main=paste(gn,"Human Log2 TPM"))
  }
}  , file.path(datapath,"neural_late_tpm.pdf"))



gn="NEUROG2"
s=1
plot.model( fits[[1]][[gn]],type = "model",showConfidenceBands=TRUE,showIndividual=F,main=paste(gn,"Mouse Log2 TPM"))
plot.model( fits[[2]][[gn]],type = "model",showConfidenceBands=TRUE,showIndividual=F,main=paste(gn,"Human Log2 TPM"))
fd <- dist.fit(fits[[s]][[gn]],fits[[pairs[s]]][[gn]])
rfd <- dist.fit(fits[[pairs[s]]][[gn]],fits[[s]][[gn]])
exp1 <- fits[[s]][[gn]][["coefficients"]][1,]
exp2 <- fits[[pairs[s]]][[gn]][["coefficients"]][1,]
dtq1 <-  dtw(fd[[1]], open.end = T , open.begin = T, step.pattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T),keep.internals = T)
rdtq1 <-  dtw(rfd[[1]], open.end = T , open.begin = T, step.pattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=T),keep.internals = T)
sumErrors <-mean(sapply(1:length(dtq1$index1),function(i) fd[[2]][dtq1$index1[i],dtq1$index2[i]])/mean(exp1[dtq1$index1[i]],exp2[dtq1$index2[i]]))
dtwPlotThreeWay(rdtq1,exp2,exp1)
warp(rdtq1,index.reference = T)
dtwPlotThreeWay(dtq1,exp1,exp2)
dtwPlotDensity(dtq1)
warp(dtq1,index.reference = T)
plot(rabinerJuangStepPattern(type=1,slope.weighting="c",smoothed=F))
#H Fwd
sapply(1:length(warp(dtq1,index.reference = T)),function(i){
  w <- warp(dtq1,index.reference = T)
  as.numeric(colnames(fits[[pairs[s]]][[1]][["coefficients"]])[w[i]])-as.numeric(colnames(fits[[s]][[1]][["coefficients"]])[i])
  #as.numeric(colnames(fits[[s]][[1]][["coefficients"]])[w[i]])- as.numeric(colnames(fits[[s]][[1]][["coefficients"]])[i])
})
#M
sapply(1:length(warp(rdtq1,index.reference = T)),function(i){
  w <- warp(rdtq1,index.reference = T)
  print(w)
  print(as.numeric(colnames(fits[[pairs[s]]][[1]][["coefficients"]])))
  print(as.numeric(colnames(fits[[s]][[1]][["coefficients"]])))
  as.numeric(colnames(fits[[pairs[s]]][[1]][["coefficients"]])[w[i]])-as.numeric(colnames(fits[[s]][[1]][["coefficients"]])[i])
  #as.numeric(colnames(fits[[s]][[1]][["coefficients"]])[w[i]])- as.numeric(colnames(fits[[s]][[1]][["coefficients"]])[i])
}) 


#nl=neural_list[!(neural_list%in%rownames(q.b.f.e)) & neural_list%in%names(fits[[1]])] 
gn <- 1

dt[[1]][gn,]
met1[gn,]
q.b.f.e <- rbind( q.b.f.e, c(met1[gn,], T))
rownames(q.b.f.e)[nrow(q.b.f.e)] <- gn
i=i+1
#q.b.f.e[gn,ncol(q.b.f.e)] <- F
q.b.f.e <- q.b.f.e[-nrow(q.b.f.e),]
q.b.f
#q.b.f <- cbind(t(met1[gn,]),T)
#rownames(q.b.f)[106:112] <- gene.list[(i-6):i]
#saveRDS(q.b.f.e,file = file.path(datapath,"qbf.RDS"))
q.b.f <- readRDS(file.path(datapath,"qbf.RDS"))
ind<-rownames(q.b.f)[(rownames(q.b.f) %in%names(fits[[1]]))]
q.b.f.e <- cbind(met1[ind,],q.b.f[ind,7])

metrics <-  list(p.error, p.lik, p.dist,range01(rancors[[1]]),xcors[[1]],dt.fc)

#   rownames(q.b.f) <- sapply(1:nrow(q.b.f),function(i){
#     names(which(unlist(sumErrors)==q.b.f[i,1]))
#   })
q.b.f
met1
f.values[[2]]

met1 <- opti[rownames(q.b.f),-ncol(opti)]
colnames(met1) <- nmet
met3 <- opti[rownames(q.b.f),-ncol(opti)]
colnames(met3) <- nmet

q.b.f.e <- cbind(met1[ind,],q.b.f[ind,7])
q.b.f.e
fff <-  randomForest(x = q.b.f.e[,-c(ncol(q.b.f.e))],y = as.factor(unlist(q.b.f.e[,ncol(q.b.f.e)])),mtry = 2,ntree = 1000,keep.forest = T)
fff$importance
fff
met.data <- q.b.f.e#[(q.b.f.e[,8]<.0025 &q.b.f.e[,6]>1 ),-c(6,8)]
met.data <- matrix(unlist(met.data),nrow = nrow(met.data),ncol=ncol(met.data),dimnames = dimnames(met.data))
rownames(met.data) <- NULL
nb= naiveBayes(x=met.data[,-c(ncol(met.data))],as.factor(unlist(met.data[,7])))
nb

mypar(3,4)
for(gn in names(p.error[100:200])[order(p.error[100:200],decreasing = T)]){
  plot.model(fits[[1]][[gn]],main=logP[[1]][gn])
  plot.model(fits[[2]][[gn]],main=logP[[3]][gn])
}

mypar(3,4)
for(gn in names(p.lik[100:300])[order(p.lik[100:300],decreasing = T)]){
  plot.model(fits[[1]][[gn]],main=dt.likelihood[[1]][gn])
  plot.model(fits[[2]][[gn]],main=p.lik[gn])
}

mypar(3,4)
for(gn in names(p.xc)[order(p.xc,decreasing = T)]){
  plot.model(fits[[1]][[gn]])
  plot.model(fits[[2]][[gn]])
}


metrics <-  list(sumErrors[[1]], dt.likelihood[[1]],dt.dist[[1]],rancors[[1]],xcors[[1]],dt.fc)
metrics <-  list(p.error, p.lik, rowMins(cbind(p.dist,p.dist2)),range01(rancors[[1]]),xcors[[1]],dt.fc)
metrics <-  list(rowMins(cbind(sumErrors[[1]],sumErrors[[2]])), rowMeans(cbind(dt.likelihood[[1]], dt.likelihood[[2]])),rowMins(cbind(dt.dist[[1]],dt.dist[[2]])),rancors[[1]],xcors[[1]],dt.fc,sapply( 1:length(f.values[[1]]),function(i) mean( f.values[[1]][i],f.values[[2]][i])))
metrics <-  list(rowMins(cbind(p.error,p.error2)), rowMeans(cbind(p.lik, p.lik2)),rowMins(cbind(p.dist,p.dist2)),range01(rancors[[1]]),xcors[[1]],dt.fc)
nmet <-  c("sumErrors", "dt.likelihood","dt.dist","cor.s","cor.x","dt.fc")
met1 <- metrics
met1 <- Reduce(cbind, met1)
colnames(met1) <- nmet
q.b.f.e
#<- colnames(opti.n)
head(met1)
met1 <- opti[rownames(q.b.f),-ncol(opti)]
q.b.f.e <- cbind(met1[ind,],q.b.f[ind,7])

fff <-  randomForest(x = q.b.f.e[,-c(ncol(q.b.f.e))],y = as.factor(unlist(q.b.f.e[,ncol(q.b.f.e)])),mtry = 2,ntree = 1000,keep.forest = T)
fff$importance
fff
q.b.f.e
met.data <- q.b.f.e
met.data <- matrix(unlist(met.data),nrow = nrow(met.data),ncol=ncol(met.data),dimnames = dimnames(met.data))
rownames(met.data) <- NULL
nb= naiveBayes(x=met.data[,-c(ncol(met.data))],as.factor(unlist(met.data[,7])),laplace = .1)
nb
xval(met.data)

hist(sumErrors[[1]])
met.data
q.b.f





fff <-  randomForest(x = q.b.f[,-c(7,ncol(q.b.f))],y = as.factor(unlist(q.b.f[,ncol(q.b.f)])),mtry = 2,ntree = 2000,keep.forest = T)
fff <-  randomForest(x = q.b.f.e[,-c(ncol(q.b.f.e))],y = as.factor(unlist(q.b.f.e[,ncol(q.b.f.e)])),mtry = 2,ntree = 5000,keep.forest = T)
fff$importance
fff
q.b.f.e
mean( (unlist(q.b.f.e[,3]) < .5 ) == q.b.f.e[,8])
mean((unlist(q.b.f.e[,3]) < .55 ) | abs((unlist(q.b.f.e[,2])) > .22 ) == q.b.f.e[,8])

q.b.f.e[(q.b.f.e[,8]<.0025 &q.b.f.e[,6]>1 ),-c(6,8)]
plot.model(fits[[2]][["ACAD11"]])
hist(dt.likelihood[[3]])
met.data
#saveRDS(fff,file = file.path(datapath,"randomforestnoz.RDS"))
#saveRDS(nb,file = file.path(datapath,"nbayes.RDS"))

xval <- function(q.b.f.e){
  best1 <- list(0)
  for(xx in 1:10){
    print(xx)
    met.data <- q.b.f.e#[(q.b.f.e[,8]<.0025 &q.b.f.e[,6]>1 ),-c(6,8)]
    met.data <- matrix(unlist(met.data),nrow = nrow(met.data),ncol=ncol(met.data),dimnames = dimnames(met.data))
    rownames(met.data) <- NULL
    #Randomly shuffle the datasetsnoz[[2]][y,]
    met.data<-met.data[sample(nrow(met.data)),]
    #Create 10 equally size folds
    folds <- cut(seq(1,nrow(met.data)),breaks=10,labels=FALSE)
    accs <- list()
    for(q in 1:4){
      accs[[q]] <- vector(mode="list",length = 0)
    }
    #Perform 10 fold cross validation
    for(i in 1:10){
      #Segement your data by fold using the which() function 
      testIndexes <- which(folds==i,arr.ind=TRUE)
      testData <- met.data[testIndexes, ]
      trainData <- met.data[-testIndexes, ]
      acc <- array()
      www <- ncol(testData)
      fac <- as.factor(unlist(trainData[,www]))
      fff <-  randomForest(x = trainData[,-www],y = fac,mtry = 2,ntree = 1000)
      ff <-  C5.0(x=trainData[,-www],fac,trials =10,costs=matrix( c(0,1,1,0),2,2,dimnames =list(list(0,1),list(0,1))),rules = T)
      acc[1]=mean(predict(object=fff, testData[,-www]) == testData[,www])
      acc[2]=mean(levels(predict.C5.0(ff, testData[,-www]))[predict.C5.0(ff, testData[,-www])]==testData[,www])
      nb <- naiveBayes(x = trainData[,-www],y = fac,laplace = .1)
      acc[3] <- mean(levels(predict(nb,testData[,-www] ))[predict(nb,testData[,-www])]==testData[,www])
      fit4 <- PART( fac~. ,data = as.data.frame(trainData[,-www]) )
      acc[4]<- mean((predict(fit4,  as.data.frame(testData[,-www]))) == testData[,www])
      for(q in 1:length(acc)){
        accs[[q]] <- c(accs[[q]],acc[q])
      }
      best1[[xx]] <- lapply(lapply(accs,unlist), mean)
    }
  }
  print(mean(unlist(lapply(best1,"[[",1))))
  print(mean(unlist(lapply(best1,"[[",2))))
  print(mean(unlist(lapply(best1,"[[",3))))
  print(mean(unlist(lapply(best1,"[[",4))))
}
min(unlist(lapply(best1,"[[",1)))
min(unlist(lapply(best1,"[[",2)))
min(unlist(lapply(best1,"[[",3)))
min(unlist(lapply(best1,"[[",4)))



frame
a=nle(y~tme,frame[[1]][[1]])
var.test(c(5,5.25,5.5), c(1,2,1.2))
pf(0.8928571,2,1 )
t.test(c(5,4), c(1.2,2))
frame[[1]][[1]]

library(nlme)
library(lme4)
library(lattice)
data(Ovary)
attach(Ovary)
names(Ovary)
plot(Ovary)
Ovary
model<-lme(follicles~sin(2*pi*Time)+cos(2*pi*Time),
           data=Ovary,random=~ 1| Mare)
summary(model)
plot(ACF(model),alpha=0.05)
model2<-update(model,correlation=corARMA(q=2))
anova(model,model2)
model3<-update(model2,correlation=corAR1())
anova(model2,model3)
plot(model3,resid(.,type="p")~fitted(.)|Mare)
qqnorm(model3,~resid(.)|Mare)

data("Rail")
?lmer
r2.lme <- lmer(travel~1+(1 | Rail),data = Rail,REML = F )
summary(r2.lme)
ranef(r2.lme)

lme4::nl
?lme
plot(rabinerJuangStepPattern(type=3,slope.weighting="d",smoothed=T))

ind <- intersect(rownames(exp_1_M), rownames(exp_2_M))
dtt <-  dtw((1-cor( exp_1_M[ind,],exp_2_M[ind,] ,method = "spearman")), open.end = T,open.begin = F,step.pattern = rabinerJuangStepPattern(type=3,slope.weighting="c",smoothed=T),keep.internals = T)
dtwPlot(dtt)
dtwPlotDensity(dtt)
tp_2[as.integer(warp(dtt, index.reference = T))]
tp_2[apply(1-cor( exp_1_M[ind,],exp_2_M[ind,] ,method = "spearman"),1, which.min)]
tp_1[1:9]
tp_2
tp.h.control

y <- "NEUROD4"
y = "SOX2"
y= 12
y=i
a <- fits[[1]][[y]]$mu$y
b <- fits[[2]][[y]]$mu$y
b = datasets$`~/code/data/cb/human_ALONE.csv`[y,]
plot(a,type="l")
plot(b, type ="l")
cor(a,b ,method = "spe")
ccf(a,b,plot = F)
gene.list

#plot(rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))
dtt <-  dtw( b,a, open.end = T , open.begin = F ,dist.method = "euclidean",step.pattern = rabinerJuangStepPattern(type=1,slope.weighting="d",smoothed=T))
dtwPlotThreeWay(dtt,xts = b,a)
dtwPlotTwoWay(dtt,xts = b,a)
dtwPlotDensity(dtt)
dtt$normalizedDistance
dtt$index2
m <- dtt$localCostMatrix
rm <- rdtt$localCostMatrix
pm <- exp(-m)
#pm[is.na(pm)] <- 0
pm <- sweep( pm,1,STATS = rowSums(pm,na.rm = T),FUN="/")
sapply(1:length(dtt$index1),function(i) pm[dtt$index1[i],dtt$index2[i]])
heatmap.2(pm,Rowv = F,Colv = F,trace = "none")
heatmap.2(m,Rowv = F,Colv = F,trace = "none")
rdtt$distance
sum(sapply(1:length(dtt$index1),function(i) m[dtt$index1[i],dtt$index2[i]]),na.rm = T)
sum(sapply(1:length(rdtt$index1),function(i) rm[rdtt$index1[i],rdtt$index2[i]]),na.rm = T)
rdtt$distance

(-1/length(dtt$index1))*log(prod(sapply(1:length(dtt$index1),function(i) pm[dtt$index1[i],dtt$index2[i]])))
mean(sapply(1:length(dtt$index1),function(i) m[dtt$index1[i],dtt$index2[i]]))

w1 <-warp(dtt,index.reference = T)-seq(1,17)
rdtt <-  dtw(b, a, open.end = T , open.begin = T ,step.pattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F),keep.internals = T)
dtwPlotThreeWay(rdtt,xts = b, a)

w2 <- warp(rdtt,index.reference = T)-seq(1,26)

check.align(w1,w2)
check.align(w2,w1)


w1
w2

dtwPlotDensity(dtt)
warp(dtt,index.reference = T)-seq(1,26)
dtt$normalizedDistance
?dist
?dtw
?dist
dtwDist(a,b)
?dtwDist()
dtt$index2
dtt$index1
as.matrix( dtw::warp(dtt,index.reference = T)-seq(1,17))
?wilcox.test()
?wilcoxGST()
typesOfsteps[[14]]

dt.adj[x,]
shift(t(a),t(b),tp_2,tp_2, seq(1,33), seq(1,33))
plot(exp_2[x,],type="l")
which.max(exp_2_H[x,])
