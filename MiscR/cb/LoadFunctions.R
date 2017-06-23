library(RColorBrewer)
library(matrixStats)
library(ggplot2)
library(gplots)
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

hist.normalized <- function(x){
  nsamp <- dim(x)[2]
  h <- hist(x[,1], plot=FALSE)
  plot(h$mids, h$density, type="l", col=rainbow(nsamp)[1], main="",
       xlab="expression value", ylab="Proportion of molecules")
  for(i in 2:nsamp){
    h <- hist(x[,i], plot=FALSE)
    lines(h$mids, h$density, col=rainbow(nsamp)[i])
  }
}

rn.merge <- function(x,y,fill=0){
  rn <- intersect(rownames(x),rownames(y))
  zerosx <- setdiff(rownames(x),rownames(y))
  zerosy <- setdiff(rownames(y),rownames(x))
  out <- cbind(x[rn,],y[rn,])
  if(length(zerosx)!=1  &length(zerosy)!=1){
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

rn.merge.normalize <- function(x,y,fill=0){
  rn <- intersect(rownames(x),rownames(y))
  zerosx <- setdiff(rownames(x),rownames(y))
  zerosy <- setdiff(rownames(y),rownames(x))
  setI <- c(rep(1,ncol(x)) , rep(2,ncol(y)))
  out <- cbind(x[rn,],y[rn,])
  if(length(zerosx)!=1  &length(zerosy)!=1){
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
  out <- GetNormalizedMat(out,MedianNorm(out))
  return(list(out[,setI==1],out[,setI==2]))
}

median.normalize <- function(x){
  GetNormalizedMat(x,MedianNorm(x))
}



rn.compare <- function(x,y,fill=0){
  rn <- as.character(intersect(rownames(x),rownames(y)))
  zerosx <- setdiff(rownames(x),rownames(y))
  zerosy <- setdiff(rownames(y),rownames(x))
  zx <- matrix(fill, nrow=length(zerosx), ncol =ncol(y), dimnames = list(zerosx,NULL))
  zy <- matrix(fill, nrow=length(zerosy), ncol =ncol(x), dimnames = list(zerosy,NULL))
  nx <- rbind(x,zy)
  ny <- rbind(y,zx)
  return(list(nx[rownames(ny),],ny[rownames(nx),]))
}

match.subtract <- function(x,y,abs=F){
  i = intersect(names(x),names(y))
  noni <- union(names(which(table(names(x))>1)),names(which(table(names(y))>1)))
  i <- setdiff(i,noni)
  a <- x[i]-y[i]
  a <- c(a,unlist(sapply(noni,function(n){
    a=x[names(x)==n] 
    b=y[names(y)==n]  
    inds=expand.grid(1:length(a),1:length(b))
    a[inds[,1]]-b[inds[,2]]
  })))
  return(a)
}

std.heatmap <- function(M,...){
  heatmap.2(M,Rowv = F,Colv = F,trace="none",col = cols,...)
}

round.log <- function(s,base=2){
  round(log(s+1, base),digits = 1)
}

go.Fisher.Test <- function(set, spec = "hs",ont = "BP"){
  if(spec=="mm"){
    libspec = "mouse"
    xx <- as.list(org.Mm.egALIAS2EG)
    names(xx) <- toupper(names(xx))
  }else{
    libspec = "human"
    xx <- as.list(org.Hs.egALIAS2EG)
    names(xx) <- toupper(names(xx))
  }
  xx
  summary(enrichGO(gene= sapply(xx[set],"[",1), organism= libspec,universe = unlist(sapply(xx[all.genes],"[",1)),ont = ont))
}

cor.compare <- function(x,y,min=0, varX=NULL ,interest.set =NULL, ...){
  d <- rn.compare(x,y)
  x <- d[[1]]
  y <- d[[2]]
  i = intersect(rownames(x), rownames(y))
  i = i[rowMaxs(x[i,],na.rm = T)>=min | rowMaxs(y[i,],na.rm = )>=min]
  if(!is.null(interest.set)){
    i = interest.set[interest.set%in%i]
  }
  if(!is.null(varX)){
    i = i[order(rowMeans(cbind(rowSds(x[i,]), rowSds(y[i,]))),decreasing = T)]
    i = i[1:ifelse(varX>length(i),length(i),varX)]
  }
  print("Num Genes:")
  print(length(i))
  return(cor(as.matrix(x[i,]),as.matrix(y[i,]), ...))
}

md <- function(s,n)s[n[n%in%rownames(s)],]


mypar <- function (a = 1, b = 1, brewer.n = 8, brewer.name = "Dark2", 
          cex.lab = 1, cex.main = 1.2, cex.axis = 1, mar = c(2.5, 2.5, 
                                                             1.6, 1.1), mgp = c(1.5, 0.5, 0), ...) 
{
  par(mar = mar, mgp = mgp, cex.lab = cex.lab, cex.main = cex.main, 
      cex.axis = cex.axis)
  par(mfrow = c(a, b), ...)
  palette(RColorBrewer::brewer.pal(brewer.n, brewer.name))
}

getMeanExprs <- function(x){
  coef(x)[1,]
}

write.pdf <- function(plot, filename){
  pdf(filename,width = 20, height = 20)
  plot
  dev.off()
}

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

get.ts.musd <- function(set,tp){
  m <- lapply(1:nrow(set),function(j){
    data.frame("y"=set[j,], "tme"=tp.i)#,"ind"=i)
  })
  #if(z) f[,"y"] <- f[,"y"]/max(f[,"y"])
  names(m) <- rownames(set)
  lapply(m,ts.mean.sd,tp )
}

ts.mean.sd <- function(m,tp.i){
  a=anova(lm(formula=y~tme, data=m))
  mu = sapply(sort(unique(m$tme)),function(x){
    list("x"=x,"y"= mean(as.numeric(m[m$tme==x,]$y)))
  })
  names(mu) <- sort(unique(m$tme))[sort(unique(m$tme))%in% tp.i]
  sd = sapply(sort(unique(m$tme)),function(x){
    sd(as.numeric(m[m$tme==x,]$y))
  })
  names(sd) <- sort(unique(m$tme))[sort(unique(m$tme))%in% tp.i]
  sd = sd[!is.na(names(sd))]
  ind = tp.i%in%names(sd)
  names(mu$y) <- mu$x
  list("mu"=mu , "anova"=a,"data"= m,"sd"= sd)
}

range01 <- function(x){ (x-min(x))/(max(x)-min(x))}

range01Max <- function(x){ x/max(x)}


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

DTW.dist <- function(my.mat){
  n <- nrow(my.mat)
  mat <- matrix(0, ncol = n, nrow = n)
  colnames(mat) <- rownames(mat) <- rownames(my.mat)
  auc.mat <- matrix(0, ncol = n, nrow = n)
  colnames(auc.mat) <- rownames(auc.mat) <- rownames(my.mat)
  o <- mclapply( 1:(nrow(my.mat)-1),function(i){
    mclapply((1+i):nrow(my.mat), function(j) {
      a = dtw(my.mat[i,],my.mat[j,])
      list(a$distance, warpArea2(a))
    })
  })
  for(i in 1:length(o)) {
    for(j in 1:(length(o[[i]]))) {
      mat[i, i+j] <-  o[[i]][[j]][[1]]
      auc.mat[i,i+j] <- o[[i]][[j]][[2]]
    }
  }
  return(list( "dist"=as.dist(t(mat)),"AUC" =as.dist(t(auc.mat))))
}

#system.time(DTW.dist(my.mat = datasets$`~/code/data/cb/human_ALONE.csv`[1:100,]))
#system.time(system.time( dist(datasets$`~/code/data/cb/human_ALONE.csv`[1:100,],method =  pr_DB$get_entry("DTW"))))

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

#######################################################


library(scales)
#which function to use to rescale genes 0-1
if(!min0){
  range01 <- function(x){
    z <- x/max(x,na.rm = T)
    sapply( z, function(i) ifelse(is.nan(i),0,i) )
  }
}else{
  range01 <- function(x){
    z <- (x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))
    sapply( z, function(i) ifelse(is.nan(i),0,i) )
    }
}

mz <- lapply(geneList, function(gn){
  if(gn %in% Reduce(intersect, lapply(sets,rownames))){
    xs <- lapply(1:length(sets), function(i){
      xi <- zoo(sets[[i]][gn,], tps[[i]])
      xi <- merge(xi,zoo(NA, tps[[i]]))[,1]
      xi
    })
    na.approx(Reduce(merge,xs))
  }
})

mzn <- geneList
mzn=mzn[!sapply(mz, is.null)] 
mz=mz[!sapply(mz, is.null)] 
if(!is.null(sortMax))mz <-  mz[order(sapply(mz,function(x) which.max(rescale(x[,sortMax]))))]

mz=Reduce(cbind, mz )
colnames(mz) <- paste0(species ,rep(mzn,1,each=length(sets)))
mz=apply(mz,2,range01)
#mz=sweep(mz,MARGIN = 2,colMaxs(mz,na.rm = T),FUN = "/")
for(s in rev(species)){
  mz <- cbind(c(seq(0,1,length.out = 11), rep(NA,nrow(mz)-11 )),mz)
  colnames(mz)[1] <-  paste("scale", s)
}

unified.tp <- Reduce(union,tps)
blocksize.vector.x <- diff(unified.tp)*(1/max(unified.tp))
blocksize.vector.x <- c(blocksize.vector.x[1],blocksize.vector.x,blocksize.vector.x[length(blocksize.vector.x)])
b.v.x <- sapply(2:length(blocksize.vector.x),function(i) blocksize.vector.x[i-1]/2 + blocksize.vector.x[i]/2 )
b.v.x <- b.v.x * (max(unified.tp)*sum(b.v.x)/sum(b.v.x))
tpx.adj <- unified.tp
tpx.n <<- tpx.adj
tpx.adj[2:(length(tpx.adj))] <- as.vector(sapply(2:(length(tpx.adj)),function(i){
  q <- tpx.n[i-1]+ .5*b.v.x[i-1] + .5*b.v.x[i]
  tpx.n[i] <<- q
  q
} ))

rownames(mz) <- tpx.adj

for(i in seq(ncol(mz)- length(species),1,by = -length(species))){
  df <- data.frame(rep(NA,nrow(mz)))
  colnames(df) <- paste0(rep(" ",i),collapse = "")
  mz <- cbind(mz[,1:i], df , mz[,(i+1):ncol(mz)] )
}

a = melt(as.matrix(mz),na.rm = F)
groupVals <- sapply( a$Var2,function(i) which( sapply(species,function(s)grepl(s,i))))
groupVals[sapply(groupVals,length)==0] <- 0
a=cbind(a, "group"= unlist(groupVals))
a$Var2 <- factor(a$Var2, levels=unique(a$Var2[order(a$Var2,decreasing = T)]),exclude = NULL)
a$rescaleoffset <- a$value + 100*(a$group)

scalerange <- range(a$value,na.rm = T)
gradientends <- scalerange + rep( seq(0, (length(sets)-1)*100, by=100), each=2)
colorends <- unlist(lapply(colors[1:length(sets)], function(co) c("black", co)))


g=ggplot(a, aes(Var1, Var2)) + 
  ggtitle(paste( paste(species,collapse = ' vs ') ,"\n( Scaled Expression)"))+
  geom_tile(aes(x = Var1,fill = rescaleoffset,width=rep(b.v.x ,nrow(a)/max(sapply(tps, length) )))) + 
  scale_fill_gradientn(colours = colorends, values = rescale(gradientends),na.value = "grey55",guide = F,name="Expression",breaks=gradientends) + 
  #scale_x_discrete("Days", expand = c(0.05, .05)) + 
  scale_y_discrete("Gene", expand = c(.05, .05)) +
  theme_bw(base_size = 9) +  
  theme(legend.position = "right",
        panel.background = element_rect(fill = "grey55", size = 2),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        #panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 330, hjust = 0))
g
