source("http://bioconductor.org/biocLite.R")
library(ggplot2)
library(RColorBrewer)
library(gplots)
#biocLite("EACI")
#library(EACI)
library(matrixStats)
library(dtw)
 library("amap")
# library(genefilter)
# library("ReporteRs")
# library(randomForest)
# library(C50)
# library(RWeka)
# library(e1071)
options("mc.cores"=3)
library(parallel)
library(doParallel)
registerDoParallel()
library(sme)
library(calibrate)
source("~/code/cb/LoadData.R")
datapath ="~/code/data/cb"

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
      #d[i,j] <-max(dist(cf1[1,i]+cv1[i], cf2[1,j]-cv2[j]),dist(cf1[1,i]-cv1[i], cf2[1,j]+cv2[j]) )
      v[i,j] <- cv1[i]+cv2[j]
    }
  }
  return(list(d,v))
}

generateCorrelatedNulls <- function(d1, d2, target.cor= 0.0,length.multiplier=1){
  nd1 <-  data.frame(matrix(nrow=0 ,ncol = ncol(d1)))
  nd2 <- data.frame(matrix(nrow=0 ,ncol = ncol(d2)))
  while(nrow(nd2)< as.integer(nrow(d2)*length.multiplier)){
    nd1r <- d1[sample(length(d1),size = ncol(d1))]
    nd2r <- d2[sample(length(d2),size = ncol(d2))]
    #print(sum(is.na(d1)))
    #print(sum(is.na(d2)))
    #if((abs(cor(nd1r,nd2r , method = "spe")) > target.cor) & (aad(nd1r)>=aad.Threshold.1) & (aad(nd2r)>=aad.Threshold.2)){
    if(max(ccf(nd1r,nd2r,plot = F,na.action = na.omit)$acf) > target.cor & ){
      nd1[(nrow(nd1)+1),] <- as.vector(nd1r)
      nd2[(nrow(nd2)+1),] <- as.vector(nd2r)
    }
  }
  return(list(as.matrix(nd1),as.matrix( nd2)))
}


range01 <- function(x){ (x-min(x))/(max(x)-min(x))}


########
#set1 = list(exp_2,exp_1[,17:24])
set1 = list(exp_1_M, exp_2_M, datasets$`~/code/data/cb/mouse_ALONE.csv`)
set2 = list( exp_1[,1:9],exp_2_H, datasets$`~/code/data/cb/human_ALONE.csv`)
set1=test.setM
set2=test.setH
tp1 =list(tp_1_M,tp_2,tp.m.control)
tp2 =list(tp_1[1:9],tp_2,tp.h.control)
output.dir = "~/code/data/cb/shiftFiles/MvsHInvitro"
gene.list = neural_list
null.length = 1
cor.Threshold = 0
fc.threshold = 2
tpm.threshold = 2
threshold.num =3
# #tp2=tp.h.control
# tp1.i =seq(0,max(tp_2))
# tp2.i =seq(0,max(tp_2))
#tp2.i=seq(0,max(tp.h.control))
open.start = F
stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F)
z = F
#plot(colSums(sapply(seq(0,1,.01),function(i){abs(rancors)>i} )[gene.list,] ))
######
####DT output from warp is the direction and magnitude of shift for set2 at each tp
#Larger aad.Thresh means more genes get in
#Set1 should be the longer of the two sequences. Called Experimental, query by dtw. "Forward shift" is set1 shift positively
run.DTW.Genes <-  function(set1, set2, output.dir , tp1, tp2, gene.list=neural_list,  open.start=F, threshold.num=3,fc.threshold= 2, tpm.threshold=2,ppt=T , stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F), z =F){
  dir.create(output.dir)
  all.rn <- as.character(Reduce(intersect, c(lapply(set1,rownames),lapply(set2, rownames))))
  set1 <- lapply(set1, function(n) n[all.rn,])
  set2 <- lapply(set2, function(n) n[all.rn,])
  s1.fc <-  lapply(set1, function(x) rowMaxs(x)-rowMins(x)+1 )
  s2.fc <-  lapply(set2, function(x) rowMaxs(x)-rowMins(x)+1 )
  pass.thresh <-  Reduce( "+" ,lapply(c(s1.fc,s2.fc), ">",fc.threshold)) >threshold.num
  set1 <- lapply(set1, function(n) n[pass.thresh,])
  set2 <- lapply(set2, function(n) n[pass.thresh,])
  s1.mx <-  lapply(set1, function(x) rowMaxs(x) )
  s2.mx <-  lapply(set2, function(x) rowMaxs(x) )
  pass.thresh <-  Reduce( "+" ,lapply(c(s1.mx,s2.mx), ">",tpm.threshold)) >=6
  set1 <- lapply(set1, function(n) n[pass.thresh,])
  set2 <- lapply(set2, function(n) n[pass.thresh,])
  if(z){
    set1 <- lapply(set1, function(n) n/rowMaxs(n))
    set2 <- lapply(set2, function(n) n/rowMaxs(n))
  }
  
  tp <- list(tp1,tp2)
  set <- list(set1,set2)
  frame <- lapply(1:length(set), function(x){
    q <-  lapply(1:nrow(set[[x]][[1]]), function(k){
      as.data.frame(Reduce(rbind,lapply(1:length(set[[x]]),function(i){
        t(sapply(1:ncol(set[[x]][[i]]),function(j){
          c("y"=set[[x]][[i]][k,j], "tme"=tp[[x]][[i]][j],"ind"=i)
        }))
      })))
    })
    names(q) <- rownames(set[[1]][[1]])
    q
  })
  
  null.frame <- lapply(1:length(set), function(x){
    q <-  lapply(1:nrow(set[[x]][[1]]), function(k){
      as.data.frame(Reduce(rbind,lapply(1:length(set[[x]]),function(i){
        t(sapply(1:ncol(set[[x]][[i]]),function(j){
          c("y"=sample( set[[x]][[i]][k,],1), "tme"=tp[[x]][[i]][j],"ind"=i)
        }))
      })))
    })
    names(q) <- rownames(set[[1]][[1]])
    q
  })
  
  if(file.exists(file.path(output.dir,"fits.RDS"))){
    fits <- readRDS(file.path(output.dir,"fits.RDS"))
  }else{
  fits1 <- mclapply(frame[[1]], sme ,lambda.mu = 1,lambda.v = 1,maxIter = 500)
  fits2 <- mclapply(frame[[2]], sme ,lambda.mu = 1,lambda.v = 1,maxIter = 500)
  fitsn.1 <- mclapply( null.frame[[1]], sme ,lambda.mu = 1,lambda.v = 1,maxIter = 500)
  fitsn.2 <- mclapply(null.frame[[2]], sme ,lambda.mu = 1,lambda.v = 1,maxIter = 500)
  fits <- list(fits1, fits2,fitsn.1,fitsn.2)
  saveRDS(fits,file.path(output.dir,"fits.RDS"))
  }
  gene.list <- gene.list[(gene.list%in%names(fits[[1]]))]
  print(fits[[2]][[1]][
  rancors <- lapply(1:length(fits[[1]]), function(x){
    ind <- intersect(colnames(fits[[1]][[x]][["coefficients"]]),colnames(fits[[2]][[x]][["coefficients"]]))
    a <- fits[[1]][[x]][["coefficients"]][1,ind]
    b <- fits[[2]][[x]][["coefficients"]][1,ind]
    cor(a,b  ,method="spearman")
  })
  xcors <- lapply(1:length(fits[[1]]), function(x){
      max(ccf(fits[[1]][[x]][["coefficients"]][1,],fits[[2]][[x]][["coefficients"]][1,] ,plot=F)$acf)
  })
  names(rancors) <- names(xcors) <- names(fits[[1]])
  
  n.rancors <- lapply(1:length(fits[[1]]), function(x){
    ind <- intersect(colnames(fits[[1]][[x]][["coefficients"]]),colnames(fits[[2]][[x]][["coefficients"]]))
    a <- fits[[3]][[x]][["coefficients"]][1,ind]
    b <- fits[[4]][[x]][["coefficients"]][1,ind]
    cor(a,b  ,method="spearman")
  })
  n.xcors <- lapply(1:length(fits[[1]]), function(x){
    max(ccf(fits[[3]][[x]][["coefficients"]][1,],fits[[4]][[x]][["coefficients"]][1,] ,plot=F)$acf)
  })

  format(Sys.time(),"%D %H:%M:%S")
  if(file.exists(file.path(output.dir,"dtw.RData"))){
    load(file.exists(file.path(output.dir,"dtw.RData")))
  }else{
  print("DTWing")
  dt.dist <- list()
  dt.likelihood <- list()
  dt.aligned.dist <- list()
  sumErrors <- list()
  dt.full <- list()
  pairs <-  c(2,1,4,3,6,5)
  for(s in 1:length(fits)){
    dt.dist[[s]] <<- vector()
    sumErrors[[s]] <<- vector()
    dt.likelihood[[s]] <<- vector()
    dt.aligned.dist[[s]] <<- matrix(nrow = length(fits[[s]]), ncol=ncol(fits[[s]][[1]][["coefficients"]]),dimnames = list(names(fits[[s]]),colnames(fits[[s]][["coefficients"]])))
    dtt <-  t(sapply(1:length(fits[[s]]),function(x){
      cat("\r",x)
      gn <- names(fits[[s]])[x]
      fd <- dist.fit(fits[[s]][[gn]],fits[[pairs[s]]][[gn]])
      exp1 <- fits[[s]][[gn]][["coefficients"]][1,]
      exp2 <- fits[[pairs[s]]][[gn]][["coefficients"]][1,]
      dtq1 <-  dtw(fd[[1]], open.end = T , open.begin = open.start, step.pattern = stepPattern,keep.internals = T)
      sE <- rowSums(sapply(1:length(dtq1$index1),function(i) c(fd[[2]][dtq1$index1[i],dtq1$index2[i]],sum(exp1[dtq1$index1[i]],exp2[dtq1$index2[i]]))))
      sumErrors[[s]][[gn]] <<- sE[1]/sE[2]
      #dtqz <- dtw(setsz[[s]][x,], setsz[[pairs[s]]][x,], open.end = T , open.begin = open.start, step.pattern = stepPattern,keep.internals = T)
      dt.dist[[s]][gn] <<- dtq1$normalizedDistance
      m <- dtq1$costMatrix
      pm <- exp(-m)
      pm <- sweep( pm,1,STATS = rowSums(pm,na.rm = T),FUN="/")
      dt.likelihood[[s]][[gn]] <<- -(1/length(dtq1$index1))*log(prod(sapply(1:length(dtq1$index1),function(i) pm[dtq1$index1[i],dtq1$index2[i]])))
      #dt.aligned.dist[[s]][gn,]<<-sapply(1:length(dtq1$index1),function(i){ setsz[[s]][gn,dtq1$index1[i]]-setsz[[pairs[s]]][gn,dtq1$index2[i]]})
      ans <- sapply(1:length(warp(dtq1,index.reference = T)),function(i){
        w <- warp(dtq1,index.reference = T)
        as.numeric(colnames(fits[[pairs[s]]][[1]][["coefficients"]])[w[i]])-as.numeric(colnames(fits[[s]][[1]][["coefficients"]])[i])
      }) 
      return(ans)
      }))
    rownames(dtt) <- names(fits[[s]])
    colnames(dtt) <- colnames(fits[[s]][[1]][["coefficients"]])
    format(Sys.time(),"%D %H:%M:%S")
    dt.full[[s]]<<-dtt
  }
  format(Sys.time(),"%D %H:%M:%S")
  save(list = list(dt,dt.dist,dt.likelihood,sumErrors),file = file.path(output.dir,"dtw.RData") )
  }

    #Could be more effieicnt
#   for(i in 1:length(fits)){
#     for(j in 1:length(fits[[i]])){
#       a= fits[[i]][[j]]$coefficients[1,]
#       b <- fits[[i]][[j]]$data
#       fits[[i]][[j]][["realSD"]] <-sapply(names(a),function(x){
#         sd.centered(as.numeric(b[b$tme==x,]$y),mu = a[x])
#       })
#       fits[[i]][[j]][["realSD.adj"]]<-sapply(names(a),function(x){
#         sd.centered(as.numeric(b[b$tme==x,]$y),mu = a[x])/a[x]
#       })
#       names(fits[[i]][[j]][["realSD.adj"]]) <- names(a)
#     }
#   }
#   dt.real.sds <-  lapply(1:nrow(dt[[1]]),function(i){
#     gn <- rownames(dt[[1]])[i]
#     a=fits[[1]][[gn]][["realSD.adj"]]
#     b=fits[[2]][[gn]][["realSD.adj"]]
#     mean(c(a,b))
#   })
#   names(dt.real.sds) = names(fits[[1]])
  

  dt.fc <- lapply(1:nrow(dt[[1]]),function(i){
    gn <- rownames(dt[[1]])[i]
    a=fits[[1]][[gn]]$coefficients[1,]
    b=fits[[2]][[gn]]$coefficients[1,]
    min((max(a)-min(a)), (max(b)-min(b)))
  })
  names(dt.fc) <- names(fits[[1]])
  
  metrics <-  list(sumErrors[[1]], dt.likelihood[[1]],dt.dist[[1]],rancors,xcors,dt.fc)
  nmet <-  c("sumErrors", "dt.likelihood","dt.dist","cor.s","cor.x","dt.fc")
  met1 <- metrics
  met1 <- Reduce(cbind, met1)
  colnames(met1) <- nmet
  
  metrics <-  list(sumErrors[[3]], dt.likelihood[[3]],dt.dist[[3]],n.rancors,n.xcors,dt.fc)
  nmet <-  c("sumErrors", "dt.likelihood","dt.dist","cor.s","cor.x","dt.fc")
  met3 <- metrics
  met3 <- Reduce(cbind, met3)
  colnames(met3) <- nmet
   
  #fff <-  randomForest(x = q.b.f[,c(-6,-8)],y = as.factor(unlist(q.b.f[,8])),mtry = 2,ntree = 2000,keep.forest = T)
  #fff$importance
  #saveRDS(fff,file = file.path(datapath,"randomforestnoz.RDS"))
  fff <- readRDS(file.path(datapath,"randomforestnoz.RDS"))
  p=predict(object=fff,met1)
  comparable <- as.logical(levels(p)[p])
  dt[[1]] <- dt.full[[1]][comparable,]
  dt[[2]] <- dt.full[[2]][comparable,]
  dt[[3]] <- dt.full[[3]][n.comparable,]
  dt[[4]] <- dt.full[[4]][n.comparable,]
  
  p=predict(object=fff,met3)
  n.comparable <- as.logical(levels(p)[p])

  write.pdf(  barplot(c("Comparable"=sum(comparable),"NonComparable"=sum(!comparable)),main = "Classifications By Random Forest",
                      xlab = "Class", ylab = "# of Genes"), file.path(output.dir,"shift_comparables.pdf"))
  
  write.pdf({
    mypar(3,4)
    for(i in 1:24){
      gn=names(om)[i]
      print(paste(gn, sapply(list(a,b,c,d,e), function(x) x[gn])))
      plot.sme( fits[[1]][[gn]],type = "model",showConfidenceBands=TRUE,showIndividual=F,main=paste(gn,"Mouse Log2 TPM"))
      plot.sme( fits[[2]][[gn]],type = "model",showConfidenceBands=TRUE,showIndividual=F,main=paste(gn,"Human Log2 TPM"))
    }
  }  , file.path(datapath,"neural_late_tpm.pdf"))
  
  
  avg.shift <-  colMedians(dt[[3]])
  rev.avg.shift <- colMedians(dt[[4]])
  dt.adj <- sweep(dt[[1]],2, avg.shift)
  nd.adj <- sweep(dt[[3]],2, avg.shift)
  nd <- dt[[3]][n.comparable,]
  nd.rev <- dt[[4]][n.comparable,]
  rev.dt.adj <- sweep(dt[[2]],2, rev.avg.shift)
  rev.nd.adj <- sweep(dt[[4]],2, rev.avg.shift)
  #######
  write.csv(set1,file = file.path(output.dir,"set1exprs.csv"),row.names = T)
  write.csv(set2,file = file.path(output.dir,"set2exprs.csv"), row.names = T)
  write.csv(s1.e,file = file.path(output.dir,"s1eexprs.csv"),row.names = T)
  write.csv(s2.c,file = file.path(output.dir,"s2cexprs.csv"), row.names = T)
  write.csv(dt[[1]],file = file.path(output.dir,"shift.csv"))
  write.csv(dt[[2]],file = file.path(output.dir,"rev_shift.csv"))
  write.csv(dt.adj,file = file.path(output.dir,"adjusted_shift.csv"))
  write.csv(dt[[3]],file = file.path(output.dir,"null_shift.csv"))
  
  print(paste("Mean Null shift",mean(rowMeans(nd))))
  print(paste("Mean shift experimental:", mean(rowMeans(dt[[1]]))))
  
  write.pdf(  barplot(colMeans(dt[[1]]),main = "Experiment Mean Shift at each Timepoint, All Genes",
                      xlab = "TP (day)", ylab = "Shift (# TPs)"), file.path(output.dir,"exp_Shift_tps.pdf"))
  
  write.pdf(   hist(dt[[1]],main = "Average Experiment Shift (across TPs for each gene), All Genes",
                    xlab = "shift (# TPs)", ylab = "Number of Observations"), file.path(output.dir,"exp_Shift_Hist.pdf"))  
  
  write.pdf(  barplot(colMeans(nd),main = "Null Mean Shift at each Timepoint, All Genes",
                      xlab = "TP (day)", ylab = "Shift (# TPs)"), file.path(output.dir,"Null_Shift_tps.pdf"))
  
  write.pdf(   hist(nd,main = "Average Null  Shift (across TPs for each gene), All Genes",
                    xlab = "shift (# TPs)", ylab = "Number of Observations"), file.path(output.dir,"Null_Shift_Hist.pdf"))
  
  #####Run T-Tests on overall data
  t <-  t.test(rowMeans(nd.adj),rowMeans(dt.adj))
  cou <- 0
  #Significance for each gene
  print("W Test Forward")
  gene.w.tests.faster <-  apply(dt[[1]], 1,function(i){
    cou<<-cou+1
    cat("\r",cou)
    p <-  wilcox.test(i,dt[[3]],alternative = "g")$p.value
    c("P Faster" = p)
  })
  print("W Test Backward")
  cou <- 0
  gene.w.tests.slower <- apply(dt[[2]], 1,function(i){
    cou<<-cou+1
    cat("\r",cou)
    p <-  wilcox.test(i,dt[[4]],alternative = "g")$p.value
    c("P Slower" = p)
  })
  
  bads <- gene.w.tests.faster < .1 & gene.w.tests.slower < .1
  gene.w.tests.slower[bads] <- gene.w.tests.faster[bads] <- .999
  pMetric <-  p.adjust(gene.w.tests.faster,method = "BH")/(p.adjust(gene.w.tests.faster,method = "BH")+p.adjust(gene.w.tests.slower,method = "BH"))
  
  w.test.genes <- data.frame("P Faster"=gene.w.tests.faster, "P Slower"= gene.w.tests.slower, 'P Faster BH Adj'= p.adjust(gene.w.tests.faster,method = "BH"),"P Slower BH Adj"= p.adjust(gene.w.tests.slower, method = "BH"),"PFaster.over.PFaster.plus.PSlower"=pMetric )
  write.pdf(barplot( colSums(w.test.genes[,3:4]<.05) ,las=1,xpd = F ,main = "# Of of Experimental Genes Shifted Fwd vs Back \nMann-Whitney-Wilcoxon Test (5% FDR)",
                     ylab = "# Of Genes"), file.path(output.dir,"exp_fwd_vs_back.pdf"))
  write.csv(w.test.genes,file = file.path(output.dir,"gene_shift_wtest.csv"))
  fratio <- -log10(w.test.genes[,3]/(w.test.genes[,3]+w.test.genes[,4]))
  fratio[which(fratio<=0.0)] <- 0.0
  fwdSigs <-  fratio
  names(fwdSigs) <- rownames(w.test.genes)
  bratio <- -log10(w.test.genes[,4]/(w.test.genes[,3]+w.test.genes[,4]))
  bratio[bratio<0] <- 0.0
  bckSigs <-  bratio
  names(bckSigs) <- rownames(w.test.genes)
  
  write.pdf(  plot(w.test.genes[,3],w.test.genes[,4],main="Pvalues of being foward vs back") , file.path(output.dir,"PvalPlot.pdf"))
  
  t.gnz <- union(names(fits[[1]]), names(fits[[2]]))
  infwd <-  (t.gnz %in% names(fwdSigs))
  inbwd <-  (t.gnz %in% names(bckSigs))
  for(i in 1:length(fits[[1]])){
    if(!infwd[i]){
      fwdSigs[t.gnz[i]] = 0.0
    }
    if(!inbwd[i]){
      bckSigs[t.gnz[i]] = 0.0
    }
  }
  fwdGO <- eacitest(fwdSigs,"org.Hs.eg","SYMBOL",sets = "GO",minsetsize = 30)$setscores
  bckGO <- eacitest(bckSigs,"org.Hs.eg","SYMBOL",sets = "GO", minsetsize=30)$setscores
  write.table(fwdGO , file =file.path(output.dir,"ForwardEnrichment.tsv") ,sep = "\t",row.names = F)
  write.table(bckGO , file =file.path(output.dir,"BackwardEnrichment.tsv") ,sep = "\t",row.names = F)
  
  pos.shift.pvals <-  sapply(1:ncol(dt[[1]]) , function(i){
    print(paste("Pval tp",i))
    sapply(1:nrow(dt[[1]]),function(j){
      mean( dt[[1]][j,i] <= dt[[3]][,i])
    })
  } )
  
  neg.shift.pvals <-  sapply(1:ncol(dt[[2]]) , function(i){
    print(paste("Pval tp",i))
    sapply(1:nrow(dt[[2]]),function(j){
      mean( dt[[2]][j,i] <= dt[[4]][,i])
    })
  } )
  
  dimnames(pos.shift.pvals) <- dimnames(dt[[1]])
  dimnames(neg.shift.pvals) <- dimnames(dt[[2]])
  matchpoints <- intersect(colnames(pos.shift.pvals),colnames(neg.shift.pvals))
  if(length(tp1.i)==length(tp2.i)){
    bads <- pos.shift.pvals[,matchpoints] < .1 & neg.shift.pvals[,matchpoints] < .1
    pos.shift.pvals[bads] <- neg.shift.pvals[bads] <- .999
    x.over.y <-  pos.shift.pvals/(pos.shift.pvals+neg.shift.pvals)
    write.pdf( heatmap.2( x.over.y[gene.list,] ,col = cols,
                          trace = "none",Rowv = F,Colv = F,cexRow = .35,cexCol = .35,
                          main = "Pfaster/(Pfaster+Pslower)"), file.path(output.dir,"PMetric.pdf"))
    
  }
  
  
  write.csv(pos.shift.pvals,file = file.path(output.dir,"forward_shift_pvals.csv"))
  write.csv(neg.shift.pvals,file = file.path(output.dir,"back_shift_pvals.csv"))
  
  write.pdf( heatmap.2( dt.adj[gene.list,] ,col = cols,
                        trace = "none",Rowv = F,Colv = F,cexRow = .45,cexCol = .45,
                        main = "Shift (# of TPs)"), file.path(output.dir,"Adjusted_Shifts.pdf"))
  
  write.pdf( heatmap.2( pos.shift.pvals[gene.list,],col = rev(cols),
                        trace = "none",Rowv = F,Colv = F,cexRow = .45,cexCol = .45,main = "Forward Shift Pvals"), file.path(output.dir,"Selected_Forward_Pvals.pdf"))
  
  write.pdf( heatmap.2( neg.shift.pvals[gene.list,],col = rev(cols),
                        trace = "none",Rowv = F,Colv = F,cexRow = .45,cexCol = .45,main = "Backward Shift Pvals"), file.path(output.dir,"Selected_Backward_Pvals.pdf"))
  
  write.pdf( heatmap.2( cbind(set2[gene.list,], set1[gene.list,] ),col = cols,
                        trace = "none",Rowv = F,Colv = F,cexRow = .45,cexCol = .45,main = "Expression of Selected Genes (log2 TPM)"), file.path(output.dir,"Selected_Expression.pdf"))
  
  write.pdf(  barplot(colMeans(pos.shift.pvals < .05),main = "% of Experimental Genes Shifted Ahead \n at each TP\n(5% FDR)",
                      xlab = "TP (day)", ylab = "% of Genes"), file.path(output.dir,"exp_forward_shift_significant.pdf"))
  
  write.pdf(  barplot(colMeans(neg.shift.pvals < .05 ),main = "% of Experimental Genes Shifted Back \n at each TP\n(5% FDR)",
                      xlab = "TP (day)", ylab = "% of Genes"), file.path(output.dir,"exp_back_shift_significant.pdf"))
  
  t.tst <-  data.frame(t$estimate[1], t$estimate[2], t$p.value, row.names = "Value")
  colnames(t.tst) <- c("Control Avg Shift", "Experimental Avg Shift", "P.value")
  write.csv(t.tst,file = file.path(output.dir,"overall_shift_T_Test.csv"))
  if(ppt){
    mydoc = pptx()
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle(mydoc, "DTW Analysis")
    mydoc = addParagraph( mydoc, value = paste(output.dir,"FC Threshold",fc.threshold , "Step Pattern", paste(stepPattern,collapse = ",") ))
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle( mydoc, "Overall T Test")
    mydoc = addParagraph(mydoc, paste("Mean Null:",t$estimate[1], "Mean Experimental:",t$estimate[2], "P Value",  t$p.value))
    mydoc = addSlide( mydoc, "Two Content" )
    mydoc = addTitle( mydoc, "Shifts Observed")
    mydoc = addPlot(mydoc, function() hist(dt[[1]],main = "Average Experiment Shift (across TPs for each gene), All Genes",
                                           xlab = "shift (# TPs)", ylab = "Number of Observations"), vector.graphic = F)
    mydoc = addPlot(mydoc, function()  hist(nd,main = "Average Null  Shift (across TPs for each gene), All Genes",
                                            xlab = "shift (# TPs)", ylab = "Number of Observations"), vector.graphic = F)
    
    mydoc = addSlide( mydoc, "Two Content" )
    mydoc = addTitle( mydoc, "Shifts at each TP")
    mydoc = addPlot(mydoc, function() barplot(colMeans(dt[[1]]),main = "Experiment Mean Shift at each Timepoint, All Genes",
                                              xlab = "TP (day)", ylab = "Shift (# TPs)"), vector.graphic = F)
    mydoc = addPlot(mydoc, function() barplot(colMeans(nd),main = "Null Mean Shift at each Timepoint, All Genes",
                                              xlab = "TP (day)", ylab = "Shift (# TPs)"), vector.graphic = F)
    
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle( mydoc, "Mann-Whitney Test")
    mydoc = addPlot(mydoc, function() barplot(colSums(w.test.genes[,3:4]<.05) ,las=1,xpd = F ,main = "# Of of Experimental Genes Shifted Fwd vs Back \nMann-Whitney Test (5% FDR)",
                                              ylab = "# Of Genes"), vector.graphic = F)
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle(mydoc, "GO Enrichment Fwd Shifted")
    mydoc = addFlexTable(mydoc,flextable = FlexTable(fwdGO[1:50,c("Term","set.size","pval")]))
    mydoc = addSlide( mydoc, "Title and Content" )
    mydoc = addTitle(mydoc, "GO Enrichment Back Shifted")
    mydoc = addFlexTable(mydoc,flextable = FlexTable(bckGO[1:50,c("Term","set.size","pval")]))
    mydoc = addSlide( mydoc, "Two Content" )
    mydoc = addTitle( mydoc, "Shifts & Expression")
    mydoc = addPlot(mydoc, function()  heatmap.2( dt.adj[gene.list,] ,col = cols,
                                                  trace = "none",Rowv = F,Colv = F,cexRow = .45,cexCol = .45,
                                                  main = "Shift (# of TPs)"), vector.graphic = F)
    mydoc = addPlot(mydoc, function()  heatmap.2( cbind(set2[gene.list,], set1[gene.list,] ),col = cols,
                                                  trace = "none",Rowv = F,Colv = F,cexRow = .45,cexCol = .45,main = "Expression of Selected Genes (log2 TPM)"), vector.graphic = F)
    
    if(length(tp1.i)==length(tp2.i)){
      mydoc = addSlide( mydoc, "Title and Content" )
      mydoc = addTitle( mydoc, "Pfaster/(Pfaster+Pslower)")
      mydoc = addPlot(mydoc, function()  heatmap.2( x.over.y[gene.list,] ,col = cols,
                                                    trace = "none",Rowv = F,Colv = F,cexRow = .45,cexCol = .45,
                                                    main = "Pfaster/(Pfaster+Pslower)"), vector.graphic = F)
    }
    
    mydoc = addSlide( mydoc, "Two Content" )
    mydoc = addTitle( mydoc, "P Values")
    mydoc = addPlot(mydoc, function()  heatmap.2( pos.shift.pvals[gene.list,],col = rev(cols),
                                                  trace = "none",Rowv = F,Colv = F,cexRow = .45,cexCol = .45,main = "Forward Shift Pvals"), vector.graphic = F)
    mydoc = addPlot(mydoc, function()  heatmap.2( neg.shift.pvals[gene.list,],col = rev(cols),
                                                  trace = "none",Rowv = F,Colv = F,cexRow = .45,cexCol = .45,main = "Backward Shift Pvals"), vector.graphic = F)
    writeDoc( mydoc, file.path(output.dir,"summary.pptx") )
  }
}

si=1:50
test.setM = list(exp_1_M[si,], exp_2_M[si,], datasets$`~/code/data/cb/mouse_ALONE.csv`)
test.setH= list( exp_1[si,1:9],exp_2_H[si,], datasets$`~/code/data/cb/human_ALONE.csv`)
setM = list(exp_1_M, exp_2_M, datasets$`~/code/data/cb/mouse_ALONE.csv`)
setH = list( exp_1[,1:9],exp_2_H, datasets$`~/code/data/cb/human_ALONE.csv`)
setMix = list(exp_1[,17:24],exp_2)
tpH =list(tp_1_M,tp_2,tp.m.control)
tpM =list(tp_1[1:9],tp_2,tp.h.control)
tpMix = list(tp_1[17:24],tp_2)

run.DTW.Genes(set1 = test.setM ,set2 = test.setH, output.dir = "~/code/data/cb/shiftFiles/MousevsHInVitro",tp1 =tpM, tp2=tpH,gene.list = neural_list, open.start = F,ppt = T ,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))#typesOfsteps[[10]])

run.DTW.Genes(set1 = setM ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MousevsHInVitro",tp1 =tpM, tp2=tpH,gene.list = neural_list, open.start = F,ppt = T ,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))#typesOfsteps[[10]])
run.DTW.Genes(set1 = setMix ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MixvsHInVitro",tp1 =tpMix, tp2=tpH,gene.list = neural_list, open.start = F,ppt = T ,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))#typesOfsteps[[10]])
run.DTW.Genes(set1 = setM ,set2 = setH, output.dir = "~/code/data/cb/shiftFiles/MousevsMouseInVitro",tp1 =tpM, tp2=tpMix,gene.list = neural_list, open.start = F,ppt = T,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))#typesOfsteps[[10]])

  

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
    plot.sme( fits[[1]][[gn]],type = "model",showConfidenceBands=TRUE,showIndividual=F,main=paste(gn,"Mouse Log2 TPM"))
    plot.sme( fits[[2]][[gn]],type = "model",showConfidenceBands=TRUE,showIndividual=F,main=paste(gn,"Human Log2 TPM"))
  }
}  , file.path(datapath,"neural_late_tpm.pdf"))


  
  gn="LGALS3"
  s=1
  plot.sme( fits[[1]][[gn]],type = "model",showConfidenceBands=TRUE,showIndividual=F,main=paste(gn,"Mouse Log2 TPM"))
  plot.sme( fits[[2]][[gn]],type = "model",showConfidenceBands=TRUE,showIndividual=F,main=paste(gn,"Human Log2 TPM"))
  fd <- dist.fit(fits[[s]][[gn]],fits[[pairs[s]]][[gn]])
  rfd <- dist.fit(fits[[pairs[s]]][[gn]],fits[[s]][[gn]])
  exp1 <- fits[[s]][[gn]][["coefficients"]][1,]
  exp2 <- fits[[pairs[s]]][[gn]][["coefficients"]][1,]
  dtq1 <-  dtw(fd[[1]], open.end = T , open.begin = open.start, step.pattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F),keep.internals = T)
  rdtq1 <-  dtw(rfd[[1]], open.end = T , open.begin = open.start, step.pattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F),keep.internals = T)
  sumErrors <-mean(sapply(1:length(dtq1$index1),function(i) fd[[2]][dtq1$index1[i],dtq1$index2[i]])/mean(exp1[dtq1$index1[i]],exp2[dtq1$index2[i]]))
  dtwPlotThreeWay(rdtq1,exp2,exp1)
  warp(rdtq1)
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
  


  
  
  
  gn <- gene.list[i]
  plot.sme( fits[[1]][[gn]],type = "model",showConfidenceBands=TRUE,showIndividual=F,main=paste(gn,"Mouse Log2 TPM"))
  plot.sme( fits[[2]][[gn]],type = "model",showConfidenceBands=TRUE,showIndividual=F,main=paste(gn,"Human Log2 TPM"))
  dt[[1]][gn,]
  met1[gn,]
  q.b.f <- rbind( q.b.f, c(met1[gn,], T))
  rownames(q.b.f)[nrow(q.b.f)] <- gn
  i=i+1
  #q.b.f <- cbind(t(met1[gn,]),T)
  #rownames(q.b.f)[106:112] <- gene.list[(i-6):i]
  saveRDS(q.b.f,file = file.path(datapath,"qbf.RDS"))
  
  #q.b.f[,"dt.real.sds"] <-dt.real.sds[rownames(q.b.f)]
  #q.b.f <-  q.b.f[-nrow(q.b.f),]
  

  best1 <- list(0)
  for(xx in 1:15){
    print(xx)
    met.data <- q.b.f[,c(-6)]
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
    for(i in 1:6){
      #Segement your data by fold using the which() function 
      testIndexes <- which(folds==i,arr.ind=TRUE)
      testData <- met.data[testIndexes, ]
      trainData <- met.data[-testIndexes, ]
      acc <- array()
      www <- ncol(testData)
      fac <- as.factor(unlist(trainData[,www]))
      fff <-  randomForest(x = trainData[,-www],y = fac,mtry = 2,ntree = 1000)
      print(fff$importance)
      ff <-  C5.0(x=trainData[,-www],fac,trials =10, costs=matrix( c(0,1,1,0),2,2,dimnames =list(list(0,1),list(0,1))),rules = T)
      acc[1]=mean(predict(object=fff, testData[,-www]) == testData[,www])
      acc[2]=mean(levels(predict.C5.0(ff, testData[,-www]))[predict.C5.0(ff, testData[,-www])]==testData[,www])
      nb <- naiveBayes(x = trainData[,-www],y = fac)
      acc[3] <- mean(levels(predict(nb,testData[,-www] ))[predict(nb,testData[,-www])]==testData[,www])
      fit4 <- PART( fac~. ,data = as.data.frame(trainData[,-www]) )
      acc[4]<- mean((predict(fit4,  as.data.frame(testData[,-www]))) == testData[,www])
      for(q in 1:length(acc)){
        accs[[q]] <- c(accs[[q]],acc[q])
      }
      best1[[xx]] <- lapply(lapply(accs,unlist), mean)
    }
  }
  mean(unlist(lapply(best1,"[[",1)))
  mean(unlist(lapply(best1,"[[",2)))
  mean(unlist(lapply(best1,"[[",3)))
  mean(unlist(lapply(best1,"[[",4)))

  
    



  #Takes two model fits and calculates their bhattacharya distance, as well as a sum of variance
  bhattacharya.fit <- function(f1,f2){
    cf1 = coef(f1)
    cf2= coef(f2)
    cv1=diag(vcov(f1))
    cv2=diag(vcov(f2))
    #anova(f1,f2)
    x=ncol(cf1)
    y=ncol(cf2)
    d <- matrix(nrow = x, ncol = y)
    v <- matrix(nrow = x, ncol = y)
    for(i in 1:x){
      for(j in 1:y){
        d[i,j] <- bhattacharyya.dist(cf1[1,i], cf2[1,j], cv1[i], cv2[j])
        print(paste(cf1[1,i], cf2[1,j], cv1[i], cv2[j]))
        print( d[i,j])
        v[i,j] <- cv1[i]+cv2[j]
      }
    }
    return(list(d,v))
  }










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

#k.cor <- (1-cor( exp_1_M[ind,],datasets$`~/code/data/cb/mouse_ALONE.csv`[ind,] ,method = "kendall"))
s.cor <-  (1-cor( exp_1_M[ind,],datasets$`~/code/data/cb/mouse_ALONE.csv`[ind,] ,method = "spearman"))
p.cor <-  (1-cor( exp_1_M[ind,],datasets$`~/code/data/cb/mouse_ALONE.csv`[ind,] ,method = "pearson"))

heatmap.2( k.cor,Rowv = F, Colv = F,trace = "none",cexRow = .5)

k.cor
s.cor

ind <- intersect(rownames(exp_1_M), rownames(datasets$`~/code/data/cb/mouse_ALONE.csv`))
dtt <-  dtw(k.cor, open.end = T,open.begin = F,step.pattern = rabinerJuangStepPattern(type=3,slope.weighting="c",smoothed=T),keep.internals = T)
dtwPlot(dtt)
dtwPlotDensity(dtt)
tp.h.control[as.integer(warp(dtt, index.reference = T))]



aad <- function (x, na.rm = FALSE) 
{
  if (!is(x, "numeric") & !is(x, "integer")) {
    stop("\"x\" must be numeric")
  }
  if (!is(na.rm, "logical") | length(na.rm) != 1) {
    stop("\"na.rm\" must be a single logical value")
  }
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  y <- mean(abs(x - median(x)))
  return(y)
}
rowAads <- function(x){
  apply(x,1,aad)
}

write.pdf <- function(plot, filename){
  pdf(filename,width = 20, height = 20)
  plot
  dev.off()
}


z.Score <- function(x){
  (x-mean(x))/sd(x)
}
z.Score.aad <- function(x){
  (x-median(x))/aad(x)
}

neural_list <- sort(as.character(read.csv("~/code/data/cb/markers/marker_union.csv",header = T)$x))
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
typesOfsteps <- list(  typeIVc, typeIIIc, typeIId, typeIIc,  typeIc,typeIcs,typeId,typeIds,symmetric2,asymmetric,asymmetricP0,asymmetricP05,asymmetricP1, mori2006)

datapath <- "~/code/data/cb"
setwd(datapath)
datafiles <-  file.path(datapath, dir()[ grepl(".csv", dir())])
datasets <- sapply( datafiles, read.csv, header=T, sep=",")
for(i in 1:length(datasets)){
  rownames(datasets[[i]]) <-  make.names(toupper(datasets[[i]][,1]),unique = T)
  colnames(datasets[[i]]) <-  gsub("X","",colnames(datasets[[i]]))
  datasets[[i]] <- log(as.matrix(datasets[[i]][sort(rownames(datasets[[i]])),-1])+1,2)
  #class(datasets[[i]]) <- "numeric"
}

datapath <- "~/code/data/cb"
setwd(datapath)
mitodatafiles <-  file.path(datapath, dir()[ grepl(".tsv", dir())])
mitodatasets <- sapply( mitodatafiles, read.csv, header=T, sep="\t",row.names=1)
head(mitodatasets[[3]])
for(i in 1:length(mitodatasets)){
  rownames(mitodatasets[[i]]) <-  make.names(toupper(rownames(mitodatasets[[i]])),unique = T)
  mitodatasets[[i]] <- log(as.matrix(mitodatasets[[i]][sort(rownames(mitodatasets[[i]])),])+1,2)
  #class(mitodatasets[[i]]) <- "numeric"
}
matplot(t(mitodatasets$`~/code/data/cb/human_TERA_mito.tsv`),type="l")
matplot(t(mitodatasets$`~/code/data/cb/human_TERA_mito_renorm.tsv`[rownames(mitodatasets$`~/code/data/cb/mouse_TERA_mito.tsv`),]),type="l")
matplot(t(mitodatasets$`~/code/data/cb/mouse_TERA_mito.tsv`),type="l")
rownames(mitodatasets$`~/code/data/cb/human_TERA_mito.tsv`)
rownames(mitodatasets$`~/code/data/cb/human_TERA_mito_renorm.tsv`)[which(grepl(pattern = "MT",x = rownames(mitodatasets$`~/code/data/cb/human_TERA_mito_renorm.tsv`)))]

dfnames <- gsub(".csv", "", dir()[ grepl(".csv", dir())])
for(i in 1:length(datasets)){
  colnames(datasets[[i]]) <-  paste(dfnames[i] ,colnames(datasets[[i]]))
}
sapply(datasets,colnames)


samples_1 <- colnames(datasets$`~/code/data/cb/TC_1.csv`)
tp_1 <- as.numeric(sapply(samples_1, function(x) strsplit(x,"d")[[1]][2]) )

samples_1_M <- colnames(datasets$`~/code/data/cb/TC_1_M.csv`)
tp_1_M <- as.numeric(gsub('.1',"",sapply(samples_1_M, function(x) strsplit(x,"d")[[1]][2])))
exp_1 <- datasets$`~/code/data/cb/TC_1.csv`[,sort(tp_1,index.return=T)$ix]
type_1<- sapply(colnames(exp_1), function(x) strsplit(x,"d")[[1]][1]) 
exp_1 <- exp_1[,sort(type_1,index.return=T)$ix]
exp_1_M <- datasets$`~/code/data/cb/TC_1_M.csv`[,sort(tp_1_M,index.return=T)$ix]
exp_1_M <- datasets$`~/code/data/cb/TC_1_M.csv`[intersect(rownames(exp_1_M),rownames(exp_1)),]
type_1_M <- gsub('_',"",sapply(colnames(exp_1_M), function(x) strsplit(x,"d")[[1]][1]))
all_exp_1 <- cbind(exp_1[intersect(rownames(exp_1_M),rownames(exp_1)),], exp_1_M[intersect(rownames(exp_1_M),rownames(exp_1)),] )
nl <- neural_list[neural_list %in% rownames(exp_1)]
nlm <- nl[nl%in%rownames(exp_1_M)]

exp_2 <- datasets$`~/code/data/cb/TC_2_Mix.csv`
exp_2_M <- datasets$`~/code/data/cb/TC_2_M.csv`
exp_2_H <- datasets$`~/code/data/cb/TC_2_H.csv`
exp_2 <- exp_2[ rowMeans(exp_2)>1,]
exp_2_H<- exp_2_H[ rowMeans(exp_2_H)>1,]
exp_2_M<- exp_2_M[ rowMeans(exp_2_M)>1,]
exp_2 <- exp_2[ intersect(rownames(exp_2), rownames(exp_2_H)),]
exp_2_H<- exp_2_H[ intersect(rownames(exp_2_H), rownames(exp_2)),]
exp_2_M<- exp_2_M[ intersect(rownames(exp_2_H), rownames(exp_2_M)),]
nl2 <- neural_list[neural_list %in% rownames(exp_2)]
tp_2 <- colnames(exp_2)
colnames(exp_2) <- paste("Mix", colnames(exp_2))
colnames(exp_2_H) <- paste("Human", colnames(exp_2_H))
colnames(exp_2_M) <- paste("Mouse", colnames(exp_2_M))
tp_2 <- as.numeric(gsub("Mix TC_2_Mix ","",colnames(exp_2)))





run.DTW.Genes(set1 = exp_2 ,set2 = exp_2_H, output.dir = "~/code/data/cb/shiftFiles/MixvsHInVitro",tp1 =tp_2, tp2=tp_2, tp1.i =seq(0,max(tp_2)), tp2.i=seq(0,max(tp_2)),gene.list = neural_list, open.start = F,ppt = T,null.length = .5 ,aad.Threshold.denom = 3,cor.Threshold = .6,stepPattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))#typesOfsteps[[10]])
run.DTW.Genes(set1 = exp_2_M ,set2 = exp_2_H, output.dir = "~/code/data/cb/shiftFiles/MousevsHHInVitro",tp1 =tp_2, tp2=tp_2, tp1.i =seq(0,max(tp_2)), tp2.i=seq(0,max(tp_2)),gene.list = neural_list, open.start = F,ppt = T,null.length = .5 ,aad.Threshold.denom = 3,cor.Threshold = .5,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))#typesOfsteps[[10]])
run.DTW.Genes(set1 = exp_2_M ,set2 = exp_2, output.dir = "~/code/data/cb/shiftFiles/MousevsMixInVitro",tp1 =tp_2, tp2=tp_2, tp1.i =seq(0,max(tp_2)), tp2.i=seq(0,max(tp_2)),gene.list = neural_list, open.start = F,ppt = T,null.length = .5 ,aad.Threshold.denom = 3,cor.Threshold = .6,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))#typesOfsteps[[10]])

run.DTW.Genes(set2 = exp_2_H ,set1 = datasets$`~/code/data/cb/human_ALONE.csv`, output.dir = "~/code/data/cb/shiftFiles/LongHvsHInVitro",tp2 =tp_2, tp1=tp.h.control, tp2.i =seq(0,max(tp_2)), tp1.i=seq(0,max(tp.h.control)),gene.list = neural_list, open.start = F,ppt = T,null.length = .5 ,aad.Threshold.denom = 3,cor.Threshold = .6,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))#typesOfsteps[[10]])
run.DTW.Genes(set1 = exp_2 ,set2 = datasets$`~/code/data/cb/human_ALONE.csv`, output.dir = "~/code/data/cb/shiftFiles/LongHvsMixInVitro",tp1 =tp_2, tp2=tp.h.control, tp1.i =seq(0,max(tp_2)), tp2.i=seq(0,max(tp.h.control)),gene.list = neural_list, open.start = F,ppt = T,null.length = .5 ,aad.Threshold.denom = 3,cor.Threshold = .6,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))#typesOfsteps[[10]])
run.DTW.Genes(set2 = exp_2_M ,set1 = datasets$`~/code/data/cb/human_ALONE.csv`, output.dir = "~/code/data/cb/shiftFiles/LongHvsMouseInVitro",tp2 =tp_2, tp1=tp.h.control, tp2.i =seq(0,max(tp_2)), tp1.i=seq(0,max(tp.h.control)),gene.list = neural_list, open.start = F,ppt = T,null.length = .5 ,aad.Threshold.denom = 3,cor.Threshold = .6,stepPattern =rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))#typesOfsteps[[10]])

y <- "NEUROD4"
y = "SOX2"
y= 12
y=i
a <- fits[[1]][[y]]$coefficients[1,]
b <- fits[[2]][[y]]$coefficients[1,]
b = datasets$`~/code/data/cb/human_ALONE.csv`[y,]
plot(a,type="l")
plot(b, type ="l")
cor(a,b ,method = "spe")
ccf(a,b,plot = F)
gene.list

#plot(rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))
dtt <-  dtw(b, a, open.end = T , open.begin = F ,dist.method = "euclidean",step.pattern = rabinerJuangStepPattern(type=1,slope.weighting="c",smoothed=T),keep.internals = T)
dtwPlotThreeWay(dtt,xts = b, a)
dtt$normalizedDistance
dtt$index2
warp(dtt)
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
dtwDist(s1.e[1:200,], s2.c[1:20,], diag=F,open.end = T , open.begin = T, step.pattern = rabinerJuangStepPattern(type=7,slope.weighting="c",smoothed=F))
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



# sapply(1:length(colnames(exp_1)), function(x) colnames(exp_1)[x] <<- strsplit( colnames(exp_1)[x], "d")[[1]][2] )
# splitSets <- sapply(unique(type_1), function(x) exp_1[,x==type_1] )
# 
# GO.by.Max.TP <-  lapply(splitSets, function(i){
#   lapply(1:length(colnames(i)), FUN =  function(t){
#     maxTP <- apply(i,1, function(x) colnames(i)[which.max(x)] ) == colnames(i)[t]
#       gnV <- as.numeric(maxTP)
#       names(gnV) <- names(maxTP)
#       #print(gnV)
#       ecT <-  eacitest(gnV,"org.Hs.eg","SYMBOL")$setscores
#   })
# })
# 
# q <-  GO.by.Max.TP
# 
# GO.by.Max.TP.EASE <-  lapply(splitSets, function(i){
#   lapply(1:length(colnames(i)), FUN =  function(t){
#     maxTP <- apply(i,1, function(x) colnames(i)[which.max(x)] ) == colnames(i)[t]
#     gnV <- as.numeric(maxTP)
#     names(gnV) <- names(maxTP)
#     #print(gnV)
#     ecT <-  easetest(gnV,"org.Hs.eg","SYMBOL")$setscores
#   })
# })
# 
# i = splitSets$MgreaterthanH_H
# t = 1
# maxTP <- apply(i,1, function(x) colnames(i)[which.max(x)] ) == colnames(i)[t]
# gnV <- as.numeric(maxTP)
# names(gnV) <- names(maxTP)
# #print(gnV)
# ecT <-  eacitest(gnV,"org.Hs.eg","SYMBOL",minsetsize = 1)$setscores
# ecT$Term[1:100]
# 
# sapply(1:length(GO.by.Max.TP), function(x){
#   names(GO.by.Max.TP[[x]]) <<- colnames(splitSets[[x]])
# })
# 
# lapply(GO.by.Max.TP, function(x){
#   sapply(1:length(x), function(i){
#     #sum(grepl("neu", x[[i]]$setscores$Term[x[[i]]$setscores$pval<.001]))
#     x[[i]]$setscores$Term[x[[i]]$setscores$pval<.00001]
#   })
# })
# 

# splitSets_Z <- sapply(splitSets, function(s) (s-rowMedians(s))/rowSds(s))
# 
# length(splitSets_Z)
# 
# sapply(splitSets_Z, function(x){
#   sapply(1:ncol(x), function(i){
#     sum(rownames(x)[x[,i]>0] %in% neural_list )
#   })
# })

