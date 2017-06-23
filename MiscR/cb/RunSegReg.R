#.libPaths(.libPaths()[2])
source("http://bioconductor.org/biocLite.R")
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(EACI)
library(matrixStats)
library(dtw)
library(EBSeq)
# library(genefilter)
options("mc.cores"=64)
library(parallel)
library(rafalib)
library(calibrate)
source("~/code/cb/LoadFunctions.R")
source("~/code/cb/LoadData.R")
#datapath ="~/code/data/cb"
#remove.packages("SegReg")
#install.packages("~/code/general/SegRegPar/package/SegRegPar",repos = NULL,type = "source")
#install.packages("~/Downloads/reactome.db_1.59.1.tar.gz",repos = NULL,type = "source")

library(SegRegPar)
library(MASS)
library(ggrepel)
#setwd("/home/mt/R/x86_64-pc-linux-gnu-library/3.2")

SegRegEvents <-  function(set1,tp1, output.dir,fn="events" , tpm.threshold =2 , gene.list=neural_list ,pvalcut = .1,maxk = 10,min.num.in.seg = 3,cutdiff = .01,logT=F){
    cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
    all.rn <- as.character(Reduce(intersect, lapply(set1,rownames)))
    set1 <- lapply(set1, function(n) n[all.rn,])
    s1.mx <-  lapply(set1, function(x) rowMaxs(x) )
    pass.thresh <-  Reduce( "+" ,lapply(s1.mx, ">",tpm.threshold)) >= length(set1)
    set1 <- lapply(set1, function(n) n[pass.thresh,])
    mat1 <- Reduce(cbind, set1)
    if(logT) mat1 <- log(mat1+1,2)
    tp1cat <- Reduce(c , tp1)

    srO <- segreg(mat1,meancut = tpm.threshold,t.vect = tp1cat ,pvalcut = pvalcut,maxk =maxk,min.num.in.seg =  min.num.in.seg,cutdiff =  cutdiff,num.try = 100,keepfit = F)
    
    detectEvents <- function(o,tp){
      mclapply(o, function(x){
        #x <- o[[n]]
        EventPoints <-  c(tp[1],x$bp, tp[length(x$fitted)])
        EventPoints <- na.omit(EventPoints)
        Event.Table <- sapply(1:(length(EventPoints)-1), function(i){
          c(start= EventPoints[i],end= EventPoints[i+1] , slp= x$slp[i],slp.pval= x$slp.pval[i],fc = x$slp[i]*(EventPoints[i+1] -EventPoints[i]))
        } )
        t(Event.Table)
      })
    }
    tpU <- unique(tp1cat) 
    tpU <- sort(tpU)
    deEvents <- detectEvents(srO, tpU)
    inds <- sapply(deEvents,length)>3
    deEvents <- deEvents[inds]
    gene.names <- unlist(lapply(names(deEvents), function(n){ rep(n, nrow(deEvents[[n]]))} ))
    de <- Reduce(rbind, deEvents)
    de <- data.frame(de,stringsAsFactors = F)
    de <- cbind(gene.names, de)
    print(srO)
    print(de)
    colnames(de) <- c("gene","start","end", "slp", "slp.pval","fc")
    dir.create(output.dir)
    save("deEvents",file = file.path(output.dir,paste0(fn,".RData")))
    saveRDS(srO,file = file.path(output.dir,paste0(fn,".rds")))
    write.table(de,file = file.path(output.dir,paste0(fn,".txt")),sep = "\t")
    return( de )
}

datapath <- "~/code/data/cb/NBNogData/segregevents"
setnums <- c(29,27,28,11,14,26,25,10,12,13,15,16,3,1,7,8)
experiments <- datasets[setnums]
tps <- tpDatasets[setnums]
srE <-  lapply(names(experiments),function(n){
  print(n)
  if(!file.exists(file.path(datapath,paste0(n,"Events.txt")))){
    print(file.path(datapath,paste0(n,"Events.txt")))
    SegRegEvents(set1 = list(experiments[[n]]), tp1 = list(tps[[n]]*24),min.num.in.seg = 3,maxk = min(round(length(unique(tps[[n]]))/5,digits = 0),10),output.dir = datapath, fn=paste0(n,"Events"),logT = T) 
  }
})

a <- read.table("~/code/data/cb/NBNogData/segregevents/MouseMegaEndodermEvents.txt")
a$gene <- toupper(a$gene)

library(plyr)
setwd("~/code/data/cb/NBNogData/segregevents/")
dir()
setnums <- c(1,3,7,8,11,14,20,25,26,27,28,29)
experiments <- datasets[setnums]
tps <- tpDatasets[setnums]
outputs <-  lapply(names(experiments),function(n){
  a <- read.table(paste0("~/code/data/cb/NBNogData/segregevents/", n,"Events.txt"))
  a.sig <- a[a$slp.pval<.1& abs(a$fc)>1,]
  a.sig$gene <- toupper(a.sig$gene)
  #a.tfs <- a.sig
  a.tfs <- a.sig[a.sig$gene%in%union(MouseTFs,HumanTFs),]
  #a.tfs <- a.sig[a.sig$gene%in%signalingMolecules,]
  dat <-  data.frame(gene=a.tfs[,1],s=a.tfs[,"start"],e=a.tfs[,"end"],y=a.tfs[,"fc"],slp=a.tfs[,"slp"],UpDown=as.factor(sign(a.tfs[,"fc"])))
  cdat <- ddply(dat, "UpDown", summarise, duration.mean=mean(e-s))
  sdat <- ddply(dat, "UpDown", summarise, slp.mean=mean(slp))
  alpha=2/log(nrow(dat)+25,5)
  pdf(paste0('~/Desktop/',n,"TFs",'.pdf'),height = 10,width = 20)
  x1 <- ggplot(dat) +
    geom_segment(aes(x = s, y = y, xend = e, yend = y),alpha=alpha/2)+
    geom_point(aes(x=s, y=y,color="Start"),alpha=alpha,show.legend = T) +
    geom_point(aes(x=e, y=y,color="End"),alpha=alpha,show.legend = T) +
    #geom_text_repel(aes(s,y, label = gene),size = 1,segment.size = 0.1,box.padding = unit(0.00, 'lines'),
    #                point.padding = unit(.02, 'lines')) +
    labs(title = paste0( nrow(a.tfs)," events, Fold Change vs time, ",length(unique(a.tfs$gene))," unique genes\n",n),x="hour", y = "Log2FC") +
    theme_classic(base_size = 16)
  x2 <- ggplot(dat) +
    stat_ecdf(data =dat ,aes(s,color=UpDown)) +
    labs(title = paste0( "Cum. Density Starts vs time"),x="hour", y = "Density") 
    #theme_classic(base_size = 16)
  x3 <- ggplot(dat) +
    stat_ecdf(data =dat ,aes(e,color=UpDown)) +
    labs(title = paste0( "Cum. Density Ends vs time"),x="hour", y = "Density") 
  x4 <- ggplot(dat) +
    geom_density(data =dat,aes(e-s, colour=UpDown))+
    labs(title = paste0("Density of Event Lengths","(Mean Duration=",round(mean(dat$e-dat$s),2),"hr)"),x="hour", y = "Density")+
    geom_vline(data=cdat, aes(xintercept=duration.mean,color=UpDown),
               linetype="dashed", size=1)
  x5 <- ggplot(dat) +
    geom_histogram(data =dat,aes(slp))+
    labs(title = paste0("Histogram of Slopes"," (Median abs log2FC/hr=",median(abs(dat$slp)),"),(Max abs log2FC/hr=",max(abs(dat$slp)),")" ),x="Log2FC/hr", y = "#Events")+
    geom_vline(data=dat, aes(xintercept=mean(slp)),
               linetype="dashed", size=1)
  crz <- cbind(melt(colMeans(cor(datasets[[n]],method="spe"))),tp=tpDatasets[[n]]*24)
  x6 <- ggplot(data = crz)+
    geom_point(aes(x=tp,y=value,colour="cor"))+
    labs(title = "Mean Cor to Dataset",x="hour", y = "cor") 
  print(multiplot(x1,x2,x3,x4,x5,x6, layout=cbind(c(1,1,1,6), c(2,3,4,5)) ))
  dev.off()
  #ggsave(paste0('~/Desktop/',n,"seg",'.pdf'),multiplot(x1,x2,cols=1),device = "pdf")  
  mypar(1,1)
  print(paste(n,max(abs(a$slp))))
  dev.off()
  #hist(dat$e-dat$s,main =paste0(n))
  a.tfs
  })

listFilters(mart)[grepl("mgi",listFilters(mart)[,2]),]

inds <- order(abs(outputs[[10]]$fc),decreasing = T)[1:10]
gnz <- outputs[[10]][inds,"gene"]
gnz="ACTA1"
for(g in gnz){
  plotmarker(data = log(datasets[[25]]+1,2),seg = T,t.vect = tpDatasets[[25]]*24,listname = c(g))
  plotmarker(data = log(datasets[[27]]+1,2),seg = T,t.vect = tpDatasets[[27]]*24,listname = c(g))
}

std.heatmap(cor(datasets[[28]],method = "spearman"))

for(g in unique(TRRUST.interactions$from)){
  mypar(1,1)
  plot(tpDatasets[[27]],datasets[[27]][g,],main=g)
  mypar(3,3)
  for(regd in TRRUST.interactions[TRRUST.interactions$from==g,"to"]){
    plot(tpDatasets[[27]],datasets[[27]][regd,],main=paste(regd,TRRUST.interactions[TRRUST.interactions$from==g & TRRUST.interactions$to==regd,c("from","type")]),plot=F)
    lines(tpDatasets[[27]],datasets[[27]][g,],col="red")
  }
}

for(g in unique(TRRUST.interactions$to)[10:15]){
  mypar(1,1)
  plot(tpDatasets[[27]],datasets[[27]][g,],main=g)
  regd <- TRRUST.interactions[ TRRUST.interactions$to==g,"from"]
  print(cor(datasets[[27]][g,],t(datasets[[27]][regd,])))
  #matplot(tpDatasets[[27]],log(cbind(datasets[[27]][g,],t(datasets[[27]][regd,]))+1,2),col = c(1,rep(2,200)))
}

for(g in unique(TRRUST.interactions$to)[10:15]){
  mypar(1,1)
  for(regd in TRRUST.interactions[ TRRUST.interactions$to==g,"from"]){
    print(TRRUST.interactions[ TRRUST.interactions$from==regd & TRRUST.interactions$to==g,])
    print(a[a$gene==g,])
    print(a[a$gene%in%regd,])
  }
}

plot(datasets[[28]]["PAX6",])
dea[a$gene=="AR",]

sort(abs(outputs[[11]]$slp),decreasing = T)[1:50]

sort(rowMeans(datasets[[11]][HumanTFs[HumanTFs%in% rownames(datasets[[11]])],]),T)
plot(datasets[[27]]["HMGB1",])

a <- readRDS("~/code/data/cb/NBNogData/segregevents/MouseTControlEvents.rds")
a.events <- detectEvents(a,tpDatasets[[14]])
gene.names <- unlist(lapply(names(a.events), function(n){ rep(n, nrow(a.events[[n]]))} ))
de <- Reduce(rbind, a.events)
de <- data.frame(de,stringsAsFactors = F)
de <- cbind(gene.names, de)
colnames(de) <- c("gene","start","end", "slp", "slp.pval","fc")



rownames(assay(experimentSummary[[1]]))
sort(colData(experimentSummary[[1]])$developmental_stage)
colData(experimentSummary[[1]])
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
rnSymbol <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),filters = c("ensembl_gene_id"),values =rownames(assay(experimentSummary[[1]])) ,mart = mart) 
inds <- rnSymbol$hgnc_symbol%in%rownames(datasets[[11]])
brainMap <- assay(experimentSummary[[1]])
brainMap <-  brainMap[rnSymbol$ensembl_gene_id,]
brainMap <-  brainMap[inds,]
colnames(brainMap) <- paste(colData(experimentSummary[[1]])$organism_part,colData(experimentSummary[[1]])$developmental_stage)
brainMap <- brainMap[,order(colnames(brainMap))]
rownames(brainMap) <- rnSymbol$hgnc_symbol[inds]

std.heatmap((cor.compare(brainMap[neural_list[neural_list%in%rownames(brainMap)],],MixSetsEC[[1]][neural_list[neural_list%in%rownames(brainMap)],which(colnames(MixSetsEC[[2]])%in%c( "10_percent_d4" ,"10_percent_d10","10_percent_d18" ,"10_percent_d24" ,"10_percent_d32" ,"100_percent_d4","100_percent_d10","100_percent_d18","100_percent_d24","100_percent_d32" ))],method="spe")),cexRow=.3)

datasets[[11]]["PAX7",]

a <- read.table("~/code/data/cb/NBNogData/segregevents/MouseMegaEndodermEvents.txt")
b <- read.table("~/code/data/cb/NBNogData/segregevents/HumanMegaEndodermEvents.txt")
a <- readRDS("~/code/data/cb/NBNogData/segregevents/33DayHumanEvents.rds")
hist(a$slp)
hist(b$slp)
max(abs(a.sig$slp))
max(abs(b.sig$slp))
a.sig <- a[a$slp.pval<.05,]
b.sig <- b[b$slp.pval<.05,]

a.tfs <- a.sig[a.sig$gene%in%union(MouseTFs,HumanTFs),]
b.tfs <- b.sig[b.sig$gene%in%union(MouseTFs,HumanTFs),]
a.tfs <- a.tfs[order(a.tfs$end),]
b.tfs <- b.tfs[order(b.tfs$end),]


a[order(a$slp,decreasing = T),]
b[order(b$slp,decreasing = T),]

qqplot(a.tfs$end,b.tfs$end)
abline(1,1)
hist(a.tfs$start)
hist(b.tfs$start)
barplot2(rbind(table(cut( a.tfs$start,0:43,right=F)), table(cut( b.tfs$start,0:43,right=F))),beside = T)
barplot2(rbind(table(cut( a.tfs$end,0:43,right=F)), table(cut( b.tfs$end,0:43,right=F))),beside = T)

length(unique(a.sig$gene))
length(unique(b.sig$gene))


compareTS <- function(a,b,pval.thresh =.1, fc.thresh=2){
  sig.a <- a[(abs(a$fc)> fc.thresh & a$slp.pval < pval.thresh),] 
  sig.b <- b[(abs(b$fc)> fc.thresh & b$slp.pval < pval.thresh),] 
  d.s.a <- outer(sig.a$start,sig.a$start,'-')
  d.s.b <- outer(sig.b$start,sig.b$start,'-')
  rownames(d.s.a) <- colnames(d.s.a) <- paste(sig.a$gene, ifelse(sign(sig.a$fc)==1 ,1,0) ,sep = "_")
  rownames(d.s.b) <- colnames(d.s.b) <-  paste(sig.b$gene, ifelse(sign(sig.b$fc)==1 ,1,0) ,sep = "_")
  
  #get matching events
  match.inds <-  rbind(cbind(1:nrow(d.s.a),match(rownames(d.s.a), rownames(d.s.b))), cbind(match(rownames(d.s.b), rownames(d.s.a)),1:nrow(d.s.b)))
  match.inds <- match.inds[complete.cases(match.inds),]
  match.inds <-  unique(match.inds,MARGIN = 1)
  
  d.s.a=d.s.a[match.inds[,1],match.inds[,1]]
  d.s.b=d.s.b[match.inds[,2],match.inds[,2]]
  
  sig.a <- sig.a[rownames(sig.a)[match.inds[,1]],]
  sig.b <- sig.b[rownames(sig.b)[match.inds[,2]],]
  
  #Events without a corresponder in other set
  rownames(d.s.a)[!rownames(d.s.a) %in% rownames(d.s.a)[match.inds[,1]]]
  rownames(d.s.b)[!rownames(d.s.b) %in% rownames(d.s.b)[match.inds[,2]]]
  
  
  tab.a <- table(sig.a$gene)
  gnz <- union(sig.a$gene,sig.b$gene)
  final.matches <- sapply(gnz, function(x){
    print(x)
    a.events <- sig.a[sig.a$gene==x,]
    b.events <- sig.b[sig.b$gene==x,]
    compare.ratios <- rn.compare(t(d.s.a[sig.a$gene==x,,drop=F]),t(d.s.b[sig.b$gene==x,,drop=F]))
    obs.events <-  c(colnames(compare.ratios[[1]]),colnames(compare.ratios[[2]]))
    #If the same class of event happens in both conditions
    sapply(obs.events,function(g){
      if(g%in%colnames(compare.ratios[[1]]) & g%in%colnames(compare.ratios[[2]])){
        e1 <- which(g==colnames(compare.ratios[[1]]))
        e2 <- which(g==colnames(compare.ratios[[2]]))
        inds=expand.grid(e1,e2)
        s <- lapply(1:nrow(inds),function(ind)match.subtract(compare.ratios[[1]][,inds[ind,1]],compare.ratios[[2]][,inds[ind,2]]))
        s <- sapply(s,function(ss) sum(abs(ss))/length(ss))
        c(rownames(a.events)[inds[which.min(s),1]],rownames(b.events)[inds[which.min(s),2]],min(s))
      }else{
        e1 <- which(g==colnames(compare.ratios[[1]]))
        e2 <- which(g==colnames(compare.ratios[[2]]))
        inds=expand.grid(e1,e2)
        return(c(rownames(a.events)[inds[1,1]],rownames(b.events)[inds[1,2]],NA))
      }
    })
  })
  final.matches <- final.matches[!sapply(final.matches,function(x)is.null(nrow(x)))]
  
  f.m=Reduce(cbind,final.matches)
  f.m[2,!is.na(f.m[3,])]
  differences <-  sig.a[unlist(f.m[1,!is.na(a[3,])]),]$start-sig.b[unlist(f.m[2,!is.na(f.m[3,])]),]$start
  end.differences <-  sig.a[unlist(f.m[1,!is.na(a[3,])]),]$end-sig.b[unlist(f.m[2,!is.na(f.m[3,])]),]$end
  names(differences) <- colnames(f.m)
  names(end.differences) <- colnames(f.m)
  summary(differences)
  summary(end.differences)
  hist(differences)
  hist(end.differences)
  gene_list <-   gsub("_\\d","",colnames(f.m))
  
  
  # plotlist <- lapply( gene_list[gene_list%in% rownames(MixSets$hs) & gene_list%in% rownames(MixSets$mm)],function(gn){
  #   meanGenes <- apply(condit.combos,1,function(x){
  #     #ts <- aggregate(GetNormalizedMat( MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ], MedianNorm(MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]))[gn,], list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), mean)
  #     ts <- aggregate(MixSets[[x[1]]][gn,conditMixSets[[x[1]]]==x[2] ], list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), mean)
  #     zoo(ts[,2],ts[,1])
  #   })
  #   ds <- Reduce(merge.zoo,meanGenes)
  #   colnames(ds) <- apply(condit.combos,1,function(aaa) paste0(aaa[1],aaa[2]))
  #   ds <- data.frame(ds)
  #   ds <- cbind("tme"=rownames(ds) , ds)
  #   melted <-  melt(ds)
  #   sdGenes <- apply(condit.combos,1,function(x){
  #     #ts <- aggregate(GetNormalizedMat( MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ], MedianNorm(MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]))[gn,], list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), function(x)sd(x)/sqrt(length(x)))
  #     ts <- aggregate(MixSets[[x[1]]][gn,conditMixSets[[x[1]]]==x[2] ], list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), function(x)sd(x)/sqrt(length(x)))
  #     zoo(ts[,2],ts[,1])
  #   })
  #   sdGenes <- Reduce(merge.zoo,sdGenes)
  #   colnames(sdGenes) <- apply(condit.combos,1,function(aaa) paste0(aaa[1],aaa[2]))
  #   sdGenes <- as.data.frame(sdGenes)
  #   sdGenes <- cbind("tme"=rownames(sdGenes) , sdGenes)
  #   meltSD <-  melt(sdGenes)
  #   melted <-  cbind(melted,"sd"= meltSD$value)
  #   #melted$tme <- factor(melted$tme,levels = melted$tme)
  #   melted$tme <- as.numeric(as.character(melted$tme))
  #   #Plotting Error Bars
  #   pd <- position_dodge(0.5) # move them .05 to the left and right
  #   sig <- sig.a[sig.a$gene==gn,,drop=F]
  #   sig <- sig[1,]
  #   sig2 <- sig.b[sig.b$gene==gn,,drop=F]
  #   sig2 <- sig2[1,]
  #   gg <-  ggplot(melted, aes(x=tme, y=value, group=variable,color=variable)) + 
  #     geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.05, position=pd) +
  #     geom_line(position=pd,size=1.2) +
  #     geom_point(position=pd,size=3) +
  #     theme_bw() +
  #     annotate("rect", xmin = sig$start, xmax = sig$end, ymin = 0, ymax = 50,
  #              alpha = .2,col="blue")+
  #     annotate("rect", xmin = sig2$start, xmax = sig2$end, ymin = 0, ymax = 50,
  #              alpha = .2,col="red")+
  #     #scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  #     labs(x = "Day", y = "TPM", 
  #          title = paste(gn,differences[grepl(gn,names(differences))],end.differences[grepl(gn,names(end.differences))]))
  #   gg
  #})
  #pdf("~/Desktop/checkSegReg.pdf", onefile = TRUE)
  #marrangeGrob(grobs=plotlist , nrow=1, ncol=1)
  #dev.off()
}




# bps <- read.table('~/code/common_up_breakpoints.txt',header = T,row.names = 1)
# amps <- read.table('~/code/common_amplitudes.txt',header=T,row.names = 1)
bps <- read.table('~/code/vitro_bp_short.txt',header = T,row.names = 1)
amps <- read.table('~/code/vitro_pk_short.txt',header=T,row.names = 1)
colnames(bps) <- c("mouse","human")
colnames(amps)[2] <- "mouse"
colnames(amps)[4] <- "human"
# tera.bps <- read.table('~/code/tera_breakpoints.txt',header = T,row.names = 1)
# tera.amps <- read.table('~/code/tera_peaks.txt',header=T,row.names = 1)
tera.bps <- read.table('~/code/tera_bp_short.txt',header = T,row.names = 1)
tera.amps <- read.table('~/code/tera_pk_short.txt',header=T,row.names = 1)
colnames(tera.bps) <- c("mouse","human")
colnames(tera.amps)[2] <- "mouse"
colnames(tera.amps)[4] <- "human"
type = "longadjusted"
colnames(amps) <- c("mouse.amp","mouse","human.amp","human")  
colnames(tera.amps) <- c("mouse.amp","mouse","human.amp","human")  

selAmps <- amps[order(amps$human.amp,decreasing = T),]
#selAmps <- selAmps[rownames(selAmps)%in%neural_related_genes,]
selAmps <- selAmps[(selAmps$mouse.amp>3) & (selAmps$human.amp>3),]
#selAmps <- selAmps[1:50,]
#selAmps$mouse <- selAmps$mouse+1
#selAmps$human <- selAmps$human

c.e <- carnegie.equivalents
c.e <- c.e[c.e$carnegie.stage %in% 6:21,]
c.e$human <- c.e$human-14 
c.e$mouse <- c.e$mouse-6
plot(c.e$human,c.e$mouse,col="red",type = "l")
text(c.e$human,c.e$mouse,c.e$carnegie.stage,col="red")
size <- sapply(1:length(selAmps$mouse), function(i) { sum(selAmps$mouse==selAmps$mouse[i] & selAmps$human==selAmps$human[i]) })
points(selAmps$human,selAmps$mouse,cex=log(size+1))
mod=lm(formula = mouse ~ human+0,data = selAmps)
abline(mod,lty=1)
text(2,1,paste0(mod$coef["human"],"x+",mod$coef["(Intercept)"]," Rsq= ",summary(mod)$r.squared), cex = .6)
mod = lm(formula = mouse ~ human+0,data = c.e)
text(2,8,paste0(mod$coef["human"],"x+",mod$coef["(Intercept)"]," Rsq= ",summary(mod)$r.squared),col="red", cex = .6)
abline(mod,col="red",lty=2)


#write.table(bps[rownames(bps)%in%neural_related_genes,],'~/code/vitro_breakpoints_neural.txt',sep = "\t")
#write.table(amps[rownames(amps)%in%neural_related_genes,],'~/code/vitro_peaks_neural.txt',sep = "\t")
#write.table(tera.bps[rownames(tera.bps)%in%neural_related_genes,],'~/code/tera_breakpoints_neural.txt',sep = "\t")
#write.table(tera.amps[rownames(tera.amps)%in%neural_related_genes,],'~/code/tera_peaks_neural.txt',sep = "\t")

pdf('~/Desktop/short_hists.pdf')
hist(bps[,2]-bps[,1],xlab = "Human BP - Mouse BP",main = "Vitro BP Day Differences")
abline(v=mean(bps[,2]-bps[,1]),col="red")
hist(amps[,4]-amps[,2],xlab = "Human Pk - Mouse Pk",main = "Vitro Peak Day Differences")
abline(v=mean(amps[,4]-amps[,2]),col="red")

gnz <- rownames(tera.bps)[rownames(tera.bps)%in%neural_related_genes]
hist(tera.bps[,2]-tera.bps[,1],xlab = "Human BP - Mouse BP",main = "Tera BP Day Differences")
abline(v=mean(tera.bps[,2]-tera.bps[,1]),col="red")
hist(tera.amps[,4]-tera.amps[,2],xlab = "Human Pk - Mouse Pk",main = "Tera Peak Day Differences")
abline(v=mean(tera.amps[,4]-tera.amps[,2]),col="red")

gnz <- rownames(bps)[rownames(bps)%in%neural_related_genes]
hist(bps[gnz,2]-bps[gnz,1],xlab = "Human BP - Mouse BP",main = "Vitro Neural BP Day Differences")
abline(v=mean(bps[gnz,2]-bps[gnz,1]),col="red")
gnz <- rownames(amps)[rownames(amps)%in%neural_related_genes]
hist(amps[gnz,4]-amps[gnz,2],xlab = "Human Pk - Mouse Pk",main = "Vitro Tera Neural Peak Day Differences")
abline(v=mean(amps[gnz,4]-amps[gnz,2]),col="red")

gnz <- rownames(tera.bps)[rownames(tera.bps)%in%neural_related_genes]
hist(tera.bps[gnz,2]-tera.bps[gnz,1],xlab = "Human BP - Mouse BP",main = "Tera Neural BP Day Differences")
abline(v=mean(tera.bps[gnz,2]-tera.bps[gnz,1]),col="red")
gnz <- rownames(tera.amps)[rownames(tera.amps)%in%neural_related_genes]
hist(tera.amps[gnz,4]-tera.amps[gnz,2],xlab = "Human Pk - Mouse Pk",main = "Tera Neural Peak Day Differences")
abline(v=mean(tera.amps[gnz,4]-tera.amps[gnz,2]),col="red")
dev.off()


bps$human <- approx(carnegie.equivalents$human,carnegie.equivalents$mouse,bps$human+14)$y
bps$mouse <- bps$mouse+6
bps <- bps[bps$mouse<=16,]

blk.gnz <- c("LIX1","ACTA2","DLK1","PAX6","LHX2","NR2F1","POU3F2","ST18","EPHA3","ZFHX3")
neural_related_genes <- unique(c(neural_related_genes,blk.gnz))
colla <- rownames(bps)[rownames(bps)%in%neural_related_genes]
colla <- ifelse(colla%in% blk.gnz,"black","grey")
nmz <- rownames(bps)
nmz <- ifelse(nmz%in% blk.gnz,nmz,"")


dat <-  data.frame(human=bps[,2],humany=jitter(rep(0,nrow(bps)),factor = 4), mouse=bps[,1],mousey=jitter(rep(1,nrow(bps)),factor = 4)) 
x <- ggplot(dat) +
  geom_point(aes(x=human, y=humany), color = 'red',show.legend = T) +
  geom_point(aes(x=mouse, y=mousey), color = 'blue') +
  geom_text_repel(aes(mouse,mousey, label = nmz),size = 4,segment.size = 0.1,box.padding = unit(0.00, 'lines'),
                  point.padding = unit(.02, 'lines')) +
  geom_text_repel(aes(human,humany, label = nmz),size = 4,segment.size = 0.1,box.padding = unit(0.00, 'lines'),
                  point.padding = unit(.02, 'lines')) +
  geom_segment(aes(x = human, y = humany, xend = mouse, yend = mousey),data = dat[rownames(bps)%in%neural_related_genes,],color=colla)+
  labs(title = "BP Compare",x="peak day", y = "spec") +
  theme_classic(base_size = 16)
ggsave(paste0('~/Desktop/VitroBP',type,'.pdf'),x,device = "pdf")   


amps$human <- approx(carnegie.equivalents$human,carnegie.equivalents$mouse,amps$human+14)$y
amps$mouse <- amps$mouse+6
amps <- amps[amps$mouse<=16,]

blk.gnz <- c("LIX1","ACTA2","DLK1","PAX6","LHX2","NR2F1","POU3F2","ST18","EPHA3","ZFHX3")
colla <- rownames(amps)[rownames(amps)%in%neural_related_genes]
colla <- ifelse(colla%in% blk.gnz,"black","grey")
nmz <- rownames(amps)
nmz <- ifelse(nmz%in% blk.gnz,nmz,"")


dat <-  data.frame(human=amps[,4],humany=jitter(rep(0,nrow(amps)),factor = 4), mouse=amps[,2],mousey=jitter(rep(1,nrow(amps)),factor = 4)) 
x <- ggplot(dat) +
  geom_point(aes(x=human, y=humany), color = 'red',show.legend = T) +
  geom_point(aes(x=mouse, y=mousey), color = 'blue') +
   geom_text_repel(aes(mouse,mousey, label = nmz),size = 4,segment.size = 0.1,box.padding = unit(0.00, 'lines'),
                   point.padding = unit(.02, 'lines')) +
   geom_text_repel(aes(human,humany, label = nmz),size = 4,segment.size = 0.1,box.padding = unit(0.00, 'lines'),
                   point.padding = unit(.02, 'lines')) +
  geom_segment(aes(x = human, y = humany, xend = mouse, yend = mousey),data = dat[rownames(amps)%in%neural_related_genes,],colour=colla)+
  labs(title = "Peak Compare",x="peak day", y = "spec") +
  theme_classic(base_size = 16)
ggsave(paste0('~/Desktop/VitroPK',type,'.pdf'),x,device = "pdf")


hist(tera.bps[,2]-tera.bps[,1])
abline(v=mean(tera.bps[,2]-tera.bps[,1]),col="red")

tera.bps$human <- approx(carnegie.equivalents$human,carnegie.equivalents$mouse,tera.bps$human+14)$y
tera.bps$mouse <- tera.bps$mouse+6
tera.bps <- tera.bps[tera.bps$mouse<=25,]
colla <- rownames(tera.bps)[rownames(tera.bps)%in%neural_related_genes]
colla <- ifelse(colla%in% blk.gnz,"black","grey")
nmz <- rownames(tera.bps)
nmz <- ifelse(nmz%in% blk.gnz,nmz,"")


dat <-  data.frame(human=tera.bps[,2],humany=jitter(rep(0,nrow(tera.bps)),factor = 4), mouse=tera.bps[,1],mousey=jitter(rep(1,nrow(tera.bps)),factor = 4)) 
x <- ggplot(dat) +
  geom_point(aes(x=human, y=humany), color = 'red',show.legend = T) +
  geom_point(aes(x=mouse, y=mousey), color = 'blue') +
     geom_text_repel(aes(mouse,mousey, label = nmz),size = 4,segment.size = 0.1,box.padding = unit(0.00, 'lines'),
                     point.padding = unit(.02, 'lines')) +
     geom_text_repel(aes(human,humany, label = nmz),size = 4,segment.size = 0.1,box.padding = unit(0.00, 'lines'),
                     point.padding = unit(.02, 'lines')) +
  geom_segment(aes(x = human, y = humany, xend = mouse, yend = mousey),data = dat[rownames(tera.bps)%in%neural_related_genes,],color=colla)+
  labs(title = "BP Compare",x="peak day", y = "spec") +
  theme_classic(base_size = 16)
ggsave(paste0('~/Desktop/TeraBP',type,'.pdf'),x,device = "pdf")


tera.amps$human <- approx(carnegie.equivalents$human,carnegie.equivalents$mouse,tera.amps$human+14)$y
tera.amps$mouse <- tera.amps$mouse+6
tera.amps <- tera.amps[tera.amps$mouse<=100,]
nmz <- rownames(tera.amps)
nmz <- ifelse(nmz%in% blk.gnz,nmz,"")
colla <- rownames(tera.amps)[rownames(tera.amps)%in%neural_related_genes]
colla <- ifelse(colla%in% blk.gnz,"black","grey")


dat <-  data.frame(human=tera.amps[,4],humany=jitter(rep(0,nrow(tera.amps)),factor = 4), mouse=tera.amps[,2],mousey=jitter(rep(1,nrow(tera.amps)),factor = 4)) 
x <- ggplot(dat) +
  geom_point(aes(x=human, y=humany), color = 'red',show.legend = T) +
  geom_point(aes(x=mouse, y=mousey), color = 'blue') +
     geom_text_repel(aes(mouse,mousey, label = nmz),size = 4,segment.size = 0.1,box.padding = unit(0.00, 'lines'),
                     point.padding = unit(.02, 'lines')) +
     geom_text_repel(aes(human,humany, label = nmz),size = 4,segment.size = 0.1,box.padding = unit(0.00, 'lines'),
                     point.padding = unit(.02, 'lines')) +
  geom_segment(aes(x = human, y = humany, xend = mouse, yend = mousey),data = dat[rownames(tera.amps)%in%neural_related_genes,],colour=colla)+
  labs(title = "Peak Compare",x="peak day", y = "spec") +
  theme_classic(base_size = 16)
ggsave(paste0('~/Desktop/TeraPK',type,'.pdf'),x,device = "pdf")







internal.dist <- function(x){
  as.dist(outer (x, x, `-`))
}

access.diag.mat <- function(i,j,mat){
  n=nrow(mat)
  if(i<j){
    mat[n*(i-1)-(i*(i+1)/2)+j]
  }else{
    -(mat[n*(j-1)-(j*(j+1)/2)+i ])
  }
}

tfl <- HumanTFs[HumanTFs%in%union(rownames(datasets[[1]]),rownames(datasets[[11]]))]
check <- datasets[[11]][tfl,]
ds1Mat <- lapply(1:ncol(datasets[[1]]),function(i){
  print(i)
  internal.dist(log(datasets[[1]]+1,2)[tfl,i])
})
#lapply if keeping whole matrices (very large)
gene10 <- sapply(1:length(tfl), function(q){
  print(q)
  colMeans(sapply(ds2Mat,function(x) unlist(sapply(1:ncol(x), access.diag.mat ,q,x))))
  })
newMat <- sapply(gene10,function(g) colMeans(g))

std.heatmap(cor(t(gene10)))
std.heatmap(cor(log(check+1,2)))
std.heatmap(cor(datasets[[11]]))

pc <- princomp(gene10[[1]])
plot(check[30,])
plot(newMat[,9])

matplot(gene10[,100:112])
matplot(t(check[100:112,]))


ds2Mat <- lapply(1:ncol(datasets[[11]]),function(i){
  print(i)
  internal.dist(log(datasets[[11]]+1,2)[tfl,i])
})


a <- internal.dist(datasets[[1]][tfl,1])
b <- internal.dist(datasets[[11]][tfl,1])

min(a-b)

a=kmeans(datasets[[1]],2000,iter.max = 40)
a$tot.withinss
order(a$withinss)
matplot(t(a$centers))
sort(table(a$cluster))
a$tot.withinss/a$totss
matplot(t(datasets[[1]][names(which(a$cluster==80)),]))

dif <- internal.dist(a$centers[,1])-internal.dist(a$centers[,2])
dif2 <- internal.dist(a$centers[,2])-internal.dist(a$centers[,3])

dif <- as.matrix(dif)
lapply(which(abs(dif)>5),function(x){
  ind=arrayInd(x, dim(dif))
  list(rownames(dif)[ind[1,1]],colnames(dif)[ind[1,2]])
})
q=rbind(dif,dif2)
dim(q)
length(dif2)
oscopeOut <- OscopeSine(a$centers)



a=Oscope::OscopeSine()
?Oscope
condit.combos <- expand.grid(c("mm","hs"),unique(conditMixSets[[1]]),stringsAsFactors = F)
condit.combos <- condit.combos[-c(2,7),]
apply(condit.combos,1,function(x){
  SegRegEvents(set1 = list(log(GetNormalizedMat( MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ], MedianNorm(MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]))+1,2)), tp1 = list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]),output.dir = file.path(datapath,"segregevents"), fn=paste0(x[1],x[2],"Events")) 
  })

apply(condit.combos,1,function(x){
  SegRegEvents(set1 = list(log(GetNormalizedMat( MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ], MedianNorm(MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]))+1,2)), tp1 = list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]),output.dir = file.path(datapath,"segregevents"), fn=paste0(x[1],x[2],"Events")) 
})

tpm.threshold =3 
pvalcut = .01
maxk = 4
min.num.in.seg = 3
cutdiff = .001
x <- unlist(condit.combos[6,])
set <- log(GetNormalizedMat( MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ], MedianNorm(MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]))+1,2)
tp <- tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]
genes <-  a[a$start>2 & a$start<10 & a$fc>4 & a$end-a$start<10,]$gene
genes <- genes[genes %in% rownames(set)[rowMeans(set[,1:9])<1]]
plotmarker(set,t.vect = tp, listname =genes,pvalcut = pvalcut,maxk =maxk,min.num.in.seg =  min.num.in.seg,cutdiff=cutdiff,num.try = 100,keepfit = F,filename = "~/code/markers.pdf",pdf = T)

b= sigSetsEvents$mm0
tpm.threshold =3 
pvalcut = .01
maxk = 4
min.num.in.seg = 3
cutdiff = .001
x <- unlist(condit.combos[1,])
set <- log(GetNormalizedMat( MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ], MedianNorm(MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]))+1,2)
tp <- tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]
m.genes <-  b[b$start>1 & b$start<10 & b$fc>4 & b$end-b$start<10,]$gene
m.genes <- m.genes[m.genes %in% rownames(set)[rowMeans(set[,1:9])<1] & m.genes %in% genes]
plotmarker(set,t.vect = tp, listname =m.genes,pvalcut = pvalcut,maxk =maxk,min.num.in.seg =  min.num.in.seg,cutdiff=cutdiff,num.try = 100,keepfit = F,filename = "~/code/both_markers.pdf",pdf = T)

plotmarker(set,t.vect = tp, listname ="PAX6",pvalcut = pvalcut,maxk =maxk,min.num.in.seg =  min.num.in.seg,cutdiff=cutdiff,num.try = 100)


setwd("~/code/data/cb/segregevents")
setsEvents <- lapply(dir(), read.table,header = T,stringsAsFactors = F)
names(setsEvents) <- gsub("Events.txt","",dir())
sigSetsEvents <- lapply(setsEvents, function(x)x[(abs(x$fc)>3 & x$slp.pval < .1),] )




lapply(names(sigSetsEvents),function(x) hist(sigSetsEvents[[x]]$end,main = x))

head(sigSetsEvents$hs100)
lapply(names(sigSetsEvents), function(x)plot(table(.bincode(sigSetsEvents[[x]]$start,seq(0,42,1),include.lowest = T)),ylab = "#Genes",xlab="day",main = x ))
lapply(names(sigSetsEvents), function(x)plot(table(.bincode(sigSetsEvents[[x]]$end,seq(0,42,1),include.lowest = T)),ylab = "#Genes",xlab="day",main = x ))



x <- ggplot(data.frame(start=de$start, fc=de$fc)) +
  geom_point(aes(start, fc), color = 'red') +
  geom_text_repel(aes(start, fc, label = de$gene),size = 1.8,segment.size = 0.1,box.padding = unit(0.02, 'lines'),
                  point.padding = unit(.02, 'lines')) +
  labs(title = "100% Human DE Genes",x="Start (Day)", y = "log2 FC") +
  theme_classic(base_size = 16)
  ggsave(filename = "~/Desktop/hs100.pdf", plot = x,width = 45, height = 45,device = "pdf")


x <- ggplot(data.frame(start=de2$start, fc=de2$fc)) +
  geom_point(aes(start, fc), color = 'red') +
  geom_text_repel(aes(start, fc, label = de2$gene),size = 1.8,segment.size = 0.1,box.padding = unit(0.02, 'lines'),
                  point.padding = unit(0.02, 'lines')) +
  labs(title = "10% Human DE Genes",x="Start (Day)", y = "log2 FC") +
  theme_classic(base_size = 16)
ggsave(filename = "~/Desktop/hs10.pdf", plot = x,width = 45, height = 45,device = "pdf")
  
  


segreg(MixSets$hs[c("NEUROD4","NEUROG2"),c(76,27:51)],meancut = 0,t.vect = tps$mm0,cutdiff = .1,min.num.in.seg = 3,maxk = 4)
nrow(de[NULL,])
join <- data.frame("e1"=integer(),"e2"=integer())

  for(gn in intersect(de$gene,de2$gene)){
    degn1 <- de[de$gene==gn,]
    degn2 <- de2[de2$gene==gn,]
    s1 <- sign(degn1$slp)*(degn1$slp.pval<.1)
    s2 <- sign(degn1$slp)*(degn1$slp.pval<.1)
    if(nrow(degn1)>0 & nrow(degn2)>0 ){
      d=dtw(s1 , s2,step.pattern = symmetricP05,keep.internals = T)
      d$normalizedDistance
      d$costMatrix
      names(sig1) <- sig1
      names(sig2) <- sig2
      print(sig1)
      print(sig2)
      corr1=sapply(sig1,function(s1){
        c(s1 ,which.min(abs(degn1$fc[s1]-degn2$fc[d$index2[sig1]])))
      })
      corr2=sapply(sig2,function(s2){
        c(which.min(abs(degn1$fc[d$index1[s2]]-degn2$fc[s2])),s2 )
      })
    print(corr1)
    print(corr2)
    }
  }

deSignificant <-  de[(de$slp.pval<=.05 & abs(de$slp*(de$end -de$start)) > 3),]
deSignificant <-  deSignificant[order(deSignificant$end),]
deSignificant <- deSignificant[order(rank(deSignificant$start,ties.method = "first" )),]


  d$normalizedDistance
  d$index1s
  d$index2s
  dim(deSignificant)
  

  
  
  a=SegReg::topsegreg(srO)
  srO
  de2<- detectEvents(srO2,tpMixSets[77:102])
  gene.names <- unlist(lapply(names(de2), function(n){ rep(n, nrow(de2[[n]]))} ))
  de2 <- Reduce(rbind, de2)
  de2 <- data.frame(de2,stringsAsFactors = F)
  de2 <- cbind(gene.names, de2)
  colnames(de2) <- c("gene","start","end", "slp", "slp.pval")
  de2Significant$start
  

  de2Significant <-  de2[(de2$slp.pval<=.05 & abs(de2$slp*(de2$end -de2$start)) > 2),]
  
  deSignificant
  de2Significant
  
  #rownames(deSignificant) <- deSignificant$gene
  #rownames(de2Significant) <- de2Significant$gene
  deSignificant <-  deSignificant[order(deSignificant$end),]
  de2Significant <-  de2Significant[order(de2Significant$end),]
  deSignificant[order(rank(deSignificant$start,ties.method = "first" )),]
  
  
  
  order(rank(de2Significant$start,ties.method = "first"))
  
  
  
  
  
  tp <- list(tp1)
  set <- list(set1)
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
  tp.i=NULL
  dir.create(output.dir)
  if(file.exists(file.path(output.dir,"fits.RData"))){
    load(file.path(output.dir, "fits.RData"))
  }else{
    fits <- mclapply(frame[[1]], smooth.model, tp.i)
    save(fits,name1,name2,file=file.path(output.dir,"fits.RData"))
  }
  gene.list <- gene.list[(gene.list%in%names(fits[[1]]))]

#   filter.fits <- lapply(frame[[1]],function(x){
#     smoothed <-  tryCatch({   tsSmooth(StructTS(zoo(x$y,x$tme)))}, error = function(e){tsSmooth(StructTS(zoo(x$y,x$tme),type = "level"))}  )
#     
#   })
#   mypar(3,3)
#   for(i in 1:99){
#     plot(filter.fits[[i]][,1],main= cor(fits[[i]]$mu$y,filter.fits[[i]][,1]))
#     plot(fits[[i]]$mu$y)
#     plot(frame[[1]][[i]]$y)
#   }
  
  f.values <-  sapply(fits, function(x){
      x$anova$`Pr(>F)`[1]
  })
  
  dt.fc <- lapply(fits,function(gn){
    (max(gn$mu$y)-min(gn$mu$y))
  })

  f.ind <-  f.values < sqrt(.05)
  f.ind[is.na(f.ind)] <- T
  fitind <- dt.fc > fc.threshold &  f.ind
  fits <- fits[fitind]

  exprsmat <- t(sapply(fits,function(x) x$mu$y))
  exprsmat <- (exprsmat-rowMedians(exprsmat))/rowSds(exprsmat)
  dtwrtn <- DTW.dist(exprsmat)
  dtw.distances <- dtwrtn[[1]]
  dtw.warps <- dtwrtn[[2]]
  
  
  
  Hto10HCor <- sapply(1:nrow(MixSets$hs), function(x) cor( cbind(MixSets$hs[,77], MixSets$hs[,27:51])[x,], MixSets$hs[x,77:102]))
  Hto85HCor <-  sapply(1:nrow(MixSets$hs), function(x) cor( cbind(MixSets$hs[,77], MixSets$hs[,52:76])[x,], MixSets$hs[x,77:102]))
  
  H10 <- cbind(MixSets$hs[,77], MixSets$hs[,27:51])
  H85 <- cbind(MixSets$hs[,77], MixSets$hs[,52:76])
  Hto10HXCor <- t(sapply(1:nrow(MixSets$hs), function(x) ccf( H10[x,], MixSets$hs[x,77:102],plot = F)$acf))
  Hto85HXCor <-  t(sapply(1:nrow(MixSets$hs), function(x) ccf( H85[x,], MixSets$hs[x,77:102],plot = F)$acf))
  apply(Hto85HXCor,1, mean, na.rm=T)
  plot(apply(Hto85HXCor,1, mean, na.rm=T))
  
  mean(Hto85HCor,na.rm = T)
  summary(Hto10HCor)
  
  
  library(STRINGdb)
  spec=STRINGdb::get_STRING_species(version = 10, "Homo sapiens")# "Mus musculus"
  string_db <- STRINGdb$new( version="10", species=9606,score_threshold=0, input_directory="" )
  maps=string_db$map(my_data_frame = data.frame("names"=neural_list),my_data_frame_id_col_names = "names")
  interactions <-  string_db$get_interactions(string_ids =  maps$STRING_id)
  rownames(maps) <- maps$names
  interactionsHS$from <- sapply(interactionsHS$from, function(x) maps$names[which(x == maps$STRING_id)])
  interactionsHS$to <- sapply(interactionsHS$to, function(x) maps$names[which(x == maps$STRING_id)])
  
  spec=STRINGdb::get_STRING_species(version = 10, "Mus musculus")# "Mus musculus"
  string_db <- STRINGdb$new( version="10", species=9606,score_threshold=0, input_directory="" )
  maps=string_db$map(my_data_frame = data.frame("names"=neural_list),my_data_frame_id_col_names = "names")
  interactionsMM <-  string_db$get_interactions(string_ids =  maps$STRING_id)
  rownames(maps) <- maps$names
  interactionsMM $from <- sapply(interactionsMM$from, function(x) maps$names[which(x == maps$STRING_id)])
  interactionsMM $to <- sapply(interactionsMM$to, function(x) maps$names[which(x == maps$STRING_id)])
  
  
  hist(interactions$combined_score/1000)
  interactions[interactions$from=="DCX",]
  string_db$plot_network(maps[neural_list,"STRING_id"],required_score = 400)
  
  TRRUST.interactions = read.table(file = "~/code/data/cb/Interactions/trrust_rawdata_20160502.txt",header = F,sep = "\t",stringsAsFactors = F)
  colnames(TRRUST.interactions) <- c("from","to", "type","?")
  #TRRUST.interactions[,c(1,2)] <- make.names(TRRUST.interactions[,c(1,2)])
  BIOGRID.interactions = read.csv2(file = "~/code/data/cb/Interactions/BIOGRID-ALL-3.4.136.tab2.txt",header = T,sep = "\t")
  BIOGRID.interactions.human <- BIOGRID.interactions[ (BIOGRID.interactions$Organism.Interactor.A == "9606" &BIOGRID.interactions$Organism.Interactor.B == "9606")  , c("Official.Symbol.Interactor.A","Official.Symbol.Interactor.B")]
  BIOGRID.interactions.mouse <- BIOGRID.interactions[ (BIOGRID.interactions$Organism.Interactor.A == "10090" &BIOGRID.interactions$Organism.Interactor.B == "10090")  , c("Official.Symbol.Interactor.A","Official.Symbol.Interactor.B")]
  BIOGRID.interactions.mouse <- apply(BIOGRID.interactions.mouse,2,toupper)
  BIOGRID.interactions.human <- apply(BIOGRID.interactions.human,2,toupper)
  ORegAnno.interactions <- read.csv2(file = "~/code/data/cb/Interactions/ORegAnno_Combined_2015.12.22.tsv",header = T,sep = "\t")
  
  x = "PHOX2B"
  TRRUST.interactions[TRRUST.interactions$to==x |  TRRUST.interactions$from==x,]
  BIOGRID.interactions.human[BIOGRID.interactions.human[,2]==x | BIOGRID.interactions.human[,1]==x, ]
  
 


DBNScoreStep1(t(set1[[1]]), predPosition = which(abs(cors[i,]) > .9), targetPosition = i)
?DBNScoreStep1
net = ebdbNet::ebdbn( set1,K$dim)
nl1 = neural_list[neural_list%in% Reduce(intersect,sapply(set1,rownames))]
K <- hankel(lapply(set1, function(x) x[nl1,]), lag = 1, cutoff = 0.99, type = "median")
?hankel
ebdneu= ebdbn(lapply(set1[-c(2,3)], function(x) x[nl1,]),K = K$dim)
interactions[interactions$from=="ASCL1"&interactions$to == "FEZF2","combined_score"]

ebdbNet::plot.ebdbNet(ebdneu,sig.level = .99,interactive = F)
proxy::
  dbn1 <- DBNScoreStep1(t(ngm(datasets$`~/code/data/cb/human_ALONE.csv`)) )#, predPosition = 1:10,targetPosition = 1:11)

dbn1
ngm(datasets$`~/code/data/cb/human_ALONE.csv`)
f.values[["IHH"]]

print("# null")
print(n.null)
print(format(Sys.time(),"%D %H:%M:%S"))

system.time( dist(datasets$`~/code/data/cb/human_ALONE.csv`[1:100,],method =  pr_DB$get_entry("DTW")))
system.time(DTW.dist(datasets$`~/code/data/cb/human_ALONE.csv`[1:100,]))
plot(sapply(1:100, function(i) system.time( dist(datasets$`~/code/data/cb/human_ALONE.csv`[1:i,],method =  pr_DB$get_entry("DTW")))["elapsed"] ))


si=1:150
test.setM = list(exp_1_M[si,], exp_2_M[si,], datasets$`~/code/data/cb/mouse_ALONE.csv`)
test.setH= list( exp_1[si,1:9],exp_2_H[si,], datasets$`~/code/data/cb/human_ALONE.csv`)
setM = list(exp_1_M, exp_2_M, datasets$`~/code/data/cb/mouse_ALONE.csv`,datasets$`~/code/data/cb/TC_L_150415.csv`)
setH = list( exp_1[,1:9],exp_2_H, datasets$`~/code/data/cb/human_ALONE.csv`)
setMix = list(exp_1[,17:24],exp_2)
tpM =list(tp_1_M,tp_2,tp.m.control,tp.m.long)
tpH =list(tp_1[1:9],tp_2,tp.h.control)
tpMix = list(tp_1[17:24],tp_2)
setMTera <- list(datasets$`~/code/data/cb/mouse_TERA.csv`)
setHTera <- list(datasets$`~/code/data/cb/human_TERA.csv`)
tpMTera <- list(tp.m.tera)
tpHTera <- list(tp.h.tera)

m = read.table("~/code/data/cb/shiftFiles/NewMix10vsHumanrj7zequalized/metrics.txt",header = T)
neural_list[neural_list %in% rownames(m)[m$p.fwd.adj<.05]]


write.table( MixSets$hs[neural_list[neural_list%in%rownames(MixSets$hs)],] , "~/code/hg19_bt2_neuralonly.txt", sep = "\t")
write.table( MixSets$mm[neural_list[neural_list%in%rownames(MixSets$mm)],] , "~/code/mm10_bt2_neuralonly.txt", sep = "\t")

write.pdf({
  mypar(3,4)
  for(gn in neural_list[neural_list%in% rownames(MixSets$mm)]){
    lims <- c(min(MixSets$mm[gn,1:26]), max(MixSets$mm[gn,1:26]))
    plot(tpMixSets[1:26],MixSets$mm[gn,1:26],main = paste(gn, "0%, Mouse\nlog2 TPMs"), type = "l",xlab = "days",ylab = "log2 TPM",ylim = lims)
    lims <- c(min(MixSets$hs[gn,27:102]), max(MixSets$hs[gn,27:102]))
    plot(tpMixSets[27:51],MixSets$hs[gn,27:51],main = paste(gn, "10%, Human \nlog2 TPMs"), type = "l",xlab = "days",ylab = "log2 TPM",ylim = lims)
    plot(tpMixSets[52:76],MixSets$hs[gn,52:76],main = paste(gn, "85%, Human\nlog2 TPMs"), type = "l",xlab = "days",ylab = "log2 TPM",ylim = lims,add=T)
    plot(tpMixSets[77:102],MixSets$hs[gn,77:102],main = paste(gn, "100%, Human\nlog2 TPMs"), type = "l",xlab = "days",ylab = "log2 TPM",ylim = lims)
  }}
  ,filename = "~/Desktop/newexperimentsexpressionofmarkers.pdf")

write.pdf({
  mypar(3,4)
  for(gn in neural_list[neural_list%in% rownames(MixSets$mm)]){
    m0 = zoo(MixSets$mm[gn,1:26],tpMixSets[1:26])
    h10 <- zoo(MixSets$hs[gn,27:51],tpMixSets[27:51])  
    h85 <- zoo(MixSets$hs[gn,52:76],tpMixSets[52:76])
    h100 <- zoo(MixSets$hs[gn,77:102],tpMixSets[77:102]) 
    plot.zoo(merge(m0,h10,h85,h100),plot.type = "single",col = 1:4, lty = 1,main=paste(gn),xlab = "days",ylab = "log2 TPM" )
    legend("topleft", c("mm0","hs10","hs,85%","hs100"), col = 1:4, lty = 1)
  }}
  ,filename = "~/Desktop/newexperimentsexpressionofmarkers.pdf")

write.pdf({
  mypar(3,4)
  for(gn in neural_list[neural_list%in% rownames(MixSets$mm)]){
    m0 = zoo(tsSmooth(StructTS(MixSets$mm[gn,1:26]))[,"level"],tpMixSets[1:26])
    h10 <- zoo(tsSmooth(StructTS(MixSets$hs[gn,27:51]))[,"level"],tpMixSets[27:51])  
    h85 <- zoo(tsSmooth(StructTS(MixSets$hs[gn,52:76]))[,"level"],tpMixSets[52:76])
    h100 <- zoo(tsSmooth(StructTS(MixSets$hs[gn,77:102]))[,"level"],tpMixSets[77:102]) 
    plot.zoo(merge(m0,h10,h85,h100),plot.type = "single",col = 1:4, lty = 1,main=paste(gn),xlab = "days",ylab = "log2 TPM" )
    legend("topleft", c("mm0","hs10","hs,85%","hs100"), col = 1:4, lty = 1)
  }}
  ,filename = "~/Desktop/newexperimentsexpressionofmarkersKalman.pdf")


# mz= sapply(neural_list, function(gn){
#   if(gn %in% rownames(x)& gn%in% rownames(y)){
#   xi = zoo(tsSmooth(StructTS(x[gn,]))[,"level"], tx)
#   yi = zoo(tsSmooth(StructTS(y[gn,]))[,"level"], ty)
#   merge(xi,yi)
#   }
# })



NeuralCrestList <-   sort(as.character(read.table("~/code/data/cb/markers/NeuralCrestList.txt",header = F)[,1]))

write.pdf({
  mypar(3,4)
  for(gn in NeuralCrestList[NeuralCrestList%in% rownames(MixSets$mm)]){
    m0 = zoo(MixSets$mm[gn,1:26],tpMixSets[1:26])
    h10 <- zoo(MixSets$hs[gn,27:51],tpMixSets[27:51])  
    h85 <- zoo(MixSets$hs[gn,52:76],tpMixSets[52:76])
    h100 <- zoo(MixSets$hs[gn,77:102],tpMixSets[77:102]) 
    plot.zoo(merge(m0,h10,h85,h100),plot.type = "single",col = 1:4, lty = 1,main=paste(gn),xlab = "days",ylab = "log2 TPM" )
    legend("topleft", c("mm0","hs10","hs,85%","hs,100%"), col = 1:4, lty = 1)
  }}
  ,filename = "~/Desktop/neuralcrestexpressionNew.pdf")

for(gn in neural_list[neural_list%in% rownames(MixSets$mm)]){
  lims <- c(min(MixSets$mm[gn,1:26]), max(MixSets$mm[gn,1:26]))
  plot(tpMixSets[1:26],MixSets$mm[gn,1:26],main = paste(gn, "0%, Mouse\nlog2 TPMs"), type = "l",xlab = "days",ylab = "log2 TPM",ylim = lims)
  lims <- c(min(MixSets$hs[gn,27:102]), max(MixSets$hs[gn,27:102]))
  plot(tpMixSets[27:51],MixSets$hs[gn,27:51],main = paste(gn, "10%, Human \nlog2 TPMs"), type = "l",xlab = "days",ylab = "log2 TPM",ylim = lims)
  plot(tpMixSets[52:76],MixSets$hs[gn,52:76],main = paste(gn, "85%, Human\nlog2 TPMs"), type = "l",xlab = "days",ylab = "log2 TPM",ylim = lims)
  plot(tpMixSets[77:102],MixSets$hs[gn,77:102],main = paste(gn, "100%, Human\nlog2 TPMs"), type = "l",xlab = "days",ylab = "log2 TPM",ylim = lims)
}




datapath <- "~/code/data/cb/shiftfiles/cubed"
setwd(datapath)
datadirs <-  file.path(datapath, dir())
ns <-  dir()
metricsTables <- lapply(datadirs, function(n){
  setwd(n)
  print(file.path(n,"metrics.txt"))
  read.table(file.path(n,"metrics.txt"),header = T,sep = "\t")
})
datadirs[c(3,5,7,11)]

head(metricsTables[[1]])
fwds <- lapply(metricsTables[c(4,6,8,12)],function(x){
  rownames(x)[x$p.fwd.adj < .10]
})
bcks <- lapply(metricsTables[c(4,6,8,12)],function(x){
  rownames(x)[x$p.bck.adj < .10]
})
names(fwds) <- names(bcks) <- ns[c(4,6,8,12)]
sum(neural_list %in% fwds[[4]])
sum(neural_list %in% bcks[[4]])

f=venn(fwds)
b=venn(bcks)

f[["intersections"]]
sapply( attr(f, "intersections"), function(x) sum(neural_list%in% x))


lh
arima(lh, order = c(3,0,0))

?arima
?KalmanLike
k=dlmFilter(as.ts(zoo(frame[[1]][[2]]$y[1:8],frame[[1]][[2]]$tme[1:8])), arima(zoo(frame[[1]][[2]]$y[1:8],frame[[1]][[2]]$tme[1:8])))
plot(tsSmooth(StructTS(zoo(frame[[1]][[2]]$y[1:8],frame[[1]][[2]]$tme[1:8]))))
dlm::dlmFilter
?dlmFilter
?dlm
dlmFilter(zoo(frame[[1]][[2]]$y[1:8],frame[[1]][[2]]$tme[1:8]))

s=tsSmooth(StructTS(zoo(frame[[1]][[2]]$y[1:8],frame[[1]][[2]]$tme[1:8])))
s
plot(zoo(frame[[1]][["PAX6"]]$y[1:8],frame[[1]][["PAX6"]]$tme[1:8]))
a=StructTS(zoo(frame[[1]][["PAX6"]]$y[1:8],frame[[1]][["PAX6"]]$tme[1:8]))
plot(tsSmooth(a))
plot(smooth.model(frame[[1]][["PAX6"]],NULL)$mu$y)

?StructTS


plot(sapply(seq(0,3,by = .01), function(x)bhattacharyya.dist(0,x,.001,.001) ))
plot(sapply(seq(0,100,by = .01), function(x)dist(0,x)^8 ))

?StructTS
ari= auto.arima(zoo(frame[[1]][["PAX6"]]$y,frame[[1]][["PAX6"]]$tme))
plot(KalmanSmooth(zoo(frame[[1]][["PAX6"]]$y,frame[[1]][["PAX6"]]$tme),mod = ari$model)$smooth[,3])
StructTS(zoo(frame[[1]][["PAX6"]]$y,frame[[1]][["PAX6"]]$tme))
plot(tsSmooth(StructTS(zoo(frame[[1]][["PAX6"]]$y,frame[[1]][["PAX6"]]$tme))))
plot(zoo(frame[[1]][["PAX6"]]$y,frame[[1]][["PAX6"]]$tme))
ari$
  
  ari$model

?arima
sum(abs(ari$residuals),na.rm = T)
coef(ari)
tsp(ari)
buildFun <- function(x) {
  dlmModARMA(coef(ari), ,sigma2 = ari$sigma2,dV = exp(x[1]))
}
fit <- dlmMLE(zoo(frame[[1]][[2]]$y[1:8],frame[[1]][[2]]$tme[1:8]), parm = c(0,0), build = buildFun)
fit


datasets$`~/code/data/cb/human_ALONE.csv`
DBNScoreStep1(t(set1[[1]]), predPosition = which(abs(cors[i,]) > .9), targetPosition = i)
?DBNScoreStep1
K <- hankel(set1, lag = 1, cutoff = 0.90, type = "median")
net = ebdbNet::ebdbn( set1,K$dim)

f.values[["IHH"]]
library(STRINGdb)
spec=STRINGdb::get_STRING_species(version = 10, "Homo sapiens")# "Mus musculus"
string_db <- STRINGdb$new( version="10", species=9606,score_threshold=0, input_directory="" )
maps=string_db$map(my_data_frame = data.frame("names"=union(neural_list,NeuralCrestList)),my_data_frame_id_col_names = "names")
interactions <-  string_db$get_interactions(string_ids =  maps$STRING_id)
rownames(maps) <- maps$names
interactions$from <- sapply(interactions$from, function(x) maps$names[which(x == maps$STRING_id)])
interactions$to <- sapply(interactions$to, function(x) maps$names[which(x == maps$STRING_id)])
hist(interactions$combined_score/1000)
string_db$plot_network(maps[neural_list,"STRING_id"],required_score = 800)
#ScanBMA(exprsmat, as.numeric(colnames(exprsmat)), verbose = T)

plot(rabinerJuangStepPattern(2,"c",T))

