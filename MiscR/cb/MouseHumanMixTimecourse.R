source("http://bioconductor.org/biocLite.R")
library(ggplot2)
library(Biobase)
library(RColorBrewer)
library(gplots)
library(plyr)
#biocLite("EACI")
library(EACI)
library(clusterProfiler)
library(matrixStats)
library(zoo)
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
#biocLite("maSigPro")
library(maSigPro)
library(edgeR)
library(EBSeq)
library(DESeq2)
library(devtools)
library(biomaRt)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(scales)
library(gridExtra)
#devtools::install_github("greenelab/TDM")
library(TDM)
library(scales)
library(reshape2)
source("~/code/cb/LoadFunctions.R")
source("~/code/cb/LoadData.R")

rownames(datasets[[11]])[grepl("MLL",rownames(datasets[[11]]))]
datasets[[11]]["MLL",]
MtoH.Orthologs("MLL2",T)
"PAX6"%in%MHOrthos$HSGENE 
MHOrthos$HSGENE[grepl("KAT",MHOrthos$HSGENE )]


library(org.Hs.eg.db)
# set up your query genes

epigenes <- c("EP300","KAT2A","KAT2B","SUZ12","EZH1","EZH2","EED","RBBP4","JARID2","PRC1","DNMT1","DNMT3A","DNMT3B","RING1A","RING1B")#,"MLL","MLL2","MLL3","MLL4","ASH2L","SMYD3","SETBP1","SETD1A","SETD1B","SETD2","SETD3","SETD4","SETD5","SETD6","SETD7","SETD8","SETDB1","SETDB2")
#epigenes <- c("EP300","EZH1","EZH2","EED","RBBP4","PRC1","DNMT1","DNMT3A","DNMT3B")
queryGeneNames <-epigenes
# use sql to get alias table and gene_info table (contains the symbols)
# first open the database connection
dbCon <- org.Hs.eg_dbconn()
# write your SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
aliasSymbol <- aliasSymbol[-31349,]
aliasSymbol <- aliasSymbol[-26170,]
aliasSymbol <- aliasSymbol[-7105,]
# subset to get your results
result <- sapply(queryGeneNames,function(g)aliasSymbol[which(aliasSymbol[,2] == g),5] )
epigenes <- result[epigenes]


for(x in 1:length(datasets)){
  print(x)
  for(n in names(result)){
    rownames(datasets[[x]])[which(rownames(datasets[[x]])==n)] <- result[n]
  }
}

for(x in 1:length(MixSetsECMat)){
  print(x)
  for(n in names(result)){
    rownames(MixSetsECMat[[x]])[which(rownames(MixSetsECMat[[x]])==n)] <- result[n]
  }
}

a=27
b=28
aset=datasets[[a]]
bset=datasets[[b]]
c <- rn.merge.normalize(aset,bset)
aset <- c[[1]]
bset <- c[[2]]

aset <- sweep(aset,2,colSums(aset),"/")*1E6
bset <- sweep(bset,2,colSums(bset),"/")*1E6

#aset=MixSetsECMat[[1]]
#bset=MixSetsECMat[[4]]
#c <- rn.merge.normalize(MixSetsECMat[[1]],MixSetsECMat[[4]])
#aset <- c[[1]]
#bset <- c[[2]]
#aset <- t(apply(aset,1,function(q)aggregate(q,by=list(tpDatasets[[a]]),mean)[,2]))
#bset <- t(apply(bset,1,function(q)aggregate(q,by=list(tpDatasets[[b]]),mean)[,2]))

MATAKAT <- t(Reduce(cbind, lapply(epigenes, function(g)cbind(log(aset[g,]+1,2),log(bset[g,]+1,2) ) )))
rownames(MATAKAT) <- paste0(rep(c("hs","mm"),length(epigenes)),rep(epigenes,each=2))
colnames(MATAKAT) <- unique(tpDatasets[[a]])
std.heatmap(MATAKAT, main="Log2 ECs\nNormalized Between Datasets")

usegenes <- c(epigenes)
MATAKAT <- t(Reduce(cbind, lapply(usegenes, function(g)log(aset[g,]+1,2)-log(bset[g,]+1,2) ) ))
rownames(MATAKAT) <- usegenes
std.heatmap(MATAKAT,main="Log2 Fold Change \nTPMs\nhs/mm")


cutoff = .001
setwd("~/code/data/cb/shiftFiles/pengOutputs")
peng.files <- dir()
peng.files <- peng.files[grepl(".txt",peng.files)]
peng.outputs <- lapply(peng.files, read.table,sep="\t",header=T)
names(peng.outputs) <- gsub(".txt","",peng.files)
tera <-  toupper(as.character(peng.outputs[[3]]$Gene))[peng.outputs[[3]]$Prob_Slow_Paired<= cutoff]
teraC <-  toupper(as.character(peng.outputs[[4]]$Gene))[peng.outputs[[4]]$Prob_Fast_Paired<= cutoff]
tera <- tera[tera%in%neural_related_genes]
teraC <- teraC[teraC%in%neural_related_genes]
library(VennDiagram)
venn.diagram(list(tera,teraC),filename = "~/code/venn.tiff")
venn(list(Teratoma=tera,InVitro= teraC))
tab <-  go.Fisher.Test(intersect(tera,teraC))
write.table(intersect(tera,teraC),"~/code/data/cb/shiftFiles/pengOutputs/TeraandInVitroFaster.txt",sep = "\t",quote = F,row.names = F)

lapply(names(peng.outputs),function(p){
  #tab <- go.Fisher.Test(toupper(as.character(peng.outputs[[p]]$Gene))[peng.outputs[[p]]$Prob_Fast_Paired< 1E-7])
  #if(!is.null(dim(tab))){
    write.table(toupper(as.character(peng.outputs[[p]]$Gene))[peng.outputs[[p]]$Prob_Slow_Paired< cutoff],file = paste0("~/code/data/cb/shiftFiles/pengOutputs","/enrich/",p,"SlowGenes.tsv"),sep = "\t",row.names = F ,quote = F )
    write.table(toupper(as.character(peng.outputs[[p]]$Gene))[peng.outputs[[p]]$Prob_Fast_Paired< cutoff],file = paste0("~/code/data/cb/shiftFiles/pengOutputs","/enrich/",p,"FastGenes.tsv"),sep = "\t",row.names = F ,quote = F )
    
    #write.table(tab,file = paste0("~/code/data/cb/shiftFiles/pengOutputs","/enrich/",p,"FasterEnrichment.tsv"),sep = "\t",row.names = F ,quote = F )
  #}
  #tab <- go.Fisher.Test(toupper(as.character(peng.outputs[[p]]$Gene))[peng.outputs[[p]]$Prob_Slow_Paired< 1E-7])
  #if(!is.null(dim(tab))){
    
    #write.table(tab,file = paste0("~/code/data/cb/shiftFiles/pengOutputs","/enrich/",p,"SlowerEnrichment.tsv"),sep = "\t",row.names = F ,quote = F )
  #}
})

theme_set(theme_bw(base_size = 6))
condit.combos <- expand.grid(c("mm","hs"),unique(conditMixSets[[2]]),stringsAsFactors = F)
condit.combos <- condit.combos[-c(2,7),]
condit.combos <- condit.combos[c(1,6,2,3,4,5),]
#condit.combos <- condit.combos[c(1,3,6),]

##Base dataset analysis
barplot2(lapply(MixSetsEC,colSums)[[1]],main = "Read Numbers MM",las=2,cex.names = .3)
barplot2(lapply(MixSetsEC,colSums)[[2]],main = "Read Numbers HS",las=2,cex.names = .3)

abline(h=1E6)
plot(hclust(as.dist(1-(cor(MixSets[[2]][,conditMixSets[["hs"]]=="100"],method = "spearman")+1)/2)),cex=.5)
plot(hclust(dist(cor(MixSets[[2]][,conditMixSets[["hs"]]=="100"]))),cex=.5)
heatmap.2(cor(MixSets[[2]][,conditMixSets[[2]]!= "0"],method = "spearman") ,Rowv = F,Colv = F,trace = "none",col=cols)


#Get a list of all the genes related to neural processes
#biocLite("GO.db")
library(biomaRt)
library(GO.db)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
t = read.table("~/Downloads/nervousdev.tsv",header = T,sep = "\t")
t$Symbol
#nervous system development "GO:0007399"
offspring = getAllPartOfBP("GO:0007399")
offspringDev = GOBPOFFSPRING[["GO:0007399"]]
#neuron part "GO:0097458"
offspringNeu = GOCCOFFSPRING[["GO:0097458"]]
#synapse "GO:0045202"
offspringSyn = GOCCOFFSPRING[["GO:0045202"]]
offspring = Reduce( union,list(offspringNeu, offspringDev,offspringSyn, "GO:0007399","GO:0097458","GO:0045202"))

gene.data <- getBM(attributes=c('hgnc_symbol'),
                   filters = 'go_id', values = offspring, mart = mart)

neural_list <- gene.data$hgnc_symbol
neural_list[!neural_list %in% gene.data$hgnc_symbol]
# neural_list[neural_list %in% t$Symbol]
# sort(setdiff(unique(t$Symbol) , gene.data$hgnc_symbol))
#write.table(data.frame("Gene"= neural_list), file = "~/Desktop/neuralrelatedgenes.csv",sep = ",",row.names = F,col.names = T)

#### Get premade gene sets
neural_related_genes <- make.names(toupper(as.character(read.csv("~/code/data/cb/markers/neuralrelatedgenes.csv",header = T)[,1])))
neural_related_genes[!neural_related_genes%in% rownames(MixSets$hs)]

#### Heatmaps for Correlation analysis

x=t(apply(datasets[[1]],1,rescale))
y=t(apply(datasets[[2]],1,rescale))
heatmap.2(cor.compare(x,y,min = 0,method ="pearson"),Rowv = F,Colv = F,trace="none")
heatmap.2(cor.compare(datasets[[1]],datasets[[2]],min = 3,method ="pearson"),Rowv = F,Colv = F,trace="none")

apply(combn(1:length(datasets),2),2, FUN= function(x){
  c.c <- cor.compare(datasets[[x[1]]], datasets[[x[2]]],method ="spearman",min = 3)
  rownames(c.c) <- tpDatasets[[x[1]]]
  colnames(c.c) <- tpDatasets[[x[2]]]
  print(paste(datasetNames[[x[1]]],datasetNames[[x[2]]]))
  write.pdf(heatmap.2(c.c, trace = "none", Colv = F,Rowv = F,cexRow = 1,cexCol = 1,srtCol=45,col=cols,xlab = "days",ylab = "days",main = paste( datasetNames[x[1]],"(Rows) vs",datasetNames[x[2]],"(Columns)" ,"\nSpearman Correlation")),filename = paste0("~/Desktop/AllGenesCorrelations/",datasetNames[x[1]],"vs",datasetNames[x[2]],".pdf"))
})

# apply(combn(1:length(datasets),2),2, FUN= function(x){
#   c.c <- cor.compare(datasets[[x[1]]] , datasets[[x[2]]],interest.set = neural_related_genes,method ="pearson")
#   rownames(c.c) <- tpDatasets[[x[1]]]
#   colnames(c.c) <- tpDatasets[[x[2]]]
#   write.pdf(heatmap.2(c.c, trace = "none", Colv = F,Rowv = F,cexRow = 1,cexCol = 1,srtCol=45,col=cols,xlab = "days",ylab = "days",main = paste( datasetNames[x[1]],"(Rows) vs",datasetNames[x[2]],"(Columns)" ,"\nPearson Correlation")),filename = paste0("~/Desktop/NeuralGenesCorrelations/",datasetNames[x[1]],"vs",datasetNames[x[2]],".pdf"))
# })

sapply(1:length(datasets), FUN= function(x){
  c.c <- cor.compare(datasets[[x]] , datasets[[x]],method ="spearman")
  rownames(c.c) <- tpDatasets[[x]]
  colnames(c.c) <- tpDatasets[[x]]
  write.pdf(heatmap.2(c.c, trace = "none", Colv = F,Rowv = F,cexRow = 1,cexCol = 1,srtCol=45,col=cols,xlab = "days",ylab = "days",main = paste( datasetNames[x],"(Rows) vs",datasetNames[x],"(Columns)" ,"\nSpearman Correlation")),filename = paste0("~/Desktop/AllGenesCorrelations/",datasetNames[x],"vs",datasetNames[x],".pdf"))
})


sapply(1:length(datasets), FUN= function(x){
  c.c <- cor.compare(datasets[[x]] , datasets[[x]], interest.set = neural_related_genes,method ="pearson")
  rownames(c.c) <- tpDatasets[[x]]
  colnames(c.c) <- tpDatasets[[x]]
  write.pdf(heatmap.2(c.c, trace = "none", Colv = F,Rowv = F,cexRow = 1,cexCol = 1,srtCol=45,col=cols,xlab = "days",ylab = "days",main = paste( datasetNames[x],"(Rows) vs",datasetNames[x],"(Columns)" ,"\nSpearman Correlation")),filename = paste0("~/Desktop/NeuralGenesCorrelations/",datasetNames[x],"vs",datasetNames[x],".pdf"))
})

#apply(combn(1:length(datasets),2),2, FUN= function(x){
for (i in 1:length(datasets)) {
  for(j in 1:length(datasets)){
    x=c(i,j)
    #c.c <- cor.compare(datasets[[x[1]]] , datasets[[x[2]]],interest.set = cc_list,method ="pearson")
    #rownames(c.c) <- tpDatasets[[x[1]]]
    #colnames(c.c) <- tpDatasets[[x[2]]]
    a <- datasets[[x[1]]]
    b <- datasets[[x[2]]]
    tpa <- tpDatasets[[x[1]]]
    tpb <- tpDatasets[[x[2]]]
    if(grepl("MouseTeraControl", datasetNames[x[1]])){
      a <- a[,-ncol(a)]
      tpa <- tpa[-length(tpa)]
    }
    if(grepl("MouseTeraControl", datasetNames[x[2]])){
      b <- b[,-ncol(b)]
      tpb <- tpb[-length(tpb)]
    }
    ggsave( filename = paste0("~/Desktop/NeuralGenesCorrelations/",datasetNames[x[1]],"vs",datasetNames[x[2]],".pdf"),scaleBoxHM(a,b,tpa,tpb,n1=datasetNames[x[1]],n2=datasetNames[x[2]],g_list = unlist(neural_related_genes)),device = 'pdf')
    #write.pdf(heatmap.2(c.c, trace = "none", Colv = F,Rowv = F,cexRow = 1,cexCol = 1,srtCol=45,col=cols,xlab = "days",ylab = "days",main = paste( datasetNames[x[1]],"(Rows) vs",datasetNames[x[2]],"(Columns)" ,"\nPearson Correlation")),filename = paste0("~/Desktop/CCGenesCorrelations/",datasetNames[x[1]],"vs",datasetNames[x[2]],".pdf"))

  }
}

in.interval <- function(x, interval){
  stopifnot(length(interval) == 2L)
  interval[1] <= x & x <= interval[2]
}

scaleBoxHM <- function(x,y,tpx,tpy,n1,n2,g_list=rownames(x),min = 0, method = "pearson"){
  o.tpx <- tpx
  o.tpy <- tpy
  namex= n1
  namey=n2
  ### FOR SCALED BOXES
  #get the correlation matrix
  cc=cor.compare(x, y,min=min,interest.set = g_list,method=method)
  #Calculate what the width of the tiles should be
  blocksize.vector.x <- diff(tpx)*(1/max(tpx))
  blocksize.vector.x <- c(blocksize.vector.x[1],blocksize.vector.x,blocksize.vector.x[length(blocksize.vector.x)])
  b.v.x <- sapply(2:length(blocksize.vector.x),function(i) blocksize.vector.x[i-1]/2 + blocksize.vector.x[i]/2 )
  blocksize.vector.y <- diff(tpy)*(1/max(tpy))
  blocksize.vector.y <- c(blocksize.vector.y[1],blocksize.vector.y,blocksize.vector.y[length(blocksize.vector.y)])
  b.v.y <- sapply(2:length(blocksize.vector.y),function(i) blocksize.vector.y[i-1]/2 + blocksize.vector.y[i]/2 )
  
  b.v.x <- b.v.x * (max(tpx)*sum(b.v.x)/sum(b.v.x))
  b.v.y <- b.v.y * (max(tpy)*sum(b.v.y)/sum(b.v.y))
  #Figure out where to center the tiles
  tpx.adj <- tpx
  tpy.adj <- tpy
  tpx.n <<- tpx.adj
  tpx.adj[2:(length(tpx.adj))] <- as.vector(sapply(2:(length(tpx.adj)),function(i){
    q <- tpx.n[i-1]+ .5*b.v.x[i-1] + .5*b.v.x[i]
    tpx.n[i] <<- q
    q
  } ))
  tpy.n <<- tpy.adj
  tpy.adj[2:(length(tpy.adj))] <- as.vector(sapply(2:(length(tpy.adj)),function(i){
    q <- tpy.n[i-1]+ .5*b.v.y[i-1] + .5*b.v.y[i]
    tpy.n[i] <<- q
    q
  } ))
  #make list of tile sizes for matrix
  widths=as.vector(t(sapply(b.v.x, function (x) rep(x,length(tpy)))))
  heights=as.vector(sapply(b.v.y, function (x) rep(x,length(tpx))))
  
  rownames(cc) <- tpx.adj
  colnames(cc) <- tpy.adj
  d <- melt(cc)
  g3=ggplot(d, aes(Var1, Var2)) + 
    geom_tile(aes(Var1, Var2,fill = value,width=widths,height = heights)) + 
    coord_fixed(ratio=1)+
    scale_fill_gradientn(colors = cols,guide = "colourbar",name="Correlation")+
    scale_x_continuous()+
    scale_y_continuous()+
    theme_bw(base_size = 6)+
    labs(x = paste0(namex,"Day"), y = paste0(namey,"Day"), 
         title = paste(namex,namey,method," Correlation"))
  g3
}



library(ggrepel)

ix=1
iy=3
species <- grepl( "Human",datasetNames)
gnzH <- gnzM <- rownames(datasets[[ix]])
if(!species[iy]){gnzM <- MtoH.Orthologs(rownames(datasets[[iy]]))[["mm"]]
gnzH <- MtoH.Orthologs(rownames(datasets[[ix]]))[["hs"]]}
x <- md(datasets[[ix]],gnzH)
y <- md(datasets[[iy]],gnzM)
y <- y[,-ncol(y)]
tpx <- tpDatasets[[ix]]
tpy <- tpDatasets[[iy]]
tpy <- tpy[-length(tpy)]

#x <- t(apply(x, 1, function(r)approx(tpx,r, seq(min(tpx),max(tpx)))$y))
#x <- t(apply(x, 1, function(r)sapply(1:length(seq(min(tpx),max(tpx))),function(i)ifelse(seq(min(tpx),max(tpx))[i] %in% tpx,r[which(tpx==seq(min(tpx),max(tpx))[i])],NA) )))
#y <- t(apply(y, 1, function(r)approx(tpy,r, seq(min(tpy),max(tpy)))$y))
#y <- t(apply(y, 1, function(r)sapply(1:length(seq(min(tpy),max(tpy))),function(i)ifelse(seq(min(tpy),max(tpy))[i] %in% tpy,r[which(tpy==seq(min(tpy),max(tpy))[i])],NA) )))

o.tpx <- tpx
o.tpy <- tpy
#tpx <- seq(min(tpx),max(tpx))
#tpy <- seq(min(tpy),max(tpy))
namex= datasetNames[[ix]]
namey=datasetNames[[iy]]
method = "pearson"
g_list = neural_related_genes
min = 0
#cor.compare(x,x,min = min,interest.set = g_list,method="spearman")-cor.compare(x,x,min = min,interest.set = g_list,method="pearson")

cc=cor.compare(x,x,min = min,interest.set = g_list,method=method)
rownames(cc) <- tpx
colnames(cc) <- tpx
d <- melt(cc)
d$Var1 <- factor(d$Var1,levels =d$Var1 )
d$Var2 <- factor(d$Var2,levels=d$Var2)

carne <- carnegie.equivalents[which(in.interval(carnegie.equivalents$human,range(as.numeric(rownames(cc))+15)) &in.interval(carnegie.equivalents$mouse,range(as.numeric(colnames(cc))+6.5))  ),]
hscale <- approx(as.numeric(rownames(cc))+15, 1:nrow(cc),carne$human)$y
mscale <- hscale

g1=ggplot(d, aes(Var1, Var2)) + 
  geom_tile(aes(Var1, Var2,fill = value)) + 
  coord_fixed(ratio=1)+
  #scale_fill_gradient(low = "black",high = "red" ,guide = "colourbar",name="Correlation",na.value = "white")+
  scale_fill_gradientn(colors = cols,guide = "colourbar",name="Correlation")+
  scale_x_discrete(breaks = o.tpx)+
  scale_y_discrete(breaks = o.tpx)+
  #geom_line(data = carne,aes(hscale,mscale))+
  #geom_text(data = carne,aes(hscale,mscale,label=carnegie.stage),size=6)+
#  geom_abline(intercept=0,slope=1,lty=2)+
  theme_bw(base_size = 6)+
  labs(x = paste0(namex,"Day"), y = paste0(namex,"Mouse Day"), 
       title = paste(namex,method,"Correlation"))
g1
ggsave( filename = paste0("~/Desktop/",namex,"CarnegieOverlay",".pdf" ),g1,device = "pdf",height = 4.5,width = 9,units = "cm" )


cc=cor.compare(y,y,min=min,interest.set = g_list,method=method)
rownames(cc) <- tpy
colnames(cc) <- tpy
d <- melt(cc)
d$Var1 <- factor(d$Var1,levels =d$Var1 )
d$Var2 <- factor(d$Var2,levels=d$Var2)


carne <- carnegie.equivalents[which(in.interval(carnegie.equivalents$mouse,range(as.numeric(colnames(cc))+5.5))  ),]
mscale <- approx(as.numeric(colnames(cc))+6.5, 1:ncol(cc),carne$mouse)$y
hscale <- mscale

g2=ggplot(d, aes(Var1, Var2)) + 
  geom_tile(aes(Var1, Var2,fill = value)) + 
  coord_fixed(ratio=1)+
  #scale_fill_gradient2(low = "blue",high = "red",midpoint = (max(d$value)+min(d$value))/2 ,guide = "colourbar",name="Correlation")+
  scale_fill_gradientn(colors = cols,guide = "colourbar",name="Correlation")+
  scale_x_discrete(breaks = o.tpy)+
  scale_y_discrete(breaks = o.tpy)+
  #geom_line(data = carne,aes(hscale,mscale))+
  #geom_text(data = carne,aes(hscale,mscale,label=carnegie.stage),size=6)+
#  geom_abline(intercept=0,slope=1,lty=2)+
  theme_bw(base_size = 6)+
  labs(x = paste0(namey,"Day"), y = paste0(namey,"Mouse Day"), 
       title = paste(namey,method,"Correlation"))

g2
ggsave( filename = paste0("~/Desktop/",namey,"CarnegieOverlay",".pdf" ),g2,device = "pdf",width = 9,height = 4.5,units = "cm" )


### FOR SCALED BOXES
#get the correlation matrix
cc=cor.compare(x, y,min=min,interest.set = g_list,method=method)
#Calculate what the width of the tiles should be
blocksize.vector.x <- diff(tpx)*(1/max(tpx))
blocksize.vector.x <- c(blocksize.vector.x[1],blocksize.vector.x,blocksize.vector.x[length(blocksize.vector.x)])
b.v.x <- sapply(2:length(blocksize.vector.x),function(i) blocksize.vector.x[i-1]/2 + blocksize.vector.x[i]/2 )
blocksize.vector.y <- diff(tpy)*(1/max(tpy))
blocksize.vector.y <- c(blocksize.vector.y[1],blocksize.vector.y,blocksize.vector.y[length(blocksize.vector.y)])
b.v.y <- sapply(2:length(blocksize.vector.y),function(i) blocksize.vector.y[i-1]/2 + blocksize.vector.y[i]/2 )

b.v.x <- b.v.x * (max(tpx)*sum(b.v.x)/sum(b.v.x))
b.v.y <- b.v.y * (max(tpy)*sum(b.v.y)/sum(b.v.y))
#Figure out where to center the tiles
tpx.adj <- tpx
tpy.adj <- tpy
tpx.n <<- tpx.adj
tpx.adj[2:(length(tpx.adj))] <- as.vector(sapply(2:(length(tpx.adj)),function(i){
  q <- tpx.n[i-1]+ .5*b.v.x[i-1] + .5*b.v.x[i]
  tpx.n[i] <<- q
  q
  } ))
tpy.n <<- tpy.adj
tpy.adj[2:(length(tpy.adj))] <- as.vector(sapply(2:(length(tpy.adj)),function(i){
  q <- tpy.n[i-1]+ .5*b.v.y[i-1] + .5*b.v.y[i]
  tpy.n[i] <<- q
  q
} ))
#make list of tile sizes for matrix
widths=as.vector(t(sapply(b.v.x, function (x) rep(x,length(tpy)))))
heights=as.vector(sapply(b.v.y, function (x) rep(x,length(tpx))))

rownames(cc) <- tpx.adj
colnames(cc) <- tpy.adj
d <- melt(cc)
#load raw carnegie frame
carne <- carnegie.equivalents[which(in.interval(carnegie.equivalents$human,range(as.numeric(tpx)+15)) &in.interval(carnegie.equivalents$mouse,range(as.numeric(tpy)+6.5))  ),]
carne.adj <- carne
carne.adj$human <- carne$human-15
carne.adj$mouse <- carne$mouse-6.5

g3=ggplot(d, aes(Var1, Var2)) + 
  geom_tile(aes(Var1, Var2,fill = value,width=widths,height = heights)) + 
  coord_fixed(ratio=1)+
  scale_fill_gradientn(colors = cols,guide = "colourbar",name="Correlation")+
  scale_x_continuous()+
  scale_y_continuous()+
  geom_point(data = carne.adj,aes(human,mouse,label=carnegie.stage))+
  geom_text_repel(data = carne.adj,aes(human,mouse,label=carnegie.stage),segment.size = 0.000001,segment.color = NA,size=3,nudge_y = 1.6,box.padding = unit(0.02, "lines"))+
  geom_line(data = carne.adj,aes(human,mouse))+
  theme_bw(base_size = 6)+
  labs(x = paste0(namex,"Day"), y = paste0(namey,"Mouse Day"), 
       title = paste(namex,namey,method," Correlation"))
g3

ggsave( filename = paste0("~/Desktop/",namex,namey,"CarnegieOverlayScaleBox",".pdf" ),g3,device = "pdf",width = 9,height = 4.5,units = "cm" )



#####FOR UNSCALED BOXES

cc=cor.compare(x, y,min=min,interest.set = g_list,method=method)
tpx.adj <- tpx
tpy.adj <- tpy
rownames(cc) <- tpx.adj
colnames(cc) <- tpy.adj
d <- melt(cc)
d$Var1 <- factor(d$Var1,levels =d$Var1 )
d$Var2 <- factor(d$Var2,levels=d$Var2)

carne <- carnegie.equivalents[which(in.interval(carnegie.equivalents$human,range(as.numeric(tpx)+14)) &in.interval(carnegie.equivalents$mouse,range(as.numeric(tpy)+6))  ),]
carne.adj <- carne
carne.adj$human <- carne$human-14
carne.adj$mouse <- carne$mouse-6
hscale <- approx(as.numeric(rownames(cc))+14, 1:nrow(cc),carne$human)$y
mscale <- approx(as.numeric(colnames(cc))+6, 1:ncol(cc),carne$mouse)$y
hequiv <- approx(as.numeric(rownames(cc)),1:nrow(cc),as.numeric(colnames(cc)))$y
mequiv <- approx(as.numeric(rownames(cc)),1:nrow(cc),as.numeric(rownames(cc))[1:ncol(cc)])$y
eq <- data.frame(mequiv,hequiv)


g3=ggplot(d, aes(Var1, Var2)) + 
  geom_tile(aes(Var1, Var2,fill = value)) + 
  coord_fixed(ratio=1)+
  #scale_fill_gradient(low = "black",high = "yellow" ,guide = "colourbar",name="Correlation",na.value = "grey70")+
  scale_fill_gradientn(colors = cols,guide = "colourbar",name="Correlation")+
  scale_x_discrete(breaks = o.tpx)+
  scale_y_discrete(breaks = o.tpy)+
  geom_point(data = carne,aes(hscale,mscale,label=carnegie.stage))+
  geom_line(data = carne,aes(hscale,mscale))+
  geom_text_repel(data = carne,aes(hscale,mscale,label=carnegie.stage),segment.size = 0.000001,segment.color = NA,size=3,nudge_y = 1.6,box.padding = unit(0.02, "lines"))+
  geom_line(data=eq,aes(hequiv,mequiv),lty=2)+
  theme_bw(base_size = 6)+
  labs(x = paste0(namex,"Day"), y = paste0(namey,"Mouse Day"), 
       title = paste(namex,namey,method," Correlation"))
g3


for(d in 1:length(datasets)){
  species <- grepl( "Human",datasetNames[d])
  gnzH <- gnzM <- rownames(BrainSpan)
  if(!species){gnzM <- MtoH.Orthologs(rownames(BrainSpan))[["mm"]]
  gnzH <- MtoH.Orthologs(rownames(BrainSpan))[["hs"]]}
  BSS <-  md(BrainSpan[,order(BrainSpanCols$structure_name)],gnzH)
  BSS <-BSS[,grepl( "pcw",colnames(BSS))]
  c.c <- cor.compare(BSS, md(datasets[[d]],gnzM),method ="spearman",min = 0)
  write.table(c.c,file = paste0("~/Desktop/BrainSpanCors/",datasetNames[d],".txt"),sep = "\t")
  write.pdf(heatmap.2(c.c, trace = "none", Colv = F,Rowv = F,cexRow = .25,cexCol = .8,col=cols,main = paste( datasetNames[d],"vs BrainSpanSamples\nSpearman Correlation")),filename = paste0("~/Desktop/BrainSpanCors/",datasetNames[d],".pdf"))
}

cors <- lapply(1:length(datasets),function(d){
  c.c <- cor.compare(BrainSpan[,order(BrainSpanCols$structure_name)], datasets[[d]],interest.set = neural_related_genes,method ="pearson")
  write.pdf(heatmap.2(c.c, trace = "none", Colv = F,Rowv = F,cexRow = .25,cexCol = .8,col=cols,main = paste( datasetNames[d],"vs BrainSpanSamples\nSpearman Correlation")),filename = paste0("~/Desktop/BrainSpanCors/",datasetNames[d],".pdf"))
})

ages <- factor( BrainSpanCols$age)
tu <- factor(BrainSpanCols$structure_name,levels = unique( BrainSpanCols$structure_name[order(BrainSpanCols$structure_name)]))

bs.maxcors <-  lapply(1:length(datasets),function(d){
  c.c <- cor.compare(BrainSpan, datasets[[d]],method ="spearman")
  cor.inds <-  t(sapply(levels(tu),function(t){ 
    in.t <- as.matrix(c.c[tu==t,,drop=FALSE])
    ind = which.max(in.t)
    mx=max(in.t )
    k <- arrayInd(ind, dim(in.t))
    c(tpDatasets[[d]][k[2]] ,rownames(in.t)[k[1]] , mx)
  }))
  cor.inds[order(cor.inds[,3],decreasing = T),]
})
names(bs.maxcors) <- datasetNames
for(x in 1:length(bs.maxcors))colnames(bs.maxcors[[x]]) <- c("Day","BS.sample","cor")


bs.maxcors <-  lapply(1:length(datasets),function(d){
  c.c <- cor.compare(BrainSpan, datasets[[d]],method ="spearman")
  sapply(levels(tu),function(t){ 
    in.t <- as.matrix(c.c[tu==t,,drop=FALSE])
    ind = which.max(in.t)
    mx=max(in.t )
    k <- arrayInd(ind, dim(in.t))
    print(k[[1]])
    heatmap.2(c.c[rownames(in.t)[k[1]],])
    #c(tpDatasets[[d]][k[2]] ,rownames(in.t)[k[1]] , mx)
  })
  #cor.inds[order(cor.inds[,3],decreasing = T),]
})



cors = apply(combn(11:13,2),2, FUN= function(xx){
  x <- datasets[[xx[1]]]
  y <- datasets[[xx[2]]]
  #interest.set=cc_list
  #interest.set=sample(rownames(x),length(cc_list))
  i = intersect(rownames(x), rownames(y))
  i = i[rowMeans(x[i,])>2 | rowMeans(y[i,])>2]
  if(!is.null(interest.set)){
    i = interest.set[interest.set%in%i]
  }
  c.c <- lapply(i,function(j){ 
    a <-  na.trim(na.approx(merge.zoo(zoo(x[j,],tpDatasets[[xx[1]]]) ,zoo(y[j,],tpDatasets[[xx[2]]]))))
    cor(a[,1],a[,2] ,method = "pearson")
  })
  print(c.c)
  #rownames(c.c) <- tpDatasets[[x[1]]]
  #colnames(c.c) <- tpDatasets[[x[2]]]
  #write.pdf(heatmap.2(c.c, trace = "none", Colv = F,Rowv = F,cexRow = 1,cexCol = 1,srtCol=45,col=cols,xlab = "days",ylab = "days",main = paste( datasetNames[x[1]],"(Rows) vs",datasetNames[x[2]],"(Columns)" ,"\nPearson Correlation")),filename = paste0("~/Desktop/CCGenesCorrelations/",datasetNames[x[1]],"vs",datasetNames[x[2]],".pdf"))
})

a=as.matrix(read.table("~/Desktop/BrainSpanCors/NewMouse0.txt",header = T,row.names = NULL,check.names = T,"\t"))
b=as.matrix(read.table("~/Desktop/BrainSpanCors/NewHuman100.txt",header = T,row.names = NULL,check.names = T,"\t"))
c=as.matrix(read.table("~/Desktop/BrainSpanCors/NewHuman10.txt",header = T,row.names = NULL,check.names = T,"\t"))
rownames(a) <- a[,1]
rownames(b) <- b[,1]
rownames(c) <- c[,1]
a <- a[,-1]
b <- b[,-1]
c <- c[,-1]
a <- matrix(as.numeric(a),nrow = nrow(a),ncol = ncol(a),dimnames = dimnames(a))
b <- matrix(as.numeric(b),nrow = nrow(b),ncol = ncol(b),dimnames = dimnames(b))
c <- matrix(as.numeric(c),nrow = nrow(c),ncol = ncol(c),dimnames = dimnames(c))
write.pdf(heatmap.2(rbind(a,b,c), trace = "none", Colv = F,Rowv = F,cexRow = .4,cexCol = 1,srtCol=45,col=cols,xlab = "days",ylab = "BRAINSPAN"),filename = "~/Desktop/MvsHvsMixBrainspan.pdf")

heatmap.2(c, trace = "none", Colv = F,Rowv = F,cexRow = .4,cexCol = 1,srtCol=45,col=cols,xlab = "days",ylab = "BRAINSPAN")

write.table(data.frame("Acronym"= unique(BrainSpanCols$structure_acronym), "Structure"=unique(BrainSpanCols$structure_name)), file = "~/Desktop/StructureKey.csv",sep = ",",row.names = F,col.names = T)


gs=rownames(MixSets$mm)[grepl(pattern = "HOX",rownames(MixSets$mm))]
xx <- as.list(org.Mm.egALIAS2EG)
names(xx) <- toupper(names(xx))
order(sapply(xx[gs],"[",1))


write.pdf({
  heatmap.2(MixSets$mm[gn,conditMixSets[[x[1]]]==x[2] ],trace = "none",Rowv = F,Colv = F,col = cols,srtCol = 45, main = "Mouse Hox Genes \n Log2 TPM")
},filename = "~/Desktop/MouseHOX.pdf" )

xx <- as.list(org.Hs.egALIAS2EG)
names(xx) <- toupper(names(xx))

write.pdf({
  heatmap.2(round.log(MixSets$mm[grepl(pattern = "HOX",rownames(MixSets$mm)) ,conditMixSets$mm != "100"]),trace = "none",Rowv = F,Colv = F,col = cols,srtCol = 45, main = "Human Hox Genes \n Log2 TPM")
},filename = "~/Desktop/MouseHOX.pdf" )

#Generate rescaled Log2 TPMs
for(i in 1:length(datasets)){
  #print(cor(apply(datasets[[i]][1:100,],1,rescale), t(datasets[[i]][1:100,]),method = "spearman"))
  write.table( t(apply(datasets[[i]],1,rescale)) ,file = paste0("~/Desktop/RescaledDatasets/",datasetNames[[i]],".tab"),sep = "\t")
}

#datasets[[3]] <- datasets[[3]][,-ncol(datasets[[3]])]
#write files of rescaled TPM expression
for(i in 1:length(datasets)){
  #print(cor(t(apply(datasets[[i]],1,rescale)), datasets[[i]],method = "spearman"))
  write.table( t(apply(datasets[[i]][neural_related_genes[neural_related_genes %in% rownames(datasets[[i]])],],1,rescale)) ,file = paste0("~/Desktop/RescaledDatasets/","NeuralOnly",datasetNames[[i]],".tab"),sep = "\t")
}

sapply(1:nrow(condit.combos),function(i){
  sapply(1:nrow(condit.combos),function(j){
  x <- unlist(condit.combos[i,])
  y <- unlist(condit.combos[j,])
  print(x)
  print(y)
  #tsx <- apply(MixSets[[x[1]]][,conditMixSets[[x[1]]]==x[2] ],1, aggregate,list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), mean )
  print("x done")
  #tsy <- apply(MixSets[[y[1]]][,conditMixSets[[y[1]]]==y[2] ],1, aggregate,list(tpMixSets[[y[1]]][conditMixSets[[y[1]]]==y[2] ]), mean )
  #tsx <- t(sapply(tsx,function(q)q[,2]))
  #tsy <- t(sapply(tsy,function(q)q[,2]))
  species <- x[1] == y[1]
  if(!species){
    if(x[1]=="hs"){
      gnzM <- MtoH.Orthologs(rownames(tsx))[["mm"]]
      gnzH <- MtoH.Orthologs(rownames(tsx))[["hs"]]
      tsx <-  md(tsx,gnzH)
      tsy <-  md(tsy,gnzM)
    }else{
      gnzM <- MtoH.Orthologs(rownames(tsy),MtoH = T)[["mm"]]
      gnzH <- MtoH.Orthologs(rownames(tsy),MtoH = T)[["hs"]]
      tsx <-  md(tsx,gnzM)
      tsy <-  md(tsy,gnzH)
    }
  }
  tsx <- MixSets[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]
  tsy <- MixSets[[y[1]]][,conditMixSets[[y[1]]]==y[2] ]
  #colnames(tsx) <- unique(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ])
  #colnames(tsy) <- unique(tpMixSets[[y[1]]][conditMixSets[[y[1]]]==y[2] ])
  c.c <- cor.compare(round.log(tsx),round.log(tsy), min = 3,method="spearman")
  write.pdf(heatmap.2(c.c[seq(nrow(c.c),1),], trace = "none", Colv = F,Rowv = F,cexRow = .25,cexCol = .25,col=cols,main = paste( paste0(x[1],x[2])," (Rows) vs",paste0(y[1],y[2]) ," (Columns)\nSpearman Correlation")),filename = paste0("~/Desktop/AllGenesCorrelations/nomean",paste0(x[1],x[2]),"vs",paste0(y[1],y[2]),".pdf"))
  write.table(c.c[seq(nrow(c.c),1),],file = paste0("~/Desktop/AllGenesCorrelations/nomean",paste0(x[1],x[2]),"vs",paste0(y[1],y[2]),".txt"),sep = "\t")
  
  })
})

condit.combos <- expand.grid(c("mm","hs"),unique(conditMixSets[[2]]),stringsAsFactors = F)
condit.combos <- condit.combos[-c(2,7),]
condit.combos <- condit.combos[c(1,4,6),]

fcHM <- apply(sapply(read.csv("~/code/data/cb/markers/FoldChandgeHeatMapLists.csv",header = T),as.character),2,toupper)
gene_list <- fcHM[,4]
gene_list <- c("PAX6", "BRN2","FAT4","SLC17A6","SLC1A3","NR4A2","DDC","RBFOX3","GABBR1","DCX","MAP2","SHH","NKX2.2","NKX6.1","NKX6.2","ISL1","AADC","FOXA1","FOXA2","TH","BLBP","OLIG2","SYT4","SNAP25","TUBB3","NEFL","ELAVL4","RAB6B","STMN2")
plotlist <- lapply( gene_list[gene_list%in% rownames(MixSets$hs) & gene_list%in% rownames(MixSets$mm)],function(gn){
    meanGenes <- apply(condit.combos,1,function(x){
      #ts <- aggregate(GetNormalizedMat( MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ], MedianNorm(MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]))[gn,], list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), mean)
      ts <- aggregate(MixSets[[x[1]]][gn,conditMixSets[[x[1]]]==x[2] ], list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), mean)
      zoo(ts[,2],ts[,1])
    })
    ds <- Reduce(merge.zoo,meanGenes)
    colnames(ds) <- apply(condit.combos,1,function(aaa) paste0(aaa[1],aaa[2]))
    ds <- data.frame(ds)
    ds <- cbind("tme"=rownames(ds) , ds)
    melted <-  melt(ds)
    sdGenes <- apply(condit.combos,1,function(x){
      #ts <- aggregate(GetNormalizedMat( MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ], MedianNorm(MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]))[gn,], list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), function(x)sd(x)/sqrt(length(x)))
      ts <- aggregate(MixSets[[x[1]]][gn,conditMixSets[[x[1]]]==x[2] ], list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), function(x)sd(x)/sqrt(length(x)))
      zoo(ts[,2],ts[,1])
    })
    sdGenes <- Reduce(merge.zoo,sdGenes)
    colnames(sdGenes) <- apply(condit.combos,1,function(aaa) paste0(aaa[1],aaa[2]))
    sdGenes <- as.data.frame(sdGenes)
    sdGenes <- cbind("tme"=rownames(sdGenes) , sdGenes)
    meltSD <-  melt(sdGenes)
    melted <-  cbind(melted,"sd"= meltSD$value)
    #melted$tme <- factor(melted$tme,levels = melted$tme)
    melted$tme <- as.numeric(as.character(melted$tme))
    #Plotting Error Bars
    pd <- position_dodge(0.5) # move them .05 to the left and right
    
    gg <-  ggplot(melted, aes(x=tme, y=value, group=variable,color=variable)) + 
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.03, position=pd) +
      geom_line(position=pd,size=.4) +
      geom_point(position=pd,size=.8) +
      theme_bw(base_size = 6) +
      theme(legend.position="none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
    #scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
      labs(x = "", y = "", 
           title = gn,family="arial",size=6)
    gg
  })
pdf("~/Desktop/ChosenMouseNoscale.pdf", onefile = TRUE)
marrangeGrob(grobs=plotlist , nrow=4, ncol=4,top = NULL)
dev.off()

###RESCALED PLOTS
plotlist <- lapply( gene_list[gene_list%in% rownames(MixSets$hs) & gene_list%in% rownames(MixSets$mm)],function(gn){
  meanGenes <- apply(condit.combos,1,function(x){
    #ts <- aggregate(GetNormalizedMat( MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ], MedianNorm(MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]))[gn,], list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), mean)
    ts <- aggregate(MixSets[[x[1]]][gn,conditMixSets[[x[1]]]==x[2] ], list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), mean)
    zoo(ts[,2],ts[,1])
  })
  ds <- Reduce(merge.zoo,meanGenes)
  colnames(ds) <- apply(condit.combos,1,function(aaa) paste0(aaa[1],aaa[2]))
  ds = apply(ds,2,rescale)
  ds <- data.frame(ds)
  ds <- cbind("tme"=rownames(ds) , ds)
  melted <-  melt(ds)
  pd <- position_dodge(0.5) # move them .05 to the left and right
  #melted$tme <- factor(melted$tme,levels = melted$tme)
  melted$tme <- as.numeric(as.character(melted$tme))
  gg <-  ggplot(melted, aes(x=tme, y=value, group=variable,color=variable)) + 
    geom_line(position=pd,size=1.2) +
    geom_point(position=pd,size=3) +
    theme_bw()+
    theme(legend.position="none",
             panel.grid.major = element_blank()
    ) +
    #scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
    labs(x="",y="", title = gn, family="Arial")
  gg
})

pdf("~/Desktop/ChosenHumanScaled.pdf", onefile = TRUE)
marrangeGrob(grobs=plotlist , nrow=1, ncol=1,top = NULL)
dev.off()

fcHM <- apply(sapply(read.csv("~/code/data/cb/markers/FoldChandgeHeatMapLists.csv",header = T),as.character),2,toupper)
gene_list <- fcHM[,3]
apply(combn(1:nrow(condit.combos),2),2,function(j){
  x <- unlist(condit.combos[j[1],])
  y <- unlist(condit.combos[j[2],])
  GL <-  gene_list[gene_list%in% rownames(MixSets[[1]]) & gene_list%in% rownames(MixSets[[2]])]
  print(x)
  print(y)
  xset <- MixSetsEC[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]
  yset <- MixSetsEC[[y[1]]][,conditMixSets[[y[1]]]==y[2] ]
  xyNorm <- rn.merge(xset,yset)
  xyNorm <- GetNormalizedMat(xyNorm, MedianNorm(xyNorm))
  xsetNorm <- xyNorm[GL,1:ncol(xset)]
  ysetNorm <- xyNorm[GL,(ncol(xset)+1):ncol(xyNorm)]
  tsx <- apply(xsetNorm,1, aggregate,list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]), mean )
  print("x done")
  tsy <- apply(ysetNorm,1, aggregate,list(tpMixSets[[y[1]]][conditMixSets[[y[1]]]==y[2] ]), mean )
  tsx <- t(sapply(tsx,function(q)q[,2]))
  tsy <- t(sapply(tsy,function(q)q[,2]))
  colnames(tsx) <- unique(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ])
  colnames(tsy) <- unique(tpMixSets[[y[1]]][conditMixSets[[y[1]]]==y[2] ])
  if(x[1]=="hs" & (x[2]=="85"| x[2]=="10"))tsx <- cbind(rep(NA, nrow(tsx)),tsx)
  if(x[1]=="mm" & (x[2]=="85"| x[2]=="10"))tsx <- cbind(rep(NA, nrow(tsx)) ,tsx)
  if(y[1]=="hs" & (y[2]=="85"| y[2]=="10"))tsy <- cbind(rep(NA, nrow(tsy)) ,tsy)
  if(y[1]=="mm" & (y[2]=="85"| y[2]=="10"))tsy <- cbind(rep(NA, nrow(tsy)) ,tsy)
  c.c <- sapply(GL ,function(gn){
    log(tsx[gn,]+1,2) - log(tsy[gn,]+1,2)
  })
  write.pdf(heatmap.2(as.matrix(t(c.c)), trace = "none", Colv = F,Rowv = F,cexRow = 1,cexCol = 1,col=cols,main = paste( paste0(x[1],x[2])," over",paste0(y[1],y[2]) ," \n Normalized log2 Fold Change")),filename = paste0("~/Desktop/FoldChangeHM/",paste0(x[1],x[2]),"vs",paste0(y[1],y[2]),".pdf"))
})




ages <- BrainSpanCols$age
tu <- BrainSpanCols$structure_name[order(BrainSpanCols$structure_name)]

#Brainspan correlations for 
bs.maxcors <- lapply(1:nrow(condit.combos),function(x){
  x <- unlist(condit.combos[x,])
  tps <- tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]
  #setMat <-  apply(MixSets[[x[1]]][,conditMixSets[[x[1]]]==x[2] ],1, aggregate,list(tps), mean )
  #setMat <- t(sapply(setMat,function(p)p[,2]))
  #colnames(setMat) <- unique(tps)
  setMat <- MixSets[[x[1]]][,conditMixSets[[x[1]]]==x[2] ]
  species <- "hs" == x[1]
  gnzH <- gnzM <- union(rownames(BrainSpan),rownames(x))
  if(!species){gnzM <- MtoH.Orthologs(rownames(BrainSpan))[["mm"]]
  gnzH <- MtoH.Orthologs(rownames(BrainSpan))[["hs"]]}
  BSS <-  md(BrainSpan[,order(BrainSpanCols$structure_name)],gnzH)
  #BSS <-BSS[,grepl( "pcw",colnames(BSS))]
  c.c <- cor.compare(round.log(BSS), round.log(md(setMat,gnzM)),method ="spearman",min = 3)
  cor.inds <-  t(sapply(unique(tu),function(q){ 
    in.t <- as.matrix(c.c[tu==q,,drop=FALSE])
    ind = which.max(in.t)
    mx=max(in.t )
    k <- arrayInd(ind, dim(in.t))
    c(colnames(in.t)[k[2]],rownames(in.t)[k[1]] , mx)
  }))
  write.table(c.c,file = paste0("~/Desktop/BrainSpanCors/",paste0(x[1],x[2]),".txt"),sep = "\t")
  write.pdf(heatmap.2(c.c, trace = "none", Colv = F,Rowv = F,cexRow = .25,cexCol = .8,col=cols,main = paste( paste0(x[1],x[2]),"vs BrainSpanSamples\nSpearman Correlation")),filename = paste0("~/Desktop/BrainSpanCors/roundlog",paste0(x[1],x[2]),".pdf"))
  return(cor.inds[order(cor.inds[,3],decreasing = T),])
})
names(bs.maxcors) <- apply(condit.combos,1,function(x)paste0(x[1],x[2]))
for(x in 1:length(bs.maxcors))colnames(bs.maxcors[[x]]) <- c("Day","BS.sample","cor")

lapply(names(bs.maxcors),function(n){
  write(n,"~/Desktop/BrainSpanMaxCors.txt",append = T)
  write.table(bs.maxcors[[n]] ,"~/Desktop/BrainSpanMaxCors.txt", append=TRUE,sep = "\t",quote = F)})

struct.cors <- lapply(1:nrow(condit.combos),function(x){
  x <- unlist(condit.combos[x,])
  m <- read.csv2(file = paste0("~/Desktop/BrainSpanCors/",paste0(x[1],x[2]),".txt"),sep = "\t",header = T,row.names = NULL)
  m <- as.matrix(m)
  rownames(m) <- m[,1]
  m <- m[,-1]
  colnames(m) <- gsub("X","",colnames(m))
  m <- matrix(as.numeric(m),nrow = nrow(m),ncol=ncol(m),dimnames = dimnames(m))
})
names(struct.cors) <- paste0(condit.combos[,1],condit.combos[,2])
heatmap.2(rbind(struct.cors$hs100[grepl("AMY",rownames(struct.cors$hs100)),], cbind( struct.cors$hs100[grepl("AMY",rownames(struct.cors$hs100)),1] ,struct.cors$hs10[grepl("AMY",rownames(struct.cors$hs10)),]),struct.cors$mm0[grepl("AMY",rownames(struct.cors$hs10)),] ),Rowv = F,Colv = F,trace = "none",col = cols,density.info = "none")

brainStructure <- "MGE 8 pcw"
MGE.c <- lapply(list("hs100"= struct.cors$hs100[grepl(brainStructure,rownames(struct.cors$hs100)),],"hs10"=struct.cors$hs10[grepl(brainStructure,rownames(struct.cors$hs10)),],"mm0"= struct.cors$mm0[grepl(brainStructure,rownames(struct.cors$mm0)),]),function(x)zoo(x,as.numeric(names(x))))

q <- Reduce(merge.zoo ,MGE.c)
q <- as.data.frame(q)
colnames(q) <- names(MGE.c)
q$tme <- as.numeric(rownames(q))
q <- melt(q,id.vars = 4)

gg <-  ggplot(q, aes(x=tme, y=value, group=variable,color=variable)) + 
  geom_line(size=1.2) +
  geom_point(size=3) +
  theme_bw() +
  #scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  labs(x = "Day", y = "Log2 TPM", 
       title = "Cor VS 8 pcw MGE")
gg


which.max(BrainSpan["FOXA2",])
plot(BrainSpan["TERC",])
plot(MixSets$hs["TERT",])
rn.compare(BrainSpan,MixSets$hs)
#Output spreadsheets of averaged values
apply(condit.combos,1,function(x){
  setMat <-  apply(MixSets[[x[1]]][,conditMixSets[[x[1]]]==x[2] ] ,1, aggregate,list(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ]) , function(x)sd(x)/sqrt(length(x)) )
  setMat <- t(sapply(setMat,function(p)p[,2]))
  colnames(setMat) <- unique(tpMixSets[[x[1]]][conditMixSets[[x[1]]]==x[2] ])
  print(setMat)
  write.table(setMat,file = paste0("~/Desktop/",paste0(x[1],x[2],"SE"),".txt"),sep = "\t")
})

apply(condit.combos,1,function(x){
  for(i in 1:length(datasets)){
    q = read.table(file = paste0("~/Desktop/",paste0(x[1],x[2],"mean"),".txt"),sep = "\t",header = T,row.names = 1)
    q = t(apply(q,1,rescale))
    colnames(q) <- gsub("X","",colnames(q))
    write.table( q ,file = paste0("~/Desktop/RescaledDatasets/",paste0(x[1],x[2],"RescaledMean"),".txt"),sep = "\t")
  }})



gs=rownames(MixSets$hs)[grepl(pattern = "HOX",rownames(MixSets$hs))]
write.pdf({
  mypar(3,3)
  
  for(gn in gs[gs%in% rownames(MixSets$hs)]){
    h10 <- zoo(MixSets$hs[gn,27:51],tpMixSets[27:51])  
    h85 <- zoo(MixSets$hs[gn,52:76],tpMixSets[52:76])
    h100 <- zoo(MixSets$hs[gn,77:102],tpMixSets[77:102]) 
    plot.zoo(merge(h10,h85,h100),plot.type = "single",col = 1:4, lty = 1,main=paste(gn),xlab = "days",ylab = "log2 TPM" )
    legend("topleft", c("hs,10%","hs,85%","hs,100%"), col = 1:3, lty = 1)
  }
}
,filename = "~/Desktop/HumanHOXTimesGraphed.pdf")

gs=rownames(MixSets$mm)[grepl(pattern = "HOX",rownames(MixSets$mm))]
write.pdf({
  mypar(3,3)
  for(gn in gs[gs%in% rownames(MixSets$mm)]){
    h10 <- zoo(MixSets$mm[gn,27:51],tpMixSets[27:51])  
    h85 <- zoo(MixSets$mm[gn,52:76],tpMixSets[52:76])
    h100 <- zoo(MixSets$mm[gn,1:26],tpMixSets[1:26]) 
    plot.zoo(merge(h10,h85,h100),plot.type = "single",col = 1:4, lty = 1,main=paste(gn),xlab = "days",ylab = "log2 TPM" )
    legend("topleft", c("mm,10%","mm,85%","mm,0%"), col = 1:3, lty = 1)
  }
}
,filename = "~/Desktop/MouseHOXTimesGraphed.pdf")


write.pdf({
  mypar(3,3)
  for(gn in cc_list[cc_list %in% rownames(MixSets$hs) & cc_list %in% rownames(MixSets$mm)]){
    m0 <- zoo(MixSets$mm[gn,1:26],tpMixSets[1:26])  
    m10 <- zoo(MixSets$mm[gn,27:51],tpMixSets[27:51])  
    h10 <- zoo(MixSets$hs[gn,27:51],tpMixSets[27:51])  
    h85 <- zoo(MixSets$hs[gn,52:76],tpMixSets[52:76])
    h100 <- zoo(MixSets$hs[gn,77:102],tpMixSets[77:102]) 
    plot.zoo(merge(m0,m10, h10,h85,h100),plot.type = "single",col = 1:6, lty = 1,main=paste(gn),xlab = "days",ylab = "log2 TPM" )
    legend("topleft", c("mm,0%","mm,10%","hs,10%","hs,85%","hs,100%"), col = 1:6, lty = 1)
  }
}
,filename = "~/Desktop/CellCycleGenesCorrelations/HumanMouseCC.pdf")


#Retrieve All Human genes associated with GO term: GO:0007399 nervous system development

getAllPartOfBP <- function(g){
  offspring = GOBPCHILDREN[[g]]
  terms <- c(g,offspring[names(offspring) == "part_of" | names(offspring) == "is_a"])
  q <- terms[names(terms)=="part_of"]
  offspring=NULL
  while(length(q)>0){
    offspring = GOBPCHILDREN[[q[1]]]
    offspring = offspring[!offspring%in% terms]
    offspring = offspring[names(offspring) == "part_of" | names(offspring) == "is_a" ]
    q <- unique(c(q,offspring))
    print(length(q))
    terms <- c(terms,offspring)
    q <- q[-1]
  }
  return(terms)
}

getAllPartOfCC <- function(g){
  offspring = GOCCCHILDREN[[g]]
  terms <- c(g,offspring[names(offspring) == "part_of" | names(offspring) == "is_a"])
  q <- terms[names(terms)=="part_of"]
  offspring=NULL
  while(length(q)>0){
    offspring = GOCCCHILDREN[[q[1]]]
    offspring = offspring[names(offspring) == "part_of" | names(offspring) == "is_a" ]
    terms <- c(terms,offspring)
    q <- q[-1]
  }
  return(terms)
}





#### Make Fancy two color heatmaps

sortneurallist <- toupper(as.character(read.csv("~/code/data/cb/markers/OrderNeuralGeneList.csv",header = F)[,1]))
xx=3
yy=1
x = datasets[[xx]][,-length(datasets[[xx]])]
y = datasets[[yy]]
x = 2^x - 1
y = 2 ^y -1 
tpx=tpDatasets[[xx]][-length(tpDatasets[[xx]])]
tpy= tpDatasets[[yy]]
n1 = datasetNames[[xx]]
n2 =datasetNames[[yy]]

mz= lapply(sortneurallist, function(gn){
  if(gn %in% rownames(x)& gn%in% rownames(y)){
    xi = zoo(x[gn,], tpx)
    yi = zoo(y[gn,], tpy)
    xi=merge(xi,zoo(NA, tpx))[,1]
    yi=merge(yi,zoo(NA, tpy))[,1]
    a=na.approx(merge(xi,yi))
    #merge(a, zoo(NA))
  }
})

mz=mz[!sapply(mz, is.null)] 
#mz= mz[order(sapply(mz,function(x) which(rescale(x[,1]) > .8 )[1]))]
mzn <- sortneurallist

mz=Reduce(merge.zoo, mz )
colnames(mz) <- paste0(c("mm","hs") ,rep(mzn,1,each=2))
mz=apply(mz,2,rescale)
#mz=sweep(mz,MARGIN = 2,colMaxs(mz,na.rm = T),FUN = "/")

mz <- cbind("scale mm"=c(seq(0,1,length.out = 11), rep(NA,nrow(mz)-11 )),mz)
mz <- cbind("scale hs"=c(seq(0,1,length.out = 11), rep(NA,nrow(mz)-11 )),mz)

### FOR SCALED BOXES
#get the correlation matrix
#Calculate what the width of the tiles should be
blocksize.vector.x <- diff(tpx)*(1/max(tpx))
blocksize.vector.x <- c(blocksize.vector.x[1],blocksize.vector.x,blocksize.vector.x[length(blocksize.vector.x)])
b.v.x <- sapply(2:length(blocksize.vector.x),function(i) blocksize.vector.x[i-1]/2 + blocksize.vector.x[i]/2 )
blocksize.vector.y <- diff(tpy)*(1/max(tpy))
blocksize.vector.y <- c(blocksize.vector.y[1],blocksize.vector.y,blocksize.vector.y[length(blocksize.vector.y)])
b.v.y <- sapply(2:length(blocksize.vector.y),function(i) blocksize.vector.y[i-1]/2 + blocksize.vector.y[i]/2 )

b.v.x <- b.v.x * (max(tpx)*sum(b.v.x)/sum(b.v.x))
b.v.y <- b.v.y * (max(tpy)*sum(b.v.y)/sum(b.v.y))

tpx.adj <- tpx
tpy.adj <- tpy
tpx.n <<- tpx.adj
tpx.adj[2:(length(tpx.adj))] <- as.vector(sapply(2:(length(tpx.adj)),function(i){
  q <- tpx.n[i-1]+ .5*b.v.x[i-1] + .5*b.v.x[i]
  tpx.n[i] <<- q
  q
} ))
tpy.n <<- tpy.adj
tpy.adj[2:(length(tpy.adj))] <- as.vector(sapply(2:(length(tpy.adj)),function(i){
  q <- tpy.n[i-1]+ .5*b.v.y[i-1] + .5*b.v.y[i]
  tpy.n[i] <<- q
  q
} ))

rownames(mz) <- tpy.adj

for(i in seq(ncol(mz)-2,1,by = -2)){
  df <- data.frame(rep(NA,nrow(mz)))
  colnames(df) <- paste0(rep(" ",i),collapse = "")
  mz <- cbind(mz[,1:i], df , mz[,(i+1):ncol(mz)] )
}

a = melt(as.matrix(mz),na.rm = F)
#a=a[nrow(a):1,]
a=cbind(a, "group"=as.integer(grepl("mm",a$Var2)))
#a$Var1 <- as.factor(a$Var1)
a$Var2 <- factor(a$Var2, levels=a$Var2[order(a$Var2,decreasing = T)],exclude = NULL)
a$rescaleoffset <- a$value + 100*(a$group)

scalerange <- range(a$value,na.rm = T)
gradientends <- scalerange + rep(c(0,100), each=2)
colorends <- c( "black", "yellow", "black", "yellow")

g=ggplot(a, aes(Var1, Var2)) + 
  ggtitle(paste(n1, "vs", n2,"\n( TPM)"))+
  geom_tile(aes(x = Var1,fill = rescaleoffset,width=rep(b.v.y,nrow(a)/max(length(tpy),length(tpx)) ))) + 
  scale_fill_gradientn(colours = colorends, values = rescale(gradientends),na.value = "grey55",guide = "legend",name="Expression",breaks=c(0,1,100,101),labels=c("mm 0%", "hs 100%", "mm 0%","mm 100%")) + 
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
ggsave(filename = paste0("~/Desktop/OverallFavNeuralHeatmapNoLogSizebox",n1,n2,".pdf" ),plot=g,width = 30,height = 30)


sapply(1:length(inds),function(j){
  write.pdf({
  barplot2(apply(DE.MixSets[[j]],2,function(x)sapply(x,length)),col = c("red","blue"),main = paste("# Genes with PostFC > 4\n",names(inds)[j]))
  legend("topright",legend = c("Up","Down"),fill = c("red","blue"))
  },filename = paste0("~/code/data/cb/shiftFiles/DifferentExpressions/","DEvs0",names(inds)[j],".pdf"))
})

apply(combn(1:length(inds),2),2,function(j){
  write.pdf({
  for(i in c(1,2)){
    mypar(2,2)
    dirx <- ifelse(i==1,"Up","Down")
    a=getDEMatchListMatrix(as.matrix(DE.MixSets[[j[1]]])[i,],as.matrix(DE.MixSets[[j[2]]])[i,])
    a= getDEMatchRefined(a)
    heatmapMatch( a,main = paste("DE",dirx,names(inds)[j]),cluster = F ,yLab = names(inds)[j[1]],xLab = names(inds)[j[2]])
    pvalPlot(hypergeoMatch(a, DE.MixSets[[j[1]]][i,],as.matrix(DE.MixSets[[j[2]]])[i,]),mainT = paste("Fisher's Test Pvals\n",dirx,names(inds)[j]),ylab=names(inds)[j[1]],xlab=names(inds)[j[2]])
  }
  },filename = paste0("~/code/data/cb/shiftFiles/DifferentExpressions/","DEvs0Heatmaps",names(inds)[j[1]],names(inds)[j[2]],".pdf"))
})

"NANOG"%in%DE.MixSets$hs100
sapply(DE.MixSets,function(x) sapply(x, function(y)"NANOG"%in%y) )

xx <- as.list(org.Hs.egALIAS2EG)
names(xx) <- toupper(names(xx))
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="useast.ensembl.org")
GOTermsHS <- getBM(attributes = c('hgnc_symbol', "go_id",
                                  "name_1006"), filters = "entrezgene",
                   values = sapply(xx[rownames(MixSets$hs)],"[",1), mart = mart)

xx <- as.list(org.Mm.egALIAS2EG)
names(xx) <- toupper(names(xx))
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="useast.ensembl.org")
GOTermsMM <- getBM(attributes = c('mgi_symbol', "go_id",
                                  "name_1006"), filters = "entrezgene",
                   values = unlist(sapply(xx[rownames(MixSets$mm)],"[",1)), mart = mart)


GOTermsHS$hgnc_symbol <- toupper(GOTermsHS$hgnc_symbol)
GOTermsMM$mgi_symbol <- toupper(GOTermsMM$mgi_symbol)
colnames(GOTermsMM)[1] <- colnames(GOTermsHS)[1] <- "Symbol"

GOTermsUD <- lapply(1:length(inds),function(i){
  print(i)
  if(grepl("mm",names(inds)[i])){
    dict <- GOTermsMM
  }else{
    dict <- GOTermsHS
  }
  apply(DE.MixSets[[i]],2,function(r){
    resultsUp <- dict$name_1006[dict$Symbol %in% unlist(r[1])]
    resultsDown <- dict$name_1006[dict$Symbol %in% unlist(r[2])]
    return(list("up"= resultsUp,"down"=resultsDown))
  })
})
save("GOTerms",file = "~/Desktop/GOobject.RData")


termz <- c("Neural Tube development","Synapse Organization","Neurological System Process","Axonogenesis","Synaptic Transmission","Transmission of nerve impulse","axon guidance","neuron differentiation","neurogenesis","neuron development","nervous system development","Hindbrain Development","Forebrain Development")
termz <- tolower(termz)
nameTables <-lapply(termz,function(n){
  lapply(1:length(inds),function(i){
    sapply(GOTermsUD[[i]],function(r){
      resultsUp=table(unlist(r[[1]]))[n]
      resultsUp[is.na(resultsUp)] <- 0
      resultsDown=table(unlist(r[[2]]))[n]
      resultsDown[is.na(resultsDown)] <- 0
      return(list("up"= sort(resultsUp,decreasing = T),"down"=sort(resultsDown,decreasing = T)))
    })
  })
})
names(nameTables) <- termz

lapply(names(nameTables),function(n){
  lapply(1:length(inds),function(i){
    d <- data.frame("tme"=as.integer(colnames(nameTables[[n]][[i]])),"down"=unlist(nameTables[[n]][[i]][1,]),"up"=unlist(nameTables[[n]][[i]][2,]),stringsAsFactors = F  )
    d <- rbind(c(0,0,0),d)
    g <- ggplot(d, aes(tme,group = 1)) + 
      geom_point(aes(y = up, colour = "Up"))+
      geom_line(aes(y = up, colour = "Up")) + 
      #geom_point(aes(y = down, colour = "Down"))+
      #geom_line(aes(y = down, colour = "Down"))+
      scale_y_continuous(limits = c(0, max(max(d$up),max(d$down))))+
      scale_color_manual(values=c("red", "blue"))+
      labs(x = "Day (vs Day 0)", y = "# of Genes", 
           title = paste(n,names(inds)[i]))+
      theme_bw()
      ggsave( filename = paste0("~/code/data/cb/shiftFiles/DifferentExpressions/GOfigure/",names(inds)[i],n,".pdf" ),g,device = "pdf")
  })
})

lapply(names(nameTables)[c(11,7,4,9,5)],function(n){
    d <- data.frame("tme"=as.integer(colnames(nameTables[[n]][[i]])),"hs100"=unlist(nameTables[[n]][[1]][2,]),"hs10"=unlist(nameTables[[n]][[3]][2,]),stringsAsFactors = F  )
    d <- rbind(c(0,0,0),d)
    g <- ggplot(d, aes(tme,group = 1)) + 
      geom_point(aes(y = hs100, colour = "hs100"))+
      geom_line(aes(y = hs100, colour = "hs100")) + 
      geom_point(aes(y = hs10, colour = "hs10"))+
      geom_line(aes(y = hs10, colour = "hs10"))+
      scale_y_continuous(limits = c(0, max(max(d$hs100),max(d$hs10))))+
      scale_color_manual(values=c("red", "blue"))+
      labs(x = "Day (vs Day 0)", y = "# of Genes", 
           title = paste(n,"hs10 vs hs100"))+
      theme_bw()
    ggsave( filename = paste0("~/code/data/cb/shiftFiles/DifferentExpressions/GOfigure/","HS100HS10Overlay",n,".pdf" ),g,device = "pdf")
})

hmgroups <- lapply(lapply(read.table("~/code/data/cb/markers/heatmapgroups.txt", sep = "\t",stringsAsFactors = F),toupper),make.names)
hmgroups <- lapply(hmgroups,function(x) x[x!="X"])

write.pdf({
  gns <- c("SHH","NKX2.2","NKX6.1","OLIG2","ELAVL4")
  lapply(gns,function(gn){
    N <- GetNormalizedMat(MixSetsEC$mm[,conditMixSets$mm=="85"],Sizes = MedianNorm(MixSetsEC$mm[,conditMixSets$mm=="85"]))
    mypar(2,2)
    x <- rbind(MixSets$mm[gn,conditMixSets$mm=="85"],MixSetsEC$mm[gn,conditMixSets$mm=="85"],N[gn,],
               colSums(MixSetsEC$mm[,conditMixSets$mm=="85"]),colSums(MixSetsEC$mm[,conditMixSets$mm=="85"]>5))
    plot(x[1,],main=paste(gn,"TPM"))
    
    plot(x[2,],main=paste(gn,"EC"))
    plot(x[3,],main=paste(gn,"NormalizedEC"))
    plot(x[4,],main=paste("TotalReads"))
    plot(x[5,],main=paste("# Genes >5"))
  })
},filename = "~/Desktop/WeirdReads.pdf")

