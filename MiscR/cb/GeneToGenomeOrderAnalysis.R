library(biomaRt)
multiCarne <- read.csv2("~/Desktop/multispeciesCarnegie.txt",sep="\t")
colnames(multiCarne) <- c("Species","Day",9:23)
rownames(multiCarne) <- multiCarne[,1]
multiCarne <- multiCarne[,-c(1,2)]
multiCarne <- as.matrix(multiCarne)
storage.mode(multiCarne) <- "numeric"
multiCarne <- rbind(multiCarne,"Human"=carnegie.equivalents$human[9:23])
sub1 <- multiCarne[,1]
sub1[is.na(sub1)] <- 0
ggplot(melt(multiCarne-sub1,varnames = c("species","stage"),value.name = "day"))+
  geom_point(aes(x=stage,y=day,color=species))+
  geom_line(aes(x=stage,y=day,color=species))
ggplot(melt(multiCarne,varnames = c("species","stage"),value.name = "day"))+
  geom_point(aes(x=stage,y=day,color=species))+
  geom_line(aes(x=stage,y=day,color=species))

pth <- "~/code/data/tfGenomes/"



listAttributes(mart)[grepl( "id",listAttributes(mart, what="description"))]

#ta <- getBM( attributes = c("ensembl_gene_id","entrezgene"),filters = "ensembl_gene_id", values=getBM(attributes = "ensembl_gene_id",values = "*",mart = mart ),mart=mart)
#gene_exon, transcript_exon,transcript_exon_intron, gene_exon_intron, cdna, coding,coding_transcript_flank,coding_gene_flank,transcript_flank,gene_flank,peptide, 3utr or 5utr
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
for(g in HumanTFs){
  if(!file.exists(paste0(pth,"/","hs/",g,".txt"))){
  print(g)
  useq <- NULL
  dseq <- NULL
  try_again(100,useq <- getSequence(id=g, type="hgnc_symbol", seqType=c("gene_flank","chromosome_name","start_position","end_position","strand"),upstream=10000, mart = mart))
  try_again(100,dseq <- getSequence(id=g, type="hgnc_symbol", seqType="gene_flank",downstream = 10000, mart = mart))
  print(useq)
  catseq <- paste0(useq[,1],dseq[,1])
  print(length(catseq))
  if(length(catseq)>=1){
    write.table(catseq,paste0(pth,"/","hs/",g,".txt"),quote = F,row.names = F,col.names = F)
  }
  }
}

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="useast.ensembl.org")
for(g in MouseTFs){
  if(!file.exists(paste0(pth,"/","mm/",g,".txt"))){
  print(g)
  useq <- NULL
  dseq <- NULL
  try_again(100,useq <- getSequence(id=g, type="mgi_symbol", seqType="gene_flank",upstream=10000, mart = mart))
  try_again(100,dseq <- getSequence(id=g, type="mgi_symbol", seqType="gene_flank",downstream = 10000, mart = mart))
  catseq <- paste0(useq[,1],dseq[,1])
  if(length(catseq)>=1){
    write.table(catseq,paste0(pth,"/","mm/",g,".txt"),quote = F,row.names = F,col.names = F)
  }
  }
}


mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "cjacchus_gene_ensembl")
for(g in union(MouseTFs, HumanTFs)){
  if(!file.exists(paste0(pth,"/","cj/",g,".txt"))){
  print(g)
  useq <- NULL
  dseq <- NULL
  try_again(100,useq <- getSequence(id=g, type="hgnc_symbol", seqType="gene_flank",upstream=10000, mart = mart))
  try_again(100,dseq <- getSequence(id=g, type="hgnc_symbol", seqType="gene_flank",downstream = 10000, mart = mart))
  catseq <- paste0(useq[,1],dseq[,1])
  if(length(catseq)>=1){
    write.table(catseq,paste0(pth,"/","cj/",g,".txt"),quote = F,row.names = F,col.names = F)
  }
  }
}
install.packages("~/data/general/phastCons60way.UCSC.mm10",repos = NULL,type = "source")
biocLite("phastCons60way.UCSC.mm10")
library(phastCons100way.UCSC.hg19)
?phastCons100way.UCSC.hg19
PhastConsDb()
phastCons100way.UCSC.hg19::phastCons100way.UCSC.hg19
scores(phastCons100way.UCSC.hg19,GRanges(seqnames="chr7", IRanges(start=117232380, width=5)))



getBM(attributes=c('gene_flank','start_position','end_position','chromosome_name','strand','ensembl_gene_id'),filters=c('upstream_flank'),values=list(ENSG='BRCA1', Upstream=1000), mart=mart, checkFilters=FALSE)

try_again(100,getSequence(id="BRCA1", type="mgi_symbol", seqType="gene_flank",upstream=10000, mart = mart,verbose = T))


useq = getSequence(id="BRCA1", type="hgnc_symbol", seqType="peptide", mart = mart)


library(data.table)
cmd <- "bedtools closest -a FileA -b FileB"
df <- fread(cmd)