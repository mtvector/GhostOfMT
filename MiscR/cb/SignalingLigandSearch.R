m1 <- round.log(datasets[[25]])#[,-79]
m2 <- round.log(datasets[[27]])
tfz <- union(HumanTFs,MouseTFs)
m1 <- rn.compare(m1,m2)[[1]]
m2 <- rn.compare(m1,m2)[[2]]

m1 <- m1[sort(rownames(m1)),]
m2 <- m2[sort(rownames(m2)),]

mean(rowSds(m1))
hist(rowMaxs(m2)-rowMins(m2),main = "Mouse NBNog GeneLog2FC")

library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")

offspring <- c("GO:0048384","GO:0035329","GO:0005109","GO:0005112","GO:0005160","GO:0005113","GO:0005118","GO:0005117","GO:0005122", "GO:0005104")
gene.data.h <- getBM(attributes=c('hgnc_symbol'),
                   filters = 'go', values = offspring, mart = mart)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org")
gene.data.m <- getBM(attributes=c('mgi_symbol'),
                   filters = 'go', values = offspring, mart = mart)

signalingMolecules <- union(y=toupper(gene.data.m[,1]),x=gene.data.h[,1])

signalingReceptors <-  c("SMO","PTCH1","PTCH2",paste0(rep("FZD",10),1:10),"NOTCH1","NOTCH2",
                        "NOTCH3","NOTCH4","TGFBR1","TGFBR2","TGFBR3",
                        "FGFR1","FGFR2","FGFR3","FGFR4","BMPR1A","BMPR1B","BMPR2",
                        "RARA","RARB","RARG","FLT1","FLT4")

defns <-  getBM(attributes=c("hgnc_symbol","go_id",'name_1006'),
      filters = 'go', values = offspring, mart = mart)

unique(defns[defns[,2]%in% offspring,3])




m <- log(m2+1,2)
#for thresholding
#m[m<=3] <- 0
#m[m>3] <- 1

m <- m[rowMaxs(m)>3,]
m <- m[order(unlist(apply(m,1,which.max))),]
m <- t(apply(m,1,range01Max))
gs1 <- tfz[tfz%in%rownames(m1)]#transcription factors
gs2 <- tfz[tfz%in%rownames(m2)]#transcription factors
m <- m[rownames(m)%in%gs1,]

cutoff <- 3
tfm1 <- m1[rowMaxs(m1)>3&rownames(m1)%in%gs1,]
inds1 <-rowMaxs(tfm1)<cutoff
tfm1 <- tfm1[order(unlist(apply(tfm1,1,which.max))),]
tfm1 <- t(apply(tfm1,1,range01Max))
tfm1[inds1,] <- rep(NA,ncol(tfm1))
tfm2 <- m2[rowMaxs(m2)>3&rownames(m2)%in%gs2,]
inds2 <- rowMaxs(tfm2)<cutoff
tfm2 <- tfm2[order(unlist(apply(tfm2,1,which.max))),]
tfm2 <- t(apply(tfm2,1,range01Max))
tfm2[inds2,] <- rep(NA,ncol(tfm2))
pdf(paste0("~/Desktop/MouseTFMouseOrderedCut",cutoff,".pdf"),width = 12,height = 12)
std.heatmap(tfm2[rownames(wmsame[order(wmsame[,2]),]),],main=paste("Mouse Endoderm\n Ordered by Mouse max\nTFs\ncutoff=",2^cutoff)) 
dev.off()
pdf(paste0("~/Desktop/HumanTFHumanOrderedCut",cutoff,".pdf"),width = 12,height = 12)
std.heatmap(tfm1[rownames(wmsame[order(wmsame[,1]),]),],main=paste("Human Endoderm\n Ordered by Human max\nTFs\ncutoff=",2^cutoff)) 
dev.off()


wm1 <- cbind(rownames(tfm1), apply(tfm1,1,which.max))
wm2 <- cbind(rownames(tfm2), apply(tfm2,1,which.max))
#write.table(wm1[,2], "~/Desktop/HumanWhichmax.txt",quote = F,col.names = F)
cor(as.numeric(wm1[intersect(rownames(wm1),rownames(wm2)),2]),as.numeric(wm2[intersect(rownames(wm1),rownames(wm2)),2]),method = "spearman")
wmsame <- cbind("human"=as.numeric(wm1[intersect(rownames(wm1),rownames(wm2)),2]),"mouse"=as.numeric(wm2[intersect(rownames(wm1),rownames(wm2)),2]))
rownames(wmsame) <- intersect(rownames(wm1),rownames(wm2))
wmsame <- wmsame[order(wmsame[,2]),]
write.table(cbind("human"=wm1[intersect(rownames(wm1),rownames(wm2)),2],"mouse"=wm2[intersect(rownames(wm1),rownames(wm2)),2]),file = "~/Desktop/WhichmaxSamegenes.txt",quote = F,col.names = T)
hist(wmsame[,1]-wmsame[,2],main="HumanMaxTP-MouseMaxTP")


#m <- (m-rowMeans(m))/rowSds(m)
#gs <- signalingMolecules[signalingMolecules%in%rownames(m)]
#gs <- signalingReceptors[signalingReceptors%in%rownames(m)]#signaling molecules
#std.heatmap(m,main="Mouse Endoderm\n Ordered by max\nTFs")






#looking for useful database of 
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_101.jdk/Contents/Home/jre/lib/server/libjvm.dylib')#path to your jvm dylib
require(rJava)
biocLite("rBiopaxParser")
library("rBiopaxParser")
biocLite("paxtoolsr")
library("paxtoolsr")
#browseVignettes("paxtoolsr")
sif <- toSif(system.file("extdata", "biopax3-short-metabolic-pathway.owl", package = "paxtoolsr"))

file = downloadBiopaxData("NCI", c( "pid", "biocarta","reactome","kegg","biogrid"))
biopax = readBiopax(file)

pw_list = listInstances(biopax, class="pathway")
pw_complete = selectInstances(biopax, class="pathway")
pwid1 = "pid_p_500240_FGFR2_ligand_binding_and_activation"
pwid2 = "pid_p_200228_tgfbrpathway"
getInstanceProperty(biopax, pwid1, property="NAME")
getInstanceProperty(biopax, pwid2, property="NAME")
pw_1 = selectInstances(biopax, class="pathway", id=pwid1)
pw_1_component_list = listPathwayComponents(biopax,pwid1)
pw_1_components = selectInstances(biopax,id=pw_1_component_list$id)
pw_2 = selectInstances(biopax, class="pathway", id=pwid2)
pw_2_component_list = listPathwayComponents(biopax,pwid2)
pw_2_components = selectInstances(biopax,id=pw_2_component_list$id)

pw_1_adj = pathway2AdjacancyMatrix(biopax, pwid1, expandSubpathways=TRUE,
                                      splitComplexMolecules=TRUE, verbose=TRUE)
pw_1_graph = pathway2RegulatoryGraph(biopax, pwid1,
                                        splitComplexMolecules=TRUE, verbose=TRUE)
pw_2_adj = pathway2AdjacancyMatrix(biopax, pwid2, expandSubpathways=TRUE,
                                      splitComplexMolecules=TRUE, verbose=TRUE)
pw_2_graph = pathway2RegulatoryGraph(biopax, pwid2,
                                        splitComplexMolecules=TRUE, verbose=TRUE)
pw_1_graph_laidout = layoutRegulatoryGraph(pw_1_graph)
pw_2_graph_laidout = layoutRegulatoryGraph(pw_2_graph)
plotRegulatoryGraph(pw_1_graph)
pdf("~/Desktop/wntGraph.pdf",width = 20,height = 20)
plotRegulatoryGraph(pw_2_graph)
dev.off()

stuff <- downloadPc2()
