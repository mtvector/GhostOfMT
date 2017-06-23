allAmby002Samplenames <- c(GSE67118samplenames, GSE35255samplenames, GSE37198samplenames )
allAmby002TP <-  c(GSE67118tp, GSE35255tp,GSE37198tp)

v <- c(GSE67118samplenames, GSE35255samplenames,GSE37198samplenames )
t <- factor(c(GSE67118tp, GSE35255tp,GSE37198tp))
allModel <- model.matrix(~ v*t+0, data= as.data.frame(allAmby002Data))
# imagemat(allModel, main="Model matrix for linear model with one factors")
colnames(allModel) <- make.names(colnames(allModel))
allFit <- lmFit(allAmby002Data, design=allModel)

regenSamplesMT <- c("regeneratinglimb","Innervatedlimb")
nonregenSamplesMT <- c("Aquatic","Terrestrial","Flankwound","Denervatedlimb")
blastemaSamples<- "regeneratinglimb"
nonblastemaSamples <- c( "Aquatic", "Terrestrial","Flankwound", "Denervatedlimb")
cont.matrix = makeContrasts(
  BlastemavsNot =contrastString(addPrefix("v",blastemaSamples), addPrefix("v",nonblastemaSamples)),
  TerrestrialvsAquatic = vTerrestrial-vAquatic,
  RegenvsNotMT = contrastString(addPrefix("v",regenSamplesMT),addPrefix("v",nonregenSamplesMT)),
  InnervatedvsDenervated =  vInnervatedlimb - vDenervatedlimb,
  InnervatedvsWound = vInnervatedlimb - vFlankwound,
  DenervatedvsWound = vDenervatedlimb - vFlankwound,
  InnervatedvsNot = vInnervatedlimb - (vFlankwound + vDenervatedlimb)/2,
  T1vsT0 = t1,
  T3vsT0 = t3,
  T7vsT0 = t7,
  T3vsT1 = t3-t1,
  T7vsT3 = t7-t3,
  levels=allModel)
colnames(cont.matrix) <- c("BlastemavsNot","TerrestrialvsAquatic","RegenvsNotMT","InnervatedvsDenervated","InnervatedvsWound" ,"DenervatedvsWound","InnervatedvsNot","T1vsT0","T3vsT0" ,"T7vsT0","T3vsT1","T7vsT3")
fit2 <- contrasts.fit(allFit,cont.matrix)
fit3 <- eBayes(fit2)

all.DE.genes = lapply(seq(1,ncol(fit3)) , function(i){
  table <-  topTable(fit3[,i], number=2000,sort.by = "logFC")
  idx <- table$logFC > 2 & table$adj.P.Val<.05
  up.genes = AllProbeDict[rownames(table[idx,]),]$Gene
  idx <- table$logFC < -2 & table$adj.P.Val<.05
  down.genes = AllProbeDict[rownames(table[idx,]),]$Gene
  list(upGenes = up.genes, downGenes=down.genes)
})
names(all.DE.genes) <- colnames(cont.matrix)
#mean(length(intersect(t,v)) / length(union(t,v)))

GSE67118Model <- model.matrix(~ GSE67118tp+0, data= as.data.frame(exprsGSE67118))
colnames(GSE67118Model) <- make.names(colnames(GSE67118Model))
GSE67118Fit <- lmFit(exprsGSE67118, design=GSE67118Model)
PreBud <- unique(addPrefix("GSE67118tp",GSE67118tp[as.numeric(GSE67118tp)<10 ]))
EarlyBud <- unique(addPrefix("GSE67118tp",GSE67118tp[as.numeric(GSE67118tp)<=15 & as.numeric(GSE67118tp)>=10]))
MidBud <- unique(addPrefix("GSE67118tp",GSE67118tp[as.numeric(GSE67118tp)<21 & as.numeric(GSE67118tp)>=16]))
LateBud <- unique(addPrefix("GSE67118tp",GSE67118tp[as.numeric(GSE67118tp)<23 & as.numeric(GSE67118tp)>=22]))
Pallet <- unique(addPrefix("GSE67118tp",GSE67118tp[as.numeric(GSE67118tp)<28 & as.numeric(GSE67118tp)>=24]))
GSE67118.cont.matrix = makeContrasts(
  T1to0 =  contrastString( addPrefix("GSE67118tp", c("0.5", "1")) , addPrefix("GSE67118tp", "0")),
  PreBudto0 <-  contrastString(PreBud,addPrefix("GSE67118tp", "0")),
  EarlyBudtoPreBud = contrastString(EarlyBud,PreBud),
  MidBudtoPreBud = contrastString(MidBud,EarlyBud),
  LateBudtoMidBud = contrastString(LateBud, MidBud),
  PallettoLateBud = contrastString(Pallet,LateBud),
  PallettoPreBud = contrastString(Pallet,PreBud),
  levels=GSE67118Model)
colnames(GSE67118.cont.matrix) <- c("T1to0" , "PreBudto0", "EarlyBudtoPreBud" ,"MidBudtoEarlyBud","LateBudtoMidBud","PallettoLateBud","PallettoPreBud")
GSE67118Fit2 <- contrasts.fit(GSE67118Fit,GSE67118.cont.matrix)
GSE67118Fit3 <- eBayes(GSE67118Fit2)

GSE67118.DE.genes = lapply(seq(1,ncol(GSE67118Fit3)) , function(i){
  table <-  topTable(GSE67118Fit3[,i], number=2000,sort.by = "logFC")
  idx <- table$logFC > 2 & table$adj.P.Val<.05
  up.genes = AllProbeDict[rownames(table[idx,]),]$Gene
  idx <- table$logFC < -2 & table$adj.P.Val<.05
  down.genes = AllProbeDict[rownames(table[idx,]),]$Gene
  list(upGenes = up.genes, downGenes=down.genes)
})
names(GSE67118.DE.genes) <- colnames(GSE67118.cont.matrix)

v <- GSE36451samplenames
v = gsub("-", "",v)
t <- GSE36451tp
GSE36451.contrast.names <- sapply(seq(1,length(v)), function(i){
  paste0("v", v[i],".","t",t[i])
})
GSE36451Model <- model.matrix(~ v*t+0, data= as.data.frame(log(normalizeBetweenArrays(exp(exprsGSE36451)),2)))
# imagemat(GSE36451Model, main="Model matrix for linear model with one factors")
colnames(GSE36451Model) <- make.names(colnames(GSE36451Model))
GSE36451Fit <- lmFit(exprsGSE36451, design=GSE36451Model)
colnames(GSE36451Model)
GSE36451.cont.matrix = makeContrasts(
  BlastemaVsAll = vregenerating - (vwoundhealing + vdeveloping + vmature)/3,
  RegenerationVsWound= vregenerating - vwoundhealing,
  levels=GSE36451Model)
#colnames(GSE36451.cont.matrix) <- c(")
GSE36451Fit2 <- contrasts.fit(GSE36451Fit,GSE36451.cont.matrix)
GSE36451Fit3 <- eBayes(GSE36451Fit2)
volcanoplot(GSE36451Fit3[,2])

GSE36451DE.genes = lapply(seq(1,ncol(GSE36451Fit3)) , function(i){
  table <-  topTable(GSE36451Fit3[,i], number=2000,sort.by = "logFC")
  idx <- table$logFC > 2 & table$adj.P.Val<.05
  up.genes = AllProbeDict[rownames(table[idx,]),]$Gene
  idx <- table$logFC < -2 & table$adj.P.Val<.05
  down.genes = AllProbeDict[rownames(table[idx,]),]$Gene
  list(upGenes = up.genes, downGenes=down.genes)
})
names(GSE36451DE.genes) <- colnames(GSE36451.cont.matrix)

v <- GSE35255samplenames
t <- GSE35255tp
GSE35255.contrast.names <- sapply(seq(1,length(v)), function(i){
  paste0("v", v[i],".","t",t[i])
})
GSE35255Model <- model.matrix(~ v*t+0, data= as.data.frame(allAmby002Data))
# imagemat(GSE35255Model, main="Model matrix for linear model with one factors")
colnames(GSE35255Model) <- make.names(colnames(GSE35255Model))
GSE35255Fit <- lmFit(exprsGSE35255, design=GSE35255Model)


v <- GSE37198samplenames
t <- GSE37198tp
GSE37198.contrast.names <- sapply(seq(1,length(v)), function(i){
  paste0("v", v[i],".","t",t[i])
})
GSE37198Model <- model.matrix(~ GSE37198.contrast.names+0, data= as.data.frame(allAmby002Data))
colnames(GSE37198Model) <- make.names(colnames(GSE37198Model))
GSE37198Fit <- lmFit(exprsGSE37198, design=GSE37198Model)

stageNumH <- c("1","2","3","4","5","6","7","8","9","10","11","12","14","16","19","24","40")

DEEmbryo <- lapply(seq(2, length(stageNumH)), function(i){
  U= read.table( file.path(datapath,"DE", "Up", paste0("Stage_",stageNumH[i] , "_Vs_",stageNumH[i-1],".Up.txt")), header = TRUE, sep = "\t")$Gene.ID
  D= read.table( file.path(datapath,"DE", "Down", paste0("Stage_", stageNumH[i] , "_Vs_",stageNumH[i-1],".Down.txt")), header = TRUE, sep = "\t")$Gene.ID
  list(upGenes= U, downGenes=D)
  })
names(DEEmbryo) <- lapply( 2:length(stageNumH), function(i){
  paste0("EmbryoStage",stageNumH[i],"vs",stageNumH[i-1])
})

DE.Universe <-  list( all.DE.genes, GSE67118.DE.genes,GSE36451DE.genes, DEEmbryo)
DE.Universe <-  unlist(DE.Universe,recursive=F)
