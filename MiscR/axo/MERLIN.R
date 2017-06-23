source("~/code/axo/LoadDatasets.R")


#~/lib/Merlin/merlin -d ~/code/data/axo/ExpressionSets/AllAmby002.txt -o ~/code/data/axo/AllAmby002Net -l ~/code/data/axo/ExpressionSets/RegulatorProbes.txt &
setwd("~/code/data/axo/ExpressionSets")
grepl("exprs",dir()) 

~/lib/Merlin/merlin -d ~/code/data/axo/ExpressionSet/AllAmby002.txt -o ~/code/data/axo/AllAmby002Net -l ~/code/data/axo/ExpressionSet/RegulatorProbes.txt

system(cmdstring, wait=F)
wound <- cbind(exprsGSE35255[,GSE35255samplenames=="Aquatic"],exprsGSE37198[,GSE37198samplenames=="Flankwound"])
regulatorProbes

woundAmby002Net <-  hLICORN( numericalExpression =  cbind(exprsGSE35255[,GSE35255samplenames=="Aquatic"],exprsGSE37198[,GSE37198samplenames=="Flankwound"]),TFlist=regulatorProbes,parallel="multicore", verbose=T)
woundAmby002Netinfluence <- regulatorInfluence(woundAmby002Net,cbind(exprsGSE35255[,GSE35255samplenames=="Aquatic"],exprsGSE37198[,GSE37198samplenames=="Flankwound"]))
write.table(coregnetToDataframe(woundAmby002Net),file = file.path(datapath, "woundAmby002Net.txt") )
save.image(file = file.path(datapath, "hLICORN_networks_WBA.Rdata"))

blastemaAmby002Net <-  hLICORN( numericalExpression =  exprsGSE67118[,GSE67118samplenames!="T0"],TFlist=regulatorProbes,parallel="multicore", verbose=T)
blastemaAmby002Netinfluence <- regulatorInfluence(blastemaAmby002Net,exprsGSE67118[,GSE67118samplenames!="T0"])
write.table(coregnetToDataframe(blastemaAmby002Net),file = file.path(datapath, "blastemaAmby002Net.txt") )
save.image(file = file.path(datapath, "hLICORN_networks_WBA.Rdata"))

allAmby002Net <-  hLICORN( numericalExpression =  cbind(exprsGSE67118,exprsGSE35255,exprsGSE37198),TFlist=regulatorProbes,parallel="multicore", verbose=T)
allAmby002Netinfluence <- regulatorInfluence(allAmby002Net,cbind(exprsGSE67118,exprsGSE35255,exprsGSE37198))
#display(allAmby002Net,cbind(exprsGSE67118,exprsGSE35255,exprsGSE37198),allAmby002Netinfluence)
write.table(coregnetToDataframe(allAmby002Net),file = file.path(datapath, "allAmby002Net.txt") )
save.image(file = file.path(datapath, "hLICORN_networks_WBA.Rdata"))

blastemaRNAseqNet <-  hLICORN( numericalExpression = NormECB ,TFlist=intersect(fullRegulatorList,rownames(NormECB) ),parallel="multicore", verbose=T)
allAmby002Netinfluence <- regulatorInfluence(blastemaRNAseqNet,NormECB)
write.table(coregnetToDataframe(blastemaRNAseqNet),file = file.path(datapath, "blastemaRNAseqNet.txt") )
save.image(file = file.path(datapath, "hLICORN_networks_WBA.Rdata"))