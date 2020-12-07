tableOfDat <- read.table(file = "/mnt/DATA4/kevhu/LorenaData/20171130UT_e2e_Combined.Uniq.txt", sep = "\t", header = TRUE)
tableOfDat <- tableOfDat[,c(1:6,which(grepl("total_cov",colnames(tableOfDat))))]
tableOfDat$Gene <- gsub("\\;.*","", tableOfDat$attributes)
tableOfDat$Gene <- gsub("GENE_ID=","", tableOfDat$Gene)

which(grepl("total_cov",colnames(tableOfDat)))

edgeRmat <- tableOfDat[,grepl("total_cov",colnames(tableOfDat))]
rownames(edgeRmat) <- tableOfDat$Gene

myFactors <- data.frame("preCan" = preCan2,
                        "CTBNN1" = CTBNN1.2)

myfilt = filtered.data(edgeRmat, factor = myFactors$preCan, norm = FALSE,
                       depth = NULL, method = 1, cpm = 1, p.adj = "fdr")


###Subset factors so I'm only comparing a set of techinical replicates to another - proof of concept for NoiSeq
preCan <- c("Pre","Can","Can","Can","Pre","Can","Can","Can","Can","Can","Pre","Pre")
CTNBB1 <- c("WT","Mut","WT","WT","Mut","Mut","WT","Mut","Mut","Mut","WT","WT")




preCan2 <- NULL
CTBNN1.2 <- NULL
for(i in seq_along(preCan)){
  preCan2 <- c(preCan2, rep(preCan[i],2))
  CTBNN1.2 <- c(CTBNN1.2, rep(CTNBB1[i],2))
}

###subsetting and randomly sampling of one samp of precursor and one samp of Can 

Pre.subset <- myfilt[,which(preCan2 == "Pre")]
Can.subset <- myfilt[,which(preCan2 == "Can")]
preCan.subsetGroup <- NULL
for(i in 1:(length(unique(Pre.subset))/2)){
  preCan.subsetGroup <- c(preCan.subsetGroup, rep(i,2))
}

preCan.subsetGroup.Can <- NULL
for(i in 1:(length(unique(Can.subset))/2)){
  preCan.subsetGroup.Can <- c(preCan.subsetGroup.Can, rep(i,2))
}

###iteratively go through pairwise comparisons
myFactors <- data.frame("preCan" = c("pre","pre","can","can"))
listOfGenes <- NULL

for(i in unique(preCan.subsetGroup)){
  print(i)
  a <- Pre.subset[,which(preCan.subsetGroup == unique(preCan.subsetGroup)[i])]
  ###grabbing and combining the subsets 1 of each
  for(j in unique(preCan.subsetGroup.Can)){
    b <- Can.subset[,which(preCan.subsetGroup.Can == unique(preCan.subsetGroup.Can)[j])]
    print(colnames(a))
    print(colnames(b))
    dummyDat <- cbind(a,b)
    myData <- readData(data = dummyDat, factors = myFactors)
    noiseq.obj <- noiseq(myData, k = 0.5, norm = "tmm", factor = "preCan",
                         replicates = "technical", conditions = c("pre","can"))
    mynoiseq.deg = degenes(noiseq.obj, q = 0.95, M = NULL)
    listOfGenes <- c(listOfGenes, rownames(mynoiseq.deg))
  }
}


saveRDS(listOfGenes, "/mnt/DATA4/kevhu/lorenaNoiSeqTestPairwiseQ.95")

#######below was just testing out original pipeline


myFactors <- data.frame("preCan" = preCan2,
                        "CTBNN1" = CTBNN1.2)

myfilt = filtered.data(edgeRmat, factor = myFactors$preCan, norm = FALSE,
                       depth = NULL, method = 1, cv.cutoff = 100, cpm = 5, p.adj = "fdr")

myData <- readData(data = edgeRmat, factors = myFactors)

mycountsbio = dat(myData, factor = NULL, type = "countsbio")
explo.plot(mycountsbio,toplot = 1, samples = NULL, plottype = "barplot")

noiseq.obj <- noiseq(myData, k = 0.5, norm = "tmm", factor = "preCan",
                     replicates = "technical", conditions = c("Pre","Can"))


mynoiseq.deg = degenes(noiseq.obj, q = 0.8, M = NULL)


DE.plot(noiseq.obj, q = 0.95, graphic = "expr", log.scale = TRUE)
DE.plot(noiseq.obj, q = 0.95, graphic = "MD")
