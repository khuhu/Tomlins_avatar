library(NOISeq)
tableOfDat <- read.table(file = "/mnt/DATA4/kevhu/KellyData/LNCap_Transcriptome_Compiled.countMat.csv", sep = ",", header = TRUE)
tableOfDat <- tableOfDat[,c(1:17)]

tableOfDat$Gene <- gsub("\\;.*","", tableOfDat$attributes)
tableOfDat$Gene <- gsub("GENE_ID=","", tableOfDat$Gene)

tableOfDat <- tableOfDat[-c(20813,20814),]

##filtering out unneeded samples for comparisons
#Hek293 <- grep("HEK293", colnames(tableOfDat))
#tableOfDat <- tableOfDat[,-c(Hek293)]
#tableOfDat <- tableOfDat[,-which(colnames(tableOfDat) == "H358.p53.total_cov")]
#tableOfDat <- tableOfDat[,-grep("PRPK", colnames(tableOfDat))]


countTable <- tableOfDat[,grepl("Totel.e2e",colnames(tableOfDat))]
rownames(countTable) <- tableOfDat$Gene
countTable <- countTable[,order(colnames(countTable))]

###weird ordering because of name for last few



myFactors <- data.frame("groups" = c("p53Hras","p53Hras","p53Hras","p53","p53","p53",
                                     "scrHras","scrHras","scrHras","scr","scr","scr"))

myFiltDat <- filtered.data(countTable, factor = myFactors$groups, norm = FALSE,
                           depth = NULL, method = 1, cpm = 1, p.adj = "fdr")

myData <- readData(data = countTable, factors = myFactors)
mycountsbio = dat(myData, factor = NULL, type = "countsbio")
explo.plot(mycountsbio,toplot = 1, samples = NULL, plottype = "barplot")

myData <- readData(data = myFiltDat, factors = myFactors)
#myTmm <- tmm(assayData(myData)$exprs)
#head(myTmm)



#myfilt = filtered.data(myData, factor = myfactors, norm = FALSE,
#                       depth = NULL, method = 1, cv.cutoff = 100, cpm = 5, p.adj = "fdr")

noiseq.obj <- noiseq(myData, k = 0.5, norm = "tmm", factor = "groups",
                     replicates = "technical", conditions = c("scr","scrHras"))


noiseq.obj2 <- noiseq(myData, k = 0.5, norm = "tmm", factor = "groups",
                     replicates = "technical", conditions = c("p53","p53Hras"))



mynoiseq.deg = degenes(noiseq.obj, q = 0.8, M = NULL)
mynoiseq.deg2 = degenes(noiseq.obj2, q = 0.8, M=NULL)


DE.plot(noiseq.obj, q = 0.8, graphic = "expr", log.scale = TRUE)
DE.plot(noiseq.obj, q = 0.8, graphic = "MD")


countTable2 <- countTable
countTable2$gene <- rownames(countTable2)


mynoiseq.deg$odds_ratio <- (mynoiseq.deg$prob * 5)

write.table(mynoiseq.deg,"/mnt/DATA4/kevhu/KellyData/scrHRasNoiSeq.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(mynoiseq.deg2,"/mnt/DATA4/kevhu/KellyData/p53HRasNoiSeq.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)


###Redoing so I leave out bad library, maybe why it's so highly differentially expressed

countTable2 <- countTable[,-c(10)]

myFactors.reduced <- data.frame("groups" = c("p53Hras","p53Hras","p53Hras","p53","p53","p53",
                                     "scrHras","scrHras","scrHras","scr","scr"))

myFiltDat <- filtered.data(countTable2, factor = myFactors.reduced$groups, norm = FALSE,
                           depth = NULL, method = 1, cpm = 1, p.adj = "fdr")

myData <- readData(data = countTable2, factors = myFactors.reduced)
mycountsbio = dat(myData, factor = NULL, type = "countsbio")
explo.plot(mycountsbio,toplot = 1, samples = NULL, plottype = "barplot")

myData <- readData(data = myFiltDat, factors = myFactors.reduced)

noiseq.obj.red <- noiseq(myData, k = 0.5, norm = "tmm", factor = "groups",
                     replicates = "technical", conditions = c("scr","scrHras"))


noiseq.obj.red.2 <- noiseq(myData, k = 0.5, norm = "tmm", factor = "groups",
                      replicates = "technical", conditions = c("p53","p53Hras"))



mynoiseq.deg.red = degenes(noiseq.obj.red, q = 0.8, M = NULL)
mynoiseq.deg2.red = degenes(noiseq.obj.red.2, q = 0.8, M=NULL)


###no it's even higher ... in this case by 2 thousand ... well I guess we conclude to keep the bad library

DE.plot(noiseq.obj, q = 0.8, graphic = "expr", log.scale = TRUE)
DE.plot(noiseq.obj, q = 0.8, graphic = "MD")
