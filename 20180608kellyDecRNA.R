#tableOfDat <- read.table(file = "/mnt/DATA4/kevhu/LorenaData/20171103LorenaRNAseq.csv", sep = ",", header = TRUE, skip = 7)
library(edgeR)
library(ggplot2)
library(ggrepel)
tableOfDat <- read.table(file = "/mnt/DATA4/kevhu/linuxRNApull/KellyRNADecDat.txt", sep = "\t", header = TRUE)

tableOfDat <- tableOfDat[,c(2,4:ncol(tableOfDat))]
tableOfDat <- tableOfDat[-1,]
edgeRmat <- tableOfDat[,2:ncol(tableOfDat)]
rownames(edgeRmat) <- tableOfDat$Gene

#keep <- rowSums(cpm(edgeRmat) > 5) >= 10
#edgeRmat <- edgeRmat[keep, ]

###edit this for so that it keeps 2 and 0. i. concordance between the two technical replicates - bit too conservative
#for(i in 1:(ncol(edgeRmat)/2)){
#  rem <- rowSums(cpm(edgeRmat[,c(2*i-1,2*i)]) > 3) == 2 | rowSums(cpm(edgeRmat[,c(2*i-1,2*i)]) > 3) == 0
#  print(length(rem))
#  edgeRmat <- edgeRmat[rem,]
#}


dummyMat <- NULL
edgeRmat.cpm <- cpm(edgeRmat)
for(i in 1:nrow(edgeRmat)){
  dummyRow <- NULL
  for(j in 1:(ncol(edgeRmat)/2)){
    if((sum(edgeRmat.cpm[i,c(2*j-1,2*j)] > 5) == 2 | sum(edgeRmat.cpm[i,c(2*j-1,2*j)] > 5) == 0) == TRUE){
      k <- 1
    }
    else{
      k <- 0
    }
    dummyRow <- c(dummyRow, k)
  }
  dummyMat <- rbind(dummyMat, dummyRow)
}

rownames(dummyMat) <- rownames(edgeRmat)
keepFirstFil <- which(rowSums(dummyMat) >= quantile(rowSums(dummyMat), .05))




hist(rowSums(dummyMat))
quantile(rowSums(dummyMat), 0.05)
summary(rowSums(dummyMat))

###filter using - how I separate the techinical replicates is a makeshift method and probably need to think harder about it
#listOfCors <- NULL
#set1 <- edgeRmat[,grep("L2", colnames(edgeRmat))]
#set2 <- edgeRmat[,-grep("L2", colnames(edgeRmat))]
#for(i in 1:nrow(edgeRmat)){
#  a <- cor(c(unlist(set1[i,])), c(unlist(set2[i,])), method = c("pearson"))
#  listOfCors <- c(listOfCors, a)
#}
#listOfCors.na <- which(is.na(listOfCors))
#listOfCors2 <- listOfCors[-which(is.na(listOfCors))]

#looking at cors by samp

listOfCors <- NULL
set1 <- edgeRmat[,grep("L2", colnames(edgeRmat))]
set2 <- edgeRmat[,-grep("L2", colnames(edgeRmat))]
for(i in 1:ncol(set1)){
  a <- cor(c(unlist(set1[,i])), c(unlist(set2[,i])), method = c("pearson"))
  listOfCors <- c(listOfCors, a)
}

summary(listOfCors)




#quantile(listOfCors2, 0.25)
#keepFirstFil <- which(listOfCors > quantile(listOfCors2, 0.25))

#edgeRmat <- edgeRmat[which(listOfCors > quantile(listOfCors, 0.25))



###cutoff is 1 sd away and emprically testing when sampling errors for things with cpm < 5
#order of cutoffs is pearson, spearman, and kendall
#edgeRmat <- edgeRmat[which(listOfCors > 0.6697699),]
#edgeRmat <- edgeRmat[which(listOfCors > 0.6457861),]
#edgeRmat <- edgeRmat[which(listOfCors > 0.5123475),]
#edgeRmat <- edgeRmat[which(listOfCors > quantile(listOfCors, 0.1)),]



###testing new filter
#testDf <- edgeRmat
#for(i in 1:(ncol(testDf)/2)){
#  rem <- rowSums(cpm(testDf[,c(2*i-1,2*i)]) > 5) == 2 | rowSums(cpm(testDf[,c(2*i-1,2*i)]) > 5) == 0
#  print(length(rem))
#  testDf <- testDf[rem,]
#}


###Test to see how filtering works
#test <- edgeRmat[10:15,]
#test.keep <- cpm(test) > 5
#rowSums(test.keep) >=2


###needed to change code b/c of triplicates
colnames(edgeRmat)[order(colnames(edgeRmat))]

edgeRmat <- edgeRmat[,order(colnames(edgeRmat))]
edgeRmat.pooled <- NULL
new.names <- NULL
###pooling technical replicates, because Marioni et at 2008 states that techincial variance i.e library replicates follows a Poisson distribution
indexes <- c(3,6,9,12,14,16,18)
for(i in indexes){
  pooled.name <- colnames(edgeRmat)[i]
  new.names <- c(new.names, pooled.name)
  if(i < 14)
  {
    pooled <- rowSums(edgeRmat[,c(i-2,i-1, i)])
    edgeRmat.pooled <- cbind(edgeRmat.pooled, pooled)
  }
  if(i > 12){
    pooled <- rowSums(edgeRmat[,c(i-1, i)])
    edgeRmat.pooled <- cbind(edgeRmat.pooled, pooled) 
  }
}

colnames(edgeRmat.pooled) <- new.names
colnames(edgeRmat.pooled)[1:4] <- c("LNCaP sh-scr +Lacz","LNCaP sh-scr +HRas-Q61R","LNCaP sh-p53 +Lacz","LNCaP sh-p53 +HRas-Q61R")
#rownames(edgeRmat.pooled) <- tableOfDat$Gene



#library(edgeR)

dgeObj <- DGEList(counts = edgeRmat.pooled)
dgeObj$samples
bioExpressFilt <- c(unlist(which(rowSums(cpm(dgeObj) > 5) >= 2)))

#keep <- intersect(bioExpressFilt, corrFilter)
keep <- intersect(bioExpressFilt, keepFirstFil)

dgeObj <- dgeObj[keep, , keep.lib.sizes = FALSE]

##checking what the actual count cutoff is with a cpm of 5
#which(dgeObj$samples$lib.size == min(dgeObj$samples$lib.size))



#can conclude the count cutoff for the smallest library is 50 counts with cpm 5

#keep <- rowSums(cpm(dgeObj) > 5) >= 10
#dgeObj <- dgeObj[keep, , keep.lib.sizes = FALSE]
dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples

### diff groups


groups <- factor(c("1","2","3","4","5","6","7"))
###preCan

design <- model.matrix(~0 + groups, data=dgeObj$samples)
design
#dgeObj <- estimateDisp(dgeObj, design)
dgeObj <- estimateCommonDisp(dgeObj, design)
fit <- glmFit(dgeObj, design)
lrt <- glmLRT(fit, contrast=c(0,0,-1,1,0,0,0))
topTags(lrt)
lrt.table <- lrt$table
lrtPvals <- lrt$table$PValue
lrtPVals.ad <- p.adjust(lrtPvals, method = c("BH"))
lrt$table$QVal <- lrtPVals.ad
lrt.tab <- lrt$table
lrt.tab$Gene <- rownames(lrt.tab)
lrt.tab$color <- NULL


###B/c none of Kelly's genes are significantly expressed, I'll highlight the top 30
#for(i in seq_along(lrt.tab$logFC)){
#  if(abs(lrt.tab$logFC[i]) > 1 & lrt.tab$QVal[i] < 0.05){
#    lrt.tab$color[i]  <- "red"
#  }
#  else{
#    lrt.tab$color[i] <- "blue"
#  }
#}


###below is what I used
topGeneNames <- rownames(topTags(lrt, n=30))
nameIdx <- NULL
for(i in seq_along(topGeneNames)){
   a <- which(rownames(lrt.table) == topGeneNames[i])
   nameIdx <- c(nameIdx, a)
}

lrt.tab$color <- "blue"
lrt.tab$color[nameIdx] <- "red"








plot2 <- ggplot(data = lrt.tab,aes(logFC, -log10(QVal)))
plot2 <- plot2 +  geom_point(pch= 20, color = lrt.tab$color, size = 1.25) + geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.25)
plot2 <-  plot2 + geom_text_repel(data = lrt.tab[which(lrt.tab$color == "red"),], aes(label = Gene, vjust =-1.0), size = 3, alpha = 1.0) + ggtitle("Precursor vs Cancer")
plot2 <-  plot2 + ggtitle("Cell lines 3 vs 4")


plot2 <- plot2 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot2


setwd("/mnt/DATA4/kevhu/KellyData/")
pdf(file = "20180611.3v4.pdf",onefile = TRUE, width = 3)
plot2
dev.off()

#setwd("/mnt/DATA4/kevhu/LorenaData/")
#png(file = "20180108.PreCancerConcordance.png", width = 3.5, height = 5,units="in", res = 200)
#plot2
#dev.off()

###write.table
sig.genes <- lrt.tab[which(lrt.tab$color == "red"),]
if(nrow(sig.genes) == 0){
  sig.genes <- topTags(lrt)
}
gNames <- rownames(sig.genes)
sig.genes <- cbind(gNames,sig.genes)
sig.genes <- sig.genes[,c(1,2,3,5,6)]
colnames(sig.genes)[1] <- "Genes"
write.table(sig.genes,"20180611.3v4.txt", quote = FALSE, row.names = FALSE, sep = "\t")


### 20190116 
###using data above to cluster Kelly's samples. nuRNA + specific gene + KEGG pathway + pol III complex
library(org.Hs.eg.db)
library(clusterProfiler)

kegg <- org.Hs.egPATH2EG
mapped <- mappedkeys(kegg)
kegg2 <- as.list(kegg[mapped])

tRnaKegg <- unlist(kegg2[which(names(kegg2) == "00970")])
tRnaKegg2 <- bitr(tRnaKegg, fromType = "ENTREZID", toType = "SYMBOL", annoDb = "org.Hs.eg.db")[,2]

nuRNAlist <- read.table("/mnt/DATA4/kevhu/KellyData/nuRNAgeneList.txt", sep = ",", stringsAsFactors = FALSE)
nuRNAlist <- unlist(nuRNAlist)
nuRNAlist <- str_remove_all(nuRNAlist, "\\(.*")

allPossibleTrnaSubunits <-  read.table("/mnt/DATA4/kevhu/KellyData/tRNAgeneList.txt", sep = ",", stringsAsFactors = FALSE)
allPossibleTrnaSubunits <- unlist(allPossibleTrnaSubunits)

allGenes <- c(tRnaKegg2, nuRNAlist, allPossibleTrnaSubunits, "POLR3GL")
allGenes <- unname(allGenes)
allGenes <- str_remove_all(allGenes," ")

###https://www.genome.jp/dbget-bin/www_bget?ko:K15202
allGenes <- c(allGenes, "GTF3C1","GTF3C2","GTF3C3","GTF3C4","GTF3C5","GTF3C6")

geneMatrix.cpm <- cpm(dgeObj)
geneMatrix.cpm <- geneMatrix.cpm[which(rownames(geneMatrix.cpm) %in% allGenes),]
geneMatrix.cpm <- log2(geneMatrix.cpm + 1)
### I feel like there were more samples ... yeah wrong file ... look at the R file in Kelly's folder
###