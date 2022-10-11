library(MASS)

data(galaxies)
X = galaxies / 1000
library(mclust, quietly=TRUE)


fit = Mclust(X, G=4, model="V")
summary(fit)

# testing this out on gene level estimates first .. easier
geneTable <- read.table("/mnt/DATA5/tmp/kev/newMouse2/cnMatrix_gene.txt", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)

geneTable[,2:ncol(geneTable)] <- log2(geneTable[,2:ncol(geneTable)])


geneManipRegions <- geneTable[grep("Del", geneTable$Gene), ]
geneManipReg_melt <- melt(geneManipRegions)

# cutoff of -1 assumes purity > 50%
geneManipReg_melt <- geneManipReg_melt[which(geneManipReg_melt$value < -1),]

# the question is how do I want to design the input for this 
# i.e flexilibity for amplicon or gene call input?

purityEst <- NULL
for (i in unique(geneManipReg_melt$variable)) {
  tmpDf <- geneManipReg_melt[which(geneManipReg_melt$variable == i),]
  if (nrow(tmpDf) < 1) {
    tmpVector <- c("Sample" = i, "Purity" = 0)
    purityEst <- rbind(purityEst, tmpVector)
  } else{
    tmpPurity <- 1 - round(2^mean(tmpDf$value), digits = 3)
    tmpVector <- c("Sample" = i, "Purity" = tmpPurity)
    purityEst <- rbind(purityEst, tmpVector)
  }
}
purityEst <- data.frame(purityEst, stringsAsFactors = FALSE)
purityEst$Purity <- as.numeric(purityEst$Purity)

geneMat <- geneTable[-grep("Del", geneTable$Gene), 2:ncol(geneTable)]
rownames(geneMat) <- geneTable$Gene[-grep("Del", geneTable$Gene)]



geneMat2 <- 2^geneMat * 2
normIdx <- match(purityEst$Sample, colnames(geneMat2))
for (i in seq_along(normIdx)) {
  geneMat2[,which(colnames(geneMat2) == purityEst$Sample[normIdx[i]])] <- geneMat2[,which(colnames(geneMat2) == purityEst$Sample[normIdx[i]])]/purityEst$Purity[i]
}
geneMat3 <- t(round(geneMat2))

load(file = "/mnt/DATA5/tmp/kev/misc/20210409mm10SegRes.Robj")
# convert segements into 2-copy states
# only need to multiple ratio by 2 to assume 2 copy state

segRes$twoCopy <- 2^segRes$seg.mean * 2

for (i in seq_along(purityEst$Sample)) {
  segRes$twoCopy[which(segRes$ID %in% purityEst$Sample[i])] <- segRes$twoCopy[which(segRes$ID %in% purityEst$Sample[i])]/purityEst$Purity[i]
}

segRes$twoCopy_whole <- round(segRes$twoCopy)

segRes_filt <- segRes[-which(segRes$ID == "MG_9X41"),]
segRes_filt$size <- segRes_filt$loc.end - segRes_filt$loc.start
gmmInput <- segRes_filt[,c("seg.mean", "num.mark", "twoCopy_whole", "size")]


BIC <- mclustBIC(geneMat3, G = 4)
summary(BIC)
mod1 <- Mclust(geneMat3, x = BIC)
summary(mod1, parameters = TRUE)

plot(mod1, what = "classification")

# input parameters I would want to try for gmm is the length of the segment
# seg mean and number of markers


