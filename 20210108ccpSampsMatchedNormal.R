library(pheatmap)
library(DNAcopy)

targetDir <- "/mnt/DATA5/tmp/kev/matchedCCP/"
setwd(targetDir)
geneCalls <- system("find . -name 'combinedCalls.txt'", intern = TRUE)
geneCalls <-str_remove(geneCalls, "^\\./")
ampCalls <- system("find . -name 'cnAmplicon_matrix.txt'", intern = TRUE)
ampCalls <-str_remove(ampCalls, "^\\./")


ccpBedFile <- read.table("/home/kevhu/data/bedFiles/CCP.noTrack.GC.bed", stringsAsFactors = FALSE, sep = "\t",
                         header = FALSE)
ccpGeneListOrdered <- unique(ccpBedFile$V8)

fullGene <- NULL
for (i in geneCalls) {
  tmpDf <- read.table(paste0(targetDir, i), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  tmpDf$CopyNumberRatio <- log2(tmpDf$CopyNumberRatio)
  tmpSig <- which(abs(tmpDf$CopyNumberRatio) > 0.20 & tmpDf$Log10QValue < -1.30103)
  tmpDf$CopyNumberRatio[-tmpSig] <- 0
  tmpMatrix <- dcast(tmpDf, Gene ~ Sample, value.var = "CopyNumberRatio")
  tmpMatrix <- tmpMatrix[match(ccpGeneListOrdered, tmpMatrix$Gene),]
  if (i == geneCalls[1]) {
    fullGene <- rbind(fullGene, tmpMatrix)
  } else {
    fullGene <- cbind(fullGene, tmpMatrix[,2:ncol(tmpMatrix)])
  }
}


fullGeneMat <- fullGene[, 2:ncol(fullGene)]
rownames(fullGeneMat) <- fullGene$Gene
fullGeneMat <- fullGeneMat[, -grep("Control",  colnames(fullGeneMat))]

heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-2,2,4/1000)

pheatmap(mat = t(fullGeneMat), cluster_rows = FALSE, cluster_cols = FALSE, color = heatMapCol,
         breaks = colors.breaks, fontsize = 3, cellwidth = 3, cellheight = 30)

### since i can't combine them all into one matrix because each of the amplicons have different dropped targets,
### segmentation needs to be done on each sample loop .......... or not - only need to look at a few examples
tmpDfAmp <- read.table(paste0(targetDir, ampCalls[2]), header = TRUE, stringsAsFactors = FALSE, sep = "\t")


ccpBedIdx <- as.numeric(str_remove(tmpDfAmp$AmpliconId, "AMP_"))
ccpAmplicons2 <- cbind(ccpBedFile$V2[ccpBedIdx], tmpDfAmp)
#ccpAmplicons2 <- cbind(ccpBedFile$V1[ccpBedIdx], ccpBedFile$V2[ccpBedIdx], tmpDfAmp)

### planning to make the output of the amplicon data to be log2 + some type of median centering - global helps, but gene level may be better
ccpAmplicons2$PR940_20_CCP <- log2(ccpAmplicons2$PR940_20_CCP)
ccpAmplicons2$PR940_20_CCP <- ccpAmplicons2$PR940_20_CCP- median(ccpAmplicons2$PR940_20_CCP)
#ccpAmplicons2$`ccpBedFile$V1[ccpBedIdx]` <- as.numeric(ccpAmplicons2$`ccpBedFile$V1[ccpBedIdx]`)


CNA.object <- CNA(cbind(ccpAmplicons2$PR940_20_CCP),
                  ccpAmplicons2$ChromNum, ccpAmplicons2$`ccpBedFile$V2[ccpBedIdx]`,
                  data.type="logratio",sampleid="PR940")

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose = 1, undo.splits = "sdundo",
                                       undo.SD = 0.2, alpha = 0.10, p.method = c("hybrid"), eta = 0.10)

plot(segment.smoothed.CNA.object, plot.type="s", ylim = c(-2,2))
tmp <- segments.p(segment.smoothed.CNA.object)
tmpAll <- tmp[which(abs(tmp$seg.mean) > .2),]



### looking at general variance of markers for gains and losses on the chr8 arm
tmpArmGains <- GRanges(seqnames = Rle(c(8,8)), 
                       ranges = IRanges(start = c(41812803, 42873524),
                                        end = c(48686723,145742872)))
ccpAmplicons2_grange <- GRanges(seqnames = Rle(ccpAmplicons2$ChromNum),
                                IRanges(start = ccpAmplicons2$StartPos,
                                        end = ccpAmplicons2$EndPos))

#for (i in seq_along(tmpArmGains)) {
#  tmpOverlap <- findOverlapPairs(tmpArmGains[i], ccpAmplicons2_grange)
#  tmpZscore <- mean()
#}


tmpOverlap <- findOverlaps(tmpArmGains[2], ccpAmplicons2_grange)
tmpCnr <- ccpAmplicons2$PR940_20_CCP[subjectHits(tmpOverlap)]

qqnorm(tmpCnr)
qqline(tmpCnr)

n <- length(tmpCnr)
z_stat <- (mean(tmpCnr) - 0)/ (sd(tmpCnr)/sqrt(n))
z_stat

### quickly showing aghCalls norm is just a persample global median norm
#data(Wilting)
#Wilting <- make_cghRaw(Wilting)
#cghdata <- preprocess(Wilting, maxmiss=30, nchrom=22)
#norm.cghdata <- normalize(cghdata, method="median", smoothOutliers=TRUE)
#tmp <- Wilting@assayData$copynumber
#tmp <- as.data.frame(Wilting@assayData$copynumber)
#tmp <- tmp[-which(is.na(tmp$AdCA10)),]
#head(tmp$AdCA10 - median(tmp$AdCA10))
#norm.cghdata <- norm.cghdata[,1:2]
#seg.cghdata <- segmentData(norm.cghdata, method="DNAcopy",undo.splits="sdundo", clen=10, relSDlong=5)


tmpPoolDensity <- data.frame("CNR" = ccpAmplicons2$PR941,
                             "AMP" = ccpAmplicons2$AmpliconId,
                             "pool" = 1,stringsAsFactors = FALSE)

tmpPoolDensity$AMP <- as.numeric(str_remove_all(tmpPoolDensity$AMP, "AMP_"))
for (i in seq_along(tmpPoolDensity$AMP)) {
  tmpPoolDensity$pool[i] <- ifelse(tmpPoolDensity$AMP[i] %% 2 == 0, 2, 1)
}

tmpPoolDensity$pool <- factor(tmpPoolDensity$pool)

ggplot(data=tmpPoolDensity, aes(x=CNR, group=pool, fill=pool)) +
  geom_histogram(alpha=.4) + facet_wrap(~pool)


summary(tmpPoolDensity$CNR[which(tmpPoolDensity$pool == 1)])
summary(tmpPoolDensity$CNR[which(tmpPoolDensity$pool == 2)])





