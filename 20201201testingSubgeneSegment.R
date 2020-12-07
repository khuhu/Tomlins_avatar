tmpCombinedCalls <- read.table("/mnt/DATA5/tmp/kev/testMouse/combinedCalls.txt", sep = "\t",
                               header = TRUE, stringsAsFactors = FALSE)

### make a table analogous to what's in the gene table .. should be easy enough, order by sample name and gene 
### values of the table are just Q-values
bedFile <- read.table("/home/kevhu/data/bedFiles/IAD124056_167_Designed.del.nopool.gc.bed",
                      sep = "\t", header = FALSE)

dfCorrectedCounts <- read.table("/mnt/DATA5/tmp/kev/testMouse/cnAmplicon_matrix.txt",
                      sep = "\t", header = TRUE)
dfCorrectedCounts2 <- cbind(bedFile$V1, bedFile$V2, dfCorrectedCounts)
dfCorrectedCounts2[,6:161] <- log2(dfCorrectedCounts2[,6:161])


dfNf1 <- dfCorrectedCounts2[which(dfCorrectedCounts2$Gene == "Nf1"),]
dfBrca1 <- dfCorrectedCounts2[which(dfCorrectedCounts2$Gene == "Brca1"),]
dfTrp53 <- dfCorrectedCounts2[which(dfCorrectedCounts2$Gene == "Trp53"),]

df3Genes <- rbind(dfNf1, dfBrca1, dfTrp53)

library(DNAcopy)

CNA.object <- CNA(cbind(dfNf1$KC_ovc23),
                  dfNf1$ChromNum, dfNf1$`bedFile$V2`,
                  data.type="logratio",sampleid="ovc23")

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
tmp <- segments.p(segment.smoothed.CNA.object)
tmp <- tmp[-which(is.na(tmp$pval)),]

plot(segment.smoothed.CNA.object, plot.type="s")


### test Nf1 (137), Brca1 (72), Trp53 (19)


tmpNames <- colnames(dfCorrectedCounts2)
tmpDf <- NULL
for (i in 6:161) {
  CNA.object <- CNA(cbind(dfNf1[ ,i]),
                    dfNf1$ChromNum, dfNf1$`bedFile$V2`,
                    data.type="logratio",sampleid = tmpNames[i])
  smoothed.CNA.object <- smooth.CNA(CNA.object);
  segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1);
  tmp <- segments.p(segment.smoothed.CNA.object);
  tmp <- tmp[-which(is.na(tmp$pval)),]
  if (nrow(tmp) == 0) {
    next()
  }
  print(i)
  tmp <- tmp[which(abs(tmp$seg.mean) == max(abs(tmp$seg.mean))),]
  check <- tmp$pval < 0.05
  if (check == FALSE) {
    next()
  }
  delInfo <- tmp[ , c("ID", "seg.mean", "pval")]
  tmpDf <- rbind(tmpDf, delInfo)
}

tmpCombinedCallsNf1 <- tmpCombinedCalls[which(tmpCombinedCalls$Gene == "Nf1"),]
tmpCombinedCallsNf1$CopyNumberRatio[match(str_remove(tmpDf$ID, "X"), tmpCombinedCallsNf1$Sample)]


### create distrubtions of all these genes deletions and then 

