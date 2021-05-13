library(ggplot2)
library(gridExtra)
library(VennDiagram)

tmpCombinedCalls <- read.table("/mnt/DATA5/tmp/kev/testMouse2/combinedCalls.txt", sep = "\t",
                               header = TRUE, stringsAsFactors = FALSE)
tmpCombinedCalls$CopyNumberRatio <- log2(tmpCombinedCalls$CopyNumberRatio)
### make a table analogous to what's in the gene table .. should be easy enough, order by sample name and gene 
### values of the table are just Q-values

brca1Ids <- readxl::read_xlsx("/home/kevhu/data/Brca1DelNames.xlsx", col_names = FALSE)
brca1Ids$...1 <- str_remove(brca1Ids$...1, "X")
brca1Ids$...1 <- tolower(brca1Ids$...1)
brca1Ids$...1 <- str_replace_all(brca1Ids$...1, "[[:punct:]]", "")

Nf1Ids <- readxl::read_xlsx("/home/kevhu/data/Nf1DelList.xlsx", col_names = FALSE)
Nf1Ids$...1 <- str_remove(Nf1Ids$...1, "X")
Nf1Ids$...1 <- tolower(Nf1Ids$...1)
Nf1Ids$...1 <- str_replace_all(Nf1Ids$...1, "[[:punct:]]", "")

bedFile <- read.table("/home/kevhu/data/bedFiles/IAD124056_167_Designed.nopool.gc.bed",
                      sep = "\t", header = FALSE)

dfCorrectedCounts <- read.table("/mnt/DATA5/tmp/kev/testMouse/cnAmplicon_matrix.txt",
                      sep = "\t", header = TRUE)
bedIdx <- as.numeric(str_remove(dfCorrectedCounts$AmpliconId, "AMP_"))
dfCorrectedCounts2 <- cbind(bedFile$V1[bedIdx], bedFile$V2[bedIdx], dfCorrectedCounts)

dfCorrectedCounts2[,6:161] <- log2(dfCorrectedCounts2[,6:161])


dfNf1 <- dfCorrectedCounts2[which(dfCorrectedCounts2$Gene == "Nf1"),]
dfBrca1 <- dfCorrectedCounts2[which(dfCorrectedCounts2$Gene == "Brca1"),]
dfTrp53 <- dfCorrectedCounts2[which(dfCorrectedCounts2$Gene == "Trp53"),]

library(DNAcopy)

CNA.object <- CNA(cbind(dfNf1$KC_ovc43),
                  dfNf1$ChromNum, dfNf1$`bedFile$V2[bedIdx]`,
                  data.type="logratio",sampleid="ovc43")

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
                    dfNf1$ChromNum, dfNf1$`bedFile$V2[bedIdx]`,
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
  delInfo <- tmp[ , c("ID", "seg.mean", "pval", "loc.start", "loc.end")]
  tmpDf <- rbind(tmpDf, delInfo)
}

tmpDf$ID <- str_remove(tmpDf$ID, "X")
tmpDf$ID <- tolower(tmpDf$ID)
tmpDf$ID <- str_replace_all(tmpDf$ID, "[[:punct:]]", "")

tmpCombinedCalls$Sample <- str_remove(tmpCombinedCalls$Sample, "X")
tmpCombinedCalls$Sample <- tolower(tmpCombinedCalls$Sample)
tmpCombinedCalls$Sample <- str_replace_all(tmpCombinedCalls$Sample, "[[:punct:]]", "")

tmpCombinedCallsNf1 <- tmpCombinedCalls[which(tmpCombinedCalls$Gene == "Nf1Del"),]

tmpNf1_reg <- tmpCombinedCallsNf1[match(tmpDf$ID, tmpCombinedCallsNf1$Sample),]
tmpNf1_reg <- tmpNf1_reg[-which(is.na(tmpNf1_reg$Sample)),]
tmpNf1_seg <- tmpDf[match(tmpNf1_reg$Sample, tmpDf$ID),]


stats::cor(tmpNf1_reg$CopyNumberRatio, tmpNf1_seg$seg.mean)
stats::cor(-tmpNf1_reg$Log10PValue, -log10(tmpNf1_seg$pval))


nf1Cnr <- ggplot(data = data.frame("reg" = tmpNf1_reg$CopyNumberRatio, "seg" = tmpNf1_seg$seg.mean)) + geom_point(aes(x = reg, y = seg)) + 
  xlab("reg log2(CNR)") + ylab("seg log2(CNR)") + ggtitle("Nf1 deletion CNR (r=0.997; n = 65, amplicons = 137)") + theme(plot.title = element_text(hjust = 0.5))

nf1Pval <- ggplot(data = data.frame("reg" = -tmpNf1_reg$Log10PValue, "seg" = -log10(tmpNf1_seg$pval))) + geom_point(aes(x = reg, y = seg)) + 
  xlab("reg -log10(pval))") + ylab("seg -log10(pval)") + ggtitle("Nf1 deletion pval (n = 65, amplicons = 137)") + theme(plot.title = element_text(hjust = 0.5))


### for pie-chart
tmpDf2 <- tmpDf[which(abs(tmpDf$seg.mean) > 0.2),]
countNf1_seg <- length(which(tmpDf2$ID %in% Nf1Ids$...1))
vennNf1Id_seg <- tmpDf2$ID

tmpCombinedCallsNf1 <- tmpCombinedCallsNf1[which(abs(tmpCombinedCallsNf1$CopyNumberRatio) > 0.2), ]
tmpCombinedCallsNf1 <- tmpCombinedCallsNf1[which(10^tmpCombinedCallsNf1$Log10PValue < 0.05), ]
vennNf1Id_reg <- tmpCombinedCallsNf1$Sample


### Brca1

tmpNames <- colnames(dfCorrectedCounts2)
tmpDf <- NULL
for (i in 6:161) {
  CNA.object <- CNA(cbind(dfBrca1[ ,i]),
                    dfBrca1$ChromNum, dfBrca1$`bedFile$V2[bedIdx]`,
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
  delInfo <- tmp[ , c("ID", "seg.mean", "pval", "loc.start", "loc.end")]
  tmpDf <- rbind(tmpDf, delInfo)
}

tmpDf$ID <- str_remove(tmpDf$ID, "X")
tmpDf$ID <- tolower(tmpDf$ID)
tmpDf$ID <- str_replace_all(tmpDf$ID, "[[:punct:]]", "")

tmpCombinedCallsBrca1 <- tmpCombinedCalls[which(tmpCombinedCalls$Gene == "Brca1Del"),]

tmpBrca1_reg <- tmpCombinedCallsBrca1[match(tmpDf$ID, tmpCombinedCallsBrca1$Sample),]
tmpBrca1_reg <- tmpBrca1_reg[-which(is.na(tmpBrca1_reg$Sample)),]
tmpBrca1_seg <- tmpDf[match(tmpBrca1_reg$Sample, tmpDf$ID),]


stats::cor(tmpBrca1_reg$CopyNumberRatio, tmpBrca1_seg$seg.mean)
stats::cor(-tmpBrca1_reg$Log10PValue, -log10(tmpBrca1_seg$pval))


Brca1Cnr <- ggplot(data = data.frame("reg" = tmpBrca1_reg$CopyNumberRatio, "seg" = tmpBrca1_seg$seg.mean)) + geom_point(aes(x = reg, y = seg)) + 
  xlab("reg log2(CNR)") + ylab("seg log2(CNR)") + ggtitle("Brca1 deletion CNR (r=0.955; n = 61, amplicons = 72)") + theme(plot.title = element_text(hjust = 0.5))

Brca1Pval <- ggplot(data = data.frame("reg" = -tmpBrca1_reg$Log10PValue, "seg" = -log10(tmpBrca1_seg$pval))) + geom_point(aes(x = reg, y = seg)) + 
  xlab("reg -log10(pval)") + ylab("seg -log10(pval)") + ggtitle("Brca1 deletion pval (r = 0.797;n = 61, amplicons = 72)") + theme(plot.title = element_text(hjust = 0.5))


tmpDf2 <- tmpDf[which(abs(tmpDf$seg.mean) > 0.2),]
countBrca1_seg <- length(which(tmpDf2$ID %in% brca1Ids$...1))
vennBrca1Id_seg <- tmpDf2$ID

tmpCombinedCallsBrca2 <- tmpCombinedCallsBrca1[which(abs(tmpCombinedCallsBrca1$CopyNumberRatio) > 0.2), ]
tmpCombinedCallsBrca2 <- tmpCombinedCallsBrca2[which(10^tmpCombinedCallsBrca2$Log10PValue < 0.05), ]
vennBrca1Id_reg <- tmpCombinedCallsBrca2$Sample

# p53


tmpNames <- colnames(dfCorrectedCounts2)
tmpDf <- NULL
for (i in 6:161) {
  CNA.object <- CNA(cbind(dfTrp53[ ,i]),
                    dfTrp53$ChromNum, dfTrp53$`bedFile$V2[bedIdx]`,
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
  #check <- tmp$pval < 0.05
  #if (check == FALSE) {
  #  next()
  #}
  delInfo <- tmp[ , c("ID", "seg.mean", "pval", "loc.start", "loc.end")]
  tmpDf <- rbind(tmpDf, delInfo)
}

tmpDf$ID <- str_remove(tmpDf$ID, "X")
tmpDf$ID <- tolower(tmpDf$ID)
tmpDf$ID <- str_replace_all(tmpDf$ID, "[[:punct:]]", "")

tmpCombinedCallsBrca1 <- tmpCombinedCalls[which(tmpCombinedCalls$Gene == "Brca1Del"),]

tmpBrca1_reg <- tmpCombinedCallsBrca1[match(tmpDf$ID, tmpCombinedCallsBrca1$Sample),]
tmpBrca1_reg <- tmpBrca1_reg[-which(is.na(tmpBrca1_reg$Sample)),]
tmpBrca1_seg <- tmpDf[match(tmpBrca1_reg$Sample, tmpDf$ID),]


stats::cor(tmpBrca1_reg$CopyNumberRatio, tmpBrca1_seg$seg.mean)
stats::cor(-tmpBrca1_reg$Log10PValue, -log10(tmpBrca1_seg$pval))

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20201207subGeneScatter.pdf", useDingbats = FALSE, height = 10, width = 12)
grid.arrange(nf1Cnr, nf1Pval, Brca1Cnr, Brca1Pval, ncol = 2)
dev.off()



venn.diagram(
  x = list(vennBrca1Id_seg, vennBrca1Id_reg,brca1Ids$...1),
  category.names = c("seg" , "reg","genotype"),
  main = "Brca1",
  filename = '/mnt/DATA5/tmp/kev/misc/20201207vennDiagBrca1.png',
  output=TRUE
)


venn.diagram(
  x = list(vennNf1Id_seg, vennNf1Id_reg, Nf1Ids$...1),
  category.names = c("seg" , "reg","genotype"),
  main = "Nf1",
  filename = '/mnt/DATA5/tmp/kev/misc/20201207vennDiagNf1.png',
  output=TRUE
)



