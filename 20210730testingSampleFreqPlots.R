
freqPlot_hg <- function(gainsDf, lossesDf, speciesType = "human", main = "no title"){ 
  require(ggplot2)
  
  if (speciesType == "human") {
    chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801hg19_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
    
  } else if(speciesType == "mouse"){
    chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  }
  
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  tmpDel <- lossesDf
  tmpDel$Freq <- tmpDel$Freq * -1
  tmpDel$Start <- tmpDel$Start/1e6
  tmpDel$End <- tmpDel$End/1e6
  tmpDel$col <- "#0000FF"
  
  tmpAmp <- gainsDf
  tmpAmp$Start <- tmpAmp$Start/1e6
  tmpAmp$End <- tmpAmp$End/1e6
  tmpAmp$col <- "#FF0000"
  
  
  tmpAmpDel_graph <- rbind(tmpAmp, tmpDel)
  for (i in unique(tmpAmpDel_graph$Chr)) {
    tmpAmpDel_graph$Start[which(tmpAmpDel_graph$Chr == i)] <- tmpAmpDel_graph$Start[which(tmpAmpDel_graph$Chr == i)] + 
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    tmpAmpDel_graph$End[which(tmpAmpDel_graph$Chr == i)] <- tmpAmpDel_graph$End[which(tmpAmpDel_graph$Chr == i)] + 
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  
  tmpAmpDel_graph2 <- NULL
  for (i in 1:nrow(tmpAmpDel_graph)) {
    tmpVec1 <- c(tmpAmpDel_graph$Chr[i], tmpAmpDel_graph$Start[i], tmpAmpDel_graph$Freq[i], tmpAmpDel_graph$col[i])
    tmpVec2 <- c(tmpAmpDel_graph$Chr[i], tmpAmpDel_graph$End[i], tmpAmpDel_graph$Freq[i], tmpAmpDel_graph$col[i])
    tmpAmpDel_graph2 <- rbind(tmpAmpDel_graph2, tmpVec1, tmpVec2)
  }
  tmpAmpDel_graph2 <- data.frame(tmpAmpDel_graph2, stringsAsFactors = FALSE)
  tmpAmpDel_graph2[,1:3] <- lapply(tmpAmpDel_graph2[,1:3], as.numeric)
  colnames(tmpAmpDel_graph2) <- c("chrom", "pos", "freq", "col")
  tmpAmpDel_graph2$freq <- tmpAmpDel_graph2$freq/100
  
  ggplot(tmpAmpDel_graph2, aes(x = pos, y = freq)) + geom_line() + geom_vline(xintercept=chromBreak) +
    geom_area(aes(x=pos, y=ifelse(freq>0, freq,0)), fill="red", position = 'identity') + 
    geom_area(aes(x=pos, y=ifelse(freq<0, freq,0)), fill="blue", position = 'identity') + theme_bw() +
    ylim(c(-1,1)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
    scale_x_continuous(breaks = NULL) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos, label = chrom)) + ggtitle(main)
  
}


freqPlot_hg(hgsc_arm_amp_bed, hgsc_arm_del_bed, speciesType = "human", main = "hgsc")





###
###
### below is how the files with locations are made

tmpTable <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genome/hg19/hg19.chrom.sizes", sep = "\t",
                       stringsAsFactors = FALSE, header = FALSE)
hg_chrom_sizes <- tmpTable[1:24,]
hg_chrom_sizes$V1 <- str_remove(hg_chrom_sizes$V1, "chr")
hg_chrom_sizes2 <- data.frame("chrom" = hg_chrom_sizes$V1, "start" = 0, "end" = hg_chrom_sizes$V2,
                              "length80" = round(hg_chrom_sizes$V2 * .80))

write.table(hg_chrom_sizes2, "/mnt/DATA5/tmp/kev/misc/20210730hg19chromSizeTable.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# create predifined limits for the graph
hg_chrom_sizes3 <- hg_chrom_sizes2
hg_chrom_sizes3 <- hg_chrom_sizes3[-which(hg_chrom_sizes3$chrom =="Y"),]
hg_chrom_sizes3$chrom <- as.numeric(str_replace(hg_chrom_sizes3$chrom, "X", "23"))
hg_chrom_sizes3 <- hg_chrom_sizes3[order(hg_chrom_sizes3$chrom),]


graphingLimits <- data.frame("chrom" = 1:23, "start" = rep(0,23), "end" = rep(0,23), stringsAsFactors = FALSE)
for (i in seq_along(unique(hg_chrom_sizes3$chrom))) {
  print(i)
  dummyDf <- hg_chrom_sizes3[hg_chrom_sizes3$chrom == i,]
  if (i == 1) {
    tmpStart <- 0
    tmpEnd <- dummyDf$end/1e6
    graphingLimits$start[i] <- tmpStart
    graphingLimits$end[i] <- tmpEnd
  } else{
    tmpStart <- graphingLimits$end[i - 1] + 10
    tmpEnd <- graphingLimits$end[i - 1] + 10 + dummyDf$end/1e6
    graphingLimits$start[i] <- tmpStart
    graphingLimits$end[i] <- tmpEnd
  }
}

graphingLimits$start <- round(graphingLimits$start)
graphingLimits$end <- round(graphingLimits$end)
chromBreak <- graphingLimits$end[1:23] + 5
chromBreak <- c(0, chromBreak)

chromTextCoords <- apply(graphingLimits[,2:3], 1, mean)
chromTextdf <- data.frame("chrom" = 1:23, "xpos" = chromTextCoords, "ypos" = 0.9,
                          "chromBreaksPos" = chromBreak[2:length(chromBreak)],
                          "graphingStart" = graphingLimits$start,
                          stringsAsFactors = FALSE)

write.table(chromTextdf, "/mnt/DATA5/tmp/kev/misc/20210801hg19_graphingLimits.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

dummyPoints <- data.frame()
for (i in chromTextdf$chrom) {
  tmpStartPos <- chromTextdf$graphingStart[which(chromTextdf$chrom == i)] + 0.5
  tmpEndPos <- chromTextdf$chromBreaksPos[which(chromTextdf$chrom == i)] - 0.05
  dummyPoints <- rbind(dummyPoints,
                       c(i, tmpStartPos, 0, "#000000"),
                       c(i, tmpEndPos, 0, "#000000"))
}
colnames(dummyPoints) <- c("chrom", "pos", "freq", "col")
dummyPoints[,1:3] <- lapply(dummyPoints[,1:3], as.numeric)

write.table(dummyPoints, "/mnt/DATA5/tmp/kev/misc/20210803hg19_dummyPoints.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# test plot

tmpDel <- hgsc_arm_del_bed
tmpDel$Freq <- tmpDel$Freq * -1
tmpDel$Start <- tmpDel$Start/1e6
tmpDel$End <- tmpDel$End/1e6

tmpAmp <- hgsc_arm_amp_bed
tmpAmp$Start <- tmpAmp$Start/1e6
tmpAmp$End <- tmpAmp$End/1e6


tmpAmp_graph <- tmpAmp
for (i in unique(tmpAmp_graph$Chr)) {
  tmpAmp_graph$Start[which(tmpAmp_graph$Chr == i)] <- tmpAmp_graph$Start[which(tmpAmp_graph$Chr == i)] + 
    graphingLimits$start[which(graphingLimits$chrom == i)]
  tmpAmp_graph$End[which(tmpAmp_graph$Chr == i)] <- tmpAmp_graph$End[which(tmpAmp_graph$Chr == i)] + 
    graphingLimits$start[which(graphingLimits$chrom == i)]
}

tmpAmp_graph2 <- NULL
for (i in 1:nrow(tmpAmp_graph)) {
  tmpVec1 <- c(tmpAmp_graph$Chr[i], tmpAmp_graph$Start[i], tmpAmp_graph$Freq[i])
  tmpVec2 <- c(tmpAmp_graph$Chr[i], tmpAmp_graph$End[i], tmpAmp_graph$Freq[i])
  tmpAmp_graph2 <- rbind(tmpAmp_graph2, tmpVec1, tmpVec2)
}
 
tmpAmp_graph2 <- data.frame(tmpAmp_graph2, stringsAsFactors = FALSE)
tmpAmp_graph2$col <- "#FF0000"
colnames(tmpAmp_graph2) <- c("chrom", "pos", "freq", "col")


tmpDel_graph <- tmpDel
for (i in unique(tmpDel_graph$Chr)) {
  tmpDel_graph$Start[which(tmpDel_graph$Chr == i)] <- tmpDel_graph$Start[which(tmpDel_graph$Chr == i)] + 
    graphingLimits$start[which(graphingLimits$chrom == i)]
  tmpDel_graph$End[which(tmpDel_graph$Chr == i)] <- tmpDel_graph$End[which(tmpDel_graph$Chr == i)] + 
    graphingLimits$start[which(graphingLimits$chrom == i)]
}

tmpDel_graph2 <- NULL
for (i in 1:nrow(tmpDel_graph)) {
  tmpVec1 <- c(tmpDel_graph$Chr[i], tmpDel_graph$Start[i], tmpDel_graph$Freq[i])
  tmpVec2 <- c(tmpDel_graph$Chr[i], tmpDel_graph$End[i], tmpDel_graph$Freq[i])
  tmpDel_graph2 <- rbind(tmpDel_graph2, tmpVec1, tmpVec2)
}

tmpDel_graph2 <- data.frame(tmpDel_graph2, stringsAsFactors = FALSE)
tmpDel_graph2$col <- "#0000FF"
colnames(tmpDel_graph2) <- c("chrom", "pos", "freq", "col")

tmpGraph_fin <- rbind(tmpAmp_graph2, tmpDel_graph2)
tmpGraph_fin$freq <- tmpGraph_fin$freq/100


ggplot(tmpGraph_fin, aes(x = pos, y = freq)) + geom_line() + geom_vline(xintercept=chromBreak) +
  geom_area(aes(x=pos, y=ifelse(freq>0, freq,0)), fill="red", position = 'identity') + 
  geom_area(aes(x=pos, y=ifelse(freq<0, freq,0)), fill="blue", position = 'identity') + theme_bw() +
  ylim(c(-1,1)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
  scale_x_continuous(breaks = NULL) +
  geom_text(data = chromTextdf, aes(x = xpos, y = ypos, label = chrom))


### creating the mouse chrom equiavlent

mm10_chromSizes <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genome/mm10/mm10.chrom.sizes",
                              sep = "\t", stringsAsFactors = FALSE, header = FALSE)
mm10_chromSizes <- mm10_chromSizes[c(1:18,20,21),]
mm10_chromSizes$V1 <- str_replace(mm10_chromSizes$V1, "chrX", "chr20")
mm10_chromSizes$V1 <- as.numeric(str_remove(mm10_chromSizes$V1, "chr"))
mm10_chromSizes <- mm10_chromSizes[order(mm10_chromSizes$V1),]
rownames(mm10_chromSizes) <- NULL

graphingLimits <- data.frame("chrom" = 1:20, "start" = rep(0,20), "end" = rep(0,20), stringsAsFactors = FALSE)
for (i in seq_along(unique(mm10_chromSizes$V1))) {
  print(i)
  dummyDf <- mm10_chromSizes[mm10_chromSizes$V1 == i,]
  if (i == 1) {
    tmpStart <- 0
    tmpEnd <- dummyDf$V2/1e6
    graphingLimits$start[i] <- tmpStart
    graphingLimits$end[i] <- tmpEnd
  } else{
    tmpStart <- graphingLimits$end[i - 1] + 10
    tmpEnd <- graphingLimits$end[i - 1] + 10 + dummyDf$V2/1e6
    graphingLimits$start[i] <- tmpStart
    graphingLimits$end[i] <- tmpEnd
  }
}

graphingLimits$start <- round(graphingLimits$start)
graphingLimits$end <- round(graphingLimits$end)
chromBreak <- graphingLimits$end[1:20] + 5
chromBreak <- c(0, chromBreak)

chromTextCoords <- apply(graphingLimits[,2:3], 1, mean)
chromTextdf <- data.frame("chrom" = 1:20, "xpos" = chromTextCoords, "ypos" = 0.9,
                          "chromBreaksPos" = chromBreak[2:length(chromBreak)],
                          "graphingStart" = graphingLimits$start,
                          stringsAsFactors = FALSE)
write.table(chromTextdf, "/mnt/DATA5/tmp/kev/misc/20210801mm10_graphingLimits.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


dummyPoints <- data.frame()
for (i in chromTextdf$chrom) {
  tmpStartPos <- chromTextdf$graphingStart[which(chromTextdf$chrom == i)] + 0.5
  tmpEndPos <- chromTextdf$chromBreaksPos[which(chromTextdf$chrom == i)] - 0.05
  dummyPoints <- rbind(dummyPoints,
                       c(i, tmpStartPos, 0, "#000000"),
                       c(i, tmpEndPos, 0, "#000000"))
}
colnames(dummyPoints) <- c("chrom", "pos", "freq", "col")
dummyPoints[,1:3] <- lapply(dummyPoints[,1:3], as.numeric)

write.table(dummyPoints, "/mnt/DATA5/tmp/kev/misc/20210803mm10_dummyPoints.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
