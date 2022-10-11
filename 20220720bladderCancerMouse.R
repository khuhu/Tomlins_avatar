### interesting gains and losses seen despite low tumor content
### running seg prune 20 for bladder samples

# source("/home/kevhu/scripts/20220427clusterSeg.R")
# mainDir <- c("/mnt/DATA6/mouseData/copynumber/")
# listOfReports <- c("Auto_user_AUS5-239-BBN_mouse_bladder_MG_493_562")
# 
# 
# allReportSegs <- NULL
# i <- listOfReports[2]
# for (i in listOfReports) {
#   amplicons <- read.table(paste0(mainDir, i, "/cnAmplicon_matrix.txt"),
#                           sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
#   zscore_gc_oe_ratios <- read.table(paste0(mainDir, i, "/gcCorrectedCounts_matrix.txt"),
#                                     sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
#   bed <- read.table(paste0(mainDir, i, "/bed.txt"),
#                     sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
#   normal <- read.table(paste0(mainDir, i, "/normals.txt"),
#                        sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
#   
#   mouseAmplicons <- amplicons
#   mouseBedFile <- bed
#   extraCols <- c("AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene", "NumGC",
#                  "Length", "GC", "TotalPool", "Weights", "MinIndex",
#                  "MaxIndex", "NumProbes", "Label", "GeneNum", "Color")
#   mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
#   mouseBed2 <- data.frame("AmpliconId" = mouseBedIdx, "ChromNum" = mouseBedFile$V1[mouseBedIdx],
#                           "StartPos" = mouseBedFile$V2[mouseBedIdx],
#                           "EndPos" = mouseBedFile$V3[mouseBedIdx])
#   mouseBed2$ChromNum <- as.numeric(str_remove(str_replace(mouseBed2$ChromNum, "chrX", "chr23"), "chr"))
#   mouseAmplicons2 <- mouseAmplicons[,-which(colnames(mouseAmplicons) %in% extraCols)]
#   mouseAmplicons2[,2:(ncol(mouseAmplicons2))] <- log2(mouseAmplicons2[,2:(ncol(mouseAmplicons2))])
#   colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)] <- nameStripper(colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)])
#   allMouseCalls <- colnames(mouseAmplicons2)[2:(ncol(mouseAmplicons2))]
#   colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)] <- nameStripper(colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)])
#   
#   all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
#                                ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
#                                                 end = zscore_gc_oe_ratios$EndPos))
#   mouseNormal <- normal$V1
#   mouseNormal <- nameStripper(mouseNormal)
#   mouseBed2$avgPos <- apply(mouseBed2[,3:4], 1, function(x) sum(x)/2)
#   
#   noiseAmps <- mouseAmplicons2[,which(colnames(mouseAmplicons2) %in% mouseNormal)]
#   if (any(grepl("\\.", colnames(noiseAmps)))) {
#     noiseAmps <- noiseAmps[, -grep("\\.", colnames(noiseAmps))]
#   }
#   noiseAmps  <- 2^noiseAmps
#   # quantile(apply(noiseAmps, 1, sd), seq(0,1,0.01))
#   
#   mouseBed2 <- mouseBed2[-which(apply(noiseAmps, 1, sd) > 0.27),]
#   mouseAmplicons2 <- mouseAmplicons2[-which(apply(noiseAmps, 1, sd) > 0.27),]
#   
#   j <- allMouseCalls[1]
#   
#   cl <- makeCluster(25)
#   registerDoParallel(cl)
#   allSeg <- foreach(j = allMouseCalls, .combine = "rbind",
#                     .packages = c("PSCBS", "copynumber")) %dopar% {
#                       source("/home/kevhu/scripts/20220427clusterSeg.R")
#                       tmp <- data.frame("chromosome"  = mouseBed2$ChromNum, "x" = mouseBed2$StartPos,
#                                         "y" = mouseAmplicons2[[j]])
#                       tmpName <- j
#                       colnames(tmp) <- c("chromosome", "x", "y")
#                       rownames(tmp) <- mouseBed2$AmpliconId
#                       tmpGaps <- findLargeGaps(tmp, minLength = 1.5e+07)
#                       tmpKnownSegs <- gapsToSegments(tmpGaps)
#                       tmpFit <- segmentByCBS(tmp, undo = 1, knownSegments = tmpKnownSegs,
#                                              p.method = c("hybrid"))
#                       tmpSegRes <- tmpFit$output
#                       tmpSegRes$sampleName <- tmpName
#                       tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$mean)),]
#                       tmpSegRes2 <- clustSeg(tmpSegRes, distVar = 0.2)
#                       tmpSegRes2
#                     }
#   stopCluster(cl)
#   colnames(allSeg) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
#   tmpZScore <- calcZscore(allSeg)
#   segZfilt <- segZscoresFilt(allSeg, tmpZScore)
#   segZfilt$ID <- nameStripper(segZfilt$ID)
#   segZfilt$ID <- str_remove(segZfilt$ID, "x.*")
#   
#   write.table(segZfilt, file = paste0(mainDir, i,"/pbsSegResPrune20.txt"), sep = "\t",
#               row.names = FALSE, col.names = TRUE, quote = FALSE)
#   
#   allReportSegs <- rbind(allReportSegs, segZfilt)
# }


freqPlot_cn <- function(df_cn, main = NULL, chromTextSpec = NULL){
  require(ggplot2)
  
  if(is.null(chromTextSpec)){
    chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "mm10"){
    chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "hg19") {
    chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801hg19_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  }
  
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  
  df_cn$loc.start <- df_cn$loc.start/1e6
  df_cn$loc.end <- df_cn$loc.end/1e6
  df_cn$col <- "#000000"
  df_cn$col[which(df_cn$seg.mean > 0.2)] <- "#FF0000"
  df_cn$col[which(df_cn$seg.mean < -0.2)] <- "#0000FF"
  
  # df_cn$col[which(df_cn$seg.mean > log2(3/2 * 0.9))] <- "#FF0000"
  # df_cn$col[which(df_cn$seg.mean < log2(1/2 / 0.9))] <- "#0000FF"
  
  for (i in unique(df_cn$chrom)) {
    df_cn$loc.start[which(df_cn$chrom == i)] <- df_cn$loc.start[which(df_cn$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    df_cn$loc.end[which(df_cn$chrom == i)] <- df_cn$loc.end[which(df_cn$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  df_cn <- df_cn[,c("chrom", "loc.start", "loc.end", "seg.mean", "col")]
  colnames(df_cn) <- c("chrom", "start", "end","cn", "col")
  
  ggplot(df_cn) + geom_hline(yintercept = seq(-3, 3, 0.5), color = "#D4D4D4") +
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn), colour = df_cn$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(-3, 3, by = 1),
                                                           limits = c(-3.2, 3.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 2.3, label = chrom), size = 2.5) + ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))
}

bladderVariants <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/Variants 7-13-22.xlsx", sheet = 2)


bladderVariants_filt <- bladderVariants[which(bladderVariants$Annotation == "Prioritized" | bladderVariants$Gene.refGene...9 == "Trp53"), c(1:21)]
colnames(bladderVariants_filt) <- str_remove(colnames(bladderVariants_filt), pattern = "\\..*")


bladderSeg <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-239-BBN_mouse_bladder_MG_493_562/pbsSegResPrune20.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)

bladderSeg_filt <- bladderSeg[grep("pp22", bladderSeg$ID), ]

stripped <- nameStripper(str_remove(str_remove(unique(bladderVariants$Sample), "\\-"), "x.*"))
nameDf <- data.frame("original" = unique(bladderVariants$Sample), "stripped" = stripped)

bladderSeg_filt$ID <- nameDf$original[match(bladderSeg_filt$ID, nameDf$stripped)]
bladderSeg_filt <- bladderSeg_filt[which(bladderSeg_filt$chrom %in% c(1:19)),]
bladderSeg_filt[, 9:12] <- lapply(bladderSeg_filt[, 9:12], function(x) signif(x, digits = 2))
bladderSeg_filt$seg.mean[which(abs(bladderSeg_filt$seg.mean) < 0.2)] <- 0

sample <- unique(bladderVariants_filt$Sample)
tc <- c(29, 14, 6, 10, 15, 25, 5, 22, 8, 13, 8, 18)
trp53muts <- c(0, 3, 2, 3, 5, 3, 3, 2, 4, 3, 3, 1)
tcDf <- data.frame("sample" = sample, "tc" = tc, "mutCount" = trp53muts)


i <- 1
for (i in seq_along(tcDf$sample)) {
  print(i)
  
  tmpSample <- tcDf$sample[i]
  tmpTc <- tcDf$tc[i]
  tmpTrp53 <- tcDf$mutCount[i]
  
  testDf_cn <- bladderSeg_filt[which(bladderSeg_filt$ID == tmpSample),]
  
  # b_2 <- freqPlot_cn(testDf_cn, main = paste(tmpSample, "tc: ", tmpTc/100, " trp53 count:", tmpTrp53))
  
  png(filename = paste0("/mnt/DATA5/tmp/kev/misc/pbsSegBafsPrune15_bladder/", tmpSample, ".png"), width = 1800, height = 500)
  print(freqPlot_cn(testDf_cn, main = paste(tmpSample, "tc: ", tmpTc/100, " trp53 count:", tmpTrp53)))
  dev.off()
  
}

bladderSeg_filt$size <- bladderSeg_filt$loc.end - bladderSeg_filt$loc.start

bladderSeg_filt2 <- bladderSeg_filt[which(bladderSeg_filt$size < 50000), ]
bladderSeg_filt2 <- bladderSeg_filt2[which(abs(bladderSeg_filt2$seg.mean) > 0), ]

write.table(bladderSeg_filt2, "/mnt/DATA5/tmp/kev/misc/20220720bladderGeneLevelSeg.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

