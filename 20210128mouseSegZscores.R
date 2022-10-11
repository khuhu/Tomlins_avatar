library(DNAcopy)
library(GenomicRanges)

mouseBedFile <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.gc.bed", stringsAsFactors = FALSE, sep = "\t",
                         header = FALSE)

mouseAmplicons <- read.table("/mnt/DATA5/tmp/kev/newMouse2/cnAmplicon_matrix.txt", sep = "\t", stringsAsFactors = FALSE,
                           header = TRUE)

zscore_gc_oe_ratios <- read.table("/mnt/DATA5/tmp/kev/newMouse2/gcCorrectedCounts_matrix.txt", sep = "\t", stringsAsFactors = FALSE,
                                    header = TRUE)

crossTableIdx <- readxl::read_xlsx("/home/kevhu/data/20201207annotations.xlsx")
crossTableIdx <- crossTableIdx[,1:3]
crossTableIdx$old_name2 <- tolower(crossTableIdx$old_name)
crossTableIdx$old_name2 <- str_remove(crossTableIdx$old_name2 , "[[:punct:]]")

newPanelCalls_names <- as.numeric(str_remove(str_remove(colnames(mouseAmplicons)[4:40], "MG_"), "X.*"))
newPanelCalls_names2 <- newPanelCalls_names
for (i in unique(newPanelCalls_names)) {
  newPanelCalls_names2[which(newPanelCalls_names %in% i)] <- crossTableIdx$old_name2[which(crossTableIdx$mg_id %in% i)]
}

colnames(mouseAmplicons)[4:40] <- newPanelCalls_names2
colnames(zscore_gc_oe_ratios)[4:40] <- newPanelCalls_names2

mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
mouseAmplicons2 <- cbind(mouseBedFile$V1[mouseBedIdx], mouseBedFile$V2[mouseBedIdx], mouseAmplicons)
mouseAmplicons2[,6:42] <- log2(mouseAmplicons2[,6:42])


mouseAmplicons2$`mouseBedFile$V1[mouseBedIdx]` <- as.numeric(mouseAmplicons2$`mouseBedFile$V1[mouseBedIdx]`)

allMouseCalls <- colnames(mouseAmplicons2)[6:42]
mouseNormal <- c("13604n", "14104t", "14154n", "14433n",
                 "2405n", "2519n", "2796n", "3867n","8234n")

BPN <- c("2027lte", "2405lot")
BPRN <- c("2163lot", "3807lot", "14085lot","14433lote")
APA <- c("14593rt", "14595rt", "14596rt", "14594rt")
PRN <- c("kcovc33", "kcovc41", "kcovc43")

allMouseCalls <- c(BPRN, BPN)

segResults <- NULL
for (i in BPN[1]) {
  tmpCNA_obj <- CNA(cbind(mouseAmplicons2[[i]]),
                    mouseAmplicons2$ChromNum, mouseAmplicons2$`mouseBedFile$V2[mouseBedIdx]`,
                    data.type="logratio",sampleid=i)
  smoothed_tmpCNA <- smooth.CNA(tmpCNA_obj)
  segment_tmpCNA <- segment(smoothed_tmpCNA, verbose = 1, undo.splits = "sdundo", undo.SD = 0.2, alpha = 0.05,
                            p.method = c("hybrid"))
  segToPlot <- segment_tmpCNA
  segToPlot$output$seg.mean[which(segToPlot$output$seg.mean < -4)] <- -4
  segToPlot$output$seg.mean[which(segToPlot$output$seg.mean > 4)] <- 4
  segToPlot$data[[3]][which(segToPlot$data[[3]] < -4)] <- -4
  segToPlot$data[[3]][which(segToPlot$data[[3]] > 4)] <- 4
  plot(segToPlot, ylim=c(-4,4))
  
  tmpCalls <- segments.p(segment_tmpCNA)
  #tmpCalls_filt <- tmpCalls[which(abs(tmpCalls$seg.mean) > .2),]
  tmpCalls_filt <- tmpCalls
  
  #tmpDf <- cbind("sample" = rep(i, nrow(tmpCalls_filt)), tmpCalls_filt)
  segResults <- rbind(segResults, tmpCalls_filt)
  ###pdf
  ###plot(segment_tmpCNA, plot.type="s", ylim = c(-2,2))
  ###dev.of()
  
}

all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                             ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                              end = zscore_gc_oe_ratios$EndPos))

segResults$ID <- str_remove(segResults$ID, "X")
z_vector <- NULL
nmean_vector <- NULL
tmean_vector <- NULL
for (i in unique(segResults$ID)) {
  tmp_seg_res <- segResults[which(segResults$ID == i), ]
  tmp_seg_ranges <- GRanges(seqnames = Rle(tmp_seg_res$chrom), 
                            ranges = IRanges(start = tmp_seg_res$loc.start,
                                             end = tmp_seg_res$loc.end))
  for (j in seq_along(tmp_seg_ranges)) {
    tmpOverlap <- findOverlaps(all_probes_grange, tmp_seg_ranges[j])
    tmpTest <- zscore_gc_oe_ratios[[i]][queryHits(tmpOverlap)]
    #tmpTest2 <- tmpTest/mean(tmpTest)
    #tmpTest2 <- tmpTest
    
    tmpNormal <- zscore_gc_oe_ratios[queryHits(tmpOverlap), mouseNormal]
    #tmpNormal2 <- apply(tmpNormal, 2, function(x) x/mean(x))
    #tmpNormal2 <- apply(tmpNormal, 2, function(x) sum(x))
    tmpNormal2 <- unlist(tmpNormal)
    
    z_stat <- (mean(tmpTest) - mean(tmpNormal2))/sd(tmpNormal2)
    normal_seg_mean <- mean(tmpNormal2)
    tumor_seg_mean <- mean(tmpTest)
    z_vector <- c(z_vector, z_stat)
    nmean_vector <- c(nmean_vector, normal_seg_mean)
    tmean_vector <- c(tmean_vector, tumor_seg_mean)
  }
}

segmental_zscores <- cbind(segResults, "z-scores" = z_vector, "normal_seg_mean" = nmean_vector,
                           "tumor_seg_mean" = tmean_vector)
segmental_zscores2 <- segmental_zscores[which(segmental_zscores$num.mark > 10),]
segmental_zscores2$length <- segmental_zscores2$loc.end - segmental_zscores2$loc.start
largeChromCutoff <- 30000000
segmental_zscores3 <- segmental_zscores2[which(segmental_zscores2$length >= largeChromCutoff),]
segmental_zscores3 <- segmental_zscores3[which(abs(segmental_zscores3$seg.mean) > .2),]
segmental_zscores3 <- segmental_zscores3[which(abs(segmental_zscores3$`z-scores`) > 1.645),]

plotDf <- data.frame("zscores" = segmental_zscores2$`z-scores`,
                     "segMean" = segmental_zscores2$seg.mean,
                     "markers" = log2(segmental_zscores2$num.mark)/max(log2(segmental_zscores2$num.mark)), stringsAsFactors = FALSE)

### the upper left quadrant calls are from a normal liver with < 80% on target, could be the reason for the phenomena
###

ggplot(data = plotDf, aes(x = segMean, y = zscores)) + geom_point(size = plotDf$markers) + 
  geom_hline(yintercept=c(-1.645, 1.645), linetype="dashed", color = "red", size=0.25) + xlab("log2(segMean)") +
  ylab("zscores(one-sided horizontal cutoffs)") +  geom_vline(xintercept=c(log2(1/2), log2(3/2)), linetype="dashed", color = "blue", size=0.25) +
  geom_vline(xintercept=0, color = "black", size=0.5) + geom_hline(yintercept=0, color = "black", size=0.5)



ggplot(data = plotDf, aes(x = segMean, y = zscores)) + geom_point(size = plotDf$markers) + 
  geom_hline(yintercept=c(-1.645, 1.645), linetype="dashed", color = "red", size=0.25) + xlab("log2(segMean)") +
  ylab("zscores(one-sided horizontal cutoffs)") +  geom_vline(xintercept=c(log2(1/2), log2(3/2)), linetype="dashed", color = "blue", size=0.25) +
  geom_vline(xintercept=0, color = "black", size=0.5) + geom_hline(yintercept=0, color = "black", size=0.5) + xlim(c(-1.5,1.0)) + ylim(c(-50,50))




### freq plot

segmental_zscores3 <- segmental_zscores2[which(segmental_zscores2$`z-scores` > 1.6 | segmental_zscores2$`z-scores` < -1.6),]
cnDf_mouse <- segmental_zscores3[,c("chrom", "loc.start", "loc.end", "num.mark", "seg.mean", "ID")]
cnDf_mouse <- cnDf_mouse[-which(cnDf_mouse$ID == "2611n"),]
colnames(cnDf_mouse) <- c("chromosome","start", "end", "probes","segmean","sample")

cnFreq(cnDf_mouse, genome="mm10", CN_low_cutoff = -0.2, CN_high_cutoff = 0.2)

