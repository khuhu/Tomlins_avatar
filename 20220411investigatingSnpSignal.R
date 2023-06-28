# install.packages("PSCBS")
library("PSCBS")

library(copynumber)
data <- PSCBS::exampleData("paired.chr01")
data <- data[, c("chromosome", "x", "CT")]
colnames(data)[3] <- "y"
data <- dropSegmentationOutliers(data)
gaps <- findLargeGaps(data, minLength = 1e+06)
gaps
knownSegments <- gapsToSegments(gaps)


### can try to use segmentByCBS from PSCBS to segment with gaps 


freqPlot_baf <- function(df, main = "no title", chromTextSpec = NULL){
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
  
  ### divide positions by megabase locations
  
  df$CHROM <- as.numeric(str_remove(df$CHROM, "chr"))
  df$point <- df$POS/1e6
  
  for (i in unique(df$CHROM)) {
    df$point[which(df$CHROM == i)] <- df$point[which(df$CHROM == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  
  df <- df[,c("CHROM", "point", "AF")]
  colnames(df) <- c("chrom", "pos", "af")
  
  
  if (nrow(df) == 0) {
    ggplot(df, aes(x = pos, y = af)) + geom_hline(yintercept = seq(0, 1, 0.1), color = "#D4D4D4")+
      geom_point(size = 0.25, color = df$color) + geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 1, by = 0.10),
                                                             limits = c(0,1.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + .2, label = chrom), size = 2.5) 
  } else{
    df$color <- "#000000"
    df$color[which(df$af < 0.3)] <- "#FF0000"
    df$color[which(df$af > 0.7)] <- "#FF0000"
    
    ggplot(df, aes(x = pos, y = af)) + geom_hline(yintercept = seq(0, 1, 0.1), color = "#D4D4D4")+
      geom_point(size = 0.25, color = df$color) + geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 1, by = 0.10),
                                                             limits = c(0,1.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + .2, label = chrom), size = 2.5) 
  }
}

nameStripper <- function(df){
  require(stringr)
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
}



calcZscore <- function(segRes){
  sampZ1 <- NULL
  sampZ2 <- NULL
  z1_vector <- NULL
  z2_vector <- NULL
  p_vector1 <- NULL
  p_vector2 <- NULL
  nmean_vector <- NULL
  tmean_vector <- NULL
  q_vector1 <- NULL
  q_vector2 <- NULL
  
  # uniqIdx <- unique(match(str_replace_all(segRes$ID, "[[:punct:]]", ""), 
  #                         str_replace_all(str_replace_all(colnames(zscore_gc_oe_ratios), "[[:punct:]]", ""), " ", "")))
  
  uniqIdx <- unique(match(segRes$ID, colnames(zscore_gc_oe_ratios)))
  zScoreColnames <- colnames(zscore_gc_oe_ratios)[uniqIdx]
  for (i in seq_along(unique(segRes$ID))) {
    tmp_seg_res <- segRes[which(segRes$ID == unique(segRes$ID)[i]), ]
    tmp_seg_ranges <- GRanges(seqnames = Rle(tmp_seg_res$chrom), 
                              ranges = IRanges(start = tmp_seg_res$loc.start,
                                               end = tmp_seg_res$loc.end))
    print(unique(segRes$ID)[i])
    print(zScoreColnames[i])
    
    sampZ1 <- NULL
    sampZ2 <- NULL
    
    for (j in seq_along(tmp_seg_ranges)) {
      tmpOverlap <- findOverlaps(all_probes_grange, tmp_seg_ranges[j])
      tmpTest <- zscore_gc_oe_ratios[[zScoreColnames[i]]][queryHits(tmpOverlap)]
      
      if (!is.numeric(tmpTest)) {
        sampZ1 <- c(sampZ1, NA)
        sampZ2 <- c(sampZ2, NA)
        nmean_vector <- c(nmean_vector, NA)
        tmean_vector <- c(tmean_vector, NA)
        next()
      } else{
        tmpNormal <- zscore_gc_oe_ratios[queryHits(tmpOverlap), mouseNormal]
        tmpNormal2 <- unlist(apply(tmpNormal, 2, mean))
        
        z1_stat <- (mean(tmpTest) - mean(tmpNormal2))/sd(tmpNormal2)
        z2_stat <- (mean(tmpTest) - mean(tmpNormal2))/(sd(tmpNormal2) * mean(tmpTest))
        pVal_z1 <- 2.0*pnorm(-abs(z1_stat))
        pVal_z2 <- 2.0*pnorm(-abs(z2_stat))
        
        
        normal_seg_mean <- mean(tmpNormal2)
        tumor_seg_mean <- mean(tmpTest)
        sampZ1 <- c(sampZ1, pVal_z1)
        sampZ2 <- c(sampZ2, pVal_z2)
        z1_vector <- c(z1_vector, z1_stat)
        z2_vector <- c(z2_vector, z2_stat)
        p_vector1 <- c(p_vector1, pVal_z1)
        p_vector2 <- c(p_vector2, pVal_z2)
        nmean_vector <- c(nmean_vector, normal_seg_mean)
        tmean_vector <- c(tmean_vector, tumor_seg_mean)
      }
    }
    qVal_z1 <- (p.adjust(sampZ1, method = "BH"))
    qVal_z2 <- (p.adjust(sampZ2, method = "BH"))
    q_vector1 <- c(q_vector1, qVal_z1)
    q_vector2 <- c(q_vector2, qVal_z2)
    
  }
  return(list("z1_vector" = z1_vector, "z2_vector" = z2_vector,
              "p1_vector" = p_vector1, "p2_vector" = p_vector2,
              "q1_vector" = q_vector1, "q2_vector" = q_vector2,
              "nmean_vector"= nmean_vector, "tmean_vector" = tmean_vector))
}


calcZscore_log2 <- function(segRes){
  sampZ1 <- NULL
  sampZ2 <- NULL
  z1_vector <- NULL
  z2_vector <- NULL
  p_vector1 <- NULL
  p_vector2 <- NULL
  nmean_vector <- NULL
  tmean_vector <- NULL
  q_vector1 <- NULL
  q_vector2 <- NULL
  
  # uniqIdx <- unique(match(str_replace_all(segRes$ID, "[[:punct:]]", ""), 
  #                         str_replace_all(str_replace_all(colnames(zscore_gc_oe_ratios), "[[:punct:]]", ""), " ", "")))
  
  uniqIdx <- unique(match(segRes$ID, colnames(zscore_gc_oe_ratios)))
  zScoreColnames <- colnames(zscore_gc_oe_ratios)[uniqIdx]
  for (i in seq_along(unique(segRes$ID))) {
    tmp_seg_res <- segRes[which(segRes$ID == unique(segRes$ID)[i]), ]
    tmp_seg_ranges <- GRanges(seqnames = Rle(tmp_seg_res$chrom), 
                              ranges = IRanges(start = tmp_seg_res$loc.start,
                                               end = tmp_seg_res$loc.end))
    print(unique(segRes$ID)[i])
    print(zScoreColnames[i])
    
    sampZ1 <- NULL
    sampZ2 <- NULL
    
    for (j in seq_along(tmp_seg_ranges)) {
      tmpOverlap <- findOverlaps(all_probes_grange, tmp_seg_ranges[j])
      tmpTest <- zscore_gc_oe_ratios[[zScoreColnames[i]]][queryHits(tmpOverlap)]
      
      if (!is.numeric(tmpTest)) {
        sampZ1 <- c(sampZ1, NA)
        sampZ2 <- c(sampZ2, NA)
        nmean_vector <- c(nmean_vector, NA)
        tmean_vector <- c(tmean_vector, NA)
        next()
      } else{
        tmpNormal <- zscore_gc_oe_ratios[queryHits(tmpOverlap), mouseNormal]
        tmpNormal2 <- unlist(apply(tmpNormal, 2, mean))
        
        tmpTest <- log2(tmpTest)
        tmpNormal2 <- log2(tmpNormal2)
        
        z1_stat <- (mean(tmpTest) - mean(tmpNormal2))/sd(tmpNormal2)
        z2_stat <- (mean(tmpTest) - mean(tmpNormal2))/(sd(tmpNormal2) * mean(tmpTest))
        pVal_z1 <- 2.0*pnorm(-abs(z1_stat))
        pVal_z2 <- 2.0*pnorm(-abs(z2_stat))
        
        
        normal_seg_mean <- mean(tmpNormal2)
        tumor_seg_mean <- mean(tmpTest)
        sampZ1 <- c(sampZ1, pVal_z1)
        sampZ2 <- c(sampZ2, pVal_z2)
        z1_vector <- c(z1_vector, z1_stat)
        z2_vector <- c(z2_vector, z2_stat)
        p_vector1 <- c(p_vector1, pVal_z1)
        p_vector2 <- c(p_vector2, pVal_z2)
        nmean_vector <- c(nmean_vector, normal_seg_mean)
        tmean_vector <- c(tmean_vector, tumor_seg_mean)
      }
    }
    qVal_z1 <- (p.adjust(sampZ1, method = "BH"))
    qVal_z2 <- (p.adjust(sampZ2, method = "BH"))
    q_vector1 <- c(q_vector1, qVal_z1)
    q_vector2 <- c(q_vector2, qVal_z2)
    
  }
  return(list("z1_vector" = z1_vector, "z2_vector" = z2_vector,
              "p1_vector" = p_vector1, "p2_vector" = p_vector2,
              "q1_vector" = q_vector1, "q2_vector" = q_vector2,
              "nmean_vector"= nmean_vector, "tmean_vector" = tmean_vector))
}


segZscoresFilt <- function(segResults, zscores){
  
  # filts for significance, cn-change and length
  segmental_zscores <- cbind(segResults, "z-scores1" = zscores[["z1_vector"]], "z-scores2" = zscores[["z2_vector"]],
                             "p-val1" = zscores[["p1_vector"]], "p-val2" = zscores[["p2_vector"]],
                             "q-val1" = zscores[["q1_vector"]], "q-val2" = zscores[["q2_vector"]],
                             "normal_seg_mean" = zscores[["nmean_vector"]], "tumor_seg_mean" = zscores[["tmean_vector"]])
  return(segmental_zscores)
}



amplicons <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/cnAmplicon_matrix.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
zscore_gc_oe_ratios <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/gcCorrectedCounts_matrix.txt",
                                  sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
bed <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/bed.txt",
                  sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
normal <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/normals.txt",
                  sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
mouseAmplicons <- amplicons
mouseBedFile <- bed
extraCols <- c("AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene", "NumGC",
               "Length", "GC", "TotalPool", "Weights", "MinIndex",
               "MaxIndex", "NumProbes", "Label", "GeneNum", "Color")
mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
mouseBed2 <- data.frame("AmpliconId" = mouseBedIdx, "ChromNum" = mouseBedFile$V1[mouseBedIdx],
                        "StartPos" = mouseBedFile$V2[mouseBedIdx],
                        "EndPos" = mouseBedFile$V3[mouseBedIdx])
mouseBed2$ChromNum <- as.numeric(str_remove(str_replace(mouseBed2$ChromNum, "chrX", "chr23"), "chr"))
mouseAmplicons2 <- mouseAmplicons[,-which(colnames(mouseAmplicons) %in% extraCols)]
mouseAmplicons2[,2:(ncol(mouseAmplicons2))] <- log2(mouseAmplicons2[,2:(ncol(mouseAmplicons2))])

colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)] <- nameStripper(colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)])
allMouseCalls <- colnames(mouseAmplicons2)[2:(ncol(mouseAmplicons2))]


colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)] <- nameStripper(colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)])

all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                             ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                              end = zscore_gc_oe_ratios$EndPos))

mouseNormal <- normal$V1
mouseNormal <- nameStripper(mouseNormal)



mouseBed2$avgPos <- apply(mouseBed2[,3:4], 1, function(x) sum(x)/2)
gapList <- c(1e+06, 2.5e+06, 5e+06, 7.5e+06, 1e+07, 1.5e+07,
             2e+07, 2.5e+07, 3e+07, 3.5e+07, 4e+07, 5e+07)


for(j in gapList){
  allSeg <- NULL
  for (i in allMouseCalls) {
    tmp <- data.frame("chromosome"  = mouseBed2$ChromNum, "x" = mouseBed2$avgPos,
                      "y" = mouseAmplicons2[[i]])
    tmpName <- i
    colnames(tmp) <- c("chromosome", "x", "y")
    rownames(tmp) <- mouseBed2$AmpliconId
    tmpGaps <- findLargeGaps(tmp, minLength = j)
    tmpKnownSegs <- gapsToSegments(tmpGaps)
    tmpFit <- segmentByCBS(tmp, undo = 1, knownSegments = tmpKnownSegs)
    tmpSegRes <- tmpFit$output
    tmpSegRes$sampleName <- tmpName
    tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$chromosome)),]
    allSeg <- rbind(allSeg, tmpSegRes)
  }
  colnames(allSeg) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
  tmpZScore <- calcZscore(allSeg)
  segZfilt <- segZscoresFilt(allSeg, tmpZScore)
  segZfilt$seg.mean[which(abs(segZfilt$seg.mean) < 0.2)] <- 0
  segZfilt$seg.mean[which(segZfilt$q.val2 > 0.05)] <- 0
  segZfilt$seg.mean[which(segZfilt$seg.mean < -3)] <- -3
  segZfilt$seg.mean[which(segZfilt$seg.mean > 3)] <- 3
  segZfilt$seg.mean[which((segZfilt$end.pos - segZfilt$start.pos) < (5 * 1e6))] <- 0
  segZfilt$ID <- nameStripper(segZfilt$ID)
  segZfilt$ID <- str_remove(segZfilt$ID, "x.*")
  write.table(segZfilt, paste0("/mnt/DATA5/tmp/kev/misc/", "20220415segGap", j, ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


### parallelized above 
library(foreach)
library(doParallel)

someSeg <- function(i, j){
  tmp <- data.frame("chromosome"  = mouseBed2$ChromNum, "x" = mouseBed2$avgPos,
                    "y" = mouseAmplicons2[[i]])
  tmpName <- i
  colnames(tmp) <- c("chromosome", "x", "y")
  rownames(tmp) <- mouseBed2$AmpliconId
  tmpGaps <- findLargeGaps(tmp, minLength = j)
  tmpKnownSegs <- gapsToSegments(tmpGaps)
  tmpFit <- segmentByCBS(tmp, undo = 1, knownSegments = tmpKnownSegs)
  tmpSegRes <- tmpFit$output
  tmpSegRes$sampleName <- tmpName
  tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$chromosome)),]
  tmpSegRes
}

cl <- makeCluster(25)
registerDoParallel(cl)
for(j in gapList){
  allSeg <- foreach(i = allMouseCalls, .combine = "rbind",
                    .packages = c("PSCBS", "copynumber")) %dopar% {
                      tmp <- data.frame("chromosome"  = mouseBed2$ChromNum, "x" = mouseBed2$avgPos,
                                        "y" = mouseAmplicons2[[i]])
                      tmpName <- i
                      colnames(tmp) <- c("chromosome", "x", "y")
                      rownames(tmp) <- mouseBed2$AmpliconId
                      tmpGaps <- findLargeGaps(tmp, minLength = j)
                      tmpKnownSegs <- gapsToSegments(tmpGaps)
                      tmpFit <- segmentByCBS(tmp, undo = 1, knownSegments = tmpKnownSegs)
                      tmpSegRes <- tmpFit$output
                      tmpSegRes$sampleName <- tmpName
                      tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$chromosome)),]
                      tmpSegRes
                    }
  colnames(allSeg) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
  tmpZScore <- calcZscore(allSeg)
  segZfilt <- segZscoresFilt(allSeg, tmpZScore)
  segZfilt$seg.mean[which(abs(segZfilt$seg.mean) < 0.2)] <- 0
  segZfilt$seg.mean[which(segZfilt$q.val2 > 0.05)] <- 0
  segZfilt$seg.mean[which(segZfilt$seg.mean < -3)] <- -3
  segZfilt$seg.mean[which(segZfilt$seg.mean > 3)] <- 3
  segZfilt$seg.mean[which((segZfilt$end.pos - segZfilt$start.pos) < (5 * 1e6))] <- 0
  segZfilt$ID <- nameStripper(segZfilt$ID)
  segZfilt$ID <- str_remove(segZfilt$ID, "x.*")
  write.table(segZfilt, paste0("/mnt/DATA5/tmp/kev/misc/", "20220415segGap", j, ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
stopCluster(cl)






logRarray <- read.table("/mnt/DATA5/tmp/kev/misc/20220405testTumorLogR3.txt", sep = "\t",
                        check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
allSeg <- NULL
for (i in 4:dim(logRarray)[2]) {
  tmp <- logRarray[,c(2:3, i)]
  tmpName <- colnames(tmp)[3]
  colnames(tmp) <- c("chromosome", "x", "y")
  rownames(tmp) <- logRarray$marker
  tmpFit <- segmentByCBS(tmp, verbose = -10)
  tmpSegRes <- tmpFit$output
  tmpSegRes$sampleName <- tmpName
  tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$chromosome)),]
  allSeg <- rbind(allSeg, tmpSegRes)
}


logRarray2 <- read.table("/mnt/DATA5/tmp/kev/misc/20220410illuminaLogRr.txt", sep = "\t",
                        check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)

allSeg2 <- NULL
for (i in 4:dim(logRarray2)[2]) {
  tmp <- logRarray2[,c(2:3, i)]
  tmpName <- colnames(tmp)[3]
  colnames(tmp) <- c("chromosome", "x", "y")
  rownames(tmp) <- logRarray2$marker
  tmpFit <- segmentByCBS(tmp, verbose = -10)
  tmpSegRes <- tmpFit$output
  tmpSegRes$sampleName <- tmpName
  tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$chromosome)),]
  allSeg2 <- rbind(allSeg2, tmpSegRes)
}

logRarray3 <- read.table("/mnt/DATA5/tmp/kev/misc/20220405testTumorLogR3.txt", sep = "\t",
                        check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)

logRarray3[,4:ncol(logRarray3)]  <- logRarray3[,4:ncol(logRarray3)] * 3
allSeg3 <- NULL
for (i in 4:dim(logRarray3)[2]) {
  tmp <- logRarray3[,c(2:3, i)]
  tmpName <- colnames(tmp)[3]
  colnames(tmp) <- c("chromosome", "x", "y")
  rownames(tmp) <- logRarray3$marker
  tmpFit <- segmentByCBS(tmp, verbose = -10, undo=1)
  tmpSegRes <- tmpFit$output
  tmpSegRes$sampleName <- tmpName
  tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$chromosome)),]
  allSeg3 <- rbind(allSeg3, tmpSegRes)
}



write.table(allSeg, "/mnt/DATA5/tmp/kev/misc/20220411segResultsLogrAdj.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(allSeg2, "/mnt/DATA5/tmp/kev/misc/20220411segResultsLogrIllu.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(allSeg3, "/mnt/DATA5/tmp/kev/misc/20220411segResultsLogrAdj3x.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)


### graphing 3 seg results side by side
### mousePanel and 2 segs of snp data

sampleMap <- read.table("/mnt/DATA5/tmp/kev/misc/Sample_Map.txt", sep = "\t",
                        check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
sampleMap$stripped <- nameStripper(sampleMap$ID)
sampleMap$stripped[21] <- "13085lt"


allSegs_mouse <- read.table("/mnt/DATA5/tmp/kev/misc/20211206all_mouse_seg.txt", sep = "\t", stringsAsFactors = FALSE,
                      header = TRUE, check.names = FALSE)

allSegs_mouse$seg.mean[which(abs(allSegs_mouse$seg.mean) < 0.2)] <- 0
allSegs_mouse$seg.mean[which(allSegs_mouse$q.val2 > 0.05)] <- 0
allSegs_mouse$seg.mean[which(allSegs_mouse$seg.mean < -3)] <- -3
allSegs_mouse$seg.mean[which(allSegs_mouse$seg.mean > 3)] <- 3
allSegs_mouse$seg.mean[which((allSegs_mouse$end.pos - allSegs_mouse$start.pos) < (10 * 1e6))] <- 0
allSegs_mouse$ID <- nameStripper(allSegs_mouse$ID)
allSegs_mouse$ID <- str_remove(allSegs_mouse$ID, "x.*")

allSegs_adj <- read.table("/mnt/DATA5/tmp/kev/misc/20220411segResultsLogrAdj.txt", sep = "\t", stringsAsFactors = FALSE,
                            header = TRUE, check.names = FALSE)
colnames(allSegs_adj) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
allSegs_adj$seg.mean[which(abs(allSegs_adj$seg.mean) < 0.2)] <- 0
allSegs_adj$seg.mean[which(allSegs_adj$seg.mean < -3)] <- -3
allSegs_adj$seg.mean[which(allSegs_adj$seg.mean > 3)] <- 3
allSegs_adj$seg.mean[which((allSegs_adj$end.pos - allSegs_adj$start.pos) < (10 * 1e6))] <- 0


allSegs_illu <- read.table("/mnt/DATA5/tmp/kev/misc/20220411segResultsLogrIllu.txt", sep = "\t", stringsAsFactors = FALSE,
                            header = TRUE, check.names = FALSE)
colnames(allSegs_illu) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
allSegs_illu$seg.mean[which(abs(allSegs_illu$seg.mean) < 0.2)] <- 0
allSegs_illu$seg.mean[which(allSegs_illu$seg.mean < -3)] <- -3
allSegs_illu$seg.mean[which(allSegs_illu$seg.mean > 3)] <- 3
allSegs_illu$seg.mean[which((allSegs_illu$end.pos - allSegs_illu$start.pos) < (10 * 1e6))] <- 0

allSegs_illu$ID[which(allSegs_illu$ID == "14656peritonealmt")] <- "14656peritnealmt"


allSegs_adj3x <- read.table("/mnt/DATA5/tmp/kev/misc/20220411segResultsLogrAdj3x.txt", sep = "\t", stringsAsFactors = FALSE,
                          header = TRUE, check.names = FALSE)
colnames(allSegs_adj3x) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
allSegs_adj$seg.mean[which(abs(allSegs_adj$seg.mean) < 0.2)] <- 0
allSegs_adj3x$seg.mean[which(allSegs_adj3x$seg.mean < -3)] <- -3
allSegs_adj3x$seg.mean[which(allSegs_adj3x$seg.mean > 3)] <- 3
allSegs_adj3x$seg.mean[which((allSegs_adj3x$end.pos - allSegs_adj3x$start.pos) < (10 * 1e6))] <- 0


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


i <- "5516-KH-36"
for (i in unique(allSegs_adj$ID)) {
  print(i)
  
  seqID <- sampleMap$stripped[which(sampleMap$Name == i)]
  if (length(seqID) == 0) {
    next()
  }
  
  matchingTc <- tcDf$tc[match(tolower(seqID), tcDf$sample)]
  if (is.na(matchingTc)) {
      matchingTc <- "NA"
  }
  
  testDf_mouse <- allSegs_mouse[which(allSegs_mouse$ID == seqID),]
  if (dim(testDf_mouse)[1] == 0) {
    next
  }
  b_mouse <- freqPlot_cn(testDf_mouse, main = paste(seqID, "MousePanel, Snp:Illu, Snp:Adj", "tc:", matchingTc))
  
  testDf_adj <- allSegs_adj[which(allSegs_adj$ID == i),]
  b_adj <- freqPlot_cn(testDf_adj)
  
  
  testDf_illu <- allSegs_illu[which(allSegs_illu$ID == seqID),]
  if (length(which(testDf_illu$chrom == 0)) > 0) {
    testDf_illu <- testDf_illu[-which(testDf_illu$chrom == 0), ]
  }
  b_illu <- freqPlot_cn(testDf_illu)
  
  
  testDf_adj3x <- allSegs_adj3x[which(allSegs_adj3x$ID == i),]
  b_adj3x <- freqPlot_cn(testDf_adj3x)
  
  
  gm <- ggplotGrob(b_mouse)
  gi <- ggplotGrob(b_illu)
  ga <- ggplotGrob(b_adj)
  g3 <- ggplotGrob(b_adj3x)
  
  dev.off()
  png(file = paste0("/mnt/DATA5/tmp/kev/comparingSegs/20220411_compSegs", seqID, ".png"), width = 1800, height = 500)
  grid::grid.newpage()
  grid::grid.draw(rbind(gm, gi, ga, g3))
  dev.off()
  
}





plot(dummyDf$size[which(dummyDf$seg.mean < -0.2)], dummyDf$seg.mean[which(dummyDf$seg.mean < -0.2)],
     xlab = "size Mb", ylab = "log2R")

tbl_1 <- read.table("/mnt/DATA5/tmp/kev/misc/20220415segGap1e+06.txt", sep = "\t", stringsAsFactors = FALSE,
                    header = TRUE)
tbl_1.5 <- read.table("/mnt/DATA5/tmp/kev/misc/20220415segGap1e+06.txt", sep = "\t", stringsAsFactors = FALSE,
                      header = TRUE)
tbl_2.5 <- read.table("/mnt/DATA5/tmp/kev/misc/20220415segGap2500000.txt", sep = "\t", stringsAsFactors = FALSE,
                      header = TRUE)
tbl_5 <- read.table("/mnt/DATA5/tmp/kev/misc/20220415segGap5e+06.txt", sep = "\t", stringsAsFactors = FALSE,
                     header = TRUE)
tbl_7.5 <- read.table("/mnt/DATA5/tmp/kev/misc/20220415segGap7500000.txt", sep = "\t", stringsAsFactors = FALSE,
                    header = TRUE)
tbl_10 <- read.table("/mnt/DATA5/tmp/kev/misc/20220415segGap1e+07.txt", sep = "\t", stringsAsFactors = FALSE,
                     header = TRUE)
tbl_15 <- read.table("/mnt/DATA5/tmp/kev/misc/20220415segGap1.5e+07.txt", sep = "\t", stringsAsFactors = FALSE,
                    header = TRUE)
tbl_20 <- read.table("/mnt/DATA5/tmp/kev/misc/20220415segGap2e+07.txt", sep = "\t", stringsAsFactors = FALSE,
                    header = TRUE)

### starting at 5Mb for graphs since everything below is too sparse i.e singletons

gapList <- c(1e+06, 2.5e+06, 5e+06, 7.5e+06, 1e+07, 1.5e+07,
             2e+07, 2.5e+07, 3e+07, 3.5e+07, 4e+07, 5e+07)

quantile(tbl_5$nbrOfLoci, seq(0,1, 0.01))
quantile(tbl_7.5$nbrOfLoci, seq(0,1, 0.01))
quantile(tbl_10$nbrOfLoci, seq(0,1, 0.01))
quantile(tbl_15$nbrOfLoci, seq(0,1, 0.01))
quantile(tbl_20$nbrOfLoci, seq(0,1, 0.01))


### running for another set of samples


amplicons <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/cnAmplicon_matrix.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
zscore_gc_oe_ratios <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/gcCorrectedCounts_matrix.txt",
                                  sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
bed <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/bed.txt",
                  sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
normal <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/normals.txt",
                     sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
mouseAmplicons <- amplicons
mouseBedFile <- bed
extraCols <- c("AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene", "NumGC",
               "Length", "GC", "TotalPool", "Weights", "MinIndex",
               "MaxIndex", "NumProbes", "Label", "GeneNum", "Color")
mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
mouseBed2 <- data.frame("AmpliconId" = mouseBedIdx, "ChromNum" = mouseBedFile$V1[mouseBedIdx],
                        "StartPos" = mouseBedFile$V2[mouseBedIdx],
                        "EndPos" = mouseBedFile$V3[mouseBedIdx])
mouseBed2$ChromNum <- as.numeric(str_remove(str_replace(mouseBed2$ChromNum, "chrX", "chr23"), "chr"))
mouseAmplicons2 <- mouseAmplicons[,-which(colnames(mouseAmplicons) %in% extraCols)]
mouseAmplicons2[,2:(ncol(mouseAmplicons2))] <- log2(mouseAmplicons2[,2:(ncol(mouseAmplicons2))])

colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)] <- nameStripper(colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)])
allMouseCalls <- colnames(mouseAmplicons2)[2:(ncol(mouseAmplicons2))]
colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)] <- nameStripper(colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)])

all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                             ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                              end = zscore_gc_oe_ratios$EndPos))
mouseNormal <- normal$V1
mouseNormal <- nameStripper(mouseNormal)
mouseBed2$avgPos <- apply(mouseBed2[,3:4], 1, function(x) sum(x)/2)

cl <- makeCluster(25)
registerDoParallel(cl)
for(j in gapList){
  allSeg <- foreach(i = allMouseCalls, .combine = "rbind",
                    .packages = c("PSCBS", "copynumber")) %dopar% {
                      tmp <- data.frame("chromosome"  = mouseBed2$ChromNum, "x" = mouseBed2$avgPos,
                                        "y" = mouseAmplicons2[[i]])
                      tmpName <- i
                      colnames(tmp) <- c("chromosome", "x", "y")
                      rownames(tmp) <- mouseBed2$AmpliconId
                      tmpGaps <- findLargeGaps(tmp, minLength = j)
                      tmpKnownSegs <- gapsToSegments(tmpGaps)
                      tmpFit <- segmentByCBS(tmp, undo = 1, knownSegments = tmpKnownSegs)
                      tmpSegRes <- tmpFit$output
                      tmpSegRes$sampleName <- tmpName
                      tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$chromosome)),]
                      tmpSegRes
                    }
  colnames(allSeg) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
  tmpZScore <- calcZscore(allSeg)
  segZfilt <- segZscoresFilt(allSeg, tmpZScore)
  segZfilt$seg.mean[which(abs(segZfilt$seg.mean) < 0.2)] <- 0
  segZfilt$seg.mean[which(segZfilt$q.val2 > 0.05)] <- 0
  segZfilt$seg.mean[which(segZfilt$seg.mean < -3)] <- -3
  segZfilt$seg.mean[which(segZfilt$seg.mean > 3)] <- 3
  segZfilt$seg.mean[which((segZfilt$end.pos - segZfilt$start.pos) < (5 * 1e6))] <- 0
  segZfilt$ID <- nameStripper(segZfilt$ID)
  segZfilt$ID <- str_remove(segZfilt$ID, "x.*")
  write.table(segZfilt, paste0("/mnt/DATA5/tmp/kev/misc/", "20220418_AUS132_segGap", j, ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
stopCluster(cl)




### plotting gapped seg results with snps so I can get more accurate comparison charts for one copy gain and loss

allPbs <- NULL
mainDir <- c("/mnt/DATA6/mouseData/copynumber/")
# listOfReports <- c("Auto_user_AUS5-120-MG_EFD4_BBN_334_304", "Auto_user_AUS5-138-MG_cho_20210621_354_343",
#                    "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
#                    "Auto_user_AUS5-156-MG_Fearon_20210809_374_382", "Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384",
#                    "Auto_user_AUS5-76-MG_test1_255_185")

listOfReports <- c("Auto_user_AUS5-120-MG_EFD4_BBN_334_304", "Auto_user_AUS5-138-MG_cho_20210621_354_343",
                   "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                   "Auto_user_AUS5-156-MG_Fearon_20210809_374_382", "Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384",
                   "Reanalysis_AUS5-76-MG_test1_217")

for (i in listOfReports) {
  tmpTable <- fread(paste0(mainDir, i, "/pbsSegRes.txt"))
  tmpTable$dir <- i 
  allPbs <- rbind(allPbs, tmpTable)
}

additionalNamesDf <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20201207annotations.xlsx") 
additionalNamesDf$mg_id <- paste0("mg", additionalNamesDf$mg_id)
additionalNamesDf$stripped <- nameStripper(additionalNamesDf$old_name)


tcDf <- read.table("/mnt/DATA5/tmp/kev/misc/20220321hgscTc.txt", sep = "\t",
                   header = TRUE, stringsAsFactors = FALSE)
non_zero_20 <- read.table("/mnt/DATA5/tmp/kev/misc/20211203mm10BafsNonzero20.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

allPbs$seg.mean[which(abs(allPbs$seg.mean) < 0.2)] <- 0
allPbs$seg.mean[which(allPbs$q.val2 > 0.05)] <- 0
allPbs$seg.mean[which(allPbs$seg.mean < -3)] <- -3
allPbs$seg.mean[which(allPbs$seg.mean > 3)] <- 3
allPbs$chrom <- str_replace_all(allPbs$chrom, "23", "20")

i <- "1628lt"
for (i in unique(allPbs$ID)) {
  print(i)
  
  
  matchingTc <- tcDf$tc[match(tolower(i), tcDf$sample)]
  
  # sampleName <- i
  # if (grepl("mg", sampleName)) {
  #   sampleName <- additionalNamesDf$stripped[which(additionalNamesDf$mg_id == i)]
  # }
  # 
  testDf_cn <- allPbs[which(allPbs$ID == i),]
  b_2 <- freqPlot_cn(testDf_cn, main = paste(i, "MousePanel, bafs tc:", matchingTc))
  
  
  testDf <- non_zero_20[which(non_zero_20$sample == i),]
  a_20 <- freqPlot_baf(testDf)
  
  g20 <- ggplotGrob(a_20)
  gB_2 <- ggplotGrob(b_2)
  
  dev.off()
  png(filename = paste0("/mnt/DATA5/tmp/kev/misc/pbsSegBafs/pbsSegBafs_", i, ".png"), width = 1800, height = 500)
  grid::grid.newpage()
  grid::grid.draw(rbind(gB_2, g20))
  dev.off()
  
  pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/pbsSegBafs/pbsSegBafs_", i, ".pdf"), width = 15, height = 7)
  grid::grid.newpage()
  grid::grid.draw(rbind(gB_2, g20))
  dev.off()
  
  
}



### trying pcsbs with h-clust pruning 0.15 & 0.25
source("/home/kevhu/scripts/20220427clusterSeg.R")
mainDir <- c("/mnt/DATA6/mouseData/copynumber/")
listOfReports <- c("Auto_user_AUS5-120-MG_EFD4_BBN_334_304", "Auto_user_AUS5-138-MG_cho_20210621_354_343",
                   "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                   "Auto_user_AUS5-156-MG_Fearon_20210809_374_382", "Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384",
                   "Reanalysis_AUS5-76-MG_test1_217", "Auto_user_AUS5-76-MG_test1_255_185")


allReportSegs <- NULL
i <- listOfReports[2]
for (i in listOfReports) {
  amplicons <- read.table(paste0(mainDir, i, "/cnAmplicon_matrix.txt"),
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  zscore_gc_oe_ratios <- read.table(paste0(mainDir, i, "/gcCorrectedCounts_matrix.txt"),
                                    sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  bed <- read.table(paste0(mainDir, i, "/bed.txt"),
                    sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  normal <- read.table(paste0(mainDir, i, "/normals.txt"),
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  
  mouseAmplicons <- amplicons
  mouseBedFile <- bed
  extraCols <- c("AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene", "NumGC",
                 "Length", "GC", "TotalPool", "Weights", "MinIndex",
                 "MaxIndex", "NumProbes", "Label", "GeneNum", "Color")
  mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
  mouseBed2 <- data.frame("AmpliconId" = mouseBedIdx, "ChromNum" = mouseBedFile$V1[mouseBedIdx],
                          "StartPos" = mouseBedFile$V2[mouseBedIdx],
                          "EndPos" = mouseBedFile$V3[mouseBedIdx])
  mouseBed2$ChromNum <- as.numeric(str_remove(str_replace(mouseBed2$ChromNum, "chrX", "chr23"), "chr"))
  mouseAmplicons2 <- mouseAmplicons[,-which(colnames(mouseAmplicons) %in% extraCols)]
  mouseAmplicons2[,2:(ncol(mouseAmplicons2))] <- log2(mouseAmplicons2[,2:(ncol(mouseAmplicons2))])
  colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)] <- nameStripper(colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)])
  allMouseCalls <- colnames(mouseAmplicons2)[2:(ncol(mouseAmplicons2))]
  colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)] <- nameStripper(colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)])
  
  all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                               ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                                end = zscore_gc_oe_ratios$EndPos))
  mouseNormal <- normal$V1
  mouseNormal <- nameStripper(mouseNormal)
  mouseBed2$avgPos <- apply(mouseBed2[,3:4], 1, function(x) sum(x)/2)
  
  noiseAmps <- mouseAmplicons2[,which(colnames(mouseAmplicons2) %in% mouseNormal)]
  if (any(grepl("\\.", colnames(noiseAmps)))) {
    noiseAmps <- noiseAmps[, -grep("\\.", colnames(noiseAmps))]
  }
  noiseAmps  <- 2^noiseAmps
  # quantile(apply(noiseAmps, 1, sd), seq(0,1,0.01))
  
  mouseBed2 <- mouseBed2[-which(apply(noiseAmps, 1, sd) > 0.27),]
  mouseAmplicons2 <- mouseAmplicons2[-which(apply(noiseAmps, 1, sd) > 0.27),]
  
  j <- allMouseCalls[1]
  
  cl <- makeCluster(25)
  registerDoParallel(cl)
  allSeg <- foreach(j = allMouseCalls, .combine = "rbind",
                    .packages = c("PSCBS", "copynumber")) %dopar% {
                      source("/home/kevhu/scripts/20220427clusterSeg.R")
                      tmp <- data.frame("chromosome"  = mouseBed2$ChromNum, "x" = mouseBed2$StartPos,
                                        "y" = mouseAmplicons2[[j]])
                      tmpName <- j
                      colnames(tmp) <- c("chromosome", "x", "y")
                      rownames(tmp) <- mouseBed2$AmpliconId
                      tmpGaps <- findLargeGaps(tmp, minLength = 1.5e+07)
                      tmpKnownSegs <- gapsToSegments(tmpGaps)
                      tmpFit <- segmentByCBS(tmp, undo = 1, knownSegments = tmpKnownSegs,
                                             p.method = c("hybrid"))
                      tmpSegRes <- tmpFit$output
                      tmpSegRes$sampleName <- tmpName
                      tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$mean)),]
                      tmpSegRes2 <- clustSeg(tmpSegRes, distVar = 0.2)
                      tmpSegRes2
                    }
  stopCluster(cl)
  colnames(allSeg) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
  tmpZScore <- calcZscore(allSeg)
  segZfilt <- segZscoresFilt(allSeg, tmpZScore)
  segZfilt$ID <- nameStripper(segZfilt$ID)
  segZfilt$ID <- str_remove(segZfilt$ID, "x.*")
  
  write.table(segZfilt, file = paste0(mainDir, i,"/pbsSegResPrune20.txt"), sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  allReportSegs <- rbind(allReportSegs, segZfilt)
}

write.table(allReportSegs, paste0("/mnt/DATA5/tmp/kev/misc/2022041915MbGapAllSegPrune20.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


### trying global clustering method instead of per chromosome

source("/home/kevhu/scripts/20220427clusterSeg.R")
mainDir <- c("/mnt/DATA6/mouseData/copynumber/")
listOfReports <- c("Auto_user_AUS5-120-MG_EFD4_BBN_334_304", "Auto_user_AUS5-138-MG_cho_20210621_354_343",
                   "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                   "Auto_user_AUS5-156-MG_Fearon_20210809_374_382", "Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384",
                   "Reanalysis_AUS5-76-MG_test1_217", "Auto_user_AUS5-76-MG_test1_255_185")


allReportSegs <- NULL
i <- listOfReports[2]
for (i in listOfReports) {
  amplicons <- read.table(paste0(mainDir, i, "/cnAmplicon_matrix.txt"),
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  zscore_gc_oe_ratios <- read.table(paste0(mainDir, i, "/gcCorrectedCounts_matrix.txt"),
                                    sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  bed <- read.table(paste0(mainDir, i, "/bed.txt"),
                    sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  normal <- read.table(paste0(mainDir, i, "/normals.txt"),
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  
  mouseAmplicons <- amplicons
  mouseBedFile <- bed
  extraCols <- c("AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene", "NumGC",
                 "Length", "GC", "TotalPool", "Weights", "MinIndex",
                 "MaxIndex", "NumProbes", "Label", "GeneNum", "Color")
  mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
  mouseBed2 <- data.frame("AmpliconId" = mouseBedIdx, "ChromNum" = mouseBedFile$V1[mouseBedIdx],
                          "StartPos" = mouseBedFile$V2[mouseBedIdx],
                          "EndPos" = mouseBedFile$V3[mouseBedIdx])
  mouseBed2$ChromNum <- as.numeric(str_remove(str_replace(mouseBed2$ChromNum, "chrX", "chr23"), "chr"))
  mouseAmplicons2 <- mouseAmplicons[,-which(colnames(mouseAmplicons) %in% extraCols)]
  mouseAmplicons2[,2:(ncol(mouseAmplicons2))] <- log2(mouseAmplicons2[,2:(ncol(mouseAmplicons2))])
  colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)] <- nameStripper(colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)])
  allMouseCalls <- colnames(mouseAmplicons2)[2:(ncol(mouseAmplicons2))]
  colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)] <- nameStripper(colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)])
  
  all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                               ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                                end = zscore_gc_oe_ratios$EndPos))
  mouseNormal <- normal$V1
  mouseNormal <- nameStripper(mouseNormal)
  mouseBed2$avgPos <- apply(mouseBed2[,3:4], 1, function(x) sum(x)/2)
  
  noiseAmps <- mouseAmplicons2[,which(colnames(mouseAmplicons2) %in% mouseNormal)]
  if (any(grepl("\\.", colnames(noiseAmps)))) {
    noiseAmps <- noiseAmps[, -grep("\\.", colnames(noiseAmps))]
  }
  noiseAmps  <- 2^noiseAmps
  # quantile(apply(noiseAmps, 1, sd), seq(0,1,0.01))
  
  mouseBed2 <- mouseBed2[-which(apply(noiseAmps, 1, sd) > 0.27),]
  mouseAmplicons2 <- mouseAmplicons2[-which(apply(noiseAmps, 1, sd) > 0.27),]
  
  j <- allMouseCalls[1]
  
  cl <- makeCluster(25)
  registerDoParallel(cl)
  allSeg <- foreach(j = allMouseCalls, .combine = "rbind",
                    .packages = c("PSCBS", "copynumber")) %dopar% {
                      source("/home/kevhu/scripts/20220427clusterSeg.R")
                      tmp <- data.frame("chromosome"  = mouseBed2$ChromNum, "x" = mouseBed2$StartPos,
                                        "y" = mouseAmplicons2[[j]])
                      tmpName <- j
                      colnames(tmp) <- c("chromosome", "x", "y")
                      rownames(tmp) <- mouseBed2$AmpliconId
                      tmpGaps <- findLargeGaps(tmp, minLength = 1.5e+07)
                      tmpKnownSegs <- gapsToSegments(tmpGaps)
                      tmpFit <- segmentByCBS(tmp, undo = 1, knownSegments = tmpKnownSegs,
                                             p.method = c("hybrid"))
                      tmpSegRes <- tmpFit$output
                      tmpSegRes$sampleName <- tmpName
                      tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$mean)),]
                      tmpSegRes2 <- clustSegG(tmpSegRes, distVar = 0.3)
                      tmpSegRes2
                    }
  stopCluster(cl)
  colnames(allSeg) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
  tmpZScore <- calcZscore(allSeg)
  segZfilt <- segZscoresFilt(allSeg, tmpZScore)
  segZfilt$ID <- nameStripper(segZfilt$ID)
  segZfilt$ID <- str_remove(segZfilt$ID, "x.*")
  
  write.table(segZfilt, file = paste0(mainDir, i,"/pbsSegResPrune20G.txt"), sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  allReportSegs <- rbind(allReportSegs, segZfilt)
}

write.table(allReportSegs, paste0("/mnt/DATA5/tmp/kev/misc/2022041915MbGapAllSegPrune20_globalCluster.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)






source("/home/kevhu/scripts/20220427clusterSeg.R")
mainDir <- c("/mnt/DATA6/mouseData/copynumber/")
listOfReports <- c("Auto_user_AUS5-120-MG_EFD4_BBN_334_304", "Auto_user_AUS5-138-MG_cho_20210621_354_343",
                   "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                   "Auto_user_AUS5-156-MG_Fearon_20210809_374_382", "Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384",
                   "Reanalysis_AUS5-76-MG_test1_217", "Auto_user_AUS5-76-MG_test1_255_185")


allReportSegs <- NULL
for (i in listOfReports) {
  amplicons <- read.table(paste0(mainDir, i, "/cnAmplicon_matrix.txt"),
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  zscore_gc_oe_ratios <- read.table(paste0(mainDir, i, "/gcCorrectedCounts_matrix.txt"),
                                    sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  bed <- read.table(paste0(mainDir, i, "/bed.txt"),
                    sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  normal <- read.table(paste0(mainDir, i, "/normals.txt"),
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  
  mouseAmplicons <- amplicons
  mouseBedFile <- bed
  extraCols <- c("AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene", "NumGC",
                 "Length", "GC", "TotalPool", "Weights", "MinIndex",
                 "MaxIndex", "NumProbes", "Label", "GeneNum", "Color")
  mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
  mouseBed2 <- data.frame("AmpliconId" = mouseBedIdx, "ChromNum" = mouseBedFile$V1[mouseBedIdx],
                          "StartPos" = mouseBedFile$V2[mouseBedIdx],
                          "EndPos" = mouseBedFile$V3[mouseBedIdx])
  mouseBed2$ChromNum <- as.numeric(str_remove(str_replace(mouseBed2$ChromNum, "chrX", "chr23"), "chr"))
  mouseAmplicons2 <- mouseAmplicons[,-which(colnames(mouseAmplicons) %in% extraCols)]
  mouseAmplicons2[,2:(ncol(mouseAmplicons2))] <- log2(mouseAmplicons2[,2:(ncol(mouseAmplicons2))])
  colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)] <- nameStripper(colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)])
  allMouseCalls <- colnames(mouseAmplicons2)[2:(ncol(mouseAmplicons2))]
  colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)] <- nameStripper(colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)])
  
  all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                               ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                                end = zscore_gc_oe_ratios$EndPos))
  mouseNormal <- normal$V1
  mouseNormal <- nameStripper(mouseNormal)
  mouseBed2$avgPos <- apply(mouseBed2[,3:4], 1, function(x) sum(x)/2)
  
  noiseAmps <- mouseAmplicons2[,which(colnames(mouseAmplicons2) %in% mouseNormal)]
  if (any(grepl("\\.", colnames(noiseAmps)))) {
    noiseAmps <- noiseAmps[, -grep("\\.", colnames(noiseAmps))]
  }
  noiseAmps  <- 2^noiseAmps
  # quantile(apply(noiseAmps, 1, sd), seq(0,1,0.01))
  
  mouseBed2 <- mouseBed2[-which(apply(noiseAmps, 1, sd) > 0.27),]
  mouseAmplicons2 <- mouseAmplicons2[-which(apply(noiseAmps, 1, sd) > 0.27),]
  
  j <- allMouseCalls[38]
  
  cl <- makeCluster(25)
  registerDoParallel(cl)
  allSeg <- foreach(j = allMouseCalls, .combine = "rbind",
                    .packages = c("PSCBS", "copynumber")) %dopar% {
                      source("/home/kevhu/scripts/20220427clusterSeg.R")
                      tmp <- data.frame("chromosome"  = mouseBed2$ChromNum, "x" = mouseBed2$StartPos,
                                        "y" = mouseAmplicons2[[j]])
                      tmpName <- j
                      colnames(tmp) <- c("chromosome", "x", "y")
                      rownames(tmp) <- mouseBed2$AmpliconId
                      tmpGaps <- findLargeGaps(tmp, minLength = 1.5e+07)
                      tmpKnownSegs <- gapsToSegments(tmpGaps)
                      tmpFit <- segmentByCBS(tmp, undo = 1, knownSegments = tmpKnownSegs,
                                             p.method = c("hybrid"))
                      tmpSegRes <- tmpFit$output
                      tmpSegRes$sampleName <- tmpName
                      tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$mean)),]
                      tmpSegRes2 <- clusterSegG_twostep(tmpSegRes)
                      tmpSegRes2
                    }
  stopCluster(cl)
  colnames(allSeg) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
  tmpZScore <- calcZscore(allSeg)
  segZfilt <- segZscoresFilt(allSeg, tmpZScore)
  segZfilt$ID <- nameStripper(segZfilt$ID)
  segZfilt$ID <- str_remove(segZfilt$ID, "x.*")
  
  write.table(segZfilt, file = paste0(mainDir, i,"/pbsSegResPrune20G_2step.txt"), sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  allReportSegs <- rbind(allReportSegs, segZfilt)
}




### pbsSeg graphs for prune 15


allPbs <- NULL
mainDir <- c("/mnt/DATA6/mouseData/copynumber/")
# listOfReports <- c("Auto_user_AUS5-120-MG_EFD4_BBN_334_304", "Auto_user_AUS5-138-MG_cho_20210621_354_343",
#                    "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
#                    "Auto_user_AUS5-156-MG_Fearon_20210809_374_382", "Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384",
#                    "Auto_user_AUS5-76-MG_test1_255_185")

listOfReports <- c("Auto_user_AUS5-120-MG_EFD4_BBN_334_304", "Auto_user_AUS5-138-MG_cho_20210621_354_343",
                   "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                   "Auto_user_AUS5-156-MG_Fearon_20210809_374_382", "Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384",
                   "Reanalysis_AUS5-76-MG_test1_217")

for (i in listOfReports) {
  tmpTable <- fread(paste0(mainDir, i, "/pbsSegResPrune15.txt"))
  tmpTable$dir <- i 
  allPbs <- rbind(allPbs, tmpTable)
}

additionalNamesDf <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20201207annotations.xlsx") 
additionalNamesDf$mg_id <- paste0("mg", additionalNamesDf$mg_id)
additionalNamesDf$stripped <- nameStripper(additionalNamesDf$old_name)


tcDf <- read.table("/mnt/DATA5/tmp/kev/misc/20220521allSampsDf.txt", sep = "\t",
                   header = TRUE, stringsAsFactors = FALSE)
non_zero_20 <- read.table("/mnt/DATA5/tmp/kev/misc/20211203mm10BafsNonzero20.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

allPbs$seg.mean[which(abs(allPbs$seg.mean) < 0.2)] <- 0
allPbs$seg.mean[which(allPbs$q.val2 > 0.05)] <- 0
allPbs$seg.mean[which(allPbs$seg.mean < -3)] <- -3
allPbs$seg.mean[which(allPbs$seg.mean > 3)] <- 3
allPbs$chrom <- str_replace_all(allPbs$chrom, "23", "20")

i <- "1628lt"
for (i in unique(allPbs$ID)) {
  print(i)
  
  
  matchingTc <- tcDf$tc[match(tolower(i), tcDf$sample)]
  
  # sampleName <- i
  # if (grepl("mg", sampleName)) {
  #   sampleName <- additionalNamesDf$stripped[which(additionalNamesDf$mg_id == i)]
  # }
  # 
  testDf_cn <- allPbs[which(allPbs$ID == i),]
  b_2 <- freqPlot_cn(testDf_cn, main = paste(i, "MousePanel, bafs tc:", matchingTc))
  
  
  testDf <- non_zero_20[which(non_zero_20$sample == i),]
  a_20 <- freqPlot_baf(testDf)
  
  g20 <- ggplotGrob(a_20)
  gB_2 <- ggplotGrob(b_2)
  
  dev.off()
  png(filename = paste0("/mnt/DATA5/tmp/kev/misc/pbsSegBafsPrune15/pbsSegBafs_", i, ".png"), width = 1800, height = 500)
  grid::grid.newpage()
  grid::grid.draw(rbind(gB_2, g20))
  dev.off()
  
  pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/pbsSegBafsPrune15/pbsSegBafs_", i, ".pdf"), width = 15, height = 7)
  grid::grid.newpage()
  grid::grid.draw(rbind(gB_2, g20))
  dev.off()
  
  
}


