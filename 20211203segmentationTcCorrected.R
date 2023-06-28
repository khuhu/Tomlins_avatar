### need to just scan and create new bed file where on things with cnAmplicon.txt are made

library(DNAcopy)
library(GenomicRanges)
library(optparse)
library(stringr)
library(foreach)
library(doParallel)

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
  
  for (i in unique(df_cn$chrom)) {
    df_cn$loc.start[which(df_cn$chrom == i)] <- df_cn$loc.start[which(df_cn$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    df_cn$loc.end[which(df_cn$chrom == i)] <- df_cn$loc.end[which(df_cn$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  df_cn <- df_cn[,c("chrom", "loc.start", "loc.end", "seg.mean", "col")]
  colnames(df_cn) <- c("chrom", "start", "end","cn", "col")
  
  
  ggplot(df_cn) + geom_hline(yintercept = c(-3, -2 , -1 , 0, 1, 2, 3), color = "#D4D4D4") +
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

option_list = list(
  make_option(c("-b", "--bed"), type="character", default=NULL, 
              help="bed file", metavar="character"),
  make_option(c("-a", "--amplicon"), type="character", default=NULL, 
              help="copy number amplicon matrix"),
  make_option(c("-g", "--gc"), type="character", default=NULL, 
              help="gc corrected count matrix"),
  make_option(c("-n", "--normal"), type="character", default=NULL, 
              help="normal file"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="normal file"),
  make_option(c("-p", "--param"), type="character", default=NULL, 
              help="segmentation parameter file")
); 



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 4){
  print_help(opt_parser)
  stop("Need all 3 arguments to run script =)", call.=FALSE)
}

# loading functions

ampSeg <- function(mouseAmplicons, mouseBedFile, outDir, minAmp){ 
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
  
  allMouseCalls <- colnames(mouseAmplicons2)[2:(ncol(mouseAmplicons2))]
  
  segResults <- NULL
  segResults <- foreach(i=allMouseCalls,
                        .combine = 'rbind', .packages = c('DNAcopy')) %dopar% {
                          tmpCNA_obj <- CNA(cbind(mouseAmplicons2[[i]]),
                                            mouseBed2$ChromNum, mouseBed2$StartPos,
                                            data.type="logratio",sampleid=i)
                          smoothed_tmpCNA <- smooth.CNA(tmpCNA_obj)
                          
                          #print(paste0(outDir,"segPlot.",i, ".pdf"))
                          # segment_tmpCNA <- segment(smoothed_tmpCNA, verbose = 1, undo.splits = "sdundo", undo.SD = 0.2, alpha = 0.05,
                          #                           p.method = c("hybrid"), min.width = minAmp)
                          # 
                          
                          
                          segment_tmpCNA <- segment(smoothed_tmpCNA, verbose = 1, undo.splits = "sdundo",
                                                    p.method = c("hybrid"), undo.SD = 1, min.width = minAmp)
                          
                          
                          # var just to plot - so no crazy scaling
                          # segToPlot <- segment_tmpCNA
                          # segToPlot$output$seg.mean[which(segToPlot$output$seg.mean < -4)] <- -4
                          # segToPlot$output$seg.mean[which(segToPlot$output$seg.mean > 4)] <- 4
                          # segToPlot$data[[3]][which(segToPlot$data[[3]] < -4)] <- -4
                          # segToPlot$data[[3]][which(segToPlot$data[[3]] > 4)] <- 4
                          # 
                          # pdf(paste0(outDir,"segPlot_tc.",i, ".pdf"), useDingbats = FALSE)
                          # plot(segToPlot, ylim=c(-4,4))
                          # dev.off()
                          
                          tmpCalls <- segments.p(segment_tmpCNA)
                          tmpCalls
                        }
  return(list(segResults, mouseAmplicons, mouseBed2))
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
  
  uniqIdx <- unique(match(str_replace_all(segRes$ID, "[[:punct:]]", ""), 
                          str_replace_all(str_replace_all(colnames(zscore_gc_oe_ratios), "[[:punct:]]", ""), " ", "")))
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

segZscoresFilt <- function(segResults, zscores){
  
  # filts for significance, cn-change and length
  segmental_zscores <- cbind(segResults, "z-scores1" = zscores[["z1_vector"]], "z-scores2" = zscores[["z2_vector"]],
                             "p-val1" = zscores[["p1_vector"]], "p-val2" = zscores[["p2_vector"]],
                             "q-val1" = zscores[["q1_vector"]], "q-val2" = zscores[["q2_vector"]],
                             "normal_seg_mean" = zscores[["nmean_vector"]], "tumor_seg_mean" = zscores[["tmean_vector"]])
  return(segmental_zscores)
}


nameStripper <- function(df){
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
}

# loading tables
mouseBedFile <- read.table(opt$bed, stringsAsFactors = FALSE, sep = "\t",
                           header = FALSE)

mouseAmplicons <- read.table(opt$amplicon, sep = "\t", stringsAsFactors = FALSE,
                             header = TRUE, check.names = FALSE)

zscore_gc_oe_ratios <- read.table(opt$gc, sep = "\t", stringsAsFactors = FALSE,
                                  header = TRUE, check.names = FALSE)

normalFile <- read.table(opt$normal, sep = "\t", stringsAsFactors = FALSE,
                         header = FALSE)

all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                             ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                              end = zscore_gc_oe_ratios$EndPos))



### 20211025 KH: add tc correction prior to segmentation on ratios 
### non_zero$sample <- nameStripper(non_zero$sample)
### non_zero$sample <- str_remove(str_remove(non_zero$sample, "x.*"), "-")

# tcDf <- read.table("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt", sep = "\t", stringsAsFactors = FALSE,
#                    header = TRUE)


###

cl <- makeCluster(5)
registerDoParallel(cl)

outDir <- str_remove(opt$out,"segResults_tc.txt")

mouseNormal <- normalFile$V1
paramFile <- read.table(opt$param)
#segRes <- ampSeg(mouseAmplicons, mouseBedFile, outDir)
segRes <- ampSeg(mouseAmplicons, mouseBedFile, outDir, 5)
segRes2 <- segRes[[1]]
segRes2$ID <- nameStripper(segRes2$ID)
segRes2$ID <- str_remove(segRes2$ID, "^X")



### 20211203: snippet for correct tcDf - this makes more sense because I'm correct aggregate probes
tcDf <- read.table("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt", header = TRUE, stringsAsFactors = FALSE)
matchingTc <- tcDf$tc[match(segRes2$ID, tcDf$sample)]
matchingTcName <- tcDf$sample[match(segRes2$ID, tcDf$sample)]

### 20211203: also need to change the names in zscore_ratio variables + normal names
nonNameVars <- c ("AmpliconId", "AmpliconIndex", "ChromNum", "StartPos", "EndPos",
                  "Gene", "NumGC", "Length", "GC", "TotalPool", "Weights")
colnames(zscore_gc_oe_ratios)[-which(colnames(zscore_gc_oe_ratios) %in% nonNameVars)] <- nameStripper(colnames(zscore_gc_oe_ratios)[-which(colnames(zscore_gc_oe_ratios) %in% nonNameVars)])
colnames(zscore_gc_oe_ratios)[-which(colnames(zscore_gc_oe_ratios) %in% nonNameVars)] <- str_remove(colnames(zscore_gc_oe_ratios)[-which(colnames(zscore_gc_oe_ratios) %in% nonNameVars)] , "^X")

colnames(mouseAmplicons)[-which(colnames(mouseAmplicons) %in% nonNameVars)] <- nameStripper(colnames(mouseAmplicons)[-which(colnames(mouseAmplicons) %in% nonNameVars)])
colnames(mouseAmplicons)[-which(colnames(mouseAmplicons) %in% nonNameVars)] <- str_remove(colnames(mouseAmplicons)[-which(colnames(mouseAmplicons) %in% nonNameVars)] , "^X")


mouseNormal <- nameStripper(mouseNormal)
mouseNormal <- str_remove(mouseNormal, "^X")


mouseAmplicons_Grange <- GRanges(seqnames = mouseAmplicons$ChromNum, 
                                 IRanges(start = mouseAmplicons$StartPos, 
                                         end = mouseAmplicons$EndPos))


segRes2$seg.mean <- 2^segRes2$seg.mean

for (i in which(!is.na(matchingTc))) {
  tmpGrange  <- GRanges(seqnames = segRes2$chrom[i], 
                        IRanges(start = segRes2$loc.start[i], 
                                end = segRes2$loc.end[i]))
  tmpVec <- segRes2$seg.mean[i]
  
  if (segRes2$seg.mean[i] > 2^0.2) {
    tmpVec <- tmpVec / matchingTc[i]
    tmpIdx <- queryHits(findOverlaps(mouseAmplicons_Grange, tmpGrange))
    mouseAmplicons[[matchingTcName[i]]][tmpIdx] <- mouseAmplicons[[matchingTcName[i]]][tmpIdx] / matchingTc[i]
  } else if(segRes2$seg.mean[i] < 2^-0.2){
    tmpVec <- tmpVec * matchingTc[i]
    tmpIdx <- queryHits(findOverlaps(mouseAmplicons_Grange, tmpGrange))
    mouseAmplicons[[matchingTcName[i]]][tmpIdx] <- mouseAmplicons[[matchingTcName[i]]][tmpIdx] * matchingTc[i]
  }
  segRes2$seg.mean[i] <- tmpVec
}

segRes2$seg.mean <- log2(segRes2$seg.mean)



if(length(grep("\\.1", segRes2$ID)) > 0){
  segRes2 <- segRes2[-grep("\\.1", segRes2$ID),]
}

segZscores <- calcZscore(segRes2)
segZfilt <- segZscoresFilt(segRes2, segZscores)

#mouseAmplicons <- segRes[[2]]
mouseBed <- segRes[[3]]

segZfilt$chrom <- str_replace_all(segZfilt$chrom, "23", "20")

# print out filtering cutoffs
write.table(segZfilt, opt$out, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
write.table(mouseAmplicons, paste0(outDir, "mouseAmplicons_tc.txt"), sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
write.table(mouseBed, paste0(outDir, "mouseProbeLoc.txt"), sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)


stopCluster(cl)

warnings()

### creating segPlots

segZfilt$seg.mean[which(segZfilt$seg.mean > 3)] <- 3
segZfilt$seg.mean[which(segZfilt$seg.mean < -3)] <- -3

for (i in unique(segZfilt$ID)) {
  
  matchingTc <- signif(tcDf$tc[match(tolower(i), tcDf$sample)], digits = 2)
  testDf_cn <- segZfilt[which(segZfilt$ID == i),]
  
  g1 <- freqPlot_cn(testDf_cn, main = paste(i, "MousePanel tc: ", matchingTc))
  grob1 <- ggplotGrob(g1)
  
  pdf(paste0(outDir,"segPlot_tc.",i, ".pdf"), useDingbats = FALSE, height = 5, width = 10)
  grid::grid.newpage()
  grid::grid.draw(rbind(grob1))
  dev.off()
  
}


