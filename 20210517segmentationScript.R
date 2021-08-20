library(DNAcopy)
library(GenomicRanges)
library(optparse)
library(stringr)
library(foreach)
library(doParallel)

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
                 "Length", "GC", "TotalPool", "Weights")
  mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
  mouseBed2 <- data.frame("AmpliconId" = mouseBedIdx, "ChromNum" = mouseBedFile$V1[mouseBedIdx],
                          "StartPos" = mouseBedFile$V2[mouseBedIdx])
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
                          segment_tmpCNA <- segment(smoothed_tmpCNA, verbose = 1, undo.splits = "sdundo", undo.SD = 0.2, alpha = 0.05,
                                                    p.method = c("hybrid"), min.width = minAmp)
                          
                          # var just to plot - so no crazy scaling
                          segToPlot <- segment_tmpCNA
                          segToPlot$output$seg.mean[which(segToPlot$output$seg.mean < -4)] <- -4
                          segToPlot$output$seg.mean[which(segToPlot$output$seg.mean > 4)] <- 4
                          segToPlot$data[[3]][which(segToPlot$data[[3]] < -4)] <- -4
                          segToPlot$data[[3]][which(segToPlot$data[[3]] > 4)] <- 4
                          
                          pdf(paste0(outDir,"segPlot.",i, ".pdf"), useDingbats = FALSE)
                          plot(segToPlot, ylim=c(-4,4))
                          dev.off()
                          
                          tmpCalls <- segments.p(segment_tmpCNA)
                          tmpCalls
                        }
  return(segResults)
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

cl <- makeCluster(5)
registerDoParallel(cl)

outDir <- str_remove(opt$out,"segResults.txt")

mouseNormal <- normalFile$V1
paramFile <- read.table(opt$param)
segRes <- ampSeg(mouseAmplicons, mouseBedFile, outDir, paramFile$V1)
segRes$ID <- str_remove(segRes$ID, "^X")
segZscores <- calcZscore(segRes)
segZfilt <- segZscoresFilt(segRes, segZscores)

# print out filtering cutoffs
write.table(segZfilt, opt$out, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

stopCluster(cl)

warnings()
