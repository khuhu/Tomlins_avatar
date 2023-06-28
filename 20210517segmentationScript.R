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
              help="normal file")
); 



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 4){
  print_help(opt_parser)
  stop("Need all 3 arguments to run script =)", call.=FALSE)
}

# loading functions

ampSeg <- function(mouseAmplicons, mouseBedFile, outDir){
  mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
  mouseAmplicons2 <- cbind(mouseBedFile$V1[mouseBedIdx], mouseBedFile$V2[mouseBedIdx], mouseAmplicons)
  mouseAmplicons2[,6:(ncol(mouseAmplicons2)-14)] <- log2(mouseAmplicons2[,6:(ncol(mouseAmplicons2)-14)])
  mouseAmplicons2$`mouseBedFile$V1[mouseBedIdx]` <- as.numeric(mouseAmplicons2$`mouseBedFile$V1[mouseBedIdx]`)
  
  
  allMouseCalls <- colnames(mouseAmplicons2)[6:(ncol(mouseAmplicons2)-14)]
  segResults <- NULL
  #for (i in allMouseCalls) {
  segResults <- foreach(i=allMouseCalls,
                        .combine = 'rbind', .packages = c('DNAcopy')) %dopar% {
    tmpCNA_obj <- CNA(cbind(mouseAmplicons2[[i]]),
                      mouseAmplicons2$ChromNum, mouseAmplicons2$`mouseBedFile$V2[mouseBedIdx]`,
                      data.type="logratio",sampleid=i)
    smoothed_tmpCNA <- smooth.CNA(tmpCNA_obj)
    
    #print(paste0(outDir,"segPlot.",i, ".pdf"))
    segment_tmpCNA <- segment(smoothed_tmpCNA, verbose = 1, undo.splits = "sdundo", undo.SD = 0.2, alpha = 0.05,
                              p.method = c("hybrid"))
    
    pdf(paste0(outDir,"segPlot.",i, ".pdf"), useDingbats = FALSE)
    plot(segment_tmpCNA)
    dev.off()
    
    tmpCalls <- segments.p(segment_tmpCNA)
    # 20210520: got rid of all filtering, can be done downstream
    
    #tmpCalls_filt <- tmpCalls[which(abs(tmpCalls$seg.mean) > .2),]
    #segResults <- rbind(segResults, tmpCalls_filt)
    
    #segResults <- rbind(segResults, tmpCalls)
    tmpCalls
  }
  return(segResults)
}


calcZscore <- function(segRes){
  z_vector <- NULL
  nmean_vector <- NULL
  tmean_vector <- NULL
  for (i in unique(segRes$ID)) {
    tmp_seg_res <- segRes[which(segRes$ID == i), ]
    tmp_seg_ranges <- GRanges(seqnames = Rle(tmp_seg_res$chrom), 
                              ranges = IRanges(start = tmp_seg_res$loc.start,
                                               end = tmp_seg_res$loc.end))
    for (j in seq_along(tmp_seg_ranges)) {
      tmpOverlap <- findOverlaps(all_probes_grange, tmp_seg_ranges[j])
      tmpTest <- zscore_gc_oe_ratios[[i]][queryHits(tmpOverlap)]
      
      if (!is.numeric(tmpTest)) {
        z_vector <- c(z_vector, NA)
        nmean_vector <- c(nmean_vector, NA)
        tmean_vector <- c(tmean_vector, NA)
        next()
      } else{
        tmpNormal <- zscore_gc_oe_ratios[queryHits(tmpOverlap), mouseNormal]
        tmpNormal2 <- unlist(tmpNormal)
        
        z_stat <- (mean(tmpTest) - mean(tmpNormal2))/sd(tmpNormal2)
        normal_seg_mean <- mean(tmpNormal2)
        tumor_seg_mean <- mean(tmpTest)
        z_vector <- c(z_vector, z_stat)
        nmean_vector <- c(nmean_vector, normal_seg_mean)
        tmean_vector <- c(tmean_vector, tumor_seg_mean)
      }
    }
  }
  return(list("z_vector" = z_vector, "nmean_vector"= nmean_vector,
              "tmean_vector" = tmean_vector))
}

segZscoresFilt <- function(segResults, zscores){
  
  # filts for significance, cn-change and length
  segmental_zscores <- cbind(segResults, "z-scores" = zscores[["z_vector"]],
                             "normal_seg_mean" = zscores[["nmean_vector"]],
                             "tumor_seg_mean" = zscores[["tmean_vector"]])
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

cl <- makeCluster(4)
registerDoParallel(cl)

outDir <- str_remove(opt$out,"segResults.txt")

mouseNormal <- normalFile$V1
segRes <- ampSeg(mouseAmplicons, mouseBedFile, outDir)
segZscores <- calcZscore(segRes)
segZfilt <- segZscoresFilt(segRes, segZscores)

# print out filtering cutoffs
write.table(segZfilt, opt$out, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

stopCluster(cl)

warnings()
