library(data.table)
library(foreach)
library(doParallel)
library(copynumber)
library(GenomicFeatures)
library(dplyr)
library(dbplyr)
library(circlize)
library(stringr)
library(DNAcopy)
library(stringr)
library(data.table)
library(philentropy)

### already do segmentaiton, but can reseg to correct for tumor content
ampSeg <- function(mouseAmplicons, mouseBedFile, outDir){
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
                                                    p.method = c("hybrid"))
                          
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
  z_vector <- NULL
  nmean_vector <- NULL
  tmean_vector <- NULL
  
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
    
    for (j in seq_along(tmp_seg_ranges)) {
      tmpOverlap <- findOverlaps(all_probes_grange, tmp_seg_ranges[j])
      tmpTest <- zscore_gc_oe_ratios[[zScoreColnames[i]]][queryHits(tmpOverlap)]
      
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

# depreceated 
# segZscoresFilt <- function(segResults, zscores){
#   
#   # filts for significance, cn-change and length
#   segmental_zscores <- cbind(segResults, "z-scores" = zscores[["z_vector"]],
#                              "normal_seg_mean" = zscores[["nmean_vector"]],
#                              "tumor_seg_mean" = zscores[["tmean_vector"]])
#   segmental_zscores2 <- segmental_zscores[which(segmental_zscores$num.mark > 10),]
#   segmental_zscores2$length <- segmental_zscores2$loc.end - segmental_zscores2$loc.start
#   largeChromCutoff <- 1000000
#   segmental_zscores3 <- segmental_zscores2[which(segmental_zscores2$length >= largeChromCutoff),]
#   segmental_zscores3 <- segmental_zscores3[which(abs(segmental_zscores3$`z-scores`) > 1.645),]
#   segInput <- segmental_zscores3[,c("ID", "chrom", "loc.start", "loc.end", "seg.mean")]
#   return(segInput)
# }

### this filt shold be used since it zeroes out instead of filtering out
### problem for copynumber getFreqData

segZscoresFilt_zeroOut <- function(segResults){
  segmental_zscores <- segResults
  # filts for significance, cn-change and length
  # segmental_zscores <- cbind(segResults, "z-scores" = zscores[["z_vector"]],
  #                            "normal_seg_mean" = zscores[["nmean_vector"]],
  #                            "tumor_seg_mean" = zscores[["tmean_vector"]])
  segmental_zscores2 <- segmental_zscores
  #segmental_zscores2$seg.mean[which(segmental_zscores2$num.mark < 10)] <- 0
  segmental_zscores2$length <- segmental_zscores2$loc.end - segmental_zscores2$loc.start
  #largeChromCutoff <- 1000000
  #segmental_zscores2$seg.mean[which(segmental_zscores2$length <= largeChromCutoff)] <- 0
  #segmental_zscores2$seg.mean[which(abs(segmental_zscores2$`z-scores`) < 1.645)] <- 0
  # this cutoff accounts for triploid
  segmental_zscores2$seg.mean[-which(segmental_zscores2$seg.mean < -0.66 | segmental_zscores2$seg.mean > 0.268)] <- 0
  segInput <- segmental_zscores2[,c("ID", "chrom", "loc.start", "loc.end", "seg.mean")]
  return(segInput)
}

#  20210705: depreceated for V2

# syntenyPlotInputs <- function(segInput){
#   segInput$chrom <- str_replace_all(segInput$chrom, "23", "X")
#   tumor_freq_ranges <- GRanges(seqnames = paste0("chr", segInput$chrom),
#                                ranges = IRanges(start = as.numeric(segInput$loc.start),
#                                                 end = as.numeric(segInput$loc.end)))
#   synteny_ranges <- GRanges(seqnames = paste0("chr", synteny_hg38mm10$comp_chr),
#                             ranges = IRanges(start = synteny_hg38mm10$comp_start_pos,
#                                              end = synteny_hg38mm10$comp_end_pos))
#   
#   
#   res_overlap <- findOverlaps(synteny_ranges, tumor_freq_ranges)
#   query_hits <- queryHits(res_overlap)
#   subject_hits <- subjectHits(res_overlap)
#   
#   
#   mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg38mm10$comp_chr[query_hits]),
#                           "start" = synteny_hg38mm10$comp_start_pos[query_hits],
#                           "end" = synteny_hg38mm10$comp_end_pos[query_hits],
#                           stringsAsFactors = FALSE)
#   
#   human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg38mm10$ref_chr[query_hits]),
#                           "start" = synteny_hg38mm10$ref_start_pos[query_hits],
#                           "end" = synteny_hg38mm10$ref_end_pos[query_hits],
#                           stringsAsFactors = FALSE)
#   
#   
#   
#   resList <- list("human_bed" = human_bed,
#                   "mouse_bed" = mouse_bed)
#   
#   return(resList)
# }
# 
# 
# syntenyPlotInputsFreq <- function(segInput){
#   segInput$chrom <- str_replace_all(segInput$chrom, "23", "X")
#   tumor_freq_ranges <- GRanges(seqnames = paste0("chr", segInput$chrom),
#                                ranges = IRanges(start = as.numeric(segInput$loc.start),
#                                                 end = as.numeric(segInput$loc.end)))
#   synteny_ranges <- GRanges(seqnames = paste0("chr", synteny_hg38mm10$comp_chr),
#                             ranges = IRanges(start = synteny_hg38mm10$comp_start_pos,
#                                              end = synteny_hg38mm10$comp_end_pos))
#   
#   
#   res_overlap <- findOverlaps(synteny_ranges, tumor_freq_ranges)
#   query_hits <- queryHits(res_overlap)
#   subject_hits <- subjectHits(res_overlap)
#   
#   
#   mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg38mm10$comp_chr[query_hits]),
#                           "start" = synteny_hg38mm10$comp_start_pos[query_hits],
#                           "end" = synteny_hg38mm10$comp_end_pos[query_hits],
#                           "freq" = segInput$seg.mean[subject_hits],
#                           stringsAsFactors = FALSE)
#   
#   human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg38mm10$ref_chr[query_hits]),
#                           "start" = synteny_hg38mm10$ref_start_pos[query_hits],
#                           "end" = synteny_hg38mm10$ref_end_pos[query_hits],
#                           stringsAsFactors = FALSE)
#   
#   
#   
#   resList <- list("human_bed" = human_bed,
#                   "mouse_bed" = mouse_bed)
#   
#   return(resList)
# }

syntenyPlotInputsFreqV2 <- function(segInput){
  # made for multi-sample freq
  tumor_freq_ranges <- GRanges(seqnames = paste0("chr", segInput$Chr),
                               ranges = IRanges(start = as.numeric(segInput$Start),
                                                end = as.numeric(segInput$End)))
  synteny_ranges <- GRanges(seqnames = paste0("chr", synteny_hg38mm10$comp_chr),
                            ranges = IRanges(start = synteny_hg38mm10$comp_start_pos,
                                             end = synteny_hg38mm10$comp_end_pos))
  
  res_overlap <- findOverlaps(synteny_ranges, tumor_freq_ranges)
  query_hits <- queryHits(res_overlap)
  subject_hits <- subjectHits(res_overlap)
  
  
  mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg38mm10$comp_chr[query_hits]),
                          "start" = synteny_hg38mm10$comp_start_pos[query_hits],
                          "end" = synteny_hg38mm10$comp_end_pos[query_hits],
                          "freq" = segInput$Freq[subject_hits],
                          stringsAsFactors = FALSE)
  
  human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg38mm10$ref_chr[query_hits]),
                          "start" = synteny_hg38mm10$ref_start_pos[query_hits],
                          "end" = synteny_hg38mm10$ref_end_pos[query_hits],
                          stringsAsFactors = FALSE)
  
  
  
  resList <- list("human_bed" = human_bed,
                  "mouse_bed" = mouse_bed)
  
  return(resList)
}

getFreqData <- function(data){
  require(copynumber)
  #Check if segments or data:
  if(colnames(data)[1]=="sampleID" || colnames(data)[2]=="arm"){
    #input is segments data frame;
    #could be on a multiseg-format -> convert to uniseg-format:
    #need to convert to appropriate format
    #first find intersection of all breakpoints:
    chr <- unique(data[,2])
    bpts <- matrix(NA,nrow=0,ncol=2)
    for(j in 1:length(chr)){
      subseg <- subsetSegments(data,chrom=chr[j])
      bpts <- rbind(bpts,data.frame(chr[j],sort(unique(c(subseg$start.pos,subseg$end.pos)))))
    } 
    colnames(bpts) <- c("chrom","pos")
    
    #Then interpolate to get pcf-val in all breakpoints
    data = interpolate.pcf(data, bpts) 
  }
  #If not segments, the input data should already be on an appropriate format (chrom, pos, estimates)
  
  return(data) 
} 


downsampleRegions <- function(df){
  finalTbl <- NULL
  for (i in unique(df$chr)) {
    tmpDf <- df[which(df$chr == i),]
    # weird artefact where the last row is nearly entire chromosome
    # for every entry - so removing it
    tmpDf <- tmpDf[1:(nrow(tmpDf) - 1),]
    
    fit <- loess(start ~ cn, degree=1, span = 0.1, data=tmpDf)
    tmpDf$loess <- fit$fitted
    tmpDf_red <- tmpDf[c(maximums(tmpDf$loess),
                         minimums(tmpDf$loess)),]
    tmpDf_red <- tmpDf_red[order(tmpDf_red$start, decreasing = FALSE),]
    redIdx <- c(1,1+which(abs(diff(tmpDf_red$cn)) >= 0.05))
    tmpDf_red2 <- NULL
    for (j in 1:length(redIdx)) {
      if (j == nrow(tmpDf_red)) {
        # if last change point is in last region
        
        tmpIdx <- redIdx[j]
        tmpVec <- c("chr" = i,
                    "start" = tmpDf_red$start[tmpIdx],
                    "end" = tmpDf_red$end[tmpIdx], 
                    "cn" = round(tmpDf_red$cn[tmpIdx],digits = 2))
      } else if (j == length(redIdx)){
        # if last change point is not in last region
        
        tmpIdx <- redIdx[j]
        tmpVec <- c("chr" = i,
                    "start" = tmpDf_red$start[tmpIdx],
                    "end" = tmpDf_red$end[nrow(tmpDf_red)], 
                    "cn" = round(mean(c(tmpDf_red$cn[tmpIdx],
                                        tmpDf_red$cn[nrow(tmpDf_red)])), digits = 2))
      } else {
        tmpIdx <- redIdx[j]
        tmpIdx2 <- redIdx[j + 1]
        tmpVec <- c("chr" = i,
                    "start" = tmpDf_red$start[tmpIdx],
                    "end" = tmpDf_red$start[tmpIdx2] - 1, 
                    "cn" = round(mean(tmpDf_red$cn[tmpIdx:tmpIdx2]),digits = 2))
      }
      tmpDf_red2 <- rbind(tmpDf_red2, tmpVec)
    }
    
    tmpDf_red2 <- data.frame(tmpDf_red2, stringsAsFactors = FALSE)
    tmpDf_red2[,2:4] <- lapply(tmpDf_red2[,2:4], as.numeric)
    
    finalTbl <- rbind(finalTbl, tmpDf_red2)
  }
  print("done")
  return(finalTbl)
}



getFreqBed <- function(amp, del){
  # modified original code for brevity
  amp2 <- amp
  del2 <- del
  
  amp2$Freq <- round(amp2$Freq/100, digits = 2)
  amp2$Chr <- paste0("h_chr", amp2$Chr)
  colnames(amp2) <- c("chr","start","end","cn")
  
  amp3 <- downsampleRegions(amp2)
  amp3$ybot <- 0
  amp3$ytop <- amp3$cn
  
  del2$Freq <- round(del2$Freq/100, digits = 2)
  del2$Chr <- paste0("h_chr", del2$Chr)
  colnames(del2) <- c("chr","start","end","cn")
  del2$cn <- -del2$cn
  
  del3 <- downsampleRegions(del2)
  del3$ytop <- 0
  del3$ybot <- del3$cn
  
  cn_track_bed <- rbind(amp3, del3)
  
  return(cn_track_bed)
}


reducingFreqBed <- function(df, idx){
  
  ### regions with zero frequency may have wrong positions, but it get's filtered out later
  ### script correctly obtains regions with same freq
  reducedDf <- NULL
  for (i in 1:(length(idx) - 1)) {
    if (length(idx) == 1) {
      reducedDf <- data.frame("Chr" = df$chr[1], "Start" = df$pos[1], "End" = df$pos[2], "Freq" = 0)
      return(reducedDf)
    } else{
      idx1 <- idx[i]
      idx2 <- idx[i+1] - 1
      tmpDf <- df[idx1:idx2,]
      tmpChr <- tmpDf$chr[1]
      tmpStart <- min(tmpDf$pos)
      tmpEnd <- max(tmpDf$pos) - 1
      tmpFreq <- tmpDf[,3][1]
      tmpVec <- c("Chr" = tmpChr, "Start" = tmpStart, "End" = tmpEnd, "Freq" = tmpFreq)
      reducedDf <- rbind(reducedDf, tmpVec)
    }
  }
  reducedDf <- data.frame(reducedDf, stringsAsFactors = FALSE)
  return(reducedDf)
}

ampsDels <- function(df){
  #getFreqOut <- getFreqData(df)
  # 20210629: for whatever reason cant find the mean of a long logical vector
  # i.e only TRUE AND FALSE
  #freqAmp <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] > 0.2)*100
  #freqDel <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] < -0.2)*100
  freqAmp <- apply(df[,-c(1:2)], 1, function(x) 
    length(which(x > 0.2))/length(x)) * 100
  freqDel <- apply(df[,-c(1:2)], 1, function(x) 
    length(which(x < -0.2))/length(x)) * 100
  
  freqDf_Amp <- cbind(df[,1:2], "amp" = freqAmp)
  freqDf_Del <- cbind(df[,1:2], "del" = freqDel)
  
  ampRedIdx <- c(1,1+which(diff(freqDf_Amp$amp)!=0), length(freqDf_Amp$amp))
  delRedIdx <- c(1,1+which(diff(freqDf_Del$del)!=0), length(freqDf_Del$del))
  
  res <- list(freqDf_Amp, ampRedIdx, freqDf_Del, delRedIdx)
}

### use to get aneuploid events and large CNAs

# separateSegments_tcga <- function(df, ploidyDf){
#   # separates the segments prior to generating frequencies - diff for tcga data and amplicon
#   # whole chromosome/arms (covers ~80% of max chr coverage from panel); changes > 1MB
#   res <- NULL
#   res2 <- NULL
#   
#   # this was genius and got rid of the loop - literally infinitely faster *pats self on back*
#   df$newTotalCN <- round2(df$Copynumber - ploidyDf$ploidy[match(df$Sample, ploidyDf$array)], 0)
#   
#   # 20210717 new method too many loops so going to parallelize it and
#   # have a specific c
#   comb <- function(x, ...) {
#     lapply(seq_along(x),
#            function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
#   }
#   
#   cl <- makeCluster(20)
#   registerDoParallel(cl)
#   
#   oper <- foreach(i=unique(df$Sample), .combine='comb', .multicombine=TRUE,
#                   .init=list(list(), list())) %dopar% {
#                     res <- NULL
#                     res2 <- NULL
#                     sampDf <- df[which(df$Sample == i),]
#                     human_arm <- read.table("/mnt/DATA5/tmp/kev/misc/20210713human_arm_syn.txt",
#                                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)
#                     for (j in unique(human_arm$chr)) {
#                       df_chr <- sampDf[sampDf$Chromosome == j,]
#                       minLength <- human_arm$length80[which(human_arm$chr == j)]
#                       armTable <- df_chr
#                       cnaTable <- df_chr
#                       
#                       # new loop below gets segments of gains or losses
#                       signChange <- c(0, diff(sign(df_chr$newTotalCN)))
#                       idxChange <- which(signChange != 0)
#                       
#                       if(length(idxChange) == 0){
#                         res <- rbind(res, armTable)
#                         res2 <- rbind(res2, cnaTable)
#                         next()
#                       }
#                       
#                       for (k in 1:(length(idxChange))){
#                         if (df_chr$newTotalCN[idxChange[k]] == 0) {
#                           next()
#                         } else if (k == length(idxChange)){
#                           idx1 <- idxChange[k]
#                           idx2 <- length(df_chr$newTotalCN)
#                         } else{
#                           idx1 <- idxChange[k]
#                           idx2 <- idxChange[k + 1] - 1
#                         }
#                         
#                         tmpLength <- sum(df_chr$length[idx1:idx2])
#                         if(tmpLength > minLength){
#                           cnaTable$newTotalCN[idx1:idx2] <- 0
#                         } else{
#                           armTable$newTotalCN[idx1:idx2] <- 0
#                         }
#                       }
#                       res <- rbind(res, armTable)
#                       res2 <- rbind(res2, cnaTable)
#                     }
#                     list(res, res2)
#                   }
#   stopCluster(cl)
#   return(oper)
#   
#   # for (i in unique(df$Sample)) {
#   #   sampDf <- df[which(df$Sample == i),]
#   #   for (j in unique(human_arm$chr)) {
#   #     df_chr <- sampDf[sampDf$Chromosome == j,]
#   #     minLength <- human_arm$length80[which(human_arm$chr == j)]
#   #     armTable <- df_chr
#   #     cnaTable <- df_chr
#   #     
#   #     # new loop below gets segments of gains or losses
#   #     signChange <- c(0, diff(sign(df_chr$newTotalCN)))
#   #     idxChange <- which(signChange != 0)
#   #     
#   #     if(length(idxChange) == 0){
#   #       res <- rbind(res, armTable)
#   #       res2 <- rbind(res2, cnaTable)
#   #       next()
#   #     }
#   #     
#   #     for (k in 1:(length(idxChange))){
#   #       if (df_chr$newTotalCN[idxChange[k]] == 0) {
#   #         next()
#   #       } else if (k == length(idxChange)){
#   #         idx1 <- idxChange[k]
#   #         idx2 <- length(df_chr$newTotalCN)
#   #       } else{
#   #         idx1 <- idxChange[k]
#   #         idx2 <- idxChange[k + 1] - 1
#   #       }
#   #       
#   #       tmpLength <- sum(df_chr$length[idx1:idx2])
#   #       if(tmpLength > minLength){
#   #         cnaTable$newTotalCN[idx1:idx2] <- 0
#   #       } else{
#   #         armTable$newTotalCN[idx1:idx2] <- 0
#   #       }
#   #     }
#   #     res <- rbind(res, armTable)
#   #     res2 <- rbind(res2, cnaTable)
#   #   }
#   # }
#   # return(list(res, res2))
#   
#   # for (i in unique(df$Chromosome)) {
#   #   tmpChrDf <- df[which(df$Chromosome == i),]
#   #   tmpChrDf$newTotalCN[which(tmpChrDf$length < 1000000)] <- 0
#   #   armTable <- tmpChrDf 
#   #   cnaTable <- tmpChrDf
#   #   tmpCyto <- human_arm[which(human_arm$chromosome == i),]
#   #   
#   #   minLength <- min(tmpCyto$length80)
#   #   
#   #   armTable$newTotalCN[which(tmpChrDf$length <= minLength)] <- 0
#   #   cnaTable$newTotalCN[which(tmpChrDf$length > minLength)] <- 0
#   #   
#   #   res <- rbind(res, armTable)
#   #   res2 <- rbind(res2, cnaTable)
#   # }
#   # return(list(res, res2))
# }

# 
# 
# separateSegments_tcgaV2 <- function(df, ploidyDf){
#   # separates the segments prior to generating frequencies - diff for tcga data and amplicon
#   # whole chromosome/arms (covers ~80% of max chr coverage from panel); changes > 1MB
#   res <- NULL
#   res2 <- NULL
# 
#   # this was genius and got rid of the loop - literally infinitely faster *pats self on back*
#   df$newTotalCN <- round2(df$Copynumber - ploidyDf$ploidy[match(df$Sample, ploidyDf$array)], 0)
# 
#   # 20210717 new method too many loops so going to parallelize it and
#   # have a specific c
#   comb <- function(x, ...) {
#     lapply(seq_along(x),
#            function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
#   }
# 
#   cl <- makeCluster(20)
#   registerDoParallel(cl)
# 
#   oper <- foreach(i=unique(df$Sample), .combine='comb', .multicombine=TRUE,
#                   .init=list(list(), list())) %dopar% {
#                     res <- NULL
#                     res2 <- NULL
#                     sampDf <- df[which(df$Sample == i),]
#                     human_arm <- read.table("/mnt/DATA5/tmp/kev/misc/20210713human_arm_syn.txt",
#                                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)
#                     for (j in unique(human_arm$chr)) {
#                       df_chr <- sampDf[sampDf$Chromosome == j,]
#                       minLength <- min(human_arm$length80[which(human_arm$chr == j)])
#                       aTable <- df_chr
#                       cTable <- df_chr
# 
#                       # new loop below gets segments of gains or losses
#                       delIdx <- which(sign(df_chr$newTotalCN) == -1)
#                       ampIdx <- which(sign(df_chr$newTotalCN) == 1)
# 
#                       if(length(delIdx) == 0 & length(ampIdx) == 0){
#                         res <- rbind(res, aTable)
#                         res2 <- rbind(res2, cTable)
#                         next()
#                       }
# 
#                       if (sum(df_chr$length[delIdx]) < minLength) {
#                         aTable$newTotalCN[delIdx] <- 0
#                       } else{
#                         cTable$newTotalCN[delIdx] <- 0
#                       }
# 
#                       if (sum(df_chr$length[ampIdx]) < minLength) {
#                         aTable$newTotalCN[ampIdx] <- 0
#                       } else{
#                         cTable$newTotalCN[ampIdx] <- 0
#                       }
# 
# 
#                       res <- rbind(res, aTable)
#                       res2 <- rbind(res2, cTable)
# 
#                     }
#                     list(res, res2)
#                   }
#   stopCluster(cl)
#   return(oper)
# }

# give this 4 outputs: arm, cna, count table and fga (do other fga method too)

separateSegments_tcgaV3 <- function(df, ploidyDf){
  
  # this was genius and got rid of the loop - literally infinitely faster *pats self on back*
  df$newTotalCN <- round2(df$Copynumber - ploidyDf$ploidy[match(df$Sample, ploidyDf$array)], 0)
  
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  
  cl <- makeCluster(20)
  registerDoParallel(cl)
  
  oper <- foreach(i=unique(df$Sample), .combine='comb', .multicombine=TRUE,
                  .init=list(list(), list(), list())) %dopar% {
                    res <- NULL
                    res2 <- NULL
                    resMat <- matrix(data = 0, nrow = 1, ncol = 4)
                    
                    sampDf <- df[which(df$Sample == i),]
                    human_arm <- read.table("/mnt/DATA5/tmp/kev/misc/20210730hg19chromSizeTable.txt",
                                            sep = "\t", stringsAsFactors = FALSE, header = TRUE)
                    human_arm <- human_arm[which(human_arm$chrom %in% c(1:22)),]
                    ### counting for anueploidy
                    for (j in unique(human_arm$chrom)) {
                      skipVar <- "no"
                      df_chr <- sampDf[sampDf$Chromosome == j,]
                      df_chr <- df_chr[order(df_chr$Start),]
                      minLength80 <- human_arm$length80[which(human_arm$chrom == j)]
                      minLength <- 15000000
                      aTable <- df_chr
                      cTable <- df_chr
                      aTable$newTotalCN <- 0
                      cTable$newTotalCN <- 0
                      
                      # new loop below gets segments of gains or losses
                      delIdx <- which(sign(df_chr$newTotalCN) == -1)
                      ampIdx <- which(sign(df_chr$newTotalCN) == 1)
                      
                      # thought experiment - this actually undercounts because of arm
                      # level events - overcount cnas, if I don't skip here
                      if(length(delIdx) == 0 & length(ampIdx) == 0){
                        res <- rbind(res, aTable)
                        res2 <- rbind(res2, cTable)
                        next()
                      }
                      
                      if (sum(df_chr$length[delIdx]) > minLength80) {
                        resMat[1,2] <- resMat[1,2] + 1
                        aTable$newTotalCN[delIdx] <- df_chr$newTotalCN[delIdx]
                        skipVar <- "yes"
                        # next()
                        # print(paste(i, j))
                      } 
                      
                      if (sum(df_chr$length[ampIdx]) > minLength80) {
                        resMat[1,1] <- resMat[1,1] + 1
                        aTable$newTotalCN[ampIdx] <- df_chr$newTotalCN[ampIdx]
                        skipVar <- "yes"
                        # next()
                        # print(paste(i, j))
                      }
                      
                      if (skipVar == "yes") {
                        res <- rbind(res, aTable)
                        res2 <- rbind(res2, cTable)
                        next() 
                      }
                      
                      
                      ### cna counting
                      
                      cn_sign <- sign(df_chr$newTotalCN)
                      k <- 1
                      idx1 <- 0
                      idx2 <- 0
                      tmpSign <- 0
                      while (k < (length(cn_sign) + 1)) {
                        if (cn_sign[k] == 0 & idx1 == 0) { #start + non end string of zeros
                          k <- k + 1
                          # print(1)
                          next()
                        } else if(cn_sign[k] != 0 & idx1 == 0 & k == length(cn_sign)){ # if there is a single segment at the end
                          idx1 <- length(cn_sign)
                          idx2 <- idx1
                          tmpLength <- sum(df_chr$length[idx1:idx2])
                          print(tmpSign <- cn_sign[k])
                          if (tmpLength > minLength & tmpLength < minLength80) {
                            if (tmpSign == 1) {
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,3] <- resMat[1,3] + 1
                            } else if(tmpSign == -1){
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,4] <- resMat[1,4] + 1
                            }
                            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
                          }
                          k <- k + 1
                          # print(2)
                          next()
                        } else if(cn_sign[k] != 0 & idx1 == 0){
                          idx1 <- k
                          tmpSign <- cn_sign[k]
                          k <- k + 1
                          # print(3)
                          next()
                        } else if(idx1 != 0 & tmpSign == cn_sign[k] & k == length(cn_sign)){
                          idx2 <- k
                          tmpLength <- sum(df_chr$length[idx1:idx2])
                          if (tmpLength > minLength & tmpLength <  minLength80) {
                            if (tmpSign == 1) {
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,3] <- resMat[1,3] + 1
                            } else if(tmpSign == -1){
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,4] <- resMat[1,4] + 1
                            }
                          }
                          k <- k + 1
                          # print("long-end")
                          print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
                          next()
                        } else if (idx1 != 0 & tmpSign == cn_sign[k] & k != length(cn_sign)){
                          idx2 <- k
                          k <- k + 1
                          # print(4)
                          next()
                        } else if(cn_sign[k] == (tmpSign * -1) & idx2 == 0){ # special case where -1 and 1 are adjacent
                          idx2 <- idx1
                          tmpLength <- sum(df_chr$length[idx1:idx2])
                          if (tmpLength > minLength & tmpLength < minLength80) {
                            if (tmpSign == 1) {
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,3] <- resMat[1,3] + 1
                            } else if(tmpSign == -1){
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,4] <- resMat[1,4] + 1
                            }
                            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
                          }
                          
                          if (k == length(cn_sign)) {
                            idx1 <- k
                            idx2 <- k
                            tmpSign <- cn_sign[k]
                            if (tmpSign == 1) {
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,3] <- resMat[1,3] + 1
                            } else if(tmpSign == -1){
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,4] <- resMat[1,4] + 1
                            }
                            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
                          }
                          
                          idx1 <- k
                          idx2 <- 0
                          tmpSign <- cn_sign[k]
                          k <- k + 1
                          # print("adjacent")
                          next()
                        } else if(tmpSign != cn_sign[k] & idx2 == 0){
                          idx2 <- idx1
                          tmpLength <- sum(df_chr$length[idx1:idx2])
                          if (tmpLength > minLength & tmpLength < minLength80) {
                            if (tmpSign == 1) {
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,3] <- resMat[1,3] + 1
                            } else if(tmpSign == -1){
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,4] <- resMat[1,4] + 1
                            }
                            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
                          }
                          idx1 <- 0
                          idx2 <- 0
                          k <- k + 1
                          # print(5)
                          next()
                        } else if(cn_sign[k] == (tmpSign * -1) & idx2 != 0){
                          tmpLength <- sum(df_chr$length[idx1:idx2])
                          if (tmpLength > minLength & tmpLength < minLength80) {
                            if (tmpSign == 1) {
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,3] <- resMat[1,3] + 1
                            } else if(tmpSign == -1){
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,4] <- resMat[1,4] + 1
                            }
                            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
                          }
                          idx1 <- k
                          tmpSign <- cn_sign[k]
                          k <- k + 1
                          idx2 <- 0
                          # print("adjacent 2")
                          next()
                        } else if(tmpSign != cn_sign[k] & idx2 != 0){
                          tmpLength <- sum(df_chr$length[idx1:idx2])
                          if (tmpLength > minLength & tmpLength < minLength80) {
                            if (tmpSign == 1) {
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,3] <- resMat[1,3] + 1
                            } else if(tmpSign == -1){
                              cTable$newTotalCN[idx1:idx2] <- df_chr$newTotalCN[idx1:idx2]
                              resMat[1,4] <- resMat[1,4] + 1
                            }
                            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
                          }
                          idx1 <- 0
                          idx2 <- 0
                          k <- k + 1
                          #print(6)
                          next()
                        }
                      }
                      res <- rbind(res, aTable)
                      res2 <- rbind(res2, cTable)
                    }
                    rownames(resMat) <- i
                    list(res, res2, resMat)
                  }
  
  stopCluster(cl)
  return(oper)
}


# 
# separateSegments_m <- function(df, chromList){
#   # whole chromosome/arms (covers ~80% of max chr coverage from panel); changes > 1MB
#   # easier with mouse since no ploidy table
#   res <- NULL
#   res2 <- NULL
#   df$length <- df$end.pos - df$start.pos
#   
#   # 20210715:the arm level calls don't match up too well with how the segmentation looks
#   # I think a lot of arm level stuff is being called as cna - new strat just do based on sign (+/-)
#   
#   for (i in unique(df$sampleID)) {
#     sampDf <- df[which(df$sampleID == i),]
#     for (j in unique(chromList$chr)) {
#       df_chr <- sampDf[sampDf$chrom == j,]
#       minLength <- chromList$length80[which(chromList$chr == j)]
#       aTable <- df_chr
#       cTable <- df_chr
# 
#       # new loop below gets segments of gains or losses
#       signChange <- c(0, diff(sign(df_chr$mean)))
#       idxChange <- which(signChange != 0)
# 
#       if(length(idxChange) == 0){
#         res <- rbind(res, aTable)
#         res2 <- rbind(res2, cTable)
#         next()
#       }
# 
#       for (k in 1:(length(idxChange))){
#         if (df_chr$mean[idxChange[k]] == 0) {
#           next()
#         } else if (k == length(idxChange)){
#           idx1 <- idxChange[k]
#           idx2 <- length(df_chr$mean)
#         } else{
#           idx1 <- idxChange[k]
#           idx2 <- idxChange[k + 1] - 1
#         }
# 
#         tmpLength <- sum(df_chr$length[idx1:idx2])
#         if(tmpLength > minLength){
#           cTable$mean[idx1:idx2] <- 0
#         } else{
#           print(paste("coord: ",aTable$sampleID[idx1], aTable$chrom[idx1],
#                        aTable$start.pos[idx1], "-",aTable$end.pos[idx1]))
#           print(paste0("before", aTable$mean[idx1:idx2]))
#           aTable$mean[idx1:idx2] <- 0
#           print(paste0("after", aTable$mean[idx1:idx2]))
#         }
#       }
#       res <- rbind(res, aTable)
#       res2 <- rbind(res2, cTable)
#     }
#   }
#   
#   # 20210714 depreceated b/c it didn't detect enough arm-level gains
#   # for (i in unique(chromList$chr)) {
#   #   df_chr <- df[which(df$chrom == i),]
#   #   minLength <- chromList$length80[which(chromList$chr ==i)]
#   #   armTable <- df_chr
#   #   cnaTable <- df_chr
#   # 
#   #   armTable$mean[which(df_chr$length <= minLength)] <- 0
#   #   cnaTable$mean[which(df_chr$length > minLength)] <- 0
#   # 
#   #   res <- rbind(res, armTable)
#   #   res2 <- rbind(res2, cnaTable)
#   # }
#   return(list(res, res2))
# }

separateSegments_mV2 <- function(df){
  res <- NULL
  res2 <- NULL
  df$length <- df$end.pos - df$start.pos
  
  for (i in unique(df$sampleID)) {
    sampDf <- df[which(df$sampleID == i),]
    for (j in unique(sampDf$chrom)) {
      df_chr <- sampDf[sampDf$chrom == j,]
      minLength <- sum(df_chr$length) * 0.8
      aTable <- df_chr
      cTable <- df_chr
      
      cn_sign <- sign(df_chr$mean)
      
      idx1 <- NULL
      idx2 <- NULL
      
      i <- 1
      idx1 <- 0
      idx2 <- 0
      tmpSign <- 0
      while (i < (length(cn_sign) + 1)) {
        if (cn_sign[i] == 0 & idx1 == 0) { #start + non end string of zeros
          i <- i + 1
          print(1)
          next()
        } else if(cn_sign[i] != 0 & idx1 == 0 & i == length(cn_sign)){ # if there is a single segment at the end
          idx1 <- length(cn_sign)
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          i <- i + 1
          print(2)
          next()
        } else if(cn_sign[i] != 0 & idx1 == 0){
          idx1 <- i
          tmpSign <- cn_sign[i]
          i <- i + 1
          print(3)
          next()
        } else if(idx1 != 0 & tmpSign == cn_sign[i] & i == length(cn_sign)){
          idx2 <- i
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          i <- i + 1
          print("long-end")
          next()
        } else if (idx1 != 0 & tmpSign == cn_sign[i] & i != length(cn_sign)){
          idx2 <- i
          i <- i + 1
          print(4)
          next()
        } else if(cn_sign[i] == (tmpSign * -1) & idx2 == 0){ # special case where -1 and 1 are adjacent
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          
          if (i ==  length(cn_sign)) {
            idx1 <- i
            idx2 <- i
            if (tmpLength > minLength) {
              cTable$mean[idx1:idx2] <- 0
            } else{
              aTable$mean[idx1:idx2] <- 0
            }
          }
          
          idx1 <- i
          tmpSign <- cn_sign[i]
          i <- i + 1
          print("adjacent")
          next()
        } else if(tmpSign != cn_sign[i] & idx2 == 0){
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          idx1 <- 0
          idx2 <- 0
          i <- i + 1
          print(5)
          next()
        } else if(cn_sign[i] == (tmpSign * -1) & idx2 != 0){
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          idx1 <- i + 1
          tmpSign <- cn_sign[i]
          i <- i + 1
          idx2 <- 0
          print("adjacent 2")
          next()
        } else if(tmpSign != cn_sign[i] & idx2 != 0){
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength) {
            cTable$mean[idx1:idx2] <- 0
          } else{
            aTable$mean[idx1:idx2] <- 0
          }
          idx1 <- 0
          i <- i + 1
          print(6)
          next()
        }
      }
      res <- rbind(res, aTable)
      res2 <- rbind(res2, cTable)
    }
  }
  return(list(res, res2))
}



# this version is a mix of total coverage >= 80% and sign
separateSegments_mV3 <- function(df){
  res <- NULL
  res2 <- NULL
  res3 <- NULL
  df$length <- df$end.pos - df$start.pos
  
  for (i in unique(df$sampleID)) {
    resMat <- matrix(data = 0, nrow = 1, ncol = 4)
    sampDf <- df[which(df$sampleID == i),]
    
    ### aneuploidy counting
    for (j in unique(sampDf$chrom)) {
      skipVar <- "no"
      df_chr <- sampDf[sampDf$chrom == j,]
      aTable <- df_chr
      cTable <- df_chr
      aTable$mean <- 0
      cTable$mean <- 0
      minLength80 <- sum(df_chr$length) * 0.8
      minLength <- 15000000
      
      # new loop below gets segments of gains or losses
      delIdx <- which(sign(df_chr$mean) == -1)
      ampIdx <- which(sign(df_chr$mean) == 1)
      
      # thought experiment - this actually undercounts because of arm
      # level events - overcount cnas, if I don't skip here
      if(length(delIdx) == 0 & length(ampIdx) == 0){
        res <- rbind(res, aTable)
        res2 <- rbind(res2, cTable)
        next()
      }
      
      if (sum(df_chr$length[delIdx]) > minLength80) {
        resMat[1,2] <- resMat[1,2] + 1
        aTable$mean[delIdx] <- df_chr$mean[delIdx]
        skipVar <- "yes"
        # next()
        # print(paste(i, j))
      } 
      
      if (sum(df_chr$length[ampIdx]) > minLength80) {
        resMat[1,1] <- resMat[1,1] + 1
        aTable$mean[ampIdx] <- df_chr$mean[ampIdx]
        skipVar <- "yes"
        # next()
        # print(paste(i, j))
      }
      
      if (skipVar == "yes") {
        res <- rbind(res, aTable)
        res2 <- rbind(res2, cTable)
        next() 
      }
      
      
      ### cna counting
      
      cn_sign <- sign(df_chr$mean)
      k <- 1
      idx1 <- 0
      idx2 <- 0
      tmpSign <- 0
      while (k < (length(cn_sign) + 1)) {
        if (cn_sign[k] == 0 & idx1 == 0) { #start + non end string of zeros
          k <- k + 1
          # print(1)
          next()
        } else if(cn_sign[k] != 0 & idx1 == 0 & k == length(cn_sign)){ # if there is a single segment at the end
          idx1 <- length(cn_sign)
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          print(tmpSign <- cn_sign[k])
          if (tmpLength > minLength & tmpLength < minLength80) {
            if (tmpSign == 1) {
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
          }
          k <- k + 1
          # print(2)
          next()
        } else if(cn_sign[k] != 0 & idx1 == 0){
          idx1 <- k
          tmpSign <- cn_sign[k]
          k <- k + 1
          # print(3)
          next()
        } else if(idx1 != 0 & tmpSign == cn_sign[k] & k == length(cn_sign)){
          idx2 <- k
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength & tmpLength <  minLength80) {
            if (tmpSign == 1) {
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,4] <- resMat[1,4] + 1
            }
          }
          k <- k + 1
          # print("long-end")
          print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
          next()
        } else if (idx1 != 0 & tmpSign == cn_sign[k] & k != length(cn_sign)){
          idx2 <- k
          k <- k + 1
          # print(4)
          next()
        } else if(cn_sign[k] == (tmpSign * -1) & idx2 == 0){ # special case where -1 and 1 are adjacent
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength & tmpLength < minLength80) {
            if (tmpSign == 1) {
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
          }
          
          if (k == length(cn_sign)) {
            idx1 <- k
            idx2 <- k
            tmpSign <- cn_sign[k]
            if (tmpSign == 1) {
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
          }
          
          idx1 <- k
          idx2 <- 0
          tmpSign <- cn_sign[k]
          k <- k + 1
          # print("adjacent")
          next()
        } else if(tmpSign != cn_sign[k] & idx2 == 0){
          idx2 <- idx1
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength & tmpLength < minLength80) {
            if (tmpSign == 1) {
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
          }
          idx1 <- 0
          idx2 <- 0
          k <- k + 1
          # print(5)
          next()
        } else if(cn_sign[k] == (tmpSign * -1) & idx2 != 0){
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength & tmpLength < minLength80) {
            if (tmpSign == 1) {
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
          }
          idx1 <- k
          tmpSign <- cn_sign[k]
          k <- k + 1
          idx2 <- 0
          # print("adjacent 2")
          next()
        } else if(tmpSign != cn_sign[k] & idx2 != 0){
          tmpLength <- sum(df_chr$length[idx1:idx2])
          if (tmpLength > minLength & tmpLength < minLength80) {
            if (tmpSign == 1) {
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,3] <- resMat[1,3] + 1
            } else if(tmpSign == -1){
              cTable$mean[idx1:idx2] <- df_chr$mean[idx1:idx2]
              resMat[1,4] <- resMat[1,4] + 1
            }
            print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
          }
          idx1 <- 0
          idx2 <- 0
          k <- k + 1
          #print(6)
          next()
        }
      }
      res <- rbind(res, aTable)
      res2 <- rbind(res2, cTable)
    }
    
    rownames(resMat) <- i
    res3 <- rbind(res3, resMat)
  }
  return(list(res, res2, res3))
}


firstUpper <- function(gene){
  firstLetter <- toupper(substr(gene, start = 1, stop = 1))
  restOfGene <- tolower(substr(gene, start = 2, stop = nchar(gene)))
  res <- paste0(firstLetter, restOfGene)
  return(res)
}

maximums <- function(x) which(x - shift(x, 10) > 0  & x - shift(x, 10, type='lead') > 0)
minimums <- function(x) which(x - shift(x, 10) < 0  & x - shift(x, 10, type='lead') < 0)


### creates final bed track for circos plot
trackBedColumns <- function(df){
  colnames(df) <- tolower(colnames(df))
  df <- cbind(df, "ybot" = NA, "ytop" = NA)
  df$freq <- df$freq/100
  df$ytop <- ifelse(df$freq > 0, df$freq, 0) 
  df$ybot <- ifelse(df$freq < 0, df$freq, 0)
  return(df)
}



### creates both the freq by freq plot and outputs precursor for pathway analysis


circosFreq <- function(m_amp, m_del, tcga_amp, tcga_del, filename = "test"){
  # all frequencies are positive numbers, made negative
  # for graphing and empirical purposes
  tcga_del$Freq  <- tcga_del$Freq * -1
  m_amp <- m_amp[which(m_amp$Freq > 0),]
  m_del <- m_del[which(m_del$Freq > 0),]
  
  mouse_bed <- NULL
  human_bed <- NULL
  
  if (nrow(m_amp) > 0) {
    ampSynteny <- syntenyPlotInputsFreqV2(m_amp)
    mouse_bed <- rbind(ampSynteny$mouse_bed)
    human_bed <- rbind(ampSynteny$human_bed)
  }
  
  if(nrow(m_del) >  0){
    delSynteny <- syntenyPlotInputsFreqV2(m_del)
    delSynteny$mouse_bed$freq <- delSynteny$mouse_bed$freq * -1
    mouse_bed <- rbind(delSynteny$mouse_bed)
    human_bed <- rbind(delSynteny$human_bed)
  }
  
  # ampSynteny <- syntenyPlotInputsFreqV2(m_amp)
  # delSynteny <- syntenyPlotInputsFreqV2(m_del)
  # delSynteny$mouse_bed$freq <- delSynteny$mouse_bed$freq * -1
  # mouse_bed <- rbind(ampSynteny$mouse_bed, delSynteny$mouse_bed)
  # human_bed <- rbind(ampSynteny$human_bed, delSynteny$human_bed)
  
  
  allSynTable <- cbind(human_bed, mouse_bed)
  colnames(allSynTable) <- c("h_chr", "h_start", "h_end", "m_chr",
                             "m_start", "m_end", "m_freq")
  
  
  cyto_interest <- unique(c(human_bed$chr, mouse_bed$chr))
  cyto_combined_red <- cyto_hg38mm10[which(cyto_hg38mm10$V1 %in% cyto_interest),]
  
  mouse_cn <- trackBedColumns(mouse_bed)
  human_cn <- trackBedColumns(rbind(tcga_amp,
                                    tcga_del))
  human_cn$chr <- paste0("h_chr", human_cn$chr)
  
  cn_track_bed <- rbind(mouse_cn, human_cn)
  cn_track_bed <- cn_track_bed[which(cn_track_bed$chr %in% cyto_interest), ]
  
  track_color <- rep("#000000", nrow(cn_track_bed))
  track_color <- ifelse(cn_track_bed$freq < 0 , "#00008B", "#8B0000") 
  cn_track_bed$col = track_color
  
  mouse_chrs <- unique(mouse_bed$chr)
  human_chrs <- unique(human_bed$chr)
  
  colorVector2 <- rep("#000000", nrow(mouse_bed))
  for (j in seq_along(mouse_chrs)) {
    colorVector2[which(mouse_bed$chr == mouse_chrs[j])] <- colorVector[j]
  }
  
  pdf(paste0("/mnt/DATA5/tmp/kev/misc/", filename, ".pdf"),
      width = 6, height = 6, useDingbats = FALSE)
  circos.par("track.height"= 0.10) +
    circos.initializeWithIdeogram(cytoband = cyto_combined_red, sort.chr = FALSE,
                                  plotType = c("ideogram", "labels")) +
    circos.genomicTrack(cn_track_bed, ylim = c(-1.00,1.00),
                        panel.fun = function(region, value, ...) {
                          circos.genomicRect(region, value, col = value$col,
                                             ytop = value$ytop,
                                             ybottom = value$ybot,
                                             border = NA,...)
                        }) +
    circos.genomicLink(mouse_bed, human_bed,
                       col = colorVector2,
                       border = NA)
  dev.off()
  
  
  ### adding h_freq for synteny tables
  freqVector <- NULL
  tmpGrange <- GRanges(seqnames = str_remove(human_cn$chr, "h_chr"),
                       IRanges(start = human_cn$start,
                               end = human_cn$end))
  for (i in 1:nrow(allSynTable)) {
    tmpGrangeH <- GRanges(seqnames = str_remove(allSynTable$h_chr[i], "h_chr"),
                          IRanges(start = allSynTable$h_start[i],
                                  end = allSynTable$h_end[i]))
    meanFreq <- mean(human_cn$freq[subjectHits(findOverlaps(tmpGrangeH, tmpGrange))] * 100)
    if (is.nan(meanFreq)) {
      meanFreq <- 0
    }
    freqVector <- c(freqVector, meanFreq)
  }
  
  allSynTable$h_freq <- freqVector
  return(allSynTable)
}

getGeneList <- function(allSynTable){
  # extracts gene list from synteny table
  colnames(allSynTable) <- c("h_chr", "h_start", "h_end","m_chr",
                             "m_start", "m_end", "m_freq","h_freq")
  geneList <- NULL
  for (i in unique(allSynTable$m_chr)) {
    
    tmpSynDf <- allSynTable[which(allSynTable$m_chr == i),]
    
    for (j in unique(tmpSynDf$m_freq)) {
      tmpSegDf <- tmpSynDf[which(tmpSynDf$m_freq == j),]
      synGrangeM <- GRanges(seqnames =  str_remove(tmpSegDf$m_chr[1],"m_chr"),
                            IRanges(start = min(tmpSegDf$m_start),
                                    end = max(tmpSegDf$m_end)))
      synGrangeH <- GRanges(seqnames = str_remove(tmpSegDf$h_chr, "h_chr"),
                            IRanges(start = tmpSegDf$h_start,
                                    end = tmpSegDf$h_end))
      
      hGenes <- h_exon_boundaries$gene[subjectHits(findOverlaps(synGrangeH, hGrange))]
      mGenes <- m_exon_boundaries$gene[subjectHits(findOverlaps(synGrangeM, mGrange))]
      hFreq <- tmpSegDf$h_freq[queryHits(findOverlaps(synGrangeH, hGrange))]
      mFreq <- tmpSegDf$m_freq[queryHits(findOverlaps(synGrangeM, mGrange))]
      
      tmpGeneLists <- c("m_chr" = tmpSegDf$m_chr[1], "m_start" = min(tmpSegDf$m_start),
                        "m_end" = max(tmpSegDf$m_end),
                        "h_gene" = paste(hGenes, collapse = ","),
                        "m_gene" = paste(mGenes, collapse = ","),
                        "h_freq" = paste(hFreq, collapse = ","),
                        "m_freq" = paste(mFreq, collapse = ","))
      geneList <- rbind(geneList, tmpGeneLists)
    }
  }
  
  
  rownames(geneList) <- NULL
  geneList <- data.frame(geneList, stringsAsFactors = FALSE)
  
  # convert mouse gene names to human for enrichment analysis
  
  geneList$convert <- "empty"
  for (i in 1:nrow(geneList)) {
    tmpMGene <- unlist(str_split(geneList$m_gene[i], ","))
    tmpConvert <- NULL
    for (j in tmpMGene) {
      tmpVar <- geneNameDf$external_gene_name[which(geneNameDf$mmusculus_homolog_associated_gene_name == j)]
      if (length(tmpVar) == 0) {
        tmpVar <- toupper(j)
      } else if(length(tmpVar) > 1){
        tmpVar <- tmpVar[1]
      }
      tmpConvert <- c(tmpConvert, tmpVar)
    }
    geneList$convert[i] <- paste(tmpConvert, collapse = ",")
  }
  return(geneList)
}


enrichmentStats <- function(geneList, pathwayList){
  # pathway list should be processed 
  pathwayRes <- NULL
  
  all_h_gene_symbol <- unlist(str_split(paste(geneList$h_gene, collapse = ","), ","))
  all_m_gene_symbol <- unlist(str_split(paste(geneList$convert, collapse = ","), ","))
  all_h_value <- as.numeric(unlist(str_split(paste(geneList$h_freq, collapse = ","), ",")))
  all_m_value <- as.numeric(unlist(str_split(paste(geneList$m_freq, collapse = ","),",")))
  
  
  
  for(i in 1:length(pathwayList)){
    
    tmpHallList <- unlist(pathwayList[[i]])
    tmpMIdx <- which(all_m_gene_symbol %in% tmpHallList)
    tmpHIdx <- which(all_h_gene_symbol %in% tmpHallList)
    
    print(paste(length(tmpMIdx), length(tmpHIdx)))
    
    tmpMGene <- all_m_value[tmpMIdx]
    names(tmpMGene) <- all_m_gene_symbol[tmpMIdx]
    tmpMGene <- tmpMGene[-which(duplicated(names(tmpMGene)))]
    
    tmpHGene <- all_h_value[tmpHIdx]
    names(tmpHGene) <- all_h_gene_symbol[tmpHIdx]
    tmpHGene <- tmpHGene[-which(duplicated(names(tmpHGene)))]
    
    # using union to help create ordered empty vectors to compare
    tmpPathway <- data.frame("gene" = union(all_m_gene_symbol[tmpMIdx], all_h_gene_symbol[tmpHIdx]),
                             "h" = 0,
                             "m" = 0)
    tmpPathway$h[match(names(tmpHGene), tmpPathway$gene)] <- tmpHGene
    tmpPathway$m[match(names(tmpMGene), tmpPathway$gene)] <- tmpMGene
    tmpPathway$h <- ifelse(abs(tmpPathway$h) > 10,  tmpPathway$h, 0)
    tmpPathway$m <- ifelse(abs(tmpPathway$m) > 10,  tmpPathway$m, 0)
    
    tmpPathway2 <- tmpPathway
    tmpPathway2$h <- ifelse(tmpPathway2$h > 0,  1, tmpPathway2$h)
    tmpPathway2$h <- ifelse(tmpPathway2$h < 0,  -1, tmpPathway2$h)
    tmpPathway2$m <- ifelse(tmpPathway2$m  > 0, 1, tmpPathway2$m)
    tmpPathway2$m <- ifelse(tmpPathway2$m  < 0, -1, tmpPathway2$m)
    
    exactMat <- matrix(ncol = 2, nrow = 2)
    exactMat[1,1] <- table(tmpPathway2$h)[3]
    exactMat[2,1] <- table(tmpPathway2$h)[1]
    exactMat[1,2] <- table(tmpPathway2$m)[3]
    exactMat[2,2] <- table(tmpPathway2$m)[1]
    
    exactMat[which(is.na(exactMat))] <- 0
    
    fisherRes <- fisher.test(exactMat)
    
    print(tmpPathway$h)
    print(tmpPathway$m)
    
    corTest <- tryCatch(cor.test(tmpPathway$h, tmpPathway$m), error = function(x) return(NULL))
    #jac <- jaccard(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
    #tani <- tanimoto(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
    cosi <- cosine_dist(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
    ruzi <- ruzicka(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
    euc <- euclidean(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
    pathw <- names(human_hallmarks_list[[i]])
    pathwayRes <- rbind(pathwayRes, c("pathway" = pathw, corTest$estimate,
                                      "pval" = corTest$p.value, "cosine" = cosi,
                                      "ruzicka" = ruzi, "euclidean" = euc,
                                      "fisher" = fisherRes$p.value))
  }
  
  pathwayRes <- data.frame(pathwayRes, stringsAsFactors = FALSE)
  pathwayRes$cor <- as.numeric(pathwayRes$cor)
  pathwayRes$pval <- as.numeric(pathwayRes$pval)
  # pathwayRes_filt <- pathwayRes[which(abs(pathwayRes$cor) > .1),]
  # pathwayRes_filt <- pathwayRes_filt[order(pathwayRes_filt$cor),]
  # pathwayRes_filt$fill <- ifelse(pathwayRes_filt$cor > 0, "firebrick1", "lightblue")
  # pathwayRes_filt$pathway <- factor(pathwayRes_filt$pathway, levels = unique(pathwayRes_filt$pathway))
  # 
  

  
  return(pathwayRes)
}


round2 = function(x, n =  0) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^n
  z*posneg
}


# 
# # added this to separateArms function
# countAneuCna_tcga <- function(df, ploidyDf){
#   
#   # this was genius and got rid of the loop - literally infinitely faster *pats self on back*
#   df$newTotalCN <- round2(df$Copynumber - ploidyDf$ploidy[match(df$Sample, ploidyDf$array)], 0)
#   
#   cl <- makeCluster(20)
#   registerDoParallel(cl)
#   
#   res <- foreach(i=unique(df$Sample), .combine= 'rbind') %dopar% {
#     
#             resMat <- matrix(data = 0, nrow = 1, ncol = 4) # aneuploidy gain, loss: CNA 15mb gain, loss
#             sampDf <- df[which(df$Sample == i),]
#             human_arm <- read.table("/mnt/DATA5/tmp/kev/misc/20210713human_arm_syn.txt",
#                                     sep = "\t", stringsAsFactors = FALSE, header = TRUE)
#             
#             
#             ### counting for anueploidy
#             for (j in unique(human_arm$chr)) {
#               skipVar <- "no"
#               df_chr <- sampDf[sampDf$Chromosome == j,]
#               df_chr <- df_chr[order(df_chr$Start),]
#               minLength80 <- min(human_arm$length[which(human_arm$chromosome == j)]) * 0.8
#               minLength <- 15000000
#               
#               # new loop below gets segments of gains or losses
#               delIdx <- which(sign(df_chr$newTotalCN) == -1)
#               ampIdx <- which(sign(df_chr$newTotalCN) == 1)
#               
#               # thought experiment - this actually undercounts because of arm
#               # level events - overcount cnas, if I don't skip here
#               if(length(delIdx) == 0 & length(ampIdx) == 0){
#                 next()
#               }
#               
#               if (sum(df_chr$length[delIdx]) > minLength80) {
#                 resMat[1,2] <- resMat[1,2] + 1
#                 skipVar <- "yes"
#                 # next()
#                 # print(paste(i, j))
#               } 
#               
#               if (sum(df_chr$length[ampIdx]) > minLength80) {
#                 resMat[1,1] <- resMat[1,1] + 1
#                 skipVar <- "yes"
#                 # next()
#                 # print(paste(i, j))
#               }
#               
#               if (skipVar == "yes") {
#                next() 
#               }
#               
#               
#               ### cna counting
#               
#               cn_sign <- sign(df_chr$newTotalCN)
#               k <- 1
#               idx1 <- 0
#               idx2 <- 0
#               tmpSign <- 0
#               while (k < (length(cn_sign) + 1)) {
#                 if (cn_sign[k] == 0 & idx1 == 0) { #start + non end string of zeros
#                   k <- k + 1
#                   # print(1)
#                   next()
#                 } else if(cn_sign[k] != 0 & idx1 == 0 & k == length(cn_sign)){ # if there is a single segment at the end
#                   idx1 <- length(cn_sign)
#                   idx2 <- idx1
#                   tmpLength <- sum(df_chr$length[idx1:idx2])
#                   print(tmpSign <- cn_sign[k])
#                   if (tmpLength > minLength & tmpLength < minLength80) {
#                     if (tmpSign == 1) {
#                       resMat[1,3] <- resMat[1,3] + 1
#                     } else if(tmpSign == -1){
#                       resMat[1,4] <- resMat[1,4] + 1
#                     }
#                     print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#                   }
#                   k <- k + 1
#                   # print(2)
#                   next()
#                 } else if(cn_sign[k] != 0 & idx1 == 0){
#                   idx1 <- k
#                   tmpSign <- cn_sign[k]
#                   k <- k + 1
#                   # print(3)
#                   next()
#                 } else if(idx1 != 0 & tmpSign == cn_sign[k] & k == length(cn_sign)){
#                   idx2 <- k
#                   tmpLength <- sum(df_chr$length[idx1:idx2])
#                   if (tmpLength > minLength & tmpLength <  minLength80) {
#                     if (tmpSign == 1) {
#                       resMat[1,3] <- resMat[1,3] + 1
#                     } else if(tmpSign == -1){
#                       resMat[1,4] <- resMat[1,4] + 1
#                     }
#                   }
#                   k <- k + 1
#                   # print("long-end")
#                   print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#                   next()
#                 } else if (idx1 != 0 & tmpSign == cn_sign[k] & k != length(cn_sign)){
#                   idx2 <- k
#                   k <- k + 1
#                   # print(4)
#                   next()
#                 } else if(cn_sign[k] == (tmpSign * -1) & idx2 == 0){ # special case where -1 and 1 are adjacent
#                   idx2 <- idx1
#                   tmpLength <- sum(df_chr$length[idx1:idx2])
#                   if (tmpLength > minLength & tmpLength < minLength80) {
#                     if (tmpSign == 1) {
#                       resMat[1,3] <- resMat[1,3] + 1
#                     } else if(tmpSign == -1){
#                       resMat[1,4] <- resMat[1,4] + 1
#                     }
#                     print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#                   }
#                   
#                   if (k == length(cn_sign)) {
#                     idx1 <- k
#                     idx2 <- k
#                     tmpSign <- cn_sign[k]
#                     if (tmpSign == 1) {
#                       resMat[1,3] <- resMat[1,3] + 1
#                     } else if(tmpSign == -1){
#                       resMat[1,4] <- resMat[1,4] + 1
#                     }
#                     print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#                   }
#                   
#                   idx1 <- k
#                   idx2 <- 0
#                   tmpSign <- cn_sign[k]
#                   k <- k + 1
#                   # print("adjacent")
#                   next()
#                 } else if(tmpSign != cn_sign[k] & idx2 == 0){
#                   idx2 <- idx1
#                   tmpLength <- sum(df_chr$length[idx1:idx2])
#                   if (tmpLength > minLength & tmpLength < minLength80) {
#                     if (tmpSign == 1) {
#                       resMat[1,3] <- resMat[1,3] + 1
#                     } else if(tmpSign == -1){
#                       resMat[1,4] <- resMat[1,4] + 1
#                     }
#                     print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#                   }
#                   idx1 <- 0
#                   idx2 <- 0
#                   k <- k + 1
#                   # print(5)
#                   next()
#                 } else if(cn_sign[k] == (tmpSign * -1) & idx2 != 0){
#                   tmpLength <- sum(df_chr$length[idx1:idx2])
#                   if (tmpLength > minLength & tmpLength < minLength80) {
#                     if (tmpSign == 1) {
#                       resMat[1,3] <- resMat[1,3] + 1
#                     } else if(tmpSign == -1){
#                       resMat[1,4] <- resMat[1,4] + 1
#                     }
#                     print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#                   }
#                   idx1 <- k
#                   tmpSign <- cn_sign[k]
#                   k <- k + 1
#                   idx2 <- 0
#                   # print("adjacent 2")
#                   next()
#                 } else if(tmpSign != cn_sign[k] & idx2 != 0){
#                   tmpLength <- sum(df_chr$length[idx1:idx2])
#                   if (tmpLength > minLength & tmpLength < minLength80) {
#                     if (tmpSign == 1) {
#                       resMat[1,3] <- resMat[1,3] + 1
#                     } else if(tmpSign == -1){
#                       resMat[1,4] <- resMat[1,4] + 1
#                     }
#                     print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#                   }
#                   idx1 <- 0
#                   idx2 <- 0
#                   k <- k + 1
#                   #print(6)
#                   next()
#                 }
#               }
#             }
#             
#             rownames(resMat) <- i
#             resMat
#           }
#   stopCluster(cl)
#   return(res)
# }


# added this to separateArms function
# countAneuCna_mouse <- function(df){
#   res <- NULL
#   df$length <- df$end.pos - df$start.pos
# 
#   for (i in unique(df$sampleID)) {
#     resMat <- matrix(data = 0, nrow = 1, ncol = 4)
#     sampDf <- df[which(df$sampleID == i),]
# 
#     ### aneuploidy counting
#     for (j in unique(sampDf$chrom)) {
#       df_chr <- sampDf[sampDf$chrom == j,]
#       minLength80 <- sum(df_chr$length) * 0.8
#       minLength <- 15000000
# 
#       # new loop below gets segments of gains or losses
#       delIdx <- which(sign(df_chr$mean) == -1)
#       ampIdx <- which(sign(df_chr$mean) == 1)
# 
#       if(length(delIdx) == 0 & length(ampIdx) == 0){
#         next()
#       }
# 
#       if (sum(df_chr$length[delIdx]) > minLength80) {
#         resMat[1,2] <- resMat[1,2] + 1
#         next()
#         #print(paste(i, j))
#       }
# 
#       if (sum(df_chr$length[ampIdx]) > minLength80) {
#         resMat[1,1] <- resMat[1,1] + 1
#         next()
#         #print(paste(i, j))
#       }
# 
# 
# 
#       ### cna counting
# 
#       cn_sign <- sign(df_chr$mean)
#       k <- 1
#       idx1 <- 0
#       idx2 <- 0
#       tmpSign <- 0
#       while (k < (length(cn_sign) + 1)) {
#         if (cn_sign[k] == 0 & idx1 == 0) { #start + non end string of zeros
#           k <- k + 1
#           # print(1)
#           next()
#         } else if(cn_sign[k] != 0 & idx1 == 0 & k == length(cn_sign)){ # if there is a single segment at the end
#           idx1 <- length(cn_sign)
#           idx2 <- idx1
#           tmpLength <- sum(df_chr$length[idx1:idx2])
#           print(tmpSign <- cn_sign[k])
#           if (tmpLength > minLength & tmpLength < minLength80) {
#             if (tmpSign == 1) {
#               resMat[1,3] <- resMat[1,3] + 1
#             } else if(tmpSign == -1){
#               resMat[1,4] <- resMat[1,4] + 1
#             }
#             print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#           }
#           k <- k + 1
#           print(2)
#           next()
#         } else if(cn_sign[k] != 0 & idx1 == 0){
#           idx1 <- k
#           tmpSign <- cn_sign[k]
#           k <- k + 1
#           print(3)
#           next()
#         } else if(idx1 != 0 & tmpSign == cn_sign[k] & k == length(cn_sign)){
#           idx2 <- k
#           tmpLength <- sum(df_chr$length[idx1:idx2])
#           if (tmpLength > minLength & tmpLength <  minLength80) {
#             if (tmpSign == 1) {
#               resMat[1,3] <- resMat[1,3] + 1
#             } else if(tmpSign == -1){
#               resMat[1,4] <- resMat[1,4] + 1
#             }
#           }
#           k <- k + 1
#           print("long-end")
#           print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#           next()
#         } else if (idx1 != 0 & tmpSign == cn_sign[k] & k != length(cn_sign)){
#           idx2 <- k
#           k <- k + 1
#           # print(4)
#           next()
#         } else if(cn_sign[k] == (tmpSign * -1) & idx2 == 0){ # special case where -1 and 1 are adjacent
#           idx2 <- idx1
#           tmpLength <- sum(df_chr$length[idx1:idx2])
#           if (tmpLength > minLength & tmpLength < minLength80) {
#             if (tmpSign == 1) {
#               resMat[1,3] <- resMat[1,3] + 1
#             } else if(tmpSign == -1){
#               resMat[1,4] <- resMat[1,4] + 1
#             }
#             print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#           }
# 
#           if (k == length(cn_sign)) {
#             idx1 <- k
#             idx2 <- k
#             tmpSign <- cn_sign[k]
#             if (tmpSign == 1) {
#               resMat[1,3] <- resMat[1,3] + 1
#             } else if(tmpSign == -1){
#               resMat[1,4] <- resMat[1,4] + 1
#             }
#             print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#           }
# 
#           idx1 <- k
#           idx2 <- 0
#           tmpSign <- cn_sign[k]
#           k <- k + 1
#           print("adjacent")
#           next()
#         } else if(tmpSign != cn_sign[k] & idx2 == 0){
#           idx2 <- idx1
#           tmpLength <- sum(df_chr$length[idx1:idx2])
#           if (tmpLength > minLength & tmpLength < minLength80) {
#             if (tmpSign == 1) {
#               resMat[1,3] <- resMat[1,3] + 1
#             } else if(tmpSign == -1){
#               resMat[1,4] <- resMat[1,4] + 1
#             }
#             print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#           }
#           idx1 <- 0
#           idx2 <- 0
#           k <- k + 1
#           print(5)
#           next()
#         } else if(cn_sign[k] == (tmpSign * -1) & idx2 != 0){
#           tmpLength <- sum(df_chr$length[idx1:idx2])
#           if (tmpLength > minLength & tmpLength < minLength80) {
#             if (tmpSign == 1) {
#               resMat[1,3] <- resMat[1,3] + 1
#             } else if(tmpSign == -1){
#               resMat[1,4] <- resMat[1,4] + 1
#             }
#             print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#           }
#           idx1 <- k
#           tmpSign <- cn_sign[k]
#           k <- k + 1
#           idx2 <- 0
#           print("adjacent 2")
#           next()
#         } else if(tmpSign != cn_sign[k] & idx2 != 0){
#           tmpLength <- sum(df_chr$length[idx1:idx2])
#           if (tmpLength > minLength & tmpLength < minLength80) {
#             if (tmpSign == 1) {
#               resMat[1,3] <- resMat[1,3] + 1
#             } else if(tmpSign == -1){
#               resMat[1,4] <- resMat[1,4] + 1
#             }
#             print(paste0(i,"/","chr", df_chr$chrom[idx1], ":",df_chr$start.pos[idx1], "-",df_chr$end.pos[idx2]))
#           }
#           idx1 <- 0
#           idx2 <- 0
#           k <- k + 1
#           print(6)
#           next()
#         }
#       }
#     }
# 
#     rownames(resMat) <- i
#     res <- rbind(res, resMat)
#   }
#   return(res)
# }


# var used by functions

colorVector <- c("#CD5C5C", "#C71585", "#FF8C00", "#F0E68C",
                 "#4B0082", "#32CD32", "#000080", "#800000",
                 "#F5F5DC", "#2F4F4F", "#808000", "#E6E6FA",
                 "#FFD700", "#FFA07A", "#7B68EE", "#008B8B",
                 "#BC8F8F", "#FAEBD7", "#FFE4E1", "#000000",
                 "#00FFFF", "#8FBC8B", "#B22222", "#D0D3D4")



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





fgaCalculator_tcga <- function(df, ploidyDf){ # should work for any absolute output
  res <- NULL
  
  df$newTotalCN <- round2(df$Copynumber - ploidyDf$ploidy[match(df$Sample, ploidyDf$array)], 0)
  for (i in unique(df$Sample)) {
    sampDf <- df[which(df$Sample == i),]
    nonZeroes <- which(sign(sampDf$newTotalCN) != 0)
    fga <- sum(sampDf$length[nonZeroes])/sum(sampDf$length)
    res <- rbind(res, c(i, fga))
  }
  return(res)
}

fgaCalculator_amp <- function(df){ # works for any amplicon panel
  res <- NULL
  df$length <- df$end.pos - df$start.pos
  for (i in unique(df$sampleID)) {
    sampDf <- df[which(df$sampleID == i),]
    nonZeroes <- which(sign(sampDf$mean) != 0)
    fga <- sum(sampDf$length[nonZeroes])/sum(sampDf$length)
    res <- rbind(res, c(i, fga))
  }
  return(res)
}
