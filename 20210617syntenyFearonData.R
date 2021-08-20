# should be a list of functions needed to be used after
# processing segmentation (input used)
library(GenomicFeatures)
library(dplyr)
library(dbplyr)
library(circlize)
library(stringr)
library(DNAcopy)
library(stringr)
#library(optparse)
library(data.table)

# four variables should be mouse amplicon files, o/e ratio,
# normal vector name for FFPE vs FF, and output path




# data and functions needed for synteny
kap_id_tab <- read.table("/mnt/DATA5/tmp/kev/misc/20210628_panCancerCoadKAPIds.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

synteny_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_syntenyDf.txt",
                               sep = "\t", stringsAsFactors = FALSE, header = TRUE)
cyto_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_cytoDf.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = FALSE)
human_gene <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210407humanGeneTable.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)
mouse_gene <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210407mouseGeneTable.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)
geneNameDf <- read.table("/home/kevhu/data/20201021geneNameBiomart.txt", header = TRUE,
                         stringsAsFactors = FALSE, sep = "\t")
mouseBedFile <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.gc.bed", stringsAsFactors = FALSE, sep = "\t",
                           header = FALSE)


tcga_coad_amp <- read.table("/mnt/DATA5/tmp/kev/misc/20210617_coad_amp.bed", sep = "\t", stringsAsFactors = FALSE,
                            header = TRUE)
tcga_coad_del <- read.table("/mnt/DATA5/tmp/kev/misc/20210617_coad_del.bed", sep = "\t", stringsAsFactors = FALSE,
                            header = TRUE)


# below are the two inputs needed for segmentation
mouseAmplicons <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/cnAmplicon_matrix.txt", sep = "\t",
                             stringsAsFactors = FALSE, header = TRUE)

zscore_gc_oe_ratios <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/gcCorrectedCounts_matrix.txt", sep = "\t",
                                  stringsAsFactors = FALSE, header = TRUE)


all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                             ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                              end = zscore_gc_oe_ratios$EndPos))

colorVector <- c("#CD5C5C", "#C71585", "#FF8C00", "#F0E68C",
                 "#4B0082", "#32CD32", "#000080", "#800000",
                 "#F5F5DC", "#2F4F4F", "#808000", "#E6E6FA",
                 "#FFD700", "#FFA07A", "#7B68EE", "#008B8B",
                 "#BC8F8F", "#FAEBD7", "#FFE4E1", "#000000",
                 "#00FFFF", "#8FBC8B", "#B22222", "#D0D3D4")

# this needs to be an input vector in the snakefile - get inputs
# mouseNormal <- c("13604n", "14104t", "14154n", "14433n", "2405n", "2519n", "2796n", "3867n","8234n")
# 2796n 13604n - have BRCA1 gene deletions
# need the MG names


mouseNormal <- c("MG_17X49", "MG_18X50", "MG_23X55", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

# used in converted human gene names to mouse
firstUpper <- function(gene){
  firstLetter <- toupper(substr(gene, start = 1, stop = 1))
  restOfGene <- tolower(substr(gene, start = 2, stop = nchar(gene)))
  res <- paste0(firstLetter, restOfGene)
  return(res)
}

maximums <- function(x) which(x - shift(x, 10) > 0  & x - shift(x, 10, type='lead') > 0)
minimums <- function(x) which(x - shift(x, 10) < 0  & x - shift(x, 10, type='lead') < 0)


# function for processing the corrected amplicons into segmentation
# when I find time: parallelize this process

ampSeg <- function(mouseAmplicons, mouseBedFile){
  mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
  mouseAmplicons2 <- cbind(mouseBedFile$V1[mouseBedIdx], mouseBedFile$V2[mouseBedIdx], mouseAmplicons)
  mouseAmplicons2[,6:42] <- log2(mouseAmplicons2[,6:42])
  mouseAmplicons2$`mouseBedFile$V1[mouseBedIdx]` <- as.numeric(mouseAmplicons2$`mouseBedFile$V1[mouseBedIdx]`)
  allMouseCalls <- colnames(mouseAmplicons2)[6:42]
  segResults <- NULL
  for (i in allMouseCalls) {
    tmpCNA_obj <- CNA(cbind(mouseAmplicons2[[i]]),
                      mouseAmplicons2$ChromNum, mouseAmplicons2$`mouseBedFile$V2[mouseBedIdx]`,
                      data.type="logratio",sampleid=i)
    smoothed_tmpCNA <- smooth.CNA(tmpCNA_obj)
    segment_tmpCNA <- segment(smoothed_tmpCNA, verbose = 1, undo.splits = "sdundo", undo.SD = 0.2, alpha = 0.05,
                              p.method = c("hybrid"))
    tmpCalls <- segments.p(segment_tmpCNA)
    tmpCalls_filt <- tmpCalls[which(abs(tmpCalls$seg.mean) > .2),]
    segResults <- rbind(segResults, tmpCalls_filt)
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
  return(list("z_vector" = z_vector, "nmean_vector"= nmean_vector,
              "tmean_vector" = tmean_vector))
}


### filt could have parameter for size of filtering
segZscoresFilt <- function(segResults, zscores){
  
  # filts for significance, cn-change and length
  segmental_zscores <- cbind(segResults, "z-scores" = zscores[["z_vector"]],
                             "normal_seg_mean" = zscores[["nmean_vector"]],
                             "tumor_seg_mean" = zscores[["tmean_vector"]])
  segmental_zscores2 <- segmental_zscores[which(segmental_zscores$num.mark > 10),]
  segmental_zscores2$length <- segmental_zscores2$loc.end - segmental_zscores2$loc.start
  largeChromCutoff <- 1000000
  segmental_zscores3 <- segmental_zscores2[which(segmental_zscores2$length >= largeChromCutoff),]
  segmental_zscores3 <- segmental_zscores3[which(abs(segmental_zscores3$`z-scores`) > 1.645),]
  segInput <- segmental_zscores3[,c("ID", "chrom", "loc.start", "loc.end", "seg.mean")]
  return(segInput)
}

segZscoresFilt_zeroOut <- function(segResults, zscores){
  
  # filts for significance, cn-change and length
  segmental_zscores <- cbind(segResults, "z-scores" = zscores[["z_vector"]],
                             "normal_seg_mean" = zscores[["nmean_vector"]],
                             "tumor_seg_mean" = zscores[["tmean_vector"]])
  segmental_zscores2 <- segmental_zscores
  segmental_zscores2$seg.mean[which(segmental_zscores2$num.mark < 10)] <- 0
  segmental_zscores2$length <- segmental_zscores2$loc.end - segmental_zscores2$loc.start
  largeChromCutoff <- 1000000
  segmental_zscores2$seg.mean[which(segmental_zscores2$length <= largeChromCutoff)] <- 0
  segmental_zscores2$seg.mean[which(abs(segmental_zscores2$`z-scores`) < 1.645)] <- 0
  segInput <- segmental_zscores2[,c("ID", "chrom", "loc.start", "loc.end", "seg.mean")]
  return(segInput)
}

# input should be segmentation results - should go from mouse to human

syntenyPlotInputs <- function(segInput){
  segInput$chrom <- str_replace_all(segInput$chrom, "23", "X")
  tumor_freq_ranges <- GRanges(seqnames = paste0("chr", segInput$chrom),
                               ranges = IRanges(start = as.numeric(segInput$loc.start),
                                                end = as.numeric(segInput$loc.end)))
  synteny_ranges <- GRanges(seqnames = paste0("chr", synteny_hg38mm10$comp_chr),
                            ranges = IRanges(start = synteny_hg38mm10$comp_start_pos,
                                             end = synteny_hg38mm10$comp_end_pos))
  
  
  res_overlap <- findOverlaps(synteny_ranges, tumor_freq_ranges)
  query_hits <- queryHits(res_overlap)
  subject_hits <- subjectHits(res_overlap)
  
  
  mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg38mm10$comp_chr[query_hits]),
                          "start" = synteny_hg38mm10$comp_start_pos[query_hits],
                          "end" = synteny_hg38mm10$comp_end_pos[query_hits],
                          stringsAsFactors = FALSE)
  
  human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg38mm10$ref_chr[query_hits]),
                          "start" = synteny_hg38mm10$ref_start_pos[query_hits],
                          "end" = synteny_hg38mm10$ref_end_pos[query_hits],
                          stringsAsFactors = FALSE)
  
  
  
  resList <- list("human_bed" = human_bed,
                  "mouse_bed" = mouse_bed)
  
  return(resList)
}


syntenyPlotInputsFreq <- function(segInput){
  segInput$chrom <- str_replace_all(segInput$chrom, "23", "X")
  tumor_freq_ranges <- GRanges(seqnames = paste0("chr", segInput$chrom),
                               ranges = IRanges(start = as.numeric(segInput$loc.start),
                                                end = as.numeric(segInput$loc.end)))
  synteny_ranges <- GRanges(seqnames = paste0("chr", synteny_hg38mm10$comp_chr),
                            ranges = IRanges(start = synteny_hg38mm10$comp_start_pos,
                                             end = synteny_hg38mm10$comp_end_pos))
  
  
  res_overlap <- findOverlaps(synteny_ranges, tumor_freq_ranges)
  query_hits <- queryHits(res_overlap)
  subject_hits <- subjectHits(res_overlap)
  
  
  mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg38mm10$comp_chr[query_hits]),
                          "start" = synteny_hg38mm10$comp_start_pos[query_hits],
                          "end" = synteny_hg38mm10$comp_end_pos[query_hits],
                          "freq" = segInput$seg.mean[subject_hits],
                          stringsAsFactors = FALSE)
  
  human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg38mm10$ref_chr[query_hits]),
                          "start" = synteny_hg38mm10$ref_start_pos[query_hits],
                          "end" = synteny_hg38mm10$ref_end_pos[query_hits],
                          stringsAsFactors = FALSE)
  
  
  
  resList <- list("human_bed" = human_bed,
                  "mouse_bed" = mouse_bed)
  
  return(resList)
}

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



# gene part incomplete b/c odd matching and mouse genes weren't not found ...
# do this later and only for genes on panel ... easiest
# not that important
syntenyGeneInputs <- function(geneList){
  mouseGeneList <- geneList
  humanGeneList <- NULL
  for (i in geneList) {
    humanOrthoName <- geneNameDf$external_gene_name[which(geneNameDf$mmusculus_homolog_associated_gene_name == i)]
    if (length(humanOrthoName)  == 0) {
      humanOrthoName <- toupper(i)
    }
    
    humanGeneList <- c(humanGeneList, humanOrthoName)
  }
  
  tmpH <- human_gene[which(human_gene$external_gene_name %in% humanGeneList), ]
  tmpM <- mouse_gene[which(mouse_gene$external_gene_name %in% geneList), ]
  return(list(tmpH, tmpM))
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
  #amp2 <- amp2[which(amp2$Chr %in% cyto_interest),]
  colnames(amp2) <- c("chr","start","end","cn")
  
  amp3 <- downsampleRegions(amp2)
  amp3$ybot <- 0
  amp3$ytop <- amp3$cn
  
  del2$Freq <- round(del2$Freq/100, digits = 2)
  del2$Chr <- paste0("h_chr", del2$Chr)
  #del2 <- del2[which(del2$Chr %in% cyto_interest),]
  colnames(del2) <- c("chr","start","end","cn")
  del2$cn <- -del2$cn
  
  del3 <- downsampleRegions(del2)
  del3$ytop <- 0
  del3$ybot <- del3$cn
  
  cn_track_bed <- rbind(amp3, del3)
  
  return(cn_track_bed)
}


reducingFreqBed <- function(df, idx){
  reducedDf <- NULL
  for (i in 1:(length(idx) - 1)) {
    idx1 <- idx[i]
    idx2 <- idx[i+1]
    tmpDf <- df[idx1:idx2,]
    tmpChr <- tmpDf$chr[1]
    tmpStart <- min(tmpDf$pos)
    tmpEnd <- max(tmpDf$pos) - 1
    tmpFreq <- tmpDf[,3][1]
    
    tmpVec <- c("Chr" = tmpChr, "Start" = tmpStart, "End" = tmpEnd, "Freq" = tmpFreq)
    reducedDf <- rbind(reducedDf, tmpVec)
  }
  reducedDf <- data.frame(reducedDf, stringsAsFactors = FALSE)
  return(reducedDf)
}

# used on df with segmentation to get frequency bed
# ampsDels <- function(df){
#   getFreqOut <- getFreqData(df)
#   freqAmp <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] > 0.2)*100
#   freqDel <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] < -0.2)*100
#   
#   freqDf_Amp <- cbind(getFreqOut[,1:2], "amp" = freqAmp)
#   freqDf_Del <- cbind(getFreqOut[,1:2], "del" = freqDel)
#   
#   ampRedIdx <- c(1,1+which(diff(freqDf_Amp$amp)!=0))
#   delRedIdx <- c(1,1+which(diff(freqDf_Del$del)!=0))
#   
#   res <- list(freqDf_Amp, ampRedIdx, freqDf_Del, delRedIdx)
# }

circosFreq <- function(m_amp, m_del, tcga_amp, tcga_del, filename = "test"){
  # all frequencies are positive numbers, made negative
  # for graphing and empirical purposes
  tcga_del$Freq  <- tcga_del$Freq * -1
  m_amp <- m_amp[which(m_amp$Freq > 0),]
  m_del <- m_del[which(m_del$Freq > 0),]
  
  
  ampSynteny <- syntenyPlotInputsFreqV2(m_amp)
  delSynteny <- syntenyPlotInputsFreqV2(m_del)
  delSynteny$mouse_bed$freq <- delSynteny$mouse_bed$freq * -1
  mouse_bed <- rbind(ampSynteny$mouse_bed, delSynteny$mouse_bed)
  human_bed <- rbind(ampSynteny$human_bed, delSynteny$human_bed)
  
  
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
  
  print(head(allSynTable))
  return(allSynTable)
}




ampsDels <- function(df){
  getFreqOut <- getFreqData(df)
  # 20210629: for whatever reason cant find the mean of a long logical vector
  # i.e only TRUE AND FALSE
  #freqAmp <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] > 0.2)*100
  #freqDel <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] < -0.2)*100
  freqAmp <- apply(df[,-c(1:2)], 1, function(x) 
    length(which(x > 0.2))/length(x)) * 100
  freqDel <- apply(df[,-c(1:2)], 1, function(x) 
    length(which(x < -0.2))/length(x)) * 100
  
  freqDf_Amp <- cbind(getFreqOut[,1:2], "amp" = freqAmp)
  freqDf_Del <- cbind(getFreqOut[,1:2], "del" = freqDel)
  
  ampRedIdx <- c(1,1+which(diff(freqDf_Amp$amp)!=0))
  delRedIdx <- c(1,1+which(diff(freqDf_Del$del)!=0))
  
  res <- list(freqDf_Amp, ampRedIdx, freqDf_Del, delRedIdx)
}

# only includes large chromosomal changes 
# the ony change is instead of genes I can use the genomic track instead 
# to represent gains or losses

#segRes <- NULL
#segRes <- ampSeg(mouseAmplicons, mouseBedFile)
#save(segRes, file = "/mnt/DATA5/tmp/kev/misc/20210514mm10SegRes.Robj")
#load(file = "/mnt/DATA5/tmp/kev/misc/20210409mm10SegRes.Robj")
#load(file = "/mnt/DATA5/tmp/kev/misc/20210514mm10SegRes.Robj")
segResults <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/segResults.txt",
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t")


mouseBedChromLim <- NULL
for (i in unique(mouseBedFile$V1)) {
  if (i == "chrX") {
    next()
  }
  tmpBed <- mouseBedFile[which(mouseBedFile$V1 == i),]
  tmpStart <- min(tmpBed$V2)
  tmpEnd <- max(tmpBed$V3)
  mouseBedChromLim <- rbind(mouseBedChromLim, c("chr" = i,
                            "start" = tmpStart,
                            "end" = tmpEnd))
}

mouseBedChromLim  <- data.frame(mouseBedChromLim, stringsAsFactors = FALSE)
mouseBedChromLim$start <- as.numeric(mouseBedChromLim$start)
mouseBedChromLim$end <- as.numeric(mouseBedChromLim$end)
mouseBedChromLim$chr <- str_remove(mouseBedChromLim$chr,"chr")
mouseBedChromLim$length <- mouseBedChromLim$end - mouseBedChromLim$start
mouseBedChromLim$length80 <- 0.8 * mouseBedChromLim$length
### need to do tc correction before segmentation


#origIds <- segResults$ID
# correct for tumor content & remove normals
#segResults$ID <- str_remove(segResults$ID, "_MG.*")
#tcDf <- read.table("/mnt/DATA5/tmp/kev/misc/20210621fearonTc.txt", sep = "\t", stringsAsFactors = FALSE,
#                   header = TRUE)
#tcDf$samples[1:9] <- paste0("EF_D0", 1:9)

# for (i in 1:nrow(tcDf)) {
#   tc <- tcDf$Tumor.content[i]
#   sa <- tcDf$samples[i]
#   segResults$seg.mean[which(segResults$ID == sa)] <- segResults$seg.mean[which(segResults$ID == sa)]/tc
# }
# segResults$ID <- origIds

wt <- c("EF_D03_MG_X14", "EF_D12_MG_X15", "EF_D13_MG_X16", "EF_D20_MG_X55")
segResults <- segResults[-which(segResults$ID %in% wt),]

segZscores <- calcZscore(segResults)
segZfilt <- segZscoresFilt_zeroOut(segResults, segZscores)

tmpSeg <- cbind(segZfilt[,1:4], NA, segZfilt[,5])
colnames(tmpSeg) <- c("sampleID","chrom", "start.pos","end.pos", "n.probes", "mean")
tmpSeg <- tmpSeg[grep("EF_D", tmpSeg$sampleID),]

# add the arm/cna data here
armRes <- separateSegments_m(tmpSeg, mouseBedChromLim)
m_arm <- armRes[[1]]
m_cna <- armRes[[2]]

m_arm_freq_out <- getFreqData(m_arm)
m_arm_coad <- ampsDels(m_arm_freq_out)

m_arm_amp_bed <- reducingFreqBed(m_arm_coad[[1]], m_arm_coad[[2]])
m_arm_del_bed <- reducingFreqBed(m_arm_coad[[3]], m_arm_coad[[4]])

write.table(m_arm_amp_bed, "/mnt/DATA5/tmp/kev/misc/20210621_fearon_arm_amp.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(m_arm_del_bed, "/mnt/DATA5/tmp/kev/misc/20210621_fearon_arm_del.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


m_cna_freq_out <- getFreqData(m_cna)
m_cna_coad <- ampsDels(m_cna_freq_out)

m_cna_amp_bed <- reducingFreqBed(m_cna_coad[[1]], m_cna_coad[[2]])
m_cna_del_bed <- reducingFreqBed(m_cna_coad[[3]], m_cna_coad[[4]])

write.table(m_cna_amp_bed, "/mnt/DATA5/tmp/kev/misc/20210621_fearon_cna_amp.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(m_cna_del_bed, "/mnt/DATA5/tmp/kev/misc/20210621_fearon_cna_del.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




# I AM DUMB, dont need to downsample each time
# just downsample once, and use that file to subset using cyto_of_interest

# downsamping method works well plots for full version using copy-nubmer package
# /mnt/DATA5/tmp/kev/misc/210420tcga_ov_brca1_loss_copynumber.pdf
# /mnt/DATA5/tmp/kev/misc/20210420tcga_ov_freq_copynumber.pdf
tcga_coad_freq <- getFreqBed(amp = tcga_coad_amp, del = tcga_coad_del)


tmpGrange <- GRanges(seqnames = str_remove(tcga_coad_freq$chr, "h_chr"),
                     IRanges(start = tcga_coad_freq$start,
                             end = tcga_coad_freq$end))

allSynTable <- NULL

for(i in unique(segZfilt$ID)){
  segZfilt_sample <- segZfilt[which(segZfilt$ID == i),]
  syntenyBed <- syntenyPlotInputs(segZfilt_sample)
  human_bed <- syntenyBed$human_bed
  mouse_bed <- syntenyBed$mouse_bed
  
  # table for gene enrichment/conservation for both human and mouse
  tmpSynTable <- cbind("Sample" = i, human_bed, mouse_bed)
  allSynTable <- rbind(allSynTable, tmpSynTable)
  

  
  #make sure to floor and ceilign to -1 and 1 for cn later
  
  mouse_cn <- segZfilt_sample[,2:5]
  colnames(mouse_cn) <- c("chr", "start", "end", "cn")
  mouse_cn$chr <- paste0("m_chr", mouse_cn$chr)
  mouse_cn$chr <- str_replace(mouse_cn$chr, "m_chr23", "m_chrX")
  
  mouse_cn$ytop <- ifelse(mouse_cn$cn > 0 , mouse_cn$cn + 0.025, mouse_cn$cn) 
  mouse_cn$ybot <- ifelse(mouse_cn$cn > 0 , mouse_cn$cn, mouse_cn$cn - 0.025) 
  
  cyto_interest <- unique(c(human_bed$chr, mouse_bed$chr))
  cyto_combined_red <- cyto_hg38mm10[which(cyto_hg38mm10$V1 %in% cyto_interest),]
  
  
  
  # get freq bed per sample
  
  cn_track_bed <- rbind(tcga_coad_freq, mouse_cn)
  cn_track_bed <- cn_track_bed[which(cn_track_bed$chr %in% cyto_interest), ]
  # subset cn_track_bed amps and dels by cyto_of_interest
  
  # some weird artefacts in the freq plot for whole-genome
  # looks like an overall arm level measurement, so filter out
  
  #assign("tmpVar", cn_track_bed)
  #print(dim(cn_track_bed))
  #print(cn_track_bed[which((cn_track_bed$start - cn_track_bed$end) > 0),])
  
  track_color <- rep("#000000", nrow(cn_track_bed))
  track_color <- ifelse(cn_track_bed$cn < 0 , "#00008B", "#8B0000") 
  cn_track_bed$col = track_color
  
  mouse_chrs <- unique(mouse_bed$chr)
  human_chrs <- unique(human_bed$chr)
  
  colorVector2 <- rep("#000000", nrow(mouse_bed))
  for (j in seq_along(mouse_chrs)) {
    colorVector2[which(mouse_bed$chr == mouse_chrs[j])] <- colorVector[j]
  }
  
  
  # plotting each of the graphs
  #pdf(paste0("/mnt/DATA5/tmp/kev/testSynteny/","circos_ov_", i, ".pdf"),
  #    width = 6, height = 6, useDingbats = FALSE)
  #circos.par("track.height"= 0.10) + 
  #  circos.initializeWithIdeogram(cytoband = cyto_combined_red, sort.chr = FALSE,
  #                                plotType = c("ideogram", "labels")) + 
  #  circos.genomicTrack(cn_track_bed, ylim = c(-1.00,1.00),
  #                      panel.fun = function(region, value, ...) {
  #                        circos.genomicRect(region, value, col = value$col, 
  #                                           ytop = value$ytop,
  #                                           ybottom = value$ybot,
  #                                           border = NA,...)
  #                      }) + 
  #  circos.genomicLink(mouse_bed, human_bed,
  #                     col = colorVector2,
  #                     border = NA)
  #dev.off()
}



###
###
### easy way to get frequency and direction i.e frequency * amp/del

freqVector <- NULL
for (i in 1:nrow(allSynTable)) {
  tmpGrangeH <- GRanges(seqnames = str_remove(allSynTable$chr[i], "h_chr"),
                        IRanges(start = allSynTable$start[i],
                                end = allSynTable$end[i]))
  meanFreq <- mean(tcga_coad_freq$cn[subjectHits(findOverlaps(tmpGrangeH, tmpGrange))])
  if (is.nan(meanFreq)) {
    meanFreq <- 0
  }
  freqVector <- c(freqVector, meanFreq)
}

allSynTable$h_freq <- freqVector

write.table(allSynTable, "/mnt/DATA5/tmp/kev/misc/20210617fearonAllSyn.bed", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)

### freq vs freq plot - still need to make correct coordinates for freq plot
fearon_amp_bed <- cbind("sample" = "all", fearon_amp_bed)
fearon_del_bed <- cbind("sample" = "all", fearon_del_bed)
colnames(fearon_amp_bed) <- c("ID", "chrom", "loc.start", "loc.end", "seg.mean")
colnames(fearon_del_bed) <- c("ID", "chrom", "loc.start", "loc.end", "seg.mean")
fearon_del_bed$seg.mean <- fearon_del_bed$seg.mean * -1

# get rid of chromosmes with no changes
fearon_amp_bed <- fearon_amp_bed[-which(fearon_amp_bed$seg.mean == 0),]
fearon_del_bed <- fearon_del_bed[-which(fearon_del_bed$seg.mean == 0),]

ampSynteny <- syntenyPlotInputsFreq(fearon_amp_bed)
delSynteny <- syntenyPlotInputsFreq(fearon_del_bed)
mouse_bed <- rbind(ampSynteny$mouse_bed, delSynteny$mouse_bed)
human_bed <- rbind(ampSynteny$human_bed, delSynteny$human_bed)

allSynTable <- cbind(human_bed, mouse_bed)
colnames(allSynTable) <- c("h_chr", "h_start", "h_end", "m_chr",
                           "m_start", "m_end", "m_freq")

cyto_interest <- unique(c(human_bed$chr, mouse_bed$chr))
cyto_combined_red <- cyto_hg38mm10[which(cyto_hg38mm10$V1 %in% cyto_interest),]

mouse_cn <- cbind(mouse_bed, "ybot" = NA, "ytop" = NA)
mouse_cn$freq <- mouse_cn$freq/100
mouse_cn$ytop <- ifelse(mouse_cn$freq > 0, mouse_cn$freq, 0) 
mouse_cn$ybot <- ifelse(mouse_cn$freq < 0, mouse_cn$freq, 0)
colnames(mouse_cn) <- colnames(tcga_coad_freq)

cn_track_bed <- rbind(tcga_coad_freq, mouse_cn)
cn_track_bed <- cn_track_bed[which(cn_track_bed$chr %in% cyto_interest), ]

track_color <- rep("#000000", nrow(cn_track_bed))
track_color <- ifelse(cn_track_bed$cn < 0 , "#00008B", "#8B0000") 
cn_track_bed$col = track_color

mouse_chrs <- unique(mouse_bed$chr)
human_chrs <- unique(human_bed$chr)

colorVector2 <- rep("#000000", nrow(mouse_bed))
for (j in seq_along(mouse_chrs)) {
  colorVector2[which(mouse_bed$chr == mouse_chrs[j])] <- colorVector[j]
}


pdf(paste0("/mnt/DATA5/tmp/kev/misc/circos_all_coad.pdf"),
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

### get m_freqs and then do the enrichment analyses with both these freq values


freqVector <- NULL
for (i in 1:nrow(allSynTable)) {
  tmpGrangeH <- GRanges(seqnames = str_remove(allSynTable$h_chr[i], "h_chr"),
                        IRanges(start = allSynTable$h_start[i],
                                end = allSynTable$h_end[i]))
  meanFreq <- mean(tcga_coad_freq$cn[subjectHits(findOverlaps(tmpGrangeH, tmpGrange))])
  if (is.nan(meanFreq)) {
    print(i)
    meanFreq <- 0
  }
  freqVector <- c(freqVector, meanFreq)
}

allSynTable$h_freq <- freqVector
allSynTable$m_freq <- allSynTable$m_freq/100

write.table(allSynTable, "/mnt/DATA5/tmp/kev/misc/20210621fearonAllSyn.bed", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)


### run through once, then make a function
###
###

arm_amp_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210621_fearon_arm_amp.bed", sep = "\t",
                             stringsAsFactors = FALSE, header = TRUE)
arm_del_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210621_fearon_arm_del.bed", sep = "\t",
                             stringsAsFactors = FALSE, header = TRUE)

arm_amp_bed <- arm_amp_bed[which(arm_amp_bed$Freq > 0),]
arm_del_bed <- arm_del_bed[which(arm_del_bed$Freq > 0),]

tcga_coad_arm_amp_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210629_coadArmAmp.bed", sep = "\t",
                                    stringsAsFactors = FALSE, header = TRUE)
tcga_coad_arm_del_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210629_coadArmDel.bed", sep = "\t",
                                    stringsAsFactors = FALSE, header = TRUE)
tcga_coad_arm_del_bed$Freq  <- tcga_coad_arm_del_bed$Freq * -1

ampSynteny <- syntenyPlotInputsFreqV2(arm_amp_bed)
delSynteny <- syntenyPlotInputsFreqV2(arm_del_bed)
delSynteny$mouse_bed$freq <- delSynteny$mouse_bed$freq * -1

mouse_bed <- rbind(ampSynteny$mouse_bed, delSynteny$mouse_bed)
human_bed <- rbind(ampSynteny$human_bed, delSynteny$human_bed)
allSynTable <- cbind(human_bed, mouse_bed)
colnames(allSynTable) <- c("h_chr", "h_start", "h_end", "m_chr",
                           "m_start", "m_end", "m_freq")

cyto_interest <- unique(c(human_bed$chr, mouse_bed$chr))
cyto_combined_red <- cyto_hg38mm10[which(cyto_hg38mm10$V1 %in% cyto_interest),]

### make function for this too
trackBedColumns <- function(df){
  colnames(df) <- tolower(colnames(df))
  df <- cbind(df, "ybot" = NA, "ytop" = NA)
  df$freq <- df$freq/100
  df$ytop <- ifelse(df$freq > 0, df$freq, 0) 
  df$ybot <- ifelse(df$freq < 0, df$freq, 0)
  return(df)
}

# mouse_cn <- cbind(mouse_bed, "ybot" = NA, "ytop" = NA)
# mouse_cn$freq <- mouse_cn$freq/100
# mouse_cn$ytop <- ifelse(mouse_cn$freq > 0, mouse_cn$freq, 0) 
# mouse_cn$ybot <- ifelse(mouse_cn$freq < 0, mouse_cn$freq, 0)
# colnames(mouse_cn) <- colnames(tcga_coad_freq)

mouse_cn <- trackBedColumns(mouse_bed)
tcga_coad <- trackBedColumns(rbind(tcga_coad_arm_amp_bed ,
                                   tcga_coad_arm_del_bed ))
tcga_coad$chr <- paste0("h_chr", tcga_coad$chr)

cn_track_bed <- rbind(mouse_cn, tcga_coad)
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


pdf(paste0("/mnt/DATA5/tmp/kev/misc/circos_all_coad_arm.pdf"),
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


### testing function

arm_amp_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210621_fearon_arm_amp.bed", sep = "\t",
                          stringsAsFactors = FALSE, header = TRUE)
arm_del_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210621_fearon_arm_del.bed", sep = "\t",
                          stringsAsFactors = FALSE, header = TRUE)

tcga_coad_arm_amp_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210629_coadArmAmp.bed", sep = "\t",
                                    stringsAsFactors = FALSE, header = TRUE)
tcga_coad_arm_del_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210629_coadArmDel.bed", sep = "\t",
                                    stringsAsFactors = FALSE, header = TRUE)

coad_arm_allSynTable <- circosFreq(arm_amp_bed, arm_del_bed,tcga_coad_arm_amp_bed, tcga_coad_arm_del_bed,
                                   filename = "20210701fearonCoad_arm")



cna_amp_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210621_fearon_cna_amp.bed", sep = "\t",
                          stringsAsFactors = FALSE, header = TRUE)
cna_del_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210621_fearon_cna_del.bed", sep = "\t",
                          stringsAsFactors = FALSE, header = TRUE)

tcga_coad_cna_amp_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210629_coadCnaAmp.bed", sep = "\t",
                                    stringsAsFactors = FALSE, header = TRUE)
tcga_coad_cna_del_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210629_coadCnaDel.bed", sep = "\t",
                                    stringsAsFactors = FALSE, header = TRUE)

coad_cna_allSynTable <- circosFreq(cna_amp_bed, cna_del_bed, tcga_coad_cna_amp_bed, tcga_coad_cna_del_bed,
                                   filename = "20210701fearonCoad_cna")

