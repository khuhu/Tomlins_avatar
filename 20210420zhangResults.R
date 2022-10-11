# synteny analysis for outside data - zhang data

library(GenomicFeatures)
library(dplyr)
library(dbplyr)
library(circlize)
library(stringr)
library(DNAcopy)
library(stringr)
library(optparse)


synteny_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_syntenyDf.txt",
                               sep = "\t", stringsAsFactors = FALSE, header = TRUE)
cyto_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_cytoDf.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = FALSE)




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
  
  # when comparing to human, change to human specific cn
  
  mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg38mm10$comp_chr[query_hits]),
                          "start" = synteny_hg38mm10$comp_start_pos[query_hits],
                          "end" = synteny_hg38mm10$comp_end_pos[query_hits],
                          "cn" = segInput$seg.mean[subject_hits],
                          stringsAsFactors = FALSE)
  
  human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg38mm10$ref_chr[query_hits]),
                          "start" = synteny_hg38mm10$ref_start_pos[query_hits],
                          "end" = synteny_hg38mm10$ref_end_pos[query_hits],
                          "cn" = segInput$seg.mean[subject_hits],
                          stringsAsFactors = FALSE)
  
  
  
  resList <- list("human_bed" = human_bed,
                  "mouse_bed" = mouse_bed)
  
  return(resList)
}


tmpDir <- "/mnt/DATA5/tmp/kev/tmpDbs/SRA/cancDiscZhang/"
setwd(tmpDir)
cnaFiles <- system('ls *CNVs*', intern = TRUE)
zhangAllCnvs <- NULL
for (i in cnaFiles) {
  tmpTable <- read.table(paste0(tmpDir, i))
  tmpTable$V1 <- paste0("chr", tmpTable$V1)
  tmpTable$diff <- tmpTable$V3 - tmpTable$V2
  tmpTable <- tmpTable[which(tmpTable$diff > 10000000),]
  tmpTable <- tmpTable[which(tmpTable$V4 > 2.5 | tmpTable$V4 < 1.5),]
  tmpTable2 <- cbind( "Sample" = str_remove_all(i, "_.*"), tmpTable)
  zhangAllCnvs <- rbind(zhangAllCnvs, tmpTable2)
}


zhang_segRes <- zhangAllCnvs[,1:5]
colnames(zhang_segRes) <- c("ID", "chrom", "loc.start", "loc.end", "seg.mean")

# let's try profiles with and without brca mutation

