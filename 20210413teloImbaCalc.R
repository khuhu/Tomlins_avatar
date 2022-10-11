# interesting, the ucsc gap file confirms what jax says about
# there being no real short arm in mice since centromeres are 
# all on the edge of the of the chromosome - all short arms are about
# 10kb long ... 
# jax link: http://www.informatics.jax.org/glossary/centromere
mouseBedFile <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.gc.bed", stringsAsFactors = FALSE, sep = "\t",
                           header = FALSE)

all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                             ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                              end = zscore_gc_oe_ratios$EndPos))


mouseNormal <- c("MG_17X49", "MG_18X50", "MG_23X55", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

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


mm10_gap_coords <- read.table("/mnt/DATA5/tmp/kev/misc/mm10_ucsc_gap.txt",
                              header = FALSE, sep = "\t", stringsAsFactors = FALSE)

centromere_gaps <- mm10_gap_coords[which(mm10_gap_coords$V8 == "centromere"),]
telomere_gaps <- mm10_gap_coords[which(mm10_gap_coords$V8 == "telomere"),]
telomere_gaps2 <- telomere_gaps[-which(telomere_gaps$V3 == 0),]


mouseBedFile <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.gc.bed", stringsAsFactors = FALSE, sep = "\t",
                           header = FALSE)

telomerePositions <- NULL
for (i in unique(mouseBedFile$V1)) {
  tmpDf <- mouseBedFile[which(mouseBedFile$V1 == i),]
  tmpGap <- telomere_gaps[which(telomere_gaps$V2 == i),]
  teloPos <- max(tmpDf$V3)
  teloGap <- abs(teloPos - max(tmpGap$V4))/100000
  telomerePositions <- rbind(telomerePositions, c(i, teloPos, teloGap))
}

telomerePositions <- data.frame(telomerePositions, stringsAsFactors = FALSE)
colnames(telomerePositions) <- c("chr", "position", "diff")

# use cutoff of ~15Mb for each TAi count that are not whole chromosome
# https://github.com/sztup/scarHRD (size maybe later)
# the change needs to extend to a telomeric end 
# so this odd because there is technically only one chromosome arm
# in mice ... i.e p arm is nonexisant meaning generally max of one TAI
# per chromosome

# need to make sure i know how far the panel sequences i.e I don't 
# target telomeric regions, but I can count them as if they hit furthest
# amplicon

zscore_gc_oe_ratios <- read.table("/mnt/DATA5/tmp/kev/newMouse2/gcCorrectedCounts_matrix.txt", sep = "\t", stringsAsFactors = FALSE,
                                  header = TRUE)

load(file = "/mnt/DATA5/tmp/kev/misc/20210409mm10SegRes.Robj")
segZscores <- calcZscore(segRes)
segZfilt <- segZscoresFilt(segRes, segZscores)

# need different filts: paper used gains/losses as 2.5/1.5 copies
# or 0.3219281 and -0.4150375 respsectively

segZfilt_2 <- segZfilt[which(segZfilt$seg.mean > log2(2.5/2) | segZfilt$seg.mean < log2(1.5/2)),]


segZfilt_3 <- NULL
for (i in unique(segZfilt_2$chrom)) {
  tmpSeg <- segZfilt_2[which(segZfilt_2$chrom == i),]
  
  if (i == 23) {
    i <- "X"
  }
  tmpDf <- mouseBedFile[which(mouseBedFile$V1 == paste0("chr",i)),]
  
  teloPos <- max(tmpDf$V3)
  tmpSeg$diff <- abs(teloPos - tmpSeg$loc.end)/1000
  segZfilt_3 <- rbind(segZfilt_3, tmpSeg)
}

segZfilt_3$diff[which(segZfilt_3$diff < 1)] <- 0
table(segZfilt_3[which(segZfilt_3$diff == 0),]$ID)

### doing this odd, I should just count BPRN, BPN and maybe PRN
### when looking at ovarian tumors
BPN <- c("MG_4X36", "MG_6X38")
BPRN <- c("MG_5X37", "MG_12X44","MG_16X48", "MG_22X54")
PRN <- c("MG_36X68","MG_37X69", "MG_27X59")


#BPN <- c("2027lte", "2405lot")
#BPRN <- c("2163lot", "3807lot", "14085lot","14433lote")
#APA <- c("14593rt", "14595rt", "14596rt", "14594rt")
#PRN <- c("kcovc33", "kcovc41", "kcovc43")

targetedSeqDf <- data.frame("Sample" = c(BPN, BPRN, PRN),
                            "Genotpye" = c(rep("BPN", length(BPN)),
                                           rep("BPRN", length(BPRN)),
                                           rep("PRN", length(PRN))),
                            "TAi"= 0)

targetedTAi <- table(segZfilt_3[which(segZfilt_3$diff == 0),]$ID)

tmpMatchIdx <- which(targetedSeqDf$Sample %in% names(targetedTAi))
targetedSeqDf$TAi[tmpMatchIdx] <- targetedTAi[match(targetedSeqDf$Sample[tmpMatchIdx], names(targetedTAi))]



# for loop to calculate per sample TAi
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




zhangTAi <- NULL
for (i in unique(zhangAllCnvs$V1)) {
  tmpSeg <- zhangAllCnvs[which(zhangAllCnvs$V1 == i),]
  
  if (i == 23) {
    i <- "X"
  }
  
  tmpDf <- telomere_gaps2[which(telomere_gaps2$V2 == i),]
  tmpSeg$diff <- abs(tmpDf$V4 - tmpSeg$V3)/1000
  zhangTAi <- rbind(zhangTAi, tmpSeg)
}

zhangSRA <- c("SRR12252465", "SRR12252466", "SRR12252467", "SRR12252468", "SRR12252469",
  "SRR12252470", "SRR12252471", "SRR12252472", "SRR12252473", "SRR12252474",
  "SRR12252475", "SRR12252476")
zhangGeno <- c("Tp53-Ccne1-Akt2-Kras-2", "Tp53-Ccne1-Akt2-Kras-1", "Tp53-Brca1-Myc-3",
  "Tp53-Brca1-Myc-1", "Tp53-Brca1-1", "Tp53-4", "Tp53-3", "Tp53-Brca1-2",
  "Tp53-Pten-Nf1-2", "Tp53-Pten-Nf1-1", "Control-4", "Control-3")
zhangGeno2 <- c("WtBrca1", "WtBrca1", "Brca1Mut",
                "Brca1Mut", "Brca1Mut", "WtBrca1", "WtBrca1", "Brca1Mut",
                "WtBrca1", "WtBrca1", "Control", "Control")

zhangTAiDf <- data.frame("SRA" = zhangSRA, "Geno1" = zhangGeno, "Geno2" = zhangGeno2,
                    "TAi" = 0, stringsAsFactors = FALSE)



zhangTAi_2 <- table(zhangTAi[which(zhangTAi$diff == 0),]$Sample)
tmpMatchIdx_z <- which(zhangTAiDf$SRA %in% names(zhangTAi_2))
zhangTAiDf$TAi[tmpMatchIdx_z] <- zhangTAi_2[match(zhangTAiDf$SRA[tmpMatchIdx], names(zhangTAi_2))]

#colnames(targetedSeqDf) <- NULL
#colnames(zhangTAiDf) <- NULL 
#rbind(targetedSeqDf, zhangTAiDf[,c(1,3,4)])

ggplot(targetedSeqDf) + geom_boxplot(aes(x = Genotpye, y = TAi))
ggplot(zhangTAiDf) + geom_boxplot(aes(x = Geno2, y = TAi))

combinedTAi <- data.frame("Geno" = c(zhangTAiDf$Geno2, as.character(targetedSeqDf$Genotpye)),
                          "TAi" = c(zhangTAiDf$TAi, targetedSeqDf$TAi),
                          stringsAsFactors = FALSE)
ggplot(combinedTAi) + geom_boxplot(aes(x = Geno, y = TAi))
