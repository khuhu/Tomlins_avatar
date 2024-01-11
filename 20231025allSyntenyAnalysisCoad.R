### (1) general copy number script from: 20230630allSampAneruploidyCna.R
### (2) permutation script from 20230725syntenyRnaAnalysis.R
### (3) rna analysis and processing from 20230830syntenyRnaAnalysis.R
### (4) pathway analysis script  is on brahm and ran only on brahm
### (5) most of the figures are from 20230823figureOrderScott.R
### if plot no on (5), then found in 1-4

### note: if rmd is made later, then can just split below code into block based on processing mouse, human, graphs, permutation etc

source("/home/kevhu/scripts/20210802syntenyFunctions.R")
source("/home/kevhu/scripts/20230713newSyntenyDfFunctions.R")
source("/home/kevhu/scripts/20230713newSyntenyDfFunctionV2.R")

### q should be ploidy
segResFilt4 <- read.table("/mnt/DATA5/tmp/kev/misc/20230907segRes2.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
segResFilt4$chrom <- paste0("chr", segResFilt4$chrom)
segResFilt4$size <- segResFilt4$loc.end - segResFilt4$loc.start

allNgsDis_all <- disjoin(GRanges(seqnames = segResFilt4$chrom,
                                 IRanges(start = segResFilt4$loc.start, end = segResFilt4$loc.end)))

### all there minus x chromosome I forgot to remove

allNgsDisDf_all <- data.frame(allNgsDis_all)
allNgsDisDf_all <- allNgsDisDf_all[-which(allNgsDisDf_all$width < 2),]

allNgsDisDf_allGr <- GRanges(seqnames = allNgsDisDf_all$seqnames,
                             IRanges(start = allNgsDisDf_all$start, end = allNgsDisDf_all$end))

ngsDisjoin2 <- NULL
for (i in unique(segResFilt4$ID)) {
  print(i)
  tmpNgs <- segResFilt4[which(segResFilt4$ID == i),]
  tmpRes <- allNgsDisDf_all[,1:3]
  tmpRes$sample <- i
  queryIdx <- 1:nrow(tmpRes)
  tmpGrange <- GRanges(seqnames = tmpNgs$chrom, IRanges(start = tmpNgs$loc.start, end = tmpNgs$loc.end))
  tmpIdx <- subjectHits(findOverlaps(query = allNgsDisDf_allGr, subject = tmpGrange))
  testQuery <- queryHits(findOverlaps(query = allNgsDisDf_allGr, subject = tmpGrange))
  
  
  if (length(which(duplicated(queryHits(findOverlaps(query = allNgsDisDf_allGr, subject = tmpGrange))))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated(queryHits(findOverlaps(query = allNgsDisDf_allGr, subject = tmpGrange))))]
    testQuery <- testQuery[-which(duplicated(queryHits(findOverlaps(query = allNgsDisDf_allGr, subject = tmpGrange))))]
  }
  ### error probably happens b/c there is no match for one of the queries
  tmpRes$cn <- NA
  tmpRes$cn[testQuery] <- tmpNgs$seg.mean[tmpIdx]
  tmpRes2 <- tmpRes[, c("sample", "seqnames", "start", "end", "cn")]
  ngsDisjoin2  <- rbind(ngsDisjoin2 , tmpRes2)
}

colnames(ngsDisjoin2)[c(2,5)] <- c("chrom", "seg.mean")
ngsDisjoin2$string <- paste(ngsDisjoin2$chrom, ngsDisjoin2$start, ngsDisjoin2$end)
ngsDisjoin2$size <- ngsDisjoin2$end - ngsDisjoin2$start


tmp_cna <- NULL
tmpDisJoinA <- read.table("/mnt/DATA5/tmp/kev/misc/20230907ngsInter_0_10Mb.bed",
                          sep = "\t", header = FALSE, stringsAsFactors = FALSE)

tmpDisJoin <- read.table(paste0("/mnt/DATA5/tmp/kev/misc/20230907ngsInter_", 3, "_", 10, "Mb.bed"),
                         sep = "\t", header = FALSE, stringsAsFactors = FALSE)

colnames(tmpDisJoin) <- c("chrom", "start", "end","sample", "seg.mean", "str", "size")
colnames(tmpDisJoinA) <- c("chrom", "start", "end","sample", "seg.mean", "str", "size")
tmpDisJoinAneuploidy <- matrix(0, nrow = length(unique(segResFilt4$ID)), ncol = length(unique(segResFilt4$chrom)))
rownames(tmpDisJoinAneuploidy) <- unique(segResFilt4$ID)
colnames(tmpDisJoinAneuploidy) <- unique(segResFilt4$chrom)
tmpCnaAneu <- NULL

k <- unique(segResFilt4$ID)[1]
for (k in rownames(tmpDisJoinAneuploidy)) {
  tmpDf <- tmpDisJoin[grep(k, tmpDisJoin$sample),]
  tmpDfA <- tmpDisJoinA[grep(k, tmpDisJoinA$sample),]
  tmpChrList <- unique(tmpDf$chrom)
  
  for (l in seq_along(tmpChrList)) {
    tmpChr <- tmpDf[which(tmpDf$chrom == tmpChrList[l]), ]
    tmpChr2 <- tmpChr
    tmpChrA <- tmpDfA[which(tmpDfA$chrom == tmpChrList[l]), ]
    tmpChr3 <- tmpChrA
    chrThres <- sum(tmpChr$size) * 0.7
    tmpAneu <- "none"
    if (sum(tmpChrA$size[which(tmpChrA$seg.mean > 0.2)]) > chrThres) {
      tmpAneu <- "gain"
      tmpChr2$seg.mean[which(tmpChr2$seg.mean > 0.2)] <- 0
      tmpChr3$seg.mean[which(tmpChr3$seg.mean < 0.2)] <- 0
    } else if(sum(tmpChrA$size[which(tmpChrA$seg.mean < -0.2)]) > chrThres){
      tmpAneu <- "loss"
      tmpChr2$seg.mean[which(tmpChr2$seg.mean < -0.2)] <- 0
      tmpChr3$seg.mean[which(tmpChr3$seg.mean > -0.2)] <- 0
    } else{
      tmpAneu <- "none"
    }
    tmpDisJoinAneuploidy[k, l] <- tmpAneu
    tmp_cna <- rbind(tmp_cna, tmpChr2)
    tmpCnaAneu <- rbind(tmpCnaAneu, tmpChr3)
  }
}

tmp_cna$seg.mean[which(abs(tmp_cna$seg.mean) < 0.2)] <- 0
tmp_cna$seg.mean[which(is.na(tmp_cna$seg.mean))] <- 0

tmpCnaAneu$seg.mean[which(abs(tmpCnaAneu$seg.mean) < 0.2)] <- 0
tmpCnaAneu$seg.mean[which(is.na(tmpCnaAneu$seg.mean))] <- 0

tmp_cna2 <- NULL
for (k in unique(tmp_cna$sample)) {
  tmpDf <- tmp_cna[which(tmp_cna$sample == k), ]
  l <- unique(tmpDf$chrom)[16]
  for (l in unique(tmpDf$chrom)) {
    tmpDfChr <- tmpDf[which(tmpDf$chrom == l),]
    if (any(tmpDfChr$seg.mean != 0)) {
      
      tmpDfChr2 <- tmpDfChr[-which(tmpDfChr$seg.mean == 0), ]
      tmpDfChr <- tmpDfChr[which(tmpDfChr$seg.mean == 0), ]
      
      for (m in unique(tmpDfChr2$seg.mean)) {
        tmpDfChr3 <- tmpDfChr2[which(tmpDfChr2$seg.mean == m),]
        tmpChrSeg <- tmpDfChr3[1, 1:6]
        tmpChrSeg$start <- min(tmpDfChr3$start)
        tmpChrSeg$end <- max(tmpDfChr3$end)
        tmpChrSeg$size <- tmpChrSeg$end - tmpChrSeg$start
        tmpDfChr <- rbind(tmpDfChr, tmpChrSeg)
      }
      tmpDfChr <- tmpDfChr[order(tmpDfChr$start), ]
      tmp_cna2 <- rbind(tmp_cna2, tmpDfChr)
    } else{
      tmp_cna2 <- rbind(tmp_cna2, tmpDfChr)
    }
  }
}

allMouseAneuploidy <- tmpDisJoinAneuploidy[,1:19]
allMouseCna <- tmp_cna2
allMouseAneu <- tmpCnaAneu


### using disjoin to find p53 deletion signal for minimum tumor content
### if assuming diploid for ploidy of all samples - less complicated b/c cutoffs get weird at different ploidy
### i.e cutoffs would be too much and we can't tell ploidy. so >30% purity abs(>30% cutoff) log2R ~ -0.5
### fearon samples a bit different since I look at both adenomas and adenocarcinomas
### only APC is required for adenoma creation while adenocar needs apc, kras and tp53
### filter adenocar by p53 and adenoma is difficult since we can't accurately detect APC
### assuming they are adenomas, we'll have a minimum an arm change

p53String <- "chr11 69588341 69588439"
apcString <- "chr18 34310771 34310875"

ngsDisjoinP53 <- ngsDisjoin2[which(ngsDisjoin2$string == p53String),]
ngsDisjoinApc <- ngsDisjoin2[which(ngsDisjoin2$string == apcString),]
goodTcHgsc <- ngsDisjoinP53$sample[which(ngsDisjoinP53$seg.mean < -0.5)]

goodTcCoadAdenoCar <- ngsDisjoinP53$sample[which(ngsDisjoinP53$seg.mean < -0.2)]
goodTcCoadAdenoCar <- goodTcCoadAdenoCar[grep("efd", goodTcCoadAdenoCar)]
goodTcCoadAdenoCar2 <- str_remove(goodTcCoadAdenoCar, "x.*")

goodTcCoadAdeno <- NULL
for (i in unique(ngsDisjoin2$sample)) {
  tmpDf <- ngsDisjoin2[which(ngsDisjoin2$sample == i), ]
  tmpDecile <- quantile(tmpDf$seg.mean, 0.05, na.rm = TRUE)
  names(tmpDecile) <- i
  goodTcCoadAdeno <- c(goodTcCoadAdeno, tmpDecile)
}
goodTcCoadAdeno <- goodTcCoadAdeno[grep("efd", names(goodTcCoadAdeno))]
goodTcCoadAdeno <- goodTcCoadAdeno[which(goodTcCoadAdeno < -0.2)]
goodTcCoadAdeno2 <- str_remove(names(goodTcCoadAdeno), "x.*")


mm <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/cytoBand.txt.gz", 
            col.names = c("chrom","chromStart","chromEnd","name","gieStain"))

mm10ChromArms <- mm[ , .(length = sum(chromEnd - chromStart)), 
                     by = .(chrom, arm = substring(name, 1, 1)) ]

tmpChrBoundaries <- data.frame("chr" = mm10ChromArms$chrom, "start" = 1,
                               "end" = mm10ChromArms$length)

tmpChrBoundaries <- tmpChrBoundaries[which(tmpChrBoundaries$chr %in% paste0("chr", 1:19)), ]
tmpChrBoundaries <- tmpChrBoundaries[match(paste0("chr", 1:19), tmpChrBoundaries$chr), ]


allMouseAneuploidyMelt <- melt(allMouseAneuploidy)
allMouseAneuploidyNum <- NULL
for (i in unique(allMouseAneuploidyMelt$Var1)) {
  tmp <- allMouseAneuploidyMelt[which(allMouseAneuploidyMelt$Var1 == i), ]
  tmpRes <- data.frame("sample" = i, tmpChrBoundaries)
  tmpRes$seg.mean <- 0
  tmpRes$seg.mean[which(tmp$value == "loss")] <- -1
  tmpRes$seg.mean[which(tmp$value == "gain")] <- 1
  allMouseAneuploidyNum <- rbind(allMouseAneuploidyNum, tmpRes)
}


allMouseAneu_coad <- allMouseAneuploidyNum[grep("efd", allMouseAneuploidyNum$sample),]
allMouseAneu_coad2 <- data.frame("sampleID" = allMouseAneu_coad$sample, "chrom" = as.numeric(str_remove(allMouseAneu_coad$chr, "chr")),
                                 "start.pos" = allMouseAneu_coad$start, "end.pos" = allMouseAneu_coad$end,
                                 "n.probes" = rep(NA , length(allMouseAneu_coad$sample)), "mean" = allMouseAneu_coad$seg.mean,
                                 "str" = paste0(allMouseAneu_coad$sample, str_remove("chr", allMouseAneu_coad$chr), allMouseAneu_coad$start, allMouseAneu_coad$end),
                                 "length" = allMouseAneu_coad$end - allMouseAneu_coad$start)
allMouseAneu_coad2$sampleID <- str_remove(allMouseAneu_coad2$sampleID, "x.*")
allMouseAneu_coad2 <- allMouseAneu_coad2[which(allMouseAneu_coad2$sampleID %in% goodTcCoadAdenoCar2), ]


### read in different colorectal annotation and fix pathology

coadHistology <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20211207yingSecondSetwithHistology.xlsx", sheet = 1)
annoTableCoad <- read.table("/mnt/DATA5/tmp/kev/misc/20220726FearonAnnoAfDelFreq.txt", sep = "\t",
                            stringsAsFactors = FALSE, header = FALSE)
annoTableCoad$histology <- NA
annoTableCoad$histology[c(14, 15, 17,18, 21, 22, 24, 25, 27:33)] <- "met"
### don't know if these primaries are adenomas or adenocarcinomas
annoTableCoad$histology[c(16, 19, 23, 26)] <- NA
annoTableCoad$histology[34:94] <- coadHistology$Histology

annoTableCoad$geno <- "wildtype"
annoTableCoad$geno[which(annoTableCoad$V4 == "CDX2P-CreERT2 Apcfl/+, KrasLSLG12D/+, p53R270H+/ex2-10 fl")] <- "AKP"
annoTableCoad$geno[grep("p53 ex2-10 fl/fl", annoTableCoad$V4)] <- "AKP"
annoTableCoad$geno[grep("p53R270H fl/ex2-10", annoTableCoad$V4)] <- "AKP"
annoTableCoad$geno[grep("p53R270H/ex2-10", annoTableCoad$V4)] <- "AKP"
annoTableCoad$geno[grep("p532-10 -/-", annoTableCoad$V4)] <- "AKP"
annoTableCoad$geno[grep("Sox9 fl/+", annoTableCoad$V4)] <- "AKPS"
annoTableCoad$geno[grep("KrasLSLG12D\\/\\+$", annoTableCoad$V4)] <- "AK"
annoTableCoad$geno[grep("CDX2 fl\\/fl Braf V600E\\/\\+$", annoTableCoad$V4)] <- "CB"
# badCoadSamples <- paste0("efd", c(14, 28, 34:36, 44, 50, 58:61, 65, 67, 69, 76:78, 85, 90, 94))
badCoadSamples <- paste0("efd", c(34:35, 44, 65, 94, 58, 88, 70, 59, 76))
annoTableCoad2 <- annoTableCoad[-which(annoTableCoad$V6 %in% badCoadSamples), ]

annoTableCoad2$histology[which(annoTableCoad2$V1 %in% paste0("EF-", c("D1", "D2", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11",
                                                                      "D16", "D19", "D23", "D26")))] <- "Adenocarcinoma"

### need to fix this later so I get all samples .. for now the 105 samples should be okay
allMouseCna_coad <- allMouseCna[grep("efd", allMouseCna$sample),]
allMouseCna_coad2 <- data.frame("sampleID" = allMouseCna_coad$sample, "chrom" = as.numeric(str_remove(allMouseCna_coad$chr, "chr")),
                                "start.pos" = allMouseCna_coad$start, "end.pos" = allMouseCna_coad$end,
                                "n.probes" = rep(NA , length(allMouseCna_coad$sample)), "mean" = allMouseCna_coad$seg.mean,
                                "str" = paste0(allMouseCna_coad$sample, str_remove("chr", allMouseCna_coad$chr), allMouseCna_coad$start, allMouseCna_coad$end),
                                "length" = allMouseCna_coad$end - allMouseCna_coad$start)
allMouseCna_coad2$sampleID <- str_remove(allMouseCna_coad2$sampleID, "x.*")


for (i in unique(allMouseCna_coad2$sampleID)) {
  tmp <- allMouseCna_coad2[which(allMouseCna_coad2$sampleID == i),]
  uchr <- unique(tmp$chrom)
  if(length(setdiff(1:19,uchr)) > 0){
    for (j in setdiff(1:19, uchr)) {
      allMouseCna_coad2 <- rbind(allMouseCna_coad2,
                                 data.frame("sampleID" = i, "chrom" = j,
                                            "start.pos" = 1, "end.pos" = 10000,
                                            "n.probes" = NA,
                                            "str" = NA, "mean" = 0,
                                            "length" = 10000))
    }
  }
}


allMouseCna_adeno <- allMouseCna_coad2[which(allMouseCna_coad2$sampleID %in% annoTableCoad2$V6[which(annoTableCoad2$histology == "Adenoma")]), ]
allMouseCna_adeno <- allMouseCna_adeno[which(allMouseCna_adeno$sampleID %in% goodTcCoadAdeno2),]
allMouseCna_adenoCar <- allMouseCna_coad2[which(allMouseCna_coad2$sampleID %in% annoTableCoad2$V6[which(annoTableCoad2$histology == "Adenocarcinoma")]), ]
allMouseCna_adenoCar <- allMouseCna_adenoCar[which(allMouseCna_adenoCar$sampleID %in% goodTcCoadAdenoCar2),]


### above is large cna, below aneuploidy
allMouseAneu_adeno <- allMouseAneu_coad2[which(allMouseAneu_coad2$sampleID %in% annoTableCoad2$V6[which(annoTableCoad2$histology == "Adenoma")]), ]
allMouseAneu_adeno <- allMouseAneu_adeno[which(allMouseAneu_adeno$sampleID %in% goodTcCoadAdeno2),]
allMouseAneu_adenoCar <- allMouseAneu_coad2[which(allMouseAneu_coad2$sampleID %in% annoTableCoad2$V6[which(annoTableCoad2$histology == "Adenocarcinoma")]), ]
allMouseAneu_adenoCar <- allMouseAneu_adenoCar[which(allMouseAneu_adenoCar$sampleID %in% goodTcCoadAdenoCar2),]


### before I got filtered out events < 10%; important to have it for permutation method

allMouseAneu_adeno_freq <- getFreqData(allMouseAneu_adeno)
allMouseAneu_adeno_freq_res <- ampsDels(allMouseAneu_adeno_freq)
allMouseAneu_adeno_amp_bed <- reducingFreqBed2(allMouseAneu_adeno_freq_res[[1]])
allMouseAneu_adeno_del_bed <- reducingFreqBed2(allMouseAneu_adeno_freq_res[[3]])


allMouseAneu_adenoCar_freq <- getFreqData(allMouseAneu_adenoCar)
allMouseAneu_adenoCar_freq_res <- ampsDels(allMouseAneu_adenoCar_freq)
allMouseAneu_adenoCar_amp_bed <- reducingFreqBed2(allMouseAneu_adenoCar_freq_res[[1]])
allMouseAneu_adenoCar_del_bed <- reducingFreqBed2(allMouseAneu_adenoCar_freq_res[[3]])



### annotation of apc and other mutations for GISTIC
coad_anno <- read.table("/mnt/DATA5/tmp/kev/misc/alterations_across_samples_coad.tsv", sep = "\t",
                        header = TRUE)
coad_anno$ApcCheck <- paste(coad_anno$APC, coad_anno$APC..MUT, coad_anno$APC..HOMDEL)
coad_anno2 <- coad_anno[grep("driver", coad_anno$ApcCheck), ]


hg19cyto <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genome/hg19/cytoBand.txt", sep = "\t",
                       col.names = c("chrom","chromStart","chromEnd","name","gieStain"))

broadBySampleGisticCoad <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/gdac.broadinstitute.org_COADREAD-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/broad_values_by_arm.txt",
                                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

hg19ArmLocations <- NULL
for (i in unique(hg19cyto$chrom)) {
  tmp <- hg19cyto[which(hg19cyto$chrom == i),]
  tmpP <- tmp[grep("p", tmp$name), ]
  tmpQ <- tmp[grep("q", tmp$name), ]
  hg19ArmLocations <- rbind(hg19ArmLocations, data.frame("chrom" = rep(i, 2), "arm" = c("p", "q"),
                                                         "start" = c(min(tmpP$chromStart), min(tmpQ$chromStart)),                                                        "end" = c(max(tmpP$chromEnd), max(tmpQ$chromEnd))))
}

hg19ArmLocations$chrStripped <- as.numeric(str_remove(hg19ArmLocations$chrom, "chr"))
hg19ArmLocations <- hg19ArmLocations[which(hg19ArmLocations$chrStripped %in% 1:22),]
hg19ArmLocations <- hg19ArmLocations[order(hg19ArmLocations$chrStripped,decreasing = FALSE),]
hg19ArmLocations$str <- paste0(hg19ArmLocations$chrStripped, hg19ArmLocations$arm)
hg19ArmLocations <- hg19ArmLocations[which(hg19ArmLocations$str %in% broadBySampleGisticCoad$`Chromosome Arm`), ]


### loading GISTIC coadread data for aneuploidy, filtering and processing for synteny table


broadBySampleGisticCoad <- broadBySampleGisticCoad[, c(1, grep(paste0(coad_anno2$Sample.ID, collapse = "|"), colnames(broadBySampleGisticCoad)))]
broadBySampleGisticCoad <- broadBySampleGisticCoad[-grep("X", broadBySampleGisticCoad$`Chromosome Arm`), ]
broadBySampleGisticMatCoad <- broadBySampleGisticCoad[, 2:ncol(broadBySampleGisticCoad)]
broadBySampleGisticCoad2 <- cbind(hg19ArmLocations[, c("chrStripped", "start", "end")], broadBySampleGisticMatCoad)
broadBySampleGisticMeltCoad <- melt(broadBySampleGisticCoad2, id.vars = c("chrStripped", "start", "end"))
broadBySampleGisticMeltCoad2 <- data.frame("sampleID" = broadBySampleGisticMeltCoad$variable, "chrom" = broadBySampleGisticMeltCoad$chrStripped,
                                           "start.pos" = broadBySampleGisticMeltCoad$start, "end.pos" = broadBySampleGisticMeltCoad$end, 
                                           "n.probes" = NA, "mean" = broadBySampleGisticMeltCoad$value)
broadBySampleGisticMeltCoad2$str <- paste0(broadBySampleGisticMeltCoad2$sampleID, broadBySampleGisticMeltCoad2$chrom, 
                                           broadBySampleGisticMeltCoad2$start.pos, broadBySampleGisticMeltCoad2$end.pos)
broadBySampleGisticMeltCoad2$length <- broadBySampleGisticMeltCoad2$end.pos - broadBySampleGisticMeltCoad2$start.pos


armGisticCoadApc_freq <- getFreqData(broadBySampleGisticMeltCoad2)
armGisticCoadApc_freq_res <- ampsDels(armGisticCoadApc_freq)
armGisticCoadApc_amp_bed <- reducingFreqBed2(armGisticCoadApc_freq_res[[1]])
armGisticCoadApc_del_bed <- reducingFreqBed2(armGisticCoadApc_freq_res[[3]])
# armGisticCoadApc_amp_bed$Freq[which(armGisticCoadApc_amp_bed$Freq < 10)] <- 0
# armGisticCoadApc_del_bed$Freq[which(armGisticCoadApc_del_bed$Freq < 10)] <- 0


coad_adeno_arm_allSynTable <- circosFreq2(allMouseAneu_adeno_amp_bed, allMouseAneu_adeno_del_bed, armGisticCoadApc_amp_bed,
                                          armGisticCoadApc_del_bed, filename = "20231025fearon_aneuAll_coad_adeno", ref = "human")

coad_adenoCar_arm_allSynTable <- circosFreq2(allMouseAneu_adenoCar_amp_bed, allMouseAneu_adenoCar_del_bed, armGisticCoadApc_amp_bed,
                                             armGisticCoadApc_del_bed, filename = "20231025fearon_aneuAll_coad_adenoCar", ref = "human")

### start of permutation

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


coad_adenoCar_arm_allSynTable$str <- paste(coad_adenoCar_arm_allSynTable$h_chr,
                                           coad_adenoCar_arm_allSynTable$h_start,
                                           coad_adenoCar_arm_allSynTable$h_end)
apcCoadArms <- broadBySampleGisticCoad2
tmpSynTable_adenoCar <- coad_adenoCar_arm_allSynTable[-which(duplicated(coad_adenoCar_arm_allSynTable$str)), ]
tmpSynTable_adenoCar$h_chr <- as.numeric(str_remove(tmpSynTable_adenoCar$h_chr, "h_chr"))
tmpSynTable_adenoCar$m_chr <- as.numeric(str_remove(tmpSynTable_adenoCar$m_chr, "m_chr"))
tmpSynTable_adenoCar <- tmpSynTable_adenoCar[order(tmpSynTable_adenoCar$h_chr, tmpSynTable_adenoCar$h_start, decreasing = FALSE),]
apcCoadArmsGr <- GRanges(seqnames = apcCoadArms$chrStripped, IRanges(start = apcCoadArms$start, end = apcCoadArms$end))
coadSyntenyGrange <- GRanges(seqnames = tmpSynTable_adenoCar$h_chr, IRanges(start = tmpSynTable_adenoCar$h_start, end = tmpSynTable_adenoCar$h_end))
coadMm10SyntenyGrange <- GRanges(seqnames = tmpSynTable_adenoCar$m_chr, IRanges(start = tmpSynTable_adenoCar$m_start, end = tmpSynTable_adenoCar$m_end))

humanCoadRegionMat <- data.frame("chr" = tmpSynTable_adenoCar$h_chr, "start" = tmpSynTable_adenoCar$h_start, "end" = tmpSynTable_adenoCar$h_end)
for (i in 4:ncol(apcCoadArms)) {
  
  tmpIdx <- subjectHits(findOverlaps(query = coadSyntenyGrange, subject = apcCoadArmsGr))
  testQuery <- queryHits(findOverlaps(query = coadSyntenyGrange, subject = apcCoadArmsGr))
  
  
  if (length(which(duplicated(testQuery))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated(testQuery))]
    testQuery <- testQuery[-which(duplicated(testQuery))]
  }
  
  ### error probably happens b/c there is no match for one of the queries
  humanCoadRegionMat[ ,i] <- NA
  humanCoadRegionMat[ ,i] <- apcCoadArms[tmpIdx, i]
}


humanCoadRegionMat2 <- humanCoadRegionMat[ , 4:ncol(humanCoadRegionMat)]
cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
synCoadPerm <- foreach(i=1:25,
                       .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                         tmpAmpAll <- NULL
                         tmpDelAll <- NULL
                         tmpAmpFreq <- NULL
                         tmpDelFreq <- NULL
                         j <- 0
                         while (j < 400) {
                           tmp <- apply(humanCoadRegionMat2, 2, function(x) sample(x, nrow(humanCoadRegionMat2)))
                           ampFreq <- apply(tmp, 1, function(x) length(which(x > 0.2))/length(x))
                           delFreq <- apply(tmp, 1, function(x) length(which(x < -0.2))/length(x))
                           tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0.2)]) * length(which(x > 0.2))/length(x))
                           tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < -0.2)])) * length(which(x < -0.2))/length(x))
                           tmpAmpAll <- cbind(tmpAmpAll, tmpAmp)
                           tmpDelAll <- cbind(tmpDelAll, tmpDel)
                           tmpAmpFreq <- cbind(tmpAmpFreq, ampFreq)
                           tmpDelFreq <- cbind(tmpDelFreq, delFreq)
                           j <- j + 1
                         }
                         return(list(tmpAmpAll, tmpDelAll, tmpAmpFreq, tmpDelFreq))
                       }


stopCluster(cl)
print( Sys.time() - start )

humanPermAmpCoadAdenoCarStat <- do.call(cbind, synCoadPerm[[1]])
humanPermDelCoadAdenoCarStat <- do.call(cbind, synCoadPerm[[2]])
humanPermAmpCoadAdenoCarFreq <- do.call(cbind, synCoadPerm[[3]])
humanPermDelCoadAdenoCarFreq <- do.call(cbind, synCoadPerm[[4]])


mouseCoadRegionMat_adenoCar <- data.frame("chr" = tmpSynTable_adenoCar$m_chr, "start" = tmpSynTable_adenoCar$m_start, "end" = tmpSynTable_adenoCar$m_end)
for (i in seq_along(unique(allMouseAneu_adenoCar$sampleID))) {
  
  tmp <- allMouseAneu_adenoCar[which(allMouseAneu_adenoCar$sampleID == unique(allMouseAneu_adenoCar$sampleID)[i]),]
  tmpGrange <- GRanges(seqnames = tmp$chrom, IRanges(start = tmp$start.pos, end = tmp$end.pos))
  
  print(unique(allMouseAneu_adenoCar$sampleID)[i])
  
  tmpIdx <- subjectHits(findOverlaps(query = coadMm10SyntenyGrange, subject = tmpGrange))
  testQuery <- queryHits(findOverlaps(query = coadMm10SyntenyGrange, subject = tmpGrange))
  
  if (length(which(duplicated(testQuery))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated( testQuery))]
    testQuery <- testQuery[-which(duplicated( testQuery))]
  }
  ### error probably happens b/c there is no match for one of the queries
  
  mouseCoadRegionMat_adenoCar[ , i + 3] <- NA
  mouseCoadRegionMat_adenoCar[testQuery, i + 3] <- tmp$mean[tmpIdx]
}

mouseCoadRegionMat_adenoCar2 <- mouseCoadRegionMat_adenoCar[ , 4:ncol(mouseCoadRegionMat_adenoCar)]
cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
synMmCoadPerm <- foreach(i=1:25,
                         .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                           tmpAmpAll <- NULL
                           tmpDelAll <- NULL
                           tmpAmpFreq <- NULL
                           tmpDelFreq <- NULL
                           j <- 0
                           while (j < 400) {
                             tmp <- apply(mouseCoadRegionMat_adenoCar2, 2, function(x) sample(x, nrow(mouseCoadRegionMat_adenoCar2)))
                             ampFreq <- apply(tmp, 1, function(x) length(which(x > 0.2))/length(x))
                             delFreq <- apply(tmp, 1, function(x) length(which(x < -0.2))/length(x))
                             tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0.2)]) * length(which(x > 0.2))/length(x))
                             tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < -0.2)])) * length(which(x < -0.2))/length(x))
                             tmpAmpAll <- cbind(tmpAmpAll, tmpAmp)
                             tmpDelAll <- cbind(tmpDelAll, tmpDel)
                             tmpAmpFreq <- cbind(tmpAmpFreq, ampFreq)
                             tmpDelFreq <- cbind(tmpDelFreq, delFreq)
                             j <- j + 1
                           }
                           return(list(tmpAmpAll, tmpDelAll, tmpAmpFreq, tmpDelFreq))
                         }


stopCluster(cl)
print( Sys.time() - start )

mousePermAmpFreq_adenoCar <- do.call(cbind, synMmCoadPerm[[3]])
mousePermDelFreq_adenoCar <- do.call(cbind, synMmCoadPerm[[4]])



coad_adenoCar_arm_allSynTable$str <- paste(coad_adenoCar_arm_allSynTable$h_chr, coad_adenoCar_arm_allSynTable$h_start, coad_adenoCar_arm_allSynTable$h_end)
tmpSynTable_adenoCar <- coad_adenoCar_arm_allSynTable[-which(duplicated(coad_adenoCar_arm_allSynTable$str)), ]
tmpSynTable_adenoCar$h_chr <- as.numeric(str_remove(tmpSynTable_adenoCar$h_chr, "h_chr"))
tmpSynTable_adenoCar$m_chr <- as.numeric(str_remove(tmpSynTable_adenoCar$m_chr, "m_chr"))
tmpSynTable_adenoCar <- tmpSynTable_adenoCar[order(tmpSynTable_adenoCar$h_chr, tmpSynTable_adenoCar$h_start, decreasing = FALSE),]

ampSynTable_adenoCar <- data.frame("chr" = tmpSynTable_adenoCar$m_chr, "start" = tmpSynTable_adenoCar$m_start, "end" = tmpSynTable_adenoCar$m_end)
delSynTable_adenoCar <- data.frame("chr" = tmpSynTable_adenoCar$m_chr, "start" = tmpSynTable_adenoCar$m_start, "end" = tmpSynTable_adenoCar$m_end)
ampSynTable_adenoCar$h_freq <- apply(humanCoadRegionMat2, 1, function(x) length(which(x > 0.2))/length(x))
ampSynTable_adenoCar$m_freq <- apply(mouseCoadRegionMat_adenoCar2, 1, function(x) length(which(x > 0.2))/length(x))
delSynTable_adenoCar$h_freq <- apply(humanCoadRegionMat2, 1, function(x) length(which(x < -0.2))/length(x))
delSynTable_adenoCar$m_freq <- apply(mouseCoadRegionMat_adenoCar2, 1, function(x) length(which(x < -0.2))/length(x))


humanPermAmp2_adenoCar <- humanPermAmpCoadAdenoCarFreq
mousePermAmp2_adenoCar <- mousePermAmpFreq_adenoCar
humanPermDel2_adenoCar <- humanPermDelCoadAdenoCarFreq
mousePermDel2_adenoCar <- mousePermDelFreq_adenoCar
# mousePermAmpStat2_adenoCar <- mousePermAmpStat_adenoCar
# mousePermDelStat2_adenoCar <- mousePermDelStat_adenoCar
# humanPermAmpStatCoad2 <- humanPermAmpCoadAdenoCarStat
# humanPermDelStatCoad2 <- humanPermDelCoadAdenoCarStat


permStatResTbl_adenoCar <- NULL
for (i in 1:396) {

  
  ### alternatively I can look at them mouse and human separately and only nominate ones where they're both significant
  tmpAmpZ_hfreq <- (ampSynTable_adenoCar$h_freq[i] - mean(humanPermAmp2_adenoCar[i, ], na.rm = TRUE))/sd(humanPermAmp2_adenoCar[i, ], na.rm = TRUE)
  tmpDelZ_hfreq <- (delSynTable_adenoCar$h_freq[i] - mean(humanPermDel2_adenoCar[i, ], na.rm = TRUE))/sd(humanPermDel2_adenoCar[i, ], na.rm = TRUE)
  
  tmpAmpZ_mfreq <- (ampSynTable_adenoCar$m_freq[i] - mean(mousePermAmp2_adenoCar[i, ], na.rm = TRUE))/sd(mousePermAmp2_adenoCar[i, ], na.rm = TRUE)
  tmpDelZ_mfreq <- (delSynTable_adenoCar$m_freq[i] - mean(mousePermDel2_adenoCar[i, ], na.rm = TRUE))/sd(mousePermDel2_adenoCar[i, ], na.rm = TRUE)
  
  
  
  permStatResTbl_adenoCar <- rbind(permStatResTbl_adenoCar,
                                   data.frame("ampZ_hfreq" = tmpAmpZ_hfreq, "delZ_hfreq" = tmpDelZ_hfreq,
                                              "ampZ_mfreq" = tmpAmpZ_mfreq, "delZ_mfreq" = tmpDelZ_mfreq))
}

permStatResTbl_adenoCar$ampZ_hfreq2 <- pnorm(q=permStatResTbl_adenoCar$ampZ_hfreq, lower.tail=FALSE)
permStatResTbl_adenoCar$delZ_hfreq2 <- pnorm(q=permStatResTbl_adenoCar$delZ_hfreq, lower.tail=FALSE)
permStatResTbl_adenoCar$ampZ_mfreq2 <- pnorm(q=permStatResTbl_adenoCar$ampZ_mfreq, lower.tail=FALSE)
permStatResTbl_adenoCar$delZ_mfreq2 <- pnorm(q=permStatResTbl_adenoCar$delZ_mfreq, lower.tail=FALSE)



adenoCar_arm_allSynTable_amp <- tmpSynTable_adenoCar
adenoCar_arm_allSynTable_amp$h_freq <- ampSynTable_adenoCar$h_freq
adenoCar_arm_allSynTable_amp$m_freq <- ampSynTable_adenoCar$m_freq
adenoCar_arm_allSynTable_amp$hfreqP <- permStatResTbl_adenoCar$ampZ_hfreq2
adenoCar_arm_allSynTable_amp$hfreqAdjP <- p.adjust(adenoCar_arm_allSynTable_amp$hfreqP, method = "BH",
                                                   n = 4 * length(adenoCar_arm_allSynTable_amp$hfreqP))
adenoCar_arm_allSynTable_amp$mfreqP <- permStatResTbl_adenoCar$ampZ_mfreq2
adenoCar_arm_allSynTable_amp$mfreqAdjP <- p.adjust(adenoCar_arm_allSynTable_amp$mfreqP, method = "BH",
                                                   n = 4 * length(adenoCar_arm_allSynTable_amp$mfreqP))
adenoCar_arm_allSynTable_amp$sigCheck <- ifelse(adenoCar_arm_allSynTable_amp$hfreqAdjP < 0.05 & adenoCar_arm_allSynTable_amp$mfreqAdjP < 0.05,
                                                "good", "bad")


adenoCar_arm_allSynTable_del <- tmpSynTable_adenoCar
adenoCar_arm_allSynTable_del$h_freq <- delSynTable_adenoCar$h_freq
adenoCar_arm_allSynTable_del$m_freq <- delSynTable_adenoCar$m_freq
adenoCar_arm_allSynTable_del$hfreqP <- permStatResTbl_adenoCar$delZ_hfreq2
adenoCar_arm_allSynTable_del$hfreqAdjP <- p.adjust(adenoCar_arm_allSynTable_del$hfreqP, method = "BH",
                                                   n = 4 * length(adenoCar_arm_allSynTable_del$hfreqP))
adenoCar_arm_allSynTable_del$mfreqP <- permStatResTbl_adenoCar$delZ_mfreq2
adenoCar_arm_allSynTable_del$mfreqAdjP <- p.adjust(adenoCar_arm_allSynTable_del$mfreqP, method = "BH",
                                                   n = 4 * length(adenoCar_arm_allSynTable_del$mfreqP))
adenoCar_arm_allSynTable_del$sigCheck <- ifelse(adenoCar_arm_allSynTable_del$hfreqAdjP < 0.05 & adenoCar_arm_allSynTable_del$mfreqAdjP < 0.05,
                                                "good", "bad")


### works, now do the RNA portion

coadRsemExpression <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/gdac.broadinstitute.org_COADREAD.Merge_rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/COADREAD.rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt",
                                 sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

coadRsemExpression2 <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/gdac.broadinstitute.org_COADREAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/COADREAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt",
                                  sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

coadRsemExpression <- coadRsemExpression[, c(1, which(coadRsemExpression[1, ] == "raw_count"))]
coadRsemExpression2 <- coadRsemExpression2[, c(1, which(coadRsemExpression2[1, ] == "raw_count"))]
colnames(coadRsemExpression)[2:ncol(coadRsemExpression)] <- gsub("(^.*?-.{3}?)-.*", "\\1",
                                                                 colnames(coadRsemExpression)[2:ncol(coadRsemExpression)])
colnames(coadRsemExpression2)[2:ncol(coadRsemExpression2)] <- gsub("(^.*?-.{3}?)-.*", "\\1",
                                                                   colnames(coadRsemExpression2)[2:ncol(coadRsemExpression2)])

coadRsemExpression <- coadRsemExpression[2:nrow(coadRsemExpression), ]
coadRsemExpression2 <- coadRsemExpression2[2:nrow(coadRsemExpression2), ]

combinedRsemRaw <- cbind(coadRsemExpression, coadRsemExpression2[, 2:ncol(coadRsemExpression2)])
combinedRsemRaw <- combinedRsemRaw[, -which(duplicated(colnames(combinedRsemRaw)))]
combinedRsemRawMat <- combinedRsemRaw[, 2:ncol(combinedRsemRaw)]
combinedRsemRawMat <- apply(combinedRsemRawMat, 2, as.numeric)
rownames(combinedRsemRawMat) <- combinedRsemRaw$`Hybridization REF`

### creating a model matrix so I can easily pull out which samples have a given aneuploidy 
coadArmStatus <- broadBySampleGisticCoad
colnames(coadArmStatus)[2:ncol(coadArmStatus)] <- gsub("(^.*?-.{3}?)-.*", "\\1", colnames(coadArmStatus)[2:ncol(coadArmStatus)])
coadArmStatusMat <- coadArmStatus[, 2:ncol(coadArmStatus)]
rownames(coadArmStatusMat) <- coadArmStatus$`Chromosome Arm`

coadArmStatusFactor <- coadArmStatus
coadArmStatusFactor$`Chromosome Arm` <- coadArmStatus$`Chromosome Arm`
coadArmStatusFactor <- coadArmStatusFactor[which(coadArmStatusFactor$`Chromosome Arm` %in% c("1p", "7p", "7q", "8p","13q",
                                                                                             "14q", "15q","18p", "18q")), ]

coadArmStatusFactor[1, 2:ncol(coadArmStatusFactor)]  <- ifelse(coadArmStatusFactor[1, 2:ncol(coadArmStatusFactor)] < -0.2, 1, 0)
coadArmStatusFactor[2, 2:ncol(coadArmStatusFactor)]  <- ifelse(coadArmStatusFactor[2, 2:ncol(coadArmStatusFactor)] > 0.2, 1, 0)
coadArmStatusFactor[3, 2:ncol(coadArmStatusFactor)]  <- ifelse(coadArmStatusFactor[3, 2:ncol(coadArmStatusFactor)] > 0.2, 1, 0)
coadArmStatusFactor[4, 2:ncol(coadArmStatusFactor)]  <- ifelse(coadArmStatusFactor[4, 2:ncol(coadArmStatusFactor)] < -0.2, 1, 0)
coadArmStatusFactor[5, 2:ncol(coadArmStatusFactor)]  <- ifelse(coadArmStatusFactor[5, 2:ncol(coadArmStatusFactor)] > 0.2, 1, 0)
coadArmStatusFactor[6, 2:ncol(coadArmStatusFactor)]  <- ifelse(coadArmStatusFactor[6, 2:ncol(coadArmStatusFactor)] < -0.2, 1, 0)
coadArmStatusFactor[7, 2:ncol(coadArmStatusFactor)]  <- ifelse(coadArmStatusFactor[7, 2:ncol(coadArmStatusFactor)] < -0.2, 1, 0)
coadArmStatusFactor[8, 2:ncol(coadArmStatusFactor)]  <- ifelse(coadArmStatusFactor[8, 2:ncol(coadArmStatusFactor)] < -0.2, 1, 0)
coadArmStatusFactor[9, 2:ncol(coadArmStatusFactor)]  <- ifelse(coadArmStatusFactor[9, 2:ncol(coadArmStatusFactor)] < -0.2, 1, 0)

status1p <- factor(unlist(coadArmStatusFactor[1, 2:ncol(coadArmStatusFactor)]))
status7p <- factor(unlist(coadArmStatusFactor[2, 2:ncol(coadArmStatusFactor)]))
status7q <- factor(unlist(coadArmStatusFactor[3, 2:ncol(coadArmStatusFactor)]))
status8p <- factor(unlist(coadArmStatusFactor[4, 2:ncol(coadArmStatusFactor)]))
status13q <- factor(unlist(coadArmStatusFactor[5, 2:ncol(coadArmStatusFactor)]))
status14q <- factor(unlist(coadArmStatusFactor[6, 2:ncol(coadArmStatusFactor)]))
status15q <- factor(unlist(coadArmStatusFactor[7, 2:ncol(coadArmStatusFactor)]))
status18p <- factor(unlist(coadArmStatusFactor[8, 2:ncol(coadArmStatusFactor)]))
status18q <- factor(unlist(coadArmStatusFactor[9, 2:ncol(coadArmStatusFactor)]))


coadArmMemebershipMat <- model.matrix(~ 0 + status1p + status7p + status7q + status8p + status13q + status14q + status15q+ status18p + status18q)

coadArmMemebershipDf <- data.frame(coadArmMemebershipMat)
rownames(coadArmMemebershipDf) <- colnames(coadArmStatus)[2:ncol(coadArmStatus)]
# write.table(coadArmMemebershipDf, "/mnt/DATA5/tmp/kev/misc/20230926coadArmChangeMembership.txt",
#             row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")


### method is a tiny bit complex but I run two DE runs each time so I can then do
### side by side comparisons of pathway analyses
### Purpose is to get rid of common altered pathways for easier viewing
### only needed to do one at the end with no reciprocal - for graphing and analytical purposes

with_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status18q1 == 1)]

no_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status18q1 == 0)]

with_1_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status1p1 == 1 &
                                                    coadArmMemebershipDf$status18q1 == 1)]
no_1_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status1p1 == 0 &
                                                  coadArmMemebershipDf$status18q1 == 1)]

with_7_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7p1 == 1 &
                                                    coadArmMemebershipDf$status18q1 == 1)]
no_7_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7p1 == 0 &
                                                  coadArmMemebershipDf$status18q1 == 1)]

with_13_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status13q1 == 1 &
                                                     coadArmMemebershipDf$status18q1 == 1)]
no_13_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status13q1 == 0 &
                                                   coadArmMemebershipDf$status18q1 == 1)]

with_14_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status14q1 == 1 &
                                                     coadArmMemebershipDf$status18q1 == 1)]
no_14_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status14q1 == 0 &
                                                   coadArmMemebershipDf$status18q1 == 1)]

with_15_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status15q1 == 1 &
                                                     coadArmMemebershipDf$status18q1 == 1)]
no_15_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status15q1 == 0 &
                                                   coadArmMemebershipDf$status18q1 == 1)]

with_7_13_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7p1 == 1 &
                                                       coadArmMemebershipDf$status13q1 == 1 &
                                                       coadArmMemebershipDf$status18p1 == 1)]


with_7_18_no_13 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7p1 == 1 &
                                                          coadArmMemebershipDf$status13q1 == 0 &
                                                          coadArmMemebershipDf$status18p1 == 1)]


with_13_18_no_7 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7p1 == 0 &
                                                          coadArmMemebershipDf$status13q1 == 1 &
                                                          coadArmMemebershipDf$status18p1 == 1)]


with_18_no_7_13 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7p1 == 0 &
                                                          coadArmMemebershipDf$status13q1 == 0 &
                                                          coadArmMemebershipDf$status18p1 == 1)]


listOfSamples2 <- list(with_18, with_1_18, with_7_18, with_13_18, with_14_18, with_15_18, with_7_13_18, with_7_18_no_13, with_13_18_no_7, with_18_no_7_13, no_18)
listOfRecip2 <- list(no_18, no_1_18, no_7_18, no_13_18, no_14_18, no_15_18, with_18_no_7_13, with_18_no_7_13, with_18_no_7_13, with_18, with_18)
listOfSamplesNames2 <- c("with_18", "with_1_18", "with_7_18", "with_13_18", "with_14_18", "with_15_18",
                         "with_7_13_18", "with_7_18_no_13", "with_13_18_no_7", "with_18_no_7_13", "no_18")

normalSamples <- colnames(combinedRsemRawMat)[grep("11A", colnames(combinedRsemRawMat))]

library(foreach)
library(doParallel)
cl <- makeCluster(11)
registerDoParallel(cl)
start <- Sys.time()
diffArmLossCombo2 <- foreach(i=seq_along(listOfSamples2),.combine = 'comb', .multicombine = TRUE,
                             .init = list(list(), list()), .packages = c("edgeR", "stringr")) %dopar% {
                               
                               tmpCoadArmNames <- unlist(listOfSamples2[i])
                               tmpRawCounts <- combinedRsemRawMat[, which(colnames(combinedRsemRawMat) %in% c(normalSamples, tmpCoadArmNames))]
                               tmpGroup <- rep(1, ncol(tmpRawCounts))
                               tmpGroup[which(colnames(tmpRawCounts) %in% tmpCoadArmNames)] <- 2
                               dgeTmp <- DGEList(counts=tmpRawCounts, group= tmpGroup)
                               keep <- filterByExpr(dgeTmp)
                               dgeTmp <- dgeTmp[keep,,keep.lib.sizes=FALSE]
                               dgeTmp <- calcNormFactors(dgeTmp)
                               designTmp <- model.matrix(~tmpGroup)
                               dgeTmp <- estimateDisp(dgeTmp, designTmp)
                               fitTmp <- glmFit(dgeTmp, designTmp)
                               lrtTmp <- glmLRT(fitTmp,coef="tmpGroup")
                               lrtTabTmp <- data.frame("genes" = rownames(lrtTmp), lrtTmp$table)
                               lrtTabTmp$qval <- p.adjust(lrtTabTmp$PValue, method = "BH", n = length(lrtTabTmp$PValue))
                               lrtTabTmp$arm <- listOfSamplesNames2[i]
                               
                               tmpCoadArmNames2 <- unlist(listOfRecip2[i])
                               tmpRawCounts2 <- combinedRsemRawMat[, which(colnames(combinedRsemRawMat) %in% c(normalSamples, tmpCoadArmNames2))]
                               tmpGroup2 <- rep(1, ncol(tmpRawCounts2))
                               tmpGroup2[which(colnames(tmpRawCounts2) %in% tmpCoadArmNames2)] <- 2
                               dgeTmp2 <- DGEList(counts=tmpRawCounts2, group= tmpGroup2)
                               keep2 <- filterByExpr(dgeTmp2)
                               dgeTmp2 <- dgeTmp2[keep2,,keep.lib.sizes=FALSE]
                               dgeTmp2 <- calcNormFactors(dgeTmp2)
                               designTmp2 <- model.matrix(~tmpGroup2)
                               dgeTmp2 <- estimateDisp(dgeTmp2, designTmp2)
                               fitTmp2 <- glmFit(dgeTmp2, designTmp2)
                               lrtTmp2 <- glmLRT(fitTmp2,coef="tmpGroup2")
                               lrtTabTmp2 <- data.frame("genes" = rownames(lrtTmp2), lrtTmp2$table)
                               lrtTabTmp2$qval <- p.adjust(lrtTabTmp2$PValue, method = "BH", n = length(lrtTabTmp2$PValue))
                               lrtTabTmp2$arm <- listOfSamplesNames2[i]
                               
                               return(list(lrtTabTmp, lrtTabTmp2))
                             }


stopCluster(cl)
print( Sys.time() - start )


diffArmLossComboDf2 <- do.call(rbind, diffArmLossCombo2[[1]])
diffArmLossComboRecipDf2 <- do.call(rbind, diffArmLossCombo2[[2]])

diffArmLossComboDf2$genes <- str_remove(diffArmLossComboDf2$genes, "\\|.*")
diffArmLossComboRecipDf2$genes <- str_remove(diffArmLossComboRecipDf2$genes, "\\|.*")


# write.table(diffArmLossComboDf2, "/mnt/DATA5/tmp/kev/misc/20231019diffArmLossComboDf2.txt", 
#             sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# 
# write.table(diffArmLossComboRecipDf2, "/mnt/DATA5/tmp/kev/misc/20231019diffArmLossRecipComboDf2.txt", 
#             sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


with_7 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7q1 == 1)]
no_7 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7q1 == 0)]
with_13 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status13q1 == 1)]
no_13 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status13q1 == 0)]
with_14 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status14q1 == 1)]
no_14 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status14q1 == 0)]
with_15 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status15q1 == 1)]
no_15 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status15q1 == 0)]
with_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status18q1 == 1)]
no_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status18q1 == 0)]
with_1 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status1p1 == 1)]
no_1 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status1p1 == 0)]


listOfSamples <- list(with_1, with_7, with_13, with_14, with_15, with_18)
listOfRecip <- list(no_1, no_7, no_13, no_14, no_15, no_18)
listOfSamplesNames <- list("with_1", "with_7", "with_13", "with_14", "with_15", "with_18")

library(foreach)
library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)
start <- Sys.time()
diffArmLossCombo <- foreach(i=seq_along(listOfSamples),.combine = 'comb', .multicombine = TRUE,
                             .init = list(list(), list()), .packages = c("edgeR", "stringr")) %dopar% {
                               
                               tmpCoadArmNames <- unlist(listOfSamples[i])
                               tmpRawCounts <- combinedRsemRawMat[, which(colnames(combinedRsemRawMat) %in% c(normalSamples, tmpCoadArmNames))]
                               tmpGroup <- rep(1, ncol(tmpRawCounts))
                               tmpGroup[which(colnames(tmpRawCounts) %in% tmpCoadArmNames)] <- 2
                               dgeTmp <- DGEList(counts=tmpRawCounts, group= tmpGroup)
                               keep <- filterByExpr(dgeTmp)
                               dgeTmp <- dgeTmp[keep,,keep.lib.sizes=FALSE]
                               dgeTmp <- calcNormFactors(dgeTmp)
                               designTmp <- model.matrix(~tmpGroup)
                               dgeTmp <- estimateDisp(dgeTmp, designTmp)
                               fitTmp <- glmFit(dgeTmp, designTmp)
                               lrtTmp <- glmLRT(fitTmp,coef="tmpGroup")
                               lrtTabTmp <- data.frame("genes" = rownames(lrtTmp), lrtTmp$table)
                               lrtTabTmp$qval <- p.adjust(lrtTabTmp$PValue, method = "BH", n = length(lrtTabTmp$PValue))
                               lrtTabTmp$arm <- listOfSamplesNames[i]
                               
                               tmpCoadArmNames2 <- unlist(listOfRecip[i])
                               tmpRawCounts2 <- combinedRsemRawMat[, which(colnames(combinedRsemRawMat) %in% c(normalSamples, tmpCoadArmNames2))]
                               tmpGroup2 <- rep(1, ncol(tmpRawCounts2))
                               tmpGroup2[which(colnames(tmpRawCounts2) %in% tmpCoadArmNames2)] <- 2
                               dgeTmp2 <- DGEList(counts=tmpRawCounts2, group= tmpGroup2)
                               keep2 <- filterByExpr(dgeTmp2)
                               dgeTmp2 <- dgeTmp2[keep2,,keep.lib.sizes=FALSE]
                               dgeTmp2 <- calcNormFactors(dgeTmp2)
                               designTmp2 <- model.matrix(~tmpGroup2)
                               dgeTmp2 <- estimateDisp(dgeTmp2, designTmp2)
                               fitTmp2 <- glmFit(dgeTmp2, designTmp2)
                               lrtTmp2 <- glmLRT(fitTmp2,coef="tmpGroup2")
                               lrtTabTmp2 <- data.frame("genes" = rownames(lrtTmp2), lrtTmp2$table)
                               lrtTabTmp2$qval <- p.adjust(lrtTabTmp2$PValue, method = "BH", n = length(lrtTabTmp2$PValue))
                               lrtTabTmp2$arm <- listOfSamplesNames[i]
                               
                               return(list(lrtTabTmp, lrtTabTmp2))
                             }


stopCluster(cl)
print( Sys.time() - start )


diffArmLossComboDf <- do.call(rbind, diffArmLossCombo[[1]])
diffArmLossComboRecipDf <- do.call(rbind, diffArmLossCombo[[2]])

diffArmLossComboDf$genes <- str_remove(diffArmLossComboDf$genes, "\\|.*")
diffArmLossComboRecipDf$genes <- str_remove(diffArmLossComboRecipDf$genes, "\\|.*")



### the pathway analysis was done in brahm since avatar's 2023102 R version was too old
###
###



library(clusterProfiler)
library(msigdbr)
packageVersion("clusterProfiler")
library(gridExtra)
library(ggplot2)
library(ggpubr)


diffArmLossComboDf <- read.table("/avatar_data5/tmp/kev/misc/20231019diffArmLossComboDf2.txt", sep = "\t",
                                 header = TRUE, stringsAsFactors = FALSE)
diffArmLossComboRecipDf <- read.table("/avatar_data5/tmp/kev/misc/20231019diffArmLossRecipComboDf2.txt", sep = "\t",
                                      header = TRUE, stringsAsFactors = FALSE)
hg19GeneLocations <- read.table("/avatar_data5/tmp/kev/misc/20231018hg19GeneLocations.txt", sep = "\t",
                                header = TRUE, stringsAsFactors = FALSE)



msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
msig_hallmark<- msig_hallmark[, c("gs_name", "gene_symbol")]



### gives me to evaluate pathways that are present and absence of pathways 
for (i in unique(diffArmLossComboDf$arm)) {
  print(i)
  
  tmpDiffDf <- diffArmLossComboDf[which(diffArmLossComboDf$arm == i), ]
  tmpDiffRecipDf <- diffArmLossComboRecipDf[which(diffArmLossComboRecipDf$arm == i), ]
  
  gene_list <- -log10(tmpDiffDf$PValue) * sign(tmpDiffDf$logFC)
  names(gene_list) <- tmpDiffDf$genes
  gene_list = sort(gene_list, decreasing = TRUE)
  
  gene_list2 <- -log10(tmpDiffRecipDf$PValue) * sign(tmpDiffRecipDf$logFC)
  names(gene_list2) <- tmpDiffRecipDf$genes
  gene_list2 = sort(gene_list2, decreasing = TRUE)
  
  if (length(which(is.infinite(gene_list))) > 0) {
    gene_list <- gene_list[-which(is.infinite(gene_list))]
  } else if (length(which(is.infinite(gene_list2))) > 0) {
    gene_list2 <- gene_list2[-which(is.infinite(gene_list2))]
  }
  
  with1 <- GSEA(gene_list, TERM2GENE = msig_hallmark, pvalueCutoff = 0.05, eps = 0)
  with1@result$IDsign <- paste(with1@result$ID, sign(with1@result$NES))
  without1 <- GSEA(gene_list2, TERM2GENE = msig_hallmark, pvalueCutoff = 0.05, eps = 0)
  without1@result$IDsign <- paste(without1@result$ID, sign(without1@result$NES))
  shared1 <- with1@result$ID[which(with1@result$IDsign %in% without1@result$IDsign)]
  
  
  with1@result <- with1@result[-which(with1@result$ID %in% shared1), ]
  without1@result <- without1@result[-which(without1@result$ID %in% shared1), ]
  
  if (nrow(with1@result) == 0) {
    with1@result <- rbind(with1@result, without1@result[1, ])
    with1@result$setSize <- 0
    with1@result$enrichmentScore <- 0
    with1@result$NES <- 0
    with1@result$pvalue <- 1
    with1@result$p.adjust <- 1
    with1@result$qvalue <- 1
  }
  
  if (nrow(without1@result) == 0) {
    without1@result <- rbind(without1@result, with1@result[1, ])
    without1@result$setSize <- 0
    without1@result$enrichmentScore <- 0
    without1@result$NES <- 0
    without1@result$pvalue <- 1
    without1@result$p.adjust <- 1
    without1@result$qvalue <- 1
  }
  
  
  
  #
  a <- tryCatch(dotplot(with1, split=".sign", font.size = 8) + facet_grid(.~.sign) + ggtitle("with arm change hallmark"),
                error = function(x) return(NULL));
  c <-tryCatch(dotplot(without1, split=".sign", font.size = 8) + facet_grid(.~.sign) + ggtitle("without arm change"),
               error = function(x) return(NULL));
  
  
  
  pdf(file = paste0("/avatar_data5/tmp/kev/misc/20231019_", i, ".pdf"),
      useDingbats = TRUE)
  print(tryCatch(grid.arrange(c, a, nrow = 2, top = text_grob(i)), error = function(x) return(NULL)))
  dev.off()
}


### will do figure based on order
###
###

### freqplots
library(gridExtra)

e <- freqPlotv2(armGisticCoadApc_amp_bed, armGisticCoadApc_del_bed,
                main = "coadread n = 349", speciesType = "human")
g <- freqPlotv2(allMouseAneu_adenoCar_amp_bed, allMouseAneu_adenoCar_del_bed,
                main = "mouse coadread n = 27", speciesType = "mouse")
f <- freqPlotv2(allMouseAneu_adeno_amp_bed, allMouseAneu_adeno_del_bed,
                main = "mouse adenoma n = 11", speciesType = "mouse")
grid.arrange(e, f, g, nrow = 3)


### creating a different heatmap for each type i.e met, adenoma and adenocarcinoma
library(pheatmap)
heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-1,1,2/1000)

allMouseAneuploidy2 <- allMouseAneuploidy
rownames(allMouseAneuploidy2) <- str_remove(rownames(allMouseAneuploidy2), "x.*")
annoCol2 <- list("Genotype" = c("AKP" = "darkred", "AKPS" = "darkorange",
                                "CB" = "purple", "AK" = "darkgreen"),
                 "Type" = c("met" = "black", "Adenoma" = "yellow",
                            "Adenocarcinoma" = "red", "hyperplastic polyp" ="darkgreen"))


allMouseAneuploidyCoadAdeno <- allMouseAneuploidy2[which(rownames(allMouseAneuploidy2) %in% annoTableCoad2$V6[which(annoTableCoad2$histology %in% c("Adenoma"))]),]
allMouseAneuploidyCoadAdeno[which(allMouseAneuploidyCoadAdeno == "none")] <- 0
allMouseAneuploidyCoadAdeno[which(allMouseAneuploidyCoadAdeno == "gain")] <- 1
allMouseAneuploidyCoadAdeno[which(allMouseAneuploidyCoadAdeno == "loss")] <- -1
allMouseAneuploidyCoadAdenoMat <- apply(allMouseAneuploidyCoadAdeno, 2, as.numeric)
rownames(allMouseAneuploidyCoadAdenoMat ) <- rownames(allMouseAneuploidyCoadAdeno)
allMouseAneuploidyCoadAdenoMat2 <- allMouseAneuploidyCoadAdenoMat[which(rownames(allMouseAneuploidyCoadAdenoMat) %in% goodTcCoadAdeno2), ]

heatmapAnnoDfCoad <- annoTableCoad2[, c("histology", "geno")]
rownames(heatmapAnnoDfCoad) <- annoTableCoad2$V6
colnames(heatmapAnnoDfCoad) <- c("Type", "Genotype")

# pdf(paste0("/mnt/DATA5/tmp/kev/misc/20231003mm10CoadAdeno", ".pdf"), useDingbats = FALSE, width = 10)
pheatmap(allMouseAneuploidyCoadAdenoMat2, cluster_rows = TRUE, cluster_cols = FALSE, cellwidth = 20,
         cellheight = 10, fontsize_row = 5, color = heatMapCol, breaks = colors.breaks,
         annotation_row = heatmapAnnoDfCoad, annotation_colors = annoCol2,
         main = "CDX2P-CreERT2 Apcfl/+, KrasLSLG12D/+")
# dev.off()

allMouseAneuploidyCoadAdenoCar <- allMouseAneuploidy2[which(rownames(allMouseAneuploidy2) %in% annoTableCoad2$V6[which(annoTableCoad2$histology %in% c("Adenocarcinoma"))]),]
allMouseAneuploidyCoadAdenoCar[which(allMouseAneuploidyCoadAdenoCar == "none")] <- 0
allMouseAneuploidyCoadAdenoCar[which(allMouseAneuploidyCoadAdenoCar == "gain")] <- 1
allMouseAneuploidyCoadAdenoCar[which(allMouseAneuploidyCoadAdenoCar == "loss")] <- -1
allMouseAneuploidyCoadAdenoCarMat <- apply(allMouseAneuploidyCoadAdenoCar, 2, as.numeric)
rownames(allMouseAneuploidyCoadAdenoCarMat ) <- rownames(allMouseAneuploidyCoadAdenoCar)
allMouseAneuploidyCoadAdenoCarMat2 <- allMouseAneuploidyCoadAdenoCarMat[which(rownames(allMouseAneuploidyCoadAdenoCarMat) %in% goodTcCoadAdenoCar2), ]

# pdf(paste0("/mnt/DATA5/tmp/kev/misc/20231003mm10CoadAdenoCar", ".pdf"), useDingbats = FALSE, width = 10)
pheatmap(allMouseAneuploidyCoadAdenoCarMat2, cluster_rows = TRUE, cluster_cols = FALSE, cellwidth = 20,
         cellheight = 10, fontsize_row = 5, color = heatMapCol, breaks = colors.breaks,
         annotation_row = heatmapAnnoDfCoad, annotation_colors = annoCol2,
         main = "CDX2P-CreERT2 Apcfl/+, KrasLSLG12D/+, p53 ex2-10fl/fl, Sox9 fl/+")
# dev.off()


allMouseAneuploidyCoadMet <- allMouseAneuploidy2[which(rownames(allMouseAneuploidy2) %in% annoTableCoad2$V6[which(annoTableCoad2$histology %in% c("met"))]),]
allMouseAneuploidyCoadMet[which(allMouseAneuploidyCoadMet == "none")] <- 0
allMouseAneuploidyCoadMet[which(allMouseAneuploidyCoadMet == "gain")] <- 1
allMouseAneuploidyCoadMet[which(allMouseAneuploidyCoadMet == "loss")] <- -1
allMouseAneuploidyCoadMetMat <- apply(allMouseAneuploidyCoadMet, 2, as.numeric)
rownames(allMouseAneuploidyCoadMetMat ) <- rownames(allMouseAneuploidyCoadMet)
# allMouseAneuploidyCoadMetMat2 <- allMouseAneuploidyCoadMetMat[which(rownames(allMouseAneuploidyCoadMetMat) %in% c(goodTcCoadAdeno2, goodTcCoadAdenoCar2)), ]

allMouseAneuploidyCoadMetMat2 <- allMouseAneuploidyCoadMetMat

# pdf(paste0("/mnt/DATA5/tmp/kev/misc/20231003mm10CoadMet", ".pdf"), useDingbats = FALSE, width = 10)
pheatmap(allMouseAneuploidyCoadMetMat2, cluster_rows = TRUE, cluster_cols = FALSE, cellwidth = 20,
         cellheight = 10, fontsize_row = 5, color = heatMapCol, breaks = colors.breaks,
         annotation_row = heatmapAnnoDfCoad[which(rownames(heatmapAnnoDfCoad) %in% rownames(allMouseAneuploidyCoadMetMat2)), ],
         annotation_colors = annoCol2,
         main = "CDX2P-CreERT2 Apcfl/+, KrasLSLG12D/+, p53 ex2-10fl/fl")
# dev.off()


### fraction of aneuploidy altered

### need to do parallel loop to calculate all fga for graph

dir <- "/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/"
setwd(dir)
allDirs <- list.dirs()
allDirs <- allDirs[grep("TP", allDirs)]

listOfFiles <- allDirs

allCans <- c("BLCA", "BRCA", "CESC", "COADREAD", "ESCA", "GBM",
             "GBMLGG", "HNSC", "KIPAN", "KIRC", "LGG", "LIHC",
             "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD",
             "SARC", "STES", "TGCT", "THYM")


allCancerFga <- NULL
for (j in seq_along(listOfFiles)) {
  setwd(dir)
  setwd(listOfFiles[j])
  
  tmpBroad <- read.table("broad_values_by_arm.txt", sep = "\t", stringsAsFactors = FALSE,
                         header = TRUE, check.names = FALSE)
  colnames(tmpBroad)[2:ncol(tmpBroad)] <- gsub("(^.*?-.{3}?)-.*", "\\1",
                                               colnames(tmpBroad)[2:ncol(tmpBroad)])
  if (allCans[j] == "COADREAD") {
    tmpBroad <- tmpBroad[, c(1, grep(paste0(coad_anno2$Sample.ID, collapse = "|"), colnames(tmpBroad)))] 
  }
  tmpBroad <- tmpBroad[-grep("X", tmpBroad$`Chromosome Arm`), ]
  tmpBroadMat <- tmpBroad[, 2:ncol(tmpBroad)]
  tmpBroad2 <- cbind(hg19ArmLocations[, c("chrStripped", "start", "end")], tmpBroadMat)
  tmpBroadMelt <- melt(tmpBroad2, id.vars = c("chrStripped", "start", "end"))
  tmpBroadMelt2 <- data.frame("sampleID" = tmpBroadMelt$variable, "chrom" = tmpBroadMelt$chrStripped,
                              "start.pos" = tmpBroadMelt$start, "end.pos" = tmpBroadMelt$end, 
                              "n.probes" = NA, "mean" = tmpBroadMelt$value)
  tmpBroadMelt2$str <- paste0(tmpBroadMelt2$sampleID, tmpBroadMelt2$chrom, 
                              tmpBroadMelt2$start.pos, tmpBroadMelt2$end.pos)
  tmpBroadMelt2$length <- tmpBroadMelt2$end.pos - tmpBroadMelt2$start.pos
  tmpRes <- fgaCalculator_tcga_aneuploidy(tmpBroadMelt2)
  tmpRes$Cancer <- allCans[j]
  allCancerFga <- rbind(allCancerFga, tmpRes)
}

adenoFga <- fgaCalculator_amp_aneuploidy(allMouseAneu_adeno)
adenoCarFga <- fgaCalculator_amp_aneuploidy(allMouseAneu_adenoCar)
adenoFga$Cancer <- "mm_adeno"
adenoCarFga$Cancer <- "mm_coadread"

allFga <- rbind(allCancerFga, adenoFga, adenoCarFga)
allFga$Cancer <- reorder(allFga$Cancer, allFga$fga, median)
allFga$Cancer <- factor(allFga$Cancer, rev(levels(allFga$Cancer)))
fgaColorVector <- rep("#000000", length(unique(allFga$Cancer)))
fgaColorVector[which(levels(allFga$Cancer) %in% c("COADREAD", "mm_coadread", "mm_adeno" ))] <- "#8B0000"

# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231005fgaComp.pdf", useDingbats = TRUE,
#     width = 7, height = 4)
ggplot(allFga, aes(x = Cancer, y = fga)) + geom_boxplot(color = fgaColorVector) + 
  ylab("fraction of genome altered (aneuploidy)") + xlab("cancer type/model") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()

### doing for just mouse adeno +  adenocarcinoma + hg coadread
allFgaRed <- allFga[which(allFga$Cancer %in% c("COADREAD", "mm_coadread", "mm_adeno")), ]

library(ggpubr)
my_comparisons <- list( c("COADREAD", "mm_coadread"), c("COADREAD", "mm_adeno"), c("mm_coadread", "mm_adeno") )
ggboxplot(allFgaRed, x = "Cancer", y = "fga") + 
  ylab("fraction of genome altered (aneuploidy)") + xlab("cancer type/model") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+ # Add pairwise comparisons p-value
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### for methodlogy figure
### 20231024: creating dummy frequency plot for methodology figure

dummyFreqAmp <- armGisticCoadApc_amp_bed
dummyFreqAmp$Freq[-which(dummyFreqAmp$Chr == "13")] <- 0
dummyFreqDel <- armGisticCoadApc_del_bed
dummyFreqDel$Freq[-which(dummyFreqDel$Chr == "18")] <- 0

# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231024exampleFreq.pdf", useDingbats = TRUE,
#     width = 10, height = 4)
freqPlotv2(dummyFreqAmp, dummyFreqDel,
           main = "example", speciesType = "human")
# dev.off()

### 20231024: single chromosome for method
circosFreqSingleChr(allMouseAneu_adenoCar_amp_bed, allMouseAneu_adenoCar_del_bed, armGisticCoadApc_amp_bed,
                    armGisticCoadApc_del_bed, filename = "20231024coadreadSingleChr13", ref = "human",
                    chromosome = "h_chr13")


### 20231024: human and mouse freq histogram for method

human_freq_10000 <- rnorm(n = 10000, mean = 0.5, sd = 0.1)
mouse_freq_10000 <- rnorm(n = 10000, mean = 0.3, sd = 0.1)

human_freq_10000 <- c(human_freq_10000, 0, 1)
mouse_freq_10000 <- c(mouse_freq_10000, 0, 1)


# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231024exampleHumanPerm.pdf", useDingbats = TRUE,
#     width = 10, height = 4)
ggplot() + geom_histogram(aes(human_freq_10000), bins = 100, color="darkgreen", fill="white") +
  geom_vline(xintercept = c(0.6), linetype = "dashed") + 
  geom_vline(xintercept = c(0.75), linetype = "dashed", color = "darkred") + scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(breaks = seq(0, 500, 100)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")) + 
  xlab("Human Frequency 10,000 Permutations") + ylab("Frequency Count")
# dev.off()

# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231024exampleMousePerm.pdf", useDingbats = TRUE,
#     width = 10, height = 4)
ggplot() + geom_histogram(aes(mouse_freq_10000), bins = 100, color="darkblue", fill="white") +
  geom_vline(xintercept = c(0.7), linetype = "dashed") +
  geom_vline(xintercept = c(0.6), linetype = "dashed", color = "darkred") +scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(breaks = seq(0, 500, 100)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")) + 
  xlab("Mouse Freqeuncy 10,000 Permutations") + ylab("Frequency Count")
# dev.off()


### barplot charts for synteny - 20231026 need to redo for first plot
### also add some quick analysis for paper stats i.e number of genes, segments 
### before and after method 

### might be able to replace this input with permutation resutls since those data already have freqs


syntenyDelIdx <- which(adenoCar_arm_allSynTable_del$m_freq >  0.1 & adenoCar_arm_allSynTable_del$h_freq > 0.1)
syntenyAmpIdx <- which(adenoCar_arm_allSynTable_amp$m_freq >  0.1 & adenoCar_arm_allSynTable_amp$h_freq > 0.1)

coad_adenoCar_arm_allSynTableMouseRef <- adenoCar_arm_allSynTable_del
coad_adenoCar_arm_allSynTableMouseRef$str <- paste( paste0("m_chr", coad_adenoCar_arm_allSynTableMouseRef$m_chr),
                                                   paste0("h_chr", coad_adenoCar_arm_allSynTableMouseRef$h_chr))

barplotFractionMouse <- NULL
i <- unique(coad_adenoCar_arm_allSynTableMouseRef$m_chr)[1]
for (i in unique(coad_adenoCar_arm_allSynTableMouseRef$m_chr)) {
  tmpDf <- coad_adenoCar_arm_allSynTableMouseRef[which(coad_adenoCar_arm_allSynTableMouseRef$m_chr == i), ]
  tmpDf$length <- tmpDf$m_end - tmpDf$m_start
  for (j in unique(tmpDf$h_chr)) {
    tmpRes <- sum(tmpDf$length[which(tmpDf$h_chr == j)])/sum(tmpDf$length)
    barplotFractionMouse <- rbind(barplotFractionMouse, data.frame("h_chr" = j, "m_chr" = i,"fraction" = tmpRes,
                                                                   "numberRegions" = length(which(tmpDf$h_chr == j))))
  }
}

barplotMChrAmp <- unique(coad_adenoCar_arm_allSynTableMouseRef$str[syntenyAmpIdx])
barplotMChrDel <- unique(coad_adenoCar_arm_allSynTableMouseRef$str[syntenyDelIdx])


### 15 in both which will be purple

barplotFractionFiltMouseCoad <- barplotFractionMouse
barplotFractionFiltMouseCoad$h_chr <- paste0("h_chr",barplotFractionFiltMouseCoad$h_chr)
barplotFractionFiltMouseCoad$m_chr <- paste0("m_chr",barplotFractionFiltMouseCoad$m_chr)
barplotFractionFiltMouseCoad$str <- paste(barplotFractionFiltMouseCoad$m_chr,
                                          barplotFractionFiltMouseCoad$h_chr)
barplotFractionFiltMouseCoad$color <- "#FFFFFF"
barplotFractionFiltMouseCoad$color[which(barplotFractionFiltMouseCoad$str %in% barplotMChrAmp)] <- "#8B0000"
barplotFractionFiltMouseCoad$color[which(barplotFractionFiltMouseCoad$str %in% barplotMChrDel)] <- "#00008B"
barplotFractionFiltMouseCoad$color[which(barplotFractionFiltMouseCoad$str %in% intersect(barplotMChrAmp, barplotMChrDel))] <- "#800080"
barplotFractionFiltMouseCoad$alpha <- ifelse(barplotFractionFiltMouseCoad$color == "#FFFFFF", 0.1, 1)

barplotFractionFiltMouseCoad$h_chr <- factor(barplotFractionFiltMouseCoad$h_chr, levels = paste0("h_chr", 1:22))
barplotFractionFiltMouseCoad$m_chr <- factor(barplotFractionFiltMouseCoad$m_chr, levels = paste0("m_chr", 1:19))


# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231027fractionSyntenyReccurentSign.pdf", useDingbats = TRUE,
#     width = 10, height = 5)
ggplot(barplotFractionFiltMouseCoad, aes(fill= h_chr, y= fraction, x = m_chr)) + 
  geom_bar(position="stack", stat="identity", color = barplotFractionFiltMouseCoad$color,
           alpha =  barplotFractionFiltMouseCoad$alpha) + 
  scale_fill_manual(values = colorVector) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + 
  ylab("human chromosome synteny fraction") + xlab("mouse chromsome")
# dev.off()



### below is the filtered version after method

barplotFractionFiltMouseCoad2 <- barplotFractionFiltMouseCoad
barplotFractionFiltMouseCoad2$str <- paste(barplotFractionFiltMouseCoad2$h_chr, barplotFractionFiltMouseCoad2$m_chr)
barplotFractionFiltMouseCoad2$alpha <- 0.1
barplotFractionFiltMouseCoad2$alpha[which(barplotFractionFiltMouseCoad2$str == "h_chr13 m_chr5" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr7 m_chr5" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr7 m_chr12" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr15 m_chr9" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr15 m_chr7" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr14 m_chr14" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr8 m_chr14" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr1 m_chr4" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr18 m_chr18")] <- 1
barplotFractionFiltMouseCoad2$color <- "#FFFFFF"
barplotFractionFiltMouseCoad2$color[which(barplotFractionFiltMouseCoad2$str == "h_chr13 m_chr5" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr7 m_chr5" | 
                                            barplotFractionFiltMouseCoad2$str == "h_chr7 m_chr12")] <- "#8B0000"

barplotFractionFiltMouseCoad2$color[which(barplotFractionFiltMouseCoad2$str == "h_chr18 m_chr18" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr15 m_chr9" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr15 m_chr7" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr14 m_chr14" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr8 m_chr14" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr1 m_chr4")] <- "#00008B"



# pdf(file = "/mnt/DATA5/tmp/kev/misc/20231011fractionSyntenyConMouseRefPerm.pdf", useDingbats = TRUE,
#     width = 10, height = 5)
ggplot(barplotFractionFiltMouseCoad2, aes(fill= h_chr, y= fraction, x = m_chr)) + 
  geom_bar(position="stack", stat="identity", color = barplotFractionFiltMouseCoad2$color,
           alpha =  barplotFractionFiltMouseCoad2$alpha) + 
  scale_fill_manual(values = colorVector) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + 
  ylab("human chromosome synteny fraction") + xlab("mouse chromsome")
# dev.off()


### the above has the right charts, going slight backwards to get the number of genes
### getting 
tmpAmpGrange <- GRanges(seqnames = coad_adenoCar_arm_allSynTableMouseRef$h_chr[syntenyAmpIdx],
                        IRanges(start = coad_adenoCar_arm_allSynTableMouseRef$h_start[syntenyAmpIdx],
                                end = coad_adenoCar_arm_allSynTableMouseRef$h_end[syntenyAmpIdx]))

tmpDelGrange <- GRanges(seqnames = coad_adenoCar_arm_allSynTableMouseRef$h_chr[syntenyDelIdx],
                        IRanges(start = coad_adenoCar_arm_allSynTableMouseRef$h_start[syntenyDelIdx],
                                end = coad_adenoCar_arm_allSynTableMouseRef$h_end[syntenyDelIdx]))


hg19biomartTable <- read.table("/home/kevhu/data/20230720hg19KnownCanBioMartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")


hg19GeneLocations <- NULL
for (i in unique(hg19biomartTable$external_gene_name)) {
  tmpDf <- hg19biomartTable[which(hg19biomartTable$external_gene_name == i), ]
  hg19GeneLocations <- rbind(hg19GeneLocations, data.frame("gene" = i, "chr"  = tmpDf$chromosome_name[1],
                                                           "start" = min(tmpDf$exon_chrom_start),                                                          "end" = max(tmpDf$exon_chrom_end)))
}

hg19ExonGr <- GRanges(seqnames = hg19GeneLocations$chr, 
                      IRanges(start = hg19GeneLocations$start, hg19GeneLocations$end))

length(subjectHits(findOverlaps(hg19ExonGr, tmpAmpGrange)))
length(subjectHits(findOverlaps(hg19ExonGr, tmpDelGrange)))

### specific gene names for nominated regions
###
###

adenoCar_arm_allSynTable_amp2 <- adenoCar_arm_allSynTable_amp[which(adenoCar_arm_allSynTable_amp$sigCheck == "good"), ]
adenoCar_arm_allSynTable_amp2$str <- paste(adenoCar_arm_allSynTable_amp2$h_chr, adenoCar_arm_allSynTable_amp2$m_chr)
i <- unique(adenoCar_arm_allSynTable_amp2$h_chr)
affectedGenesGains <- NULL
for (i in unique(adenoCar_arm_allSynTable_amp2$h_chr)) {
  tmp <- adenoCar_arm_allSynTable_amp2[which(adenoCar_arm_allSynTable_amp2$h_chr == i),]
  tmpGr <- GRanges(seqnames = tmp$h_chr, IRanges(start = tmp$h_start, end = tmp$h_end))
  tmpGene <- hg19GeneLocations$gene[queryHits(findOverlaps(hg19ExonGr, tmpGr))]
  affectedGenesGains <- rbind(affectedGenesGains, data.frame("h_chr" = tmp$h_chr[subjectHits(findOverlaps(hg19ExonGr, tmpGr))],
                                                             "mchr" = tmp$m_chr[subjectHits(findOverlaps(hg19ExonGr, tmpGr))],
                                                             "gene" = tmpGene,
                                                             "type" = rep("gains", length(tmpGene))))
}

adenoCar_arm_allSynTable_del2 <- adenoCar_arm_allSynTable_del[which(adenoCar_arm_allSynTable_del$sigCheck == "good"), ]
adenoCar_arm_allSynTable_del2$str <- paste(adenoCar_arm_allSynTable_del2$h_chr, adenoCar_arm_allSynTable_del2$m_chr)
affectedGenesLosses <- NULL
for (i in unique(adenoCar_arm_allSynTable_del2$h_chr)) {
  tmp <- adenoCar_arm_allSynTable_del2[which(adenoCar_arm_allSynTable_del2$h_chr == i),]
  tmpGr <- GRanges(seqnames = tmp$h_chr, IRanges(start = tmp$h_start, end = tmp$h_end))
  tmpGene <- hg19GeneLocations$gene[queryHits(findOverlaps(hg19ExonGr, tmpGr))]
  affectedGenesLosses <- rbind(affectedGenesLosses, data.frame("h_chr" = tmp$h_chr[subjectHits(findOverlaps(hg19ExonGr, tmpGr))],
                                                               "mchr" = tmp$m_chr[subjectHits(findOverlaps(hg19ExonGr, tmpGr))],
                                                               "gene" = tmpGene,
                                                               "type" = rep("losses", length(tmpGene))))
}

affectedGenesDf <- rbind(affectedGenesGains, affectedGenesLosses)
affectedGenesDf <- affectedGenesDf[which(affectedGenesDf$h_chr %in% c(1, 7, 8, 13, 14, 15, 18)), ]
affectedGenesDf2 <- affectedGenesDf[-which(affectedGenesDf$h_chr == 8),]
diffArmLossComboDf$arm <- unlist(diffArmLossComboDf$arm)
diffArmLossComboDf$arm2 <- str_remove(diffArmLossComboDf$arm, "with\\_")

### I need to rework this

interestedGenes <- NULL
allAffectedGeneDf <- NULL
for (i in seq_along(unique(diffArmLossComboDf$arm2))) {
  tmpDf <- diffArmLossComboDf[which(diffArmLossComboDf$arm2 == unique(diffArmLossComboDf$arm2)[i]), ]
  tmpGeneDf <- affectedGenesDf2[which(affectedGenesDf2$h_chr == unique(diffArmLossComboDf$arm2)[i]), ]
  tmpDf <- tmpDf[which(tmpDf$genes %in%  tmpGeneDf$gene), ]
  allAffectedGeneDf <- rbind(allAffectedGeneDf, tmpDf)
  if (tmpGeneDf$type[1] == "gains") {
    tmpGenes <- tmpDf$genes[which(tmpDf$logFC > (log2(3/2) * 0.9) & tmpDf$qval < 0.05)]
  } else if(tmpGeneDf$type[1] == "losses"){
    tmpGenes <- tmpDf$genes[which(tmpDf$logFC < (log2(1/2) * 0.9) & tmpDf$qval < 0.05)]
  }
  interestedGenes <- c(interestedGenes, tmpGenes)
}

interestedGenesMat <- NULL
interestedGenesDf <- NULL
for (i in seq_along(unique(diffArmLossComboDf$arm))) {
  tmpDf <- diffArmLossComboDf[which(diffArmLossComboDf$arm == unique(diffArmLossComboDf$arm)[i]), ]
  tmpDf <- tmpDf[which(tmpDf$genes %in% interestedGenes), ]
  
  ### weird but some of the genes are fitlered out in some analysis - assign it 0 if missing
  tmpDf2 <- data.frame("genes" = interestedGenes, "logFC" = 0)
  tmpDf2$logFC[match(tmpDf$genes, tmpDf2$genes)] <- tmpDf$logFC
  interestedGenesMat <- rbind(interestedGenesMat, tmpDf2$logFC)
  interestedGenesDf <- rbind(interestedGenesDf, tmpDf2)
}

rownames(interestedGenesMat) <- c("1_loss", "7_gain", "13_gain", "14_loss", "15_loss", "18_loss")
colnames(interestedGenesMat) <- interestedGenes


heatmapColGenes <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colorsBreaksGenes <- seq(-3,3, 6/1000)

interestedGenesMat2 <- interestedGenesMat
interestedGenesMat2[interestedGenesMat2 < (log2(3/2) * 0.9) & interestedGenesMat2 > (log2(1/2) * 0.9)] <- 0
pheatmap(t(interestedGenesMat2), cluster_rows = FALSE, cluster_cols = FALSE,
         color = heatmapColGenes, breaks = colorsBreaksGenes, cellwidth = 5, cellheight = 4)

### for each of these genes, look for overlaps in pathways affected

msigDf <- data.frame("genes" = interestedGenes)
msigDf$sig <- "none"
msigDf$sig <- msig_hallmark$gs_name[match(msigDf$genes, msig_hallmark$gene_symbol)]
msigDf$cgs <- "none"
msigDf$cgs <- ifelse(!is.na(match(msigDf$genes, cancerGeneCensus2$Gene.Symbol)), "yes", "no")
msigDf$cgsAll <- "none"
msigDf$cgsAll <- ifelse(!is.na(match(msigDf$genes, cancerGeneCensusAll2$Gene.Symbol)), "yes", "no")
msigDf$chr <- hg19AllGeneExonLocs$chrom[match(msigDf$genes, hg19AllGeneExonLocs$gene)]




### pathways analysis chart -made on brahm

dotplotDf <- NULL
### can't try catch the error - so adding it separately
unique(diffArmLossComboDf$arm)
tmpVector <- c("no_18", "with_18", "with_18_no_7_13", "with_7_13_18", "with_7_18_no_13", "with_13_18_no_7")
i <- tmpVector[5]
for (i in tmpVector) {
  print(i)
  
  tmpDiffDf <- diffArmLossComboDf[which(diffArmLossComboDf$arm == i), ]
  
  gene_list <- -log10(tmpDiffDf$PValue) * sign(tmpDiffDf$logFC)
  names(gene_list) <- tmpDiffDf$genes
  gene_list = sort(gene_list, decreasing = TRUE)
  
  with1 <- GSEA(gene_list, TERM2GENE = msig_hallmark, pvalueCutoff = 0.05, eps = 0)
  tmpDotplot <- data.frame("Pathways" = c("ACT. WNT_BETA_CATENIN",
                                          "ACT. MTORC1",
                                          "SUPP. APOPTOSIS"),
                           "Condition" = rep(i, 3),
                           "GeneRatio" = rep(1, 3),
                           "p.adjust" = rep(NA, 3))
  
  
  if (length(with1@result$p.adjust[which(with1@result$ID == "HALLMARK_WNT_BETA_CATENIN_SIGNALING")]) == 0) {
    tmpDotplot$p.adjust[1] <- NA
    tmpDotplot$GeneRatio1[1] <- NA
  } else{
    tmpDotplot$p.adjust[1] <- with1@result$p.adjust[which(with1@result$ID == "HALLMARK_WNT_BETA_CATENIN_SIGNALING")]
    tmpCore <- length(unlist(strsplit(with1@result$core_enrichment[which(with1@result$ID == "HALLMARK_WNT_BETA_CATENIN_SIGNALING")], "/")))
    tmpCore <- tmpCore/with1@result$setSize[which(with1@result$ID == "HALLMARK_WNT_BETA_CATENIN_SIGNALING")]
    tmpDotplot$GeneRatio[1] <- tmpCore
  }
  
  
  if (length(with1@result$p.adjust[which(with1@result$ID == "HALLMARK_MTORC1_SIGNALING")]) == 0) {
    tmpDotplot$p.adjust[2] <- NA
    tmpDotplot$GeneRatio[2] <- NA
  } else{
    tmpDotplot$p.adjust[2] <- with1@result$p.adjust[which(with1@result$ID == "HALLMARK_MTORC1_SIGNALING")]
    tmpCore <- length(unlist(strsplit(with1@result$core_enrichment[which(with1@result$ID == "HALLMARK_MTORC1_SIGNALING")], "/")))
    tmpCore <- tmpCore/with1@result$setSize[which(with1@result$ID == "HALLMARK_MTORC1_SIGNALING")]
    tmpDotplot$GeneRatio[2] <- tmpCore
  }
  
  
  if (length(with1@result$p.adjust[which(with1@result$ID == "HALLMARK_APOPTOSIS")]) == 0) {
    tmpDotplot$p.adjust[3] <- NA
    tmpDotplot$GeneRatio[3] <- NA
  } else{
    tmpDotplot$p.adjust[3] <- with1@result$p.adjust[which(with1@result$ID == "HALLMARK_APOPTOSIS")]
    tmpCore <- length(unlist(strsplit(with1@result$core_enrichment[which(with1@result$ID == "HALLMARK_APOPTOSIS")], "/")))
    tmpCore <- tmpCore/with1@result$setSize[which(with1@result$ID == "HALLMARK_APOPTOSIS")]
    tmpDotplot$GeneRatio[3] <- tmpCore
  }
  
  
  
  dotplotDf <- rbind(dotplotDf, tmpDotplot)
}


dotplotDf2 <- dotplotDf
dotplotDf2$Condition <- factor(dotplotDf2$Condition, levels = c("no_18", "with_18",
                                                                "with_18_no_7_13","with_7_13_18",
                                                                "with_7_18_no_13","with_13_18_no_7"))

pdf(file = paste0("/avatar_data5/tmp/kev/misc/20231024aneuploidyGsea", ".pdf"),
    useDingbats = FALSE, width = 10, height = 5)
ggplot(data = dotplotDf2, aes(x = Condition, y = Pathways, 
                              color = `p.adjust`, size = GeneRatio)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("Gene set enrichment")
dev.off()



### supplementary figures
coadreadMutTable <- read.table("/mnt/DATA5/tmp/kev/misc/20231012coadread_all_relevant_mut_sample_matrixV2.txt",
                               sep = "\t", header = TRUE)
coadreadMutTable2 <- coadreadMutTable[, c(1, 3:ncol(coadreadMutTable))]
coadreadMutTable2$studyID.sampleId <- str_remove(coadreadMutTable2$studyID.sampleId, ".*\\:")
coadreadMutTable2$chr1pLoss <- 0
coadreadMutTable2$chr7pGain <- 0
coadreadMutTable2$chr7qGain <- 0
coadreadMutTable2$chr8pLoss <- 0
coadreadMutTable2$chr14qLoss <- 0
coadreadMutTable2$chr15qLoss <- 0
coadreadMutTable2$chr13qGain <- 0
coadreadMutTable2$chr18pLoss <- 0
coadreadMutTable2$chr18qLoss <- 0
coadreadMutTable2$studyID.sampleId <- paste0(coadreadMutTable2$studyID.sampleId, "A")
coadreadMutTable2 <- coadreadMutTable2[which(coadreadMutTable2$studyID.sampleId %in% rownames(coadArmMemebershipDf)),]
coadArmMemebershipDf2 <- coadArmMemebershipDf[which(rownames(coadArmMemebershipDf) %in% coadreadMutTable2$studyID.sampleId), ]

coadreadMutTable2$studyID.sampleId[-which(coadreadMutTable2$studyID.sampleId %in% rownames(coadArmMemebershipDf))]
rownames(coadArmMemebershipDf2)[-which(rownames(coadArmMemebershipDf2) %in% coadreadMutTable2$studyID.sampleId)]

coadreadMutTable2$chr1pLoss[match(rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status1p1 == 1)],
                                  coadreadMutTable2$studyID.sampleId)] <- 1
coadreadMutTable2$chr7pGain[match(rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status7p1 == 1)],
                                  coadreadMutTable2$studyID.sampleId)] <- 1
coadreadMutTable2$chr7qGain[match(rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status7q1 == 1)],
                                  coadreadMutTable2$studyID.sampleId)] <- 1
coadreadMutTable2$chr8pLoss[match(rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status8p1 == 1)],
                                  coadreadMutTable2$studyID.sampleId)] <- 1
coadreadMutTable2$chr13qGain[match(rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status13q1 == 1)],
                                   coadreadMutTable2$studyID.sampleId)] <- 1
coadreadMutTable2$chr14qLoss[match(rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status14q1 == 1)],
                                   coadreadMutTable2$studyID.sampleId)] <- 1
coadreadMutTable2$chr15qLoss[match(rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status15q1 == 1)],
                                   coadreadMutTable2$studyID.sampleId)] <- 1
coadreadMutTable2$chr18pLoss[match(rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status18p1 == 1)],
                                   coadreadMutTable2$studyID.sampleId)] <- 1
coadreadMutTable2$chr18qLoss[match(rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status18q1 == 1)],
                                   coadreadMutTable2$studyID.sampleId)] <- 1

coadreadMutTable3 <- coadreadMutTable2[, -which(colnames(coadreadMutTable2) %in% c("NRAS", "HRAS"))]
colnames(coadreadMutTable3)[2] <- "KRAS/NRAS/HRAS"
coadreadMutTable3$`KRAS/NRAS/HRAS` <- 0
coadreadMutTable3$`KRAS/NRAS/HRAS`[which(apply(coadreadMutTable2[c("KRAS", "NRAS", "HRAS")],1, sum) > 0)] <- 1


coadreadMutTableMat2 <- coadreadMutTable3[, 2:ncol(coadreadMutTable3)]
rownames(coadreadMutTableMat2) <- coadreadMutTable3$studyID.sampleId
coadreadMutTableMat2 <- t(coadreadMutTableMat2)

### calculating total of with known effectors from cbioportal numbers
coadreadMutTableMainMuts <- coadreadMutTable3[, c("TP53", "chr18qLoss", "KRAS/NRAS/HRAS")]
length(which(apply(coadreadMutTableMainMuts, 1, sum) == 3))/ (378 + 155 + 61)


