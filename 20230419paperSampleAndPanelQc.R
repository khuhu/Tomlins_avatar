### tasks - can create the statistics needed for the qc portion of the paper
### i.e # of ffpe, ff tissues

### do the general qc for the panel. i.e amplicons used for segmentation
### should create on directory - for a mass run i.e tmp.txt would contain all the
### pertinenet sample names + normal. run far enough to obtain the gcCorrected amplicon ratios
### then it'd be simple to do some type of filtering for segmentation after

### can use some of the bladder sample data for qc of variants

library(stringr)

expectedLogRCalc <- function(p, na, nb, pl){
  ### p is purity, pl is ploidy and nb and na are copies of parental allele
  res <- log2( (2* (1-p) + p* (na + nb)) / pl )
}

listOfDirs <- c("Auto_user_AUS5-138-MG_cho_20210621_354_343", "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349",
                "Auto_user_AUS5-142-MG_cho_20210701_357_353", "Auto_user_AUS5-156-MG_Fearon_20210809_374_382",
                "Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384", "Auto_user_AUS5-239-BBN_mouse_bladder_MG_493_562",
                "Auto_user_AUS5-260-BBN_mouse_bladder_MG_2_514_605", "Reanalysis_AUS5-76-MG_test1_217",
                "Auto_user_AUS5-120-MG_EFD4_BBN_334_304")
seqDir <- c("/mnt/DATA3/eros_tmp/")

listOfDirs <- paste0(seqDir, listOfDirs)

i <- listOfDirs[1]
allQcSamps <- NULL
for (i in listOfDirs) {
  setwd(i)
  tmp <- system("find . -name '*\\.bc_summary.xls'", intern = TRUE)
  if (i == "/mnt/DATA3/eros_tmp/Reanalysis_AUS5-76-MG_test1_217") {
    tmp <- "./plugin_out/dummyCov/Reanalysis_AUS5-76-MG_test1_eros_217.bc_summary.xls"
  }
  tmpXls <- read.table(tmp, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  allQcSamps <- rbind(allQcSamps, tmpXls)
}

allQcSamps$UniformityNum <- as.numeric(str_remove(allQcSamps$Uniformity, "\\%"))
allQcSamps$OnTargetNum <- as.numeric(str_remove(allQcSamps$On.Target, "\\%"))

nonMouseSamps <- c("UM-MsUCAFF2_RX1", "UM-MsUCAC_1_RX2", "UM-MsUCAC_2_RX3")
allQcSamps <- allQcSamps[-which(allQcSamps$Sample.Name %in% nonMouseSamps),]

allQcSamps$type <- "FFPE"
allQcSamps$strippedName <- str_remove(allQcSamps$Sample.Name, "\\_X.*")
allQcSamps$strippedName <- str_remove(allQcSamps$strippedName, "\\_.*")
sampsPfa <- paste0("EF-D", c(14:33))
sampsPfa[which(sampsPfa %in% allQcSamps$strippedName)]

allQcSamps$type[which(sampsPfa %in% allQcSamps$strippedName)] <- "PFA"

sampsFr <- c(paste0("EF-D0", 4:9), paste0("EF-D", 10:13))
sampsFr[which(sampsFr %in% allQcSamps$strippedName)]

allQcSamps$type[which(sampsFr %in% allQcSamps$strippedName)] <- "Fr"

table(allQcSamps$type)

allQcSamps$qc <- "good"
allQcSamps$qc[which(allQcSamps$Mean.Depth < 200 | allQcSamps$UniformityNum < 80 | allQcSamps$OnTargetNum < 80)] <- "bad"

# write.table(allQcSamps, "/mnt/DATA5/tmp/kev/misc/20230614allQcSampsNewMousePanel.txt",
#             sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

base::table(allQcSamps$type, allQcSamps$qc)

### loading in all samps to look at amplicon based filtering

tmpIdx <- NULL

i <- listOfDirs[8]
for (i in listOfDirs) {
  setwd(i)
  tmp <- system("find . -name '*coverageAnalysis_out*'", intern = TRUE)
  # if (i == "/mnt/DATA3/eros_tmp/Reanalysis_AUS5-76-MG_test1_217") {
  #   tmp <- "./plugin_out/dummyCov/Reanalysis_AUS5-76-MG_test1_eros_217.bc_summary.xls"
  # }
  tmpIdx  <- c(tmpIdx , paste0(i, str_remove(tmp, "\\."), "/"))
}

writeLines(con = "/mnt/DATA5/tmp/kev/qcAllMousePaper/tmp.txt", text = tmpIdx)

ampliconAll <- read.table("/mnt/DATA5/tmp/kev/qcAllMousePaper/cnAmplicon_matrix.txt", 
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

unwantedColums <- 305:314
ampliconAll_otherCol <- ampliconAll[, unwantedColums]

ampliconAll <- ampliconAll[,-unwantedColums]



### remove duplicate columns for normal samples, can still use normal samples in analysis
### just look for consistently bad amplicons

mouseBed <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed",
                     header = FALSE, stringsAsFactors = FALSE, sep = "\t")

mouseBed$ampId <- ampliconAll$AmpliconId

normalSamples <- c("MG_8X40", "MG_11X43", "MG_13X45", "EF_D03_MG_X14",
                   "MG_18X50", "MG_6X38", "MG_21X53")
ampliconAll_norm <- ampliconAll[, c(which(colnames(ampliconAll) %in% normalSamples))]
ampliconAll_tumors <- ampliconAll[, -c(1, which(colnames(ampliconAll) %in% normalSamples))]

ampliconAll_norm <- log2(ampliconAll_norm)
ampliconAll_tumors <- log2(ampliconAll_tumors)

ampliconAll_noAmp <- ampliconAll[, -1]
ampliconAll_noAmp <- log2(ampliconAll_noAmp)
# ampliconAll_normMelt <- melt(ampliconAll_norm)
# ampliconAll_tumorMelt <- melt(ampliconAll_tumors)


# ampliconAll_normMelt$log2 <- log2(ampliconAll_normMelt$value)
# ampliconAll_tumorMelt$log2 <- log2(ampliconAll_tumorMelt$value)
# 
# ggplot(ampliconAll_normMelt, aes(x = AmpliconId, y = log2)) + 
#   geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ampliconAll_normMed <- apply(ampliconAll_norm, 1, median)
ampliconAll_tumorMed <- apply(ampliconAll_tumors, 1, median)
ampliconAll_Med <- apply(ampliconAll_noAmp, 1, median)

ampliconAll_normMean <- apply(ampliconAll_norm, 1, mean)
ampliconAll_tumorMean <- apply(ampliconAll_tumors, 1, mean)
ampliconAll_Mean <- apply(ampliconAll_noAmp, 1, mean)

res <- PMCMRplus::gesdTest(ampliconAll_tumorMed, maxr = 0.05 * length(ampliconAll_tumorMed))
res <- PMCMRplus::gesdTest(ampliconAll_tumorMed, maxr = 0.05 * length(ampliconAll_tumorMed))
res$p.value[which(res$p.value < (0.01/213))]
summary(res)


IQR.outliers <- function(x) {
  if(any(is.na(x)))
    stop("x is missing values")
  if(!is.numeric(x))
    stop("x is not numeric")
  Q3<-quantile(x,0.75)
  Q1<-quantile(x,0.25)
  IQR<-(Q3-Q1)
  left<- (Q1-(1.5*IQR))
  right<- (Q3+(1.5*IQR))
  c(x[x <left],x[x>right])
}

delGenes <- c("Rb1", "Brca1", "Trp53", "Nf1")

tumorOutliers <- IQR.outliers(ampliconAll_tumorMed)
ampliconAll$AmpliconId[which(ampliconAll_tumorMed %in% tumorOutliers)]
lqAmps <- mouseBed[which(ampliconAll_tumorMed %in% tumorOutliers),]


allOutliers <- IQR.outliers(ampliconAll_Med)
lqAmpsAll <- mouseBed[which(ampliconAll_Med %in% allOutliers),]
lqAmpsAll <- lqAmpsAll[-which(lqAmpsAll$V8 %in% delGenes),]

zscore <- function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

allMeanZ <- zscore(ampliconAll_Mean)

which(abs(allMeanZ) > 2.58)
which(abs(allMeanZ) > 1.96)
lqAmpsAllZ <- mouseBed[which(abs(allMeanZ) > 1.96),]
lqAmpsAllZ <- lqAmpsAllZ[-which(lqAmpsAllZ$V8 %in% delGenes),]
lqAmpsAllZ$idNum <- as.numeric(str_remove(lqAmpsAllZ$ampId, "AMP\\_"))

save(list = "lqAmpsAllZ", file = "/mnt/DATA5/tmp/kev/misc/20230426lqAmpsAll.RData")


### looking at the general segmentation with and without the bad amplicons and
### comparing it to the best fit 2n, 3n and 4n samples from ASCAT(snp)
### pick 10 samples to look at

### (1) no filtering whatsoever
### (2) z-score + bstat
### (3) anytype of gap filtering - regions just consistently bad

### samples 15723lt, 12641lt, 2027lte, 14397lt, 14399rt, 12669rt, 2027lts, 13268rt, 13085rot(l), kc10
### 7 real good, 3 straddling that line

segResFilt <- read.table("/mnt/DATA5/tmp/kev/qcAllMousePaper/segResultsFilt.txt", sep = "\t",
                         stringsAsFactors = FALSE, header = TRUE)

segRes <- read.table("/mnt/DATA5/tmp/kev/qcAllMousePaper/segResults.txt", sep = "\t",
                     stringsAsFactors = FALSE, header = TRUE)

load("/mnt/DATA5/tmp/kev/misc/20230424sampleQcDf.Rdata")
snpRes05 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/snpVsPanel/20220503snpRes05V3.xlsx")

# allPloidyCalls <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220517allHgscPloidySnpComp.xlsx")
allPloidyCalls <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220517allHgscPloidySnpComp_driverZero.xlsx")


panelV3 <- allPloidyCalls
colnames(panelV3)[2] <- "V3"


eros <- c("/mnt/DATA6/mouseData/copynumber/")
listOfDirectories <- c("Auto_user_AUS5-138-MG_cho_20210621_354_343", "Auto_user_AUS5-76-MG_test1_255_185",
                       "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                       "Auto_user_AUS5-120-MG_EFD4_BBN_334_304")

# 
# allNoshadSeg <- read.table("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genes/absoluteRes.txt", sep = "\t",
#                            stringsAsFactors = FALSE, header = TRUE, fill = TRUE)


panelV3$sample[which(panelV3$sample == "2027lte")] <- "mg4"
panelV3$sample[which(panelV3$sample == "13085lt")] <- "mg15"
panelV3$sample[which(panelV3$sample == "14399rt")] <- "mg20"
panelV3$sample[which(panelV3$sample == "14656peritonealmt")] <- "14656peritnealmt"
panelV3$sample[which(panelV3$sample == "133576rt")] <- "13576rt"
panelV3$sample[which(panelV3$sample == "14150rt")] <- "14150lt"
panelV3$sample[which(panelV3$sample == "14154rot")] <- "14154lt"


snpRes05$sample[which(snpRes05$sample == "2027lte")] <- "mg4"
snpRes05$sample[which(snpRes05$sample == "13085lt")] <- "mg15"
snpRes05$sample[which(snpRes05$sample == "14399rt")] <- "mg20"
snpRes05$sample[which(snpRes05$sample == "14656peritonealmt")] <- "14656peritnealmt"
snpRes05$sample[which(snpRes05$sample == "133576rt")] <- "13576rt"
snpRes05$sample[which(snpRes05$sample == "14150rt")] <- "14150lt"
snpRes05$sample[which(snpRes05$sample == "14154rot")] <- "14154lt"

panelV3 <- panelV3[match(snpRes05$sample, panelV3$sample),]

### changed to res 05 gamma 30
allAscatSegs <- read.table("/mnt/DATA5/tmp/kev/misc/20230424allSegResults_0.3.txt", sep = "\t",
                           stringsAsFactors = FALSE, header = TRUE, fill = TRUE)

load("/mnt/DATA5/tmp/kev/misc/20230424forestPlotResultSnp.Rdata")
### reduce the samples to match with noshad calls for now
allPloidyCalls$string <- paste(allPloidyCalls$sample, allPloidyCalls$ploidy_int)
allSampsCorDf2_bestV2$sampleStripped[which(allSampsCorDf2_bestV2$sampleStripped == "13085ROT(L)")] <- "13085lt"
allSampsCorDf2_bestV2$sampleStripped[which(allSampsCorDf2_bestV2$sampleStripped == "14399RT")] <- "14399rt"
allSampsCorDf2_bestV2$sampleStripped[which(allSampsCorDf2_bestV2$sampleStripped == "2027")] <- "2027lte"
allSampsCorDf2_bestV2$string4 <- paste(allSampsCorDf2_bestV2$sampleStripped, allSampsCorDf2_bestV2$ploidy)
allSampsCorDf2_bestV2$string5 <- unlist(lapply(str_split(allSampsCorDf2_bestV2$sample, " "), "[", 1))
# allSampsCorDf2_bestV3 <- allSampsCorDf2_bestV2[which(allSampsCorDf2_bestV2$string4 %in% allPloidyCalls$string), ]
# allSampsCorDf2_bestV3$string5 <- unlist(lapply(str_split(allSampsCorDf2_bestV3$sample, " "), "[", 1))

allAscatSegs_filt <- allAscatSegs[which(allAscatSegs$sample %in% allSampsCorDf2_bestV2$string5),]
allAscatSegs_filt[, 2:6] <- lapply(allAscatSegs_filt[, 2:6], as.numeric)
allAscatSegs_filt <- allAscatSegs_filt[-which(is.na(allAscatSegs_filt$chr)),]
allAscatSegs_filt$totalSc <- as.numeric(allAscatSegs_filt$nMajor) + as.numeric(allAscatSegs_filt$nMinor)
allAscatSegs_filt$chr <- paste0("chr", allAscatSegs_filt$chr)
# allAscatSegs_filt <-  allAscatSegs_filt[which(allAscatSegs_filt$chr %in% unique(allNoshadSeg_filt$chr)),]
allAscatSegs_filt$size <- allAscatSegs_filt$endpos - allAscatSegs_filt$startpos

# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "14399RT", "14399rt")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "2027_LT-E", "2027lte")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "13085ROT\\(L\\)", "13085lt")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "15774rt", "tmp")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "15774lt", "15774rt")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "tmp", "15774lt")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "2027lte",  "mg4")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "13085lt",  "mg15")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "14399rt",  "mg20")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "14656peritonealmt",  "14656peritnealmt")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "133576rt",  "13576rt")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "14150rt",  "14150lt")
# allAscatSegs_filt$sample <- str_replace_all(allAscatSegs_filt$sample, "14154rot",  "14154lt")
# 

### function taken from segment plotting, altered to have 6 graphs 3 for ngs 3 for snp
### also added predicted one copy gains and losses based on expected tumor content

### checklist
### all needed samples for snp are there ... 
### just need to iterate over the 3 samples with matching tumor contents
source("/home/kevhu/scripts/20210802syntenyFunctions.R")

rownames(allQcSamps) <- NULL
allQcSamps$strippedName2 <- str_remove(nameStripper(allQcSamps$strippedName), "\\-")
allQcSamps$strippedName2[236:272] <- tolower(str_remove(str_remove(allQcSamps$Sample.Name[236:272], "X.*"), "\\_"))
allQcSamps$strippedName2 <- str_replace_all(allQcSamps$strippedName2, "15774rt", "tmp")
allQcSamps$strippedName2 <- str_replace_all(allQcSamps$strippedName2, "15774lt", "15774rt")
allQcSamps$strippedName2 <- str_replace_all(allQcSamps$strippedName2, "tmp", "15774lt")
allQcSamps$strippedName2 <- str_replace_all(allQcSamps$strippedName2, "14656peritonealmt",  "14656peritnealmt")
allQcSamps$strippedName2 <- str_replace_all(allQcSamps$strippedName2, "133576rt",  "13576rt")
allQcSamps$strippedName2 <- str_replace_all(allQcSamps$strippedName2, "14150rt",  "14150lt")
allQcSamps$strippedName2 <- str_replace_all(allQcSamps$strippedName2, "14154rot",  "14154lt")

allQcSamps$sampleName2 <- str_remove(allQcSamps$Sample.Name, " ")
allQcSamps$sampleName2 <- str_remove_all(str_replace_all(allQcSamps$sampleName2, "\\-", "\\_"), " ")



segRes2 <- segRes
segRes2$ID <- allQcSamps$strippedName2[match(segRes2$ID, allQcSamps$sampleName2)]

segResFilt2 <- segResFilt
segResFilt2$ID <- allQcSamps$strippedName2[match(segResFilt2$ID, allQcSamps$sampleName2)]

tcDf <- read.table("/mnt/DATA5/tmp/kev/misc/20220718allSampsDf40.txt", sep = "\t",
                   stringsAsFactors = FALSE, header = TRUE)


### probably just iterate through the samples and use skeleton code of cn_graph
### the correct i.e best 2n, 3n, and 4n should be filtered in ascat data
### now just iterate per sample that I'd want to filter through


chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801mm10_graphingLimits.txt",
                          sep = "\t", stringsAsFactors = FALSE, header = TRUE)
chromTextdf <- chromTextdf[1:19, ]
chromBreak <- c(0, chromTextdf$chromBreaksPos)


### samples 15723lt, 12641lt, 2027lte, 14397lt, 14399rt, 12669rt, 2027lts, 13268rt, 13085rot(l), kc10
# sampleNames <- c("15723lt", "12641lt", "mg4", "14397lt", "mg20", "12669rt", "2027lts", "13268rt", "mg15", "kc10")

library(gridExtra)

allSampsCorDf2_bestV2$sampleStripped2 <- allSampsCorDf2_bestV2$sampleStripped
allSampsCorDf2_bestV2$sampleStripped2 <- str_replace_all(allSampsCorDf2_bestV2$sampleStripped2, "2027lte",  "mg4")
allSampsCorDf2_bestV2$sampleStripped2 <- str_replace_all(allSampsCorDf2_bestV2$sampleStripped2, "13085lt",  "mg15")
allSampsCorDf2_bestV2$sampleStripped2 <- str_replace_all(allSampsCorDf2_bestV2$sampleStripped2, "14399rt",  "mg20")

### for all samples - uncomment out when looking at all
sampleNames <- unique(allSampsCorDf2_bestV2$sampleStripped2)


y <- sampleNames[1]
for (y in sampleNames) {
  
  snpName <- allSampsCorDf2_bestV2$string5[which(allSampsCorDf2_bestV2$sampleStripped2 == y)]
  tmpSnpDf <- allAscatSegs_filt[which(allAscatSegs_filt$sample %in% snpName),]
  # colnames(tmpSnpDf)[2] <- "log2cnr"
  # tmpSnpDf <- tmpSnpDf[-which(tmpSnpDf$chr %in% xGenesLeft), ]
  
  ngsName <- y
  
  if (ngsName == "133576rt") {
    ngsName <- "13576rt"
  } else if(ngsName == "14150rt"){
    ngsName <- "14150lt"
  } else if(ngsName == "14154rot"){
    ngsName <- "14154lt"
  } else if(ngsName == "14656peritonealmt"){
    ngsName <- "14656peritnealmt"
  }
  
  ngsDf <- segRes2[which(segRes2$ID %in% ngsName), ]
  ngsDf <- ngsDf[-which(ngsDf$chrom %in% 20), ]
  
  ngsDf$seg.mean[ngsDf$seg.mean > 3] <- 3
  ngsDf$seg.mean[ngsDf$seg.mean < -3] <- -3
  
  ### need to alter code such that each transformation done to df_cn
  ### is also done to tmpSnpDf and ngsDf
  
  tmpSnpDf$chr <- str_remove(tmpSnpDf$chr, "chr")
  tmpSnpDf$startpos <- tmpSnpDf$startpos/1e6
  tmpSnpDf$endpos <- tmpSnpDf$endpos/1e6
  
  ngsDf$loc.start <- ngsDf$loc.start/1e6
  ngsDf$loc.end <- ngsDf$loc.end/1e6

  ngsDf$col <- "#000000"
  tmpSnpDf$col <- "#000000"
  
  i <- unique(tmpSnpDf$chr)[1]
  for (i in unique(tmpSnpDf$chr)) {
    tmpSnpDf$startpos[which(tmpSnpDf$chr == i)] <- tmpSnpDf$startpos[which(tmpSnpDf$chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    tmpSnpDf$endpos[which(tmpSnpDf$chr == i)] <- tmpSnpDf$endpos[which(tmpSnpDf$chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  tmpSnpDf <- tmpSnpDf[,c("sample", "chr", "startpos", "endpos", "totalSc", "col")]
  colnames(tmpSnpDf) <- c("sample", "chrom", "start", "end","cn", "col")
  
  
  for (i in unique(ngsDf$chrom)) {
    ngsDf$loc.start[which(ngsDf$chrom == i)] <- ngsDf$loc.start[which(ngsDf$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    ngsDf$loc.end[which(ngsDf$chrom == i)] <- ngsDf$loc.end[which(ngsDf$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  ngsDf <- ngsDf[,c("chrom", "loc.start", "loc.end", "seg.mean", "col")]
  colnames(ngsDf) <- c("chrom", "start", "end","cn", "col")
  
  ngsDf$cn[which(abs(ngsDf$cn) < 0.2)] <- 0
  
  ### calculating predicted tc
  
  tmpTc <- tcDf$tc[which(tcDf$sample %in% ngsName)]
  
  logRTable <- data.frame("ploidy" = c(rep(2,4), rep(3,6), rep(4,7)),
                          "tc" = rep(tmpTc, 17),
                          "na" = c(c(0,1,1,2), c(0,1,1,2,2,3), c(0,1,1,2,2,3,4)),
                          "nb" = c(c(0,0,1,1), c(0,0,1,1,2,2), c(0,0,1,1,2,2,2)))
  graphRes <- NULL
  # j <- 1
  for (z in 1:nrow(logRTable)) {
    tmpRes <- round(expectedLogRCalc(logRTable$tc[z], logRTable$na[z], logRTable$nb[z], logRTable$ploidy[z]), digits = 2)
    graphRes <- c(graphRes, tmpRes)
  }
  
  logRTable$logR <- graphRes
  logRTable$logR[which(logRTable$logR < -3)] <- -3
  logRTable$logR[which(logRTable$logR > 3)] <- 3
  
  ### getting rid of new zero lines for 3n and 4n since it's not exactly 0
  logRTable$logR[which(logRTable$ploidy == 3 & apply(logRTable[,3:4], 1, sum) == 3)] <- 0 
  logRTable$logR[which(logRTable$ploidy == 4 & apply(logRTable[,3:4], 1, sum) == 4)] <- 0 
  
  ### size 6 for png, size 3 for pdf for 10 x 6
  ### the color changes for each ngs ploidy diff - simplest way is to iterate color changes before each graph
  ### done for all except for snps - doesn't change
  
  
  for (z in 1:length(ngsDf$col)) {
  
    ### just for detecting a 1 copy change at least
    ngsExpectedGain  <- min(logRTable$logR[which(logRTable$ploidy == 2)][which(logRTable$logR[which(logRTable$ploidy == 2)] > 0)])
    ngsExpectedLoss <- max(logRTable$logR[which(logRTable$ploidy == 2)][which(logRTable$logR[which(logRTable$ploidy == 2)] < 0)])
    
    
    if(ngsDf$cn[z] >= ngsExpectedGain){
      ngsDf$col[z] <- "#FF0000"
    } else if(ngsDf$cn[z] <= ngsExpectedLoss){
      ngsDf$col[z] <- "#0000FF"
    }
  }
  
  tmpGraphGamma <- allSampsCorDf2_bestV2[which(allSampsCorDf2_bestV2$sampleStripped2 %in% y), ]
  tmpSnpDf1 <- tmpSnpDf[which(tmpSnpDf$sample %in% unique(tmpGraphGamma$string5[which(tmpGraphGamma$ploidy == 2)])), ]
  tmpSnpDf2 <- tmpSnpDf[which(tmpSnpDf$sample %in% unique(tmpGraphGamma$string5[which(tmpGraphGamma$ploidy == 3)])), ]
  tmpSnpDf3 <- tmpSnpDf[which(tmpSnpDf$sample %in% unique(tmpGraphGamma$string5[which(tmpGraphGamma$ploidy == 4)])), ]
  
  snp2 <- ggplot(tmpSnpDf1) + 
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn), colour = tmpSnpDf1$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 4, by = 1),
                                                           limits = c(0, 4.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 3.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Abs cn 2n", y, "tc:", round(tmpTc, digits = 2)))
  
  snp3 <- ggplot(tmpSnpDf2) + 
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn), colour = tmpSnpDf2$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 5, by = 1),
                                                           limits = c(0, 5.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 4.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Abs cn 3n", y, "tc:", round(tmpTc, digits = 2)))
  
  
  snp4 <- ggplot(tmpSnpDf3) +
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn), colour = tmpSnpDf3$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 6, by = 1),
                                                           limits = c(0, 6.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 5.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Abs cn 4n", y, "tc:", round(tmpTc, digits = 2)))
  
  
  ngs2 <- ggplot(ngsDf) + geom_hline(yintercept = c(logRTable$logR[which(logRTable$ploidy == 2)]), linetype = 2, alpha =  0.3) + 
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn)) + scale_colour_manual(values = ngsDf$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3.2, 3.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 2.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Log2Cnr(ngs) 2n", y, "tc:", round(tmpTc, digits = 2)))
  
  ### 3n
  ngsDf$col <- "#000000"
  
  for (z in 1:length(ngsDf$col)) {
    
    ### just for detecting a 1 copy change at least
    ngsExpectedGain  <- min(logRTable$logR[which(logRTable$ploidy == 3)][which(logRTable$logR[which(logRTable$ploidy == 3)] > 0)])
    ngsExpectedLoss <- max(logRTable$logR[which(logRTable$ploidy == 3)][which(logRTable$logR[which(logRTable$ploidy == 3)] < 0)])
    
    
    if(ngsDf$cn[z] >= ngsExpectedGain){
      ngsDf$col[z] <- "#FF0000"
    } else if(ngsDf$cn[z] <= ngsExpectedLoss){
      ngsDf$col[z] <- "#0000FF"
    }
  }
  
  ngs3 <- ggplot(ngsDf) + geom_hline(yintercept = c(logRTable$logR[which(logRTable$ploidy == 2)]), linetype = 2, alpha =  0.3) + 
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn)) + scale_colour_manual(values = ngsDf$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3.2, 3.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 2.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Log2Cnr(ngs) 3n", y, "tc:", round(tmpTc, digits = 2)))
  
  
  ### 4n
  ngsDf$col <- "#000000"
  
  for (z in 1:length(ngsDf$col)) {
    
    ### just for detecting a 1 copy change at least
    ngsExpectedGain  <- min(logRTable$logR[which(logRTable$ploidy == 4)][which(logRTable$logR[which(logRTable$ploidy == 4)] > 0)])
    ngsExpectedLoss <- max(logRTable$logR[which(logRTable$ploidy == 4)][which(logRTable$logR[which(logRTable$ploidy == 4)] < 0)])
    
    
    if(ngsDf$cn[z] >= ngsExpectedGain){
      ngsDf$col[z] <- "#FF0000"
    } else if(ngsDf$cn[z] <= ngsExpectedLoss){
      ngsDf$col[z] <- "#0000FF"
    }
  }
  
  ngs4 <- ggplot(ngsDf) + geom_hline(yintercept = c(logRTable$logR[which(logRTable$ploidy == 2)]), linetype = 2, alpha =  0.3) + 
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn)) + scale_colour_manual(values = ngsDf$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3.2, 3.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 2.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Log2Cnr(ngs) 4n", y, "tc:", round(tmpTc, digits = 2)))
  
  
  # 
  # png(paste0("/br_z1/kevin_storage/misc/20230202CompSnpNgsLogR/20230202logRComp", i, ".png"), width = 1800, height = 1100)
  # grid.arrange(ngs2, ngs3, ngs4, snp2, nrow = 2)
  # dev.off()
  
  pdf(paste0("/mnt/DATA5/tmp/kev/misc/20230501ascatSegs/", "segComp_", y, ".pdf"),
      width = 10, height = 7, useDingbats = FALSE)
  grid.arrange(ngs2, snp2, ngs3, snp3, ngs4, snp4, ncol = 2)
  dev.off()
  
}


### same but with filterd amplicons

for (y in sampleNames) {
  
  snpName <- allSampsCorDf2_bestV2$string5[which(allSampsCorDf2_bestV2$sampleStripped2 == y)]
  tmpSnpDf <- allAscatSegs_filt[which(allAscatSegs_filt$sample %in% snpName),]
  # colnames(tmpSnpDf)[2] <- "log2cnr"
  # tmpSnpDf <- tmpSnpDf[-which(tmpSnpDf$chr %in% xGenesLeft), ]
  
  ngsName <- y
  
  if (ngsName == "133576rt") {
    ngsName <- "13576rt"
  } else if(ngsName == "14150rt"){
    ngsName <- "14150lt"
  } else if(ngsName == "14154rot"){
    ngsName <- "14154lt"
  } else if(ngsName == "14656peritonealmt"){
    ngsName <- "14656peritnealmt"
  }
  
  ngsDf <- segResFilt2[which(segResFilt2$ID %in% ngsName), ]
  ngsDf <- ngsDf[-which(ngsDf$chrom %in% 20), ]
  
  if (nrow(ngsDf) == 0) {
    print(y)
    next()
  }
  
  ngsDf$seg.mean[ngsDf$seg.mean > 3] <- 3
  ngsDf$seg.mean[ngsDf$seg.mean < -3] <- -3
  
  ### need to alter code such that each transformation done to df_cn
  ### is also done to tmpSnpDf and ngsDf
  
  tmpSnpDf$chr <- str_remove(tmpSnpDf$chr, "chr")
  tmpSnpDf$startpos <- tmpSnpDf$startpos/1e6
  tmpSnpDf$endpos <- tmpSnpDf$endpos/1e6
  
  ngsDf$loc.start <- ngsDf$loc.start/1e6
  ngsDf$loc.end <- ngsDf$loc.end/1e6
  
  ngsDf$col <- "#000000"
  tmpSnpDf$col <- "#000000"
  
  i <- unique(tmpSnpDf$chr)[1]
  for (i in unique(tmpSnpDf$chr)) {
    tmpSnpDf$startpos[which(tmpSnpDf$chr == i)] <- tmpSnpDf$startpos[which(tmpSnpDf$chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    tmpSnpDf$endpos[which(tmpSnpDf$chr == i)] <- tmpSnpDf$endpos[which(tmpSnpDf$chr == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  tmpSnpDf <- tmpSnpDf[,c("sample", "chr", "startpos", "endpos", "totalSc", "col")]
  colnames(tmpSnpDf) <- c("sample", "chrom", "start", "end","cn", "col")
  
  
  for (i in unique(ngsDf$chrom)) {
    ngsDf$loc.start[which(ngsDf$chrom == i)] <- ngsDf$loc.start[which(ngsDf$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    ngsDf$loc.end[which(ngsDf$chrom == i)] <- ngsDf$loc.end[which(ngsDf$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  ngsDf <- ngsDf[,c("chrom", "loc.start", "loc.end", "seg.mean", "col")]
  colnames(ngsDf) <- c("chrom", "start", "end","cn", "col")
  
  ngsDf$cn[which(abs(ngsDf$cn) < 0.2)] <- 0
  
  ### calculating predicted tc
  
  tmpTc <- tcDf$tc[which(tcDf$sample %in% ngsName)]
  
  logRTable <- data.frame("ploidy" = c(rep(2,4), rep(3,6), rep(4,7)),
                          "tc" = rep(tmpTc, 17),
                          "na" = c(c(0,1,1,2), c(0,1,1,2,2,3), c(0,1,1,2,2,3,4)),
                          "nb" = c(c(0,0,1,1), c(0,0,1,1,2,2), c(0,0,1,1,2,2,2)))
  graphRes <- NULL
  # j <- 1
  for (z in 1:nrow(logRTable)) {
    tmpRes <- round(expectedLogRCalc(logRTable$tc[z], logRTable$na[z], logRTable$nb[z], logRTable$ploidy[z]), digits = 2)
    graphRes <- c(graphRes, tmpRes)
  }
  
  logRTable$logR <- graphRes
  logRTable$logR[which(logRTable$logR < -3)] <- -3
  logRTable$logR[which(logRTable$logR > 3)] <- 3
  
  ### getting rid of new zero lines for 3n and 4n since it's not exactly 0
  logRTable$logR[which(logRTable$ploidy == 3 & apply(logRTable[,3:4], 1, sum) == 3)] <- 0 
  logRTable$logR[which(logRTable$ploidy == 4 & apply(logRTable[,3:4], 1, sum) == 4)] <- 0 
  
  ### size 6 for png, size 3 for pdf for 10 x 6
  ### the color changes for each ngs ploidy diff - simplest way is to iterate color changes before each graph
  ### done for all except for snps - doesn't change
  
  
  for (z in 1:length(ngsDf$col)) {
    
    ### just for detecting a 1 copy change at least
    ngsExpectedGain  <- min(logRTable$logR[which(logRTable$ploidy == 2)][which(logRTable$logR[which(logRTable$ploidy == 2)] > 0)])
    ngsExpectedLoss <- max(logRTable$logR[which(logRTable$ploidy == 2)][which(logRTable$logR[which(logRTable$ploidy == 2)] < 0)])
    
    
    if(ngsDf$cn[z] >= ngsExpectedGain){
      ngsDf$col[z] <- "#FF0000"
    } else if(ngsDf$cn[z] <= ngsExpectedLoss){
      ngsDf$col[z] <- "#0000FF"
    }
  }
  
  tmpGraphGamma <- allSampsCorDf2_bestV2[which(allSampsCorDf2_bestV2$sampleStripped2 %in% y), ]
  tmpSnpDf1 <- tmpSnpDf[which(tmpSnpDf$sample %in% unique(tmpGraphGamma$string5[which(tmpGraphGamma$ploidy == 2)])), ]
  tmpSnpDf2 <- tmpSnpDf[which(tmpSnpDf$sample %in% unique(tmpGraphGamma$string5[which(tmpGraphGamma$ploidy == 3)])), ]
  tmpSnpDf3 <- tmpSnpDf[which(tmpSnpDf$sample %in% unique(tmpGraphGamma$string5[which(tmpGraphGamma$ploidy == 4)])), ]
  
  snp2 <- ggplot(tmpSnpDf1) + 
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn), colour = tmpSnpDf1$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 4, by = 1),
                                                           limits = c(0, 4.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 3.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Abs cn 2n", y, "tc:", round(tmpTc, digits = 2)))
  
  snp3 <- ggplot(tmpSnpDf2) + 
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn), colour = tmpSnpDf2$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 5, by = 1),
                                                           limits = c(0, 5.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 4.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Abs cn 3n", y, "tc:", round(tmpTc, digits = 2)))
  
  
  snp4 <- ggplot(tmpSnpDf3) +
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn), colour = tmpSnpDf3$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 6, by = 1),
                                                           limits = c(0, 6.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 5.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Abs cn 4n", y, "tc:", round(tmpTc, digits = 2)))
  
  
  ngs2 <- ggplot(ngsDf) + geom_hline(yintercept = c(logRTable$logR[which(logRTable$ploidy == 2)]), linetype = 2, alpha =  0.3) + 
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn)) + scale_colour_manual(values = ngsDf$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3.2, 3.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 2.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Log2Cnr(ngs) 2n", y, "tc:", round(tmpTc, digits = 2)))
  
  ### 3n
  ngsDf$col <- "#000000"
  
  for (z in 1:length(ngsDf$col)) {
    
    ### just for detecting a 1 copy change at least
    ngsExpectedGain  <- min(logRTable$logR[which(logRTable$ploidy == 3)][which(logRTable$logR[which(logRTable$ploidy == 3)] > 0)])
    ngsExpectedLoss <- max(logRTable$logR[which(logRTable$ploidy == 3)][which(logRTable$logR[which(logRTable$ploidy == 3)] < 0)])
    
    
    if(ngsDf$cn[z] >= ngsExpectedGain){
      ngsDf$col[z] <- "#FF0000"
    } else if(ngsDf$cn[z] <= ngsExpectedLoss){
      ngsDf$col[z] <- "#0000FF"
    }
  }
  
  ngs3 <- ggplot(ngsDf) + geom_hline(yintercept = c(logRTable$logR[which(logRTable$ploidy == 2)]), linetype = 2, alpha =  0.3) + 
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn)) + scale_colour_manual(values = ngsDf$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3.2, 3.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 2.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Log2Cnr(ngs) 3n", y, "tc:", round(tmpTc, digits = 2)))
  
  
  ### 4n
  ngsDf$col <- "#000000"
  
  for (z in 1:length(ngsDf$col)) {
    
    ### just for detecting a 1 copy change at least
    ngsExpectedGain  <- min(logRTable$logR[which(logRTable$ploidy == 4)][which(logRTable$logR[which(logRTable$ploidy == 4)] > 0)])
    ngsExpectedLoss <- max(logRTable$logR[which(logRTable$ploidy == 4)][which(logRTable$logR[which(logRTable$ploidy == 4)] < 0)])
    
    
    if(ngsDf$cn[z] >= ngsExpectedGain){
      ngsDf$col[z] <- "#FF0000"
    } else if(ngsDf$cn[z] <= ngsExpectedLoss){
      ngsDf$col[z] <- "#0000FF"
    }
  }
  
  ngs4 <- ggplot(ngsDf) + geom_hline(yintercept = c(logRTable$logR[which(logRTable$ploidy == 2)]), linetype = 2, alpha =  0.3) + 
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn)) + scale_colour_manual(values = ngsDf$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3.2, 3.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 2.3, label = chrom), size = 2.5) + 
    theme(plot.title = element_text(hjust = 0.5)) + ggtitle(paste("Log2Cnr(ngs) 4n", y, "tc:", round(tmpTc, digits = 2)))
  
  
  # 
  # png(paste0("/br_z1/kevin_storage/misc/20230202CompSnpNgsLogR/20230202logRComp", i, ".png"), width = 1800, height = 1100)
  # grid.arrange(ngs2, ngs3, ngs4, snp2, nrow = 2)
  # dev.off()
  
  pdf(paste0("/mnt/DATA5/tmp/kev/misc/20230501ascatSegs/", "segCompAmpFilt_", y, ".pdf"),
      width = 10, height = 7, useDingbats = FALSE)
  grid.arrange(ngs2, snp2, ngs3, snp3, ngs4, snp4, ncol = 2)
  dev.off()
  
}



### after I create the graphs do disjoin for all samples and can do comparison by portion of the genome
### also calculate percentage of genome altered - i.e calling method may only be good for certain samples with >= pga




### calculating fga



y <- sampleNames[1]
allFga <- NULL
for (y in sampleNames) {
  tmpFga <- data.frame()

  
  snpName <- allSampsCorDf2_bestV2$string5[which(allSampsCorDf2_bestV2$sampleStripped2 == y)]
  tmpSnpDf <- allAscatSegs_filt[which(allAscatSegs_filt$sample %in% snpName),]
  tmpGraphGamma <- allSampsCorDf2_bestV2[which(allSampsCorDf2_bestV2$sampleStripped2 %in% y), ]
  tmpSnpDf1 <- tmpSnpDf[which(tmpSnpDf$sample %in% unique(tmpGraphGamma$string5[which(tmpGraphGamma$ploidy == 2)])), ]
  tmpSnpDf2 <- tmpSnpDf[which(tmpSnpDf$sample %in% unique(tmpGraphGamma$string5[which(tmpGraphGamma$ploidy == 3)])), ]
  tmpSnpDf3 <- tmpSnpDf[which(tmpSnpDf$sample %in% unique(tmpGraphGamma$string5[which(tmpGraphGamma$ploidy == 4)])), ]
  
  
  tmpSnpDf1$cnChange <- tmpSnpDf1$totalSc - 2
  tmpSnpDf2$cnChange <- tmpSnpDf2$totalSc - 3
  tmpSnpDf3$cnChange <- tmpSnpDf3$totalSc - 4
  
  tmpSnpDf1$length <- tmpSnpDf1$endpos - tmpSnpDf1$startpos
  tmpSnpDf2$length <- tmpSnpDf2$endpos - tmpSnpDf2$startpos
  tmpSnpDf3$length <- tmpSnpDf3$endpos - tmpSnpDf3$startpos

  snpFga1 <- sum(tmpSnpDf1$length[(tmpSnpDf1$cnChange != 0)])/sum(tmpSnpDf1$length)
  snpFga2 <- sum(tmpSnpDf2$length[(tmpSnpDf2$cnChange != 0)])/sum(tmpSnpDf2$length)
  snpFga3 <- sum(tmpSnpDf3$length[(tmpSnpDf3$cnChange != 0)])/sum(tmpSnpDf3$length)
  
  tmpFga <- rbind(tmpFga, data.frame("sample" = c(tmpSnpDf1$sample[1], tmpSnpDf2$sample[2], tmpSnpDf3$sample[3]),
                                     "fga_snp" = c(round(as.numeric(snpFga1), digits = 2), round(as.numeric(snpFga2), digits = 2),
                                                      round(as.numeric(snpFga3), digits = 2)),
                                     "ploidy" = c(2:4)))
  
  
  ngsName <- y
  if (ngsName == "133576rt") {
    ngsName <- "13576rt"
  } else if(ngsName == "14150rt"){
    ngsName <- "14150lt"
  } else if(ngsName == "14154rot"){
    ngsName <- "14154lt"
  } else if(ngsName == "14656peritonealmt"){
    ngsName <- "14656peritnealmt"
  }
  
  ngsDf <- segResFilt2[which(segResFilt2$ID %in% ngsName), ]
  ngsDf <- ngsDf[-which(ngsDf$chrom %in% 20), ]
  ngsDf$seg.mean[which(abs(ngsDf$seg.mean) < 0.2)] <- 0
  ngsDf$length <- ngsDf$loc.end - ngsDf$loc.start
  ### calculating predicted tc
  
  tmpTc <- tcDf$tc[which(tcDf$sample %in% ngsName)]
  
  logRTable <- data.frame("ploidy" = c(rep(2,4), rep(3,6), rep(4,7)),
                          "tc" = rep(tmpTc, 17),
                          "na" = c(c(0,1,1,2), c(0,1,1,2,2,3), c(0,1,1,2,2,3,4)),
                          "nb" = c(c(0,0,1,1), c(0,0,1,1,2,2), c(0,0,1,1,2,2,2)))
  graphRes <- NULL
  # j <- 1
  for (z in 1:nrow(logRTable)) {
    tmpRes <- round(expectedLogRCalc(logRTable$tc[z], logRTable$na[z], logRTable$nb[z], logRTable$ploidy[z]), digits = 2)
    graphRes <- c(graphRes, tmpRes)
  }
  
  logRTable$logR <- graphRes
  logRTable$logR[which(logRTable$logR < -3)] <- -3
  logRTable$logR[which(logRTable$logR > 3)] <- 3
  
  ### getting rid of new zero lines for 3n and 4n since it's not exactly 0
  logRTable$logR[which(logRTable$ploidy == 3 & apply(logRTable[,3:4], 1, sum) == 3)] <- 0 
  logRTable$logR[which(logRTable$ploidy == 4 & apply(logRTable[,3:4], 1, sum) == 4)] <- 0 
  
  tmpNgs <- NULL
  tmpNgsLen <- NULL
  for (i in 2:4) {
    ngsExpectedGain  <- min(logRTable$logR[which(logRTable$ploidy == i)][which(logRTable$logR[which(logRTable$ploidy == i)] > 0)])
    ngsExpectedLoss <- max(logRTable$logR[which(logRTable$ploidy == i)][which(logRTable$logR[which(logRTable$ploidy == i)] < 0)])
    
    ngsExpectedGain <- ngsExpectedGain * 0.9
    ngsExpectedLoss <- ngsExpectedLoss * 0.9
    
    tryCatch({ngsExpectedGain  <- min(logRTable$logR[which(logRTable$ploidy == i)][which(logRTable$logR[which(logRTable$ploidy == i)] > 0)])},
             warning=function(w) print(y))
    
    ngsFga <- sum(ngsDf$length[which(ngsDf$seg.mean > ngsExpectedGain | ngsDf$seg.mean < ngsExpectedLoss)])/sum(ngsDf$length)
    ngsFgaLenient <- sum(ngsDf$length[which(abs(ngsDf$seg.mean) > 0.2)])/sum(ngsDf$length)
    tmpNgs <- c(tmpNgs, round(ngsFga, digits = 2))
    tmpNgsLen <- c(tmpNgsLen, round(ngsFgaLenient, digits = 2))
  }
  
  tmpFga <- cbind(tmpFga, data.frame("fga_ngs" = tmpNgs, "fga_ngs_lenient" = tmpNgsLen))
  allFga <- rbind(allFga, tmpFga)
}




allFga$logChiSq <- allSampsCorDf2_bestV2$logChiSq[match(allFga$sample, allSampsCorDf2_bestV2$string5)]
allFga$spearman <- allSampsCorDf2_bestV2$spearman[match(allFga$sample, allSampsCorDf2_bestV2$string5)]
allFga$diff <- abs(allFga$fga_snp - allFga$fga_ngs)

minChi  <- min(allFga$logChiSq)/25
maxChi  <- max(allFga$logChiSq)/25

library(ggplot2)

a <- ggplot(data = allFga, aes(x = fga_snp, y = fga_ngs, color = factor(ploidy))) + geom_point(aes(size = logChiSq)) + facet_wrap(~factor(ploidy)) +
  geom_abline(intercept =0 , slope = 1) + scale_size_continuous(range = c(minChi, maxChi)) + ggtitle("1 copy loss fraction of genome altered") + 
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1)) + scale_x_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5))


b <- ggplot(data = allFga, aes(x = fga_snp, y = fga_ngs_lenient, color = factor(ploidy))) + geom_point(aes(size = logChiSq)) + facet_wrap(~factor(ploidy)) +
  geom_abline(intercept =0 , slope = 1) + scale_size_continuous(range = c(minChi, maxChi)) + ggtitle("1 copy loss fraction of genome altered (0.2 thres)") + 
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1)) + scale_x_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5))


grid.arrange(a, b, nrow = 2)


### title the graph and then make it - shows its either lack of subclonal calling from the caller
### make the confusion matrix but with subclonal calling
### using code from comp 4 to properly get 1 of each ploidy and then look at subclonal calling


### make loop similar to previous the ones for graph above to produce 2n, 3n and 4n comparison
###
###

### make ngs matrix outside of loop since it doesn't change
sampleNames2 <- sampleNames
sampleNames2[which(sampleNames2 =="133576rt" )] <- "13576rt"
sampleNames2[which(sampleNames2 =="14150rt" )] <- "14150lt"
sampleNames2[which(sampleNames2 =="14154rot" )] <- "14154lt"
sampleNames2[which(sampleNames2 =="14656peritonealmt" )] <- "14656peritnealmt"

segResFilt3 <- segResFilt2[which(segResFilt2$ID %in% sampleNames2), ]
segResFilt3$chrom <- paste0("chr", segResFilt3$chrom)
segResFilt3$size <- segResFilt3$loc.end - segResFilt3$loc.start
noshadAneploidySub <- matrix(0, nrow = length(unique(sampleNames2)), ncol = length(unique(segResFilt3$chrom)))
rownames(noshadAneploidySub) <- unique(sampleNames2)
colnames(noshadAneploidySub) <- unique(segResFilt3$chrom)
tmpChrList <- unique(segResFilt3$chrom)[1:19]

### use vector of cutoffs
# i <- rownames(noshadAneploidySub)[39]
# for (i in rownames(noshadAneploidySub)) {
#   tmpDf <- segResFilt3[grep(i, segResFilt3$ID),]
#   # tmpPloidy <- as.numeric(panelV3$V3[which(panelV3$sample == i)])
#   # tmpDf$totalSc_n <- tmpDf$sC- tmpPloidy
#   
#   j <- 1
#   for (j in seq_along(tmpChrList)) {
#     tmpChr <- tmpDf[which(tmpDf$chrom == tmpChrList[j]), ]
#     # chrThres <- sum(tmpChr$size) * vectorOfCutoffs[j]
#     chrThres <- sum(tmpChr$size) * 0.7
#     tmpAneu <- "none"
#     if (sum(tmpChr$size[which(tmpChr$seg.mean > 0.2)]) > chrThres) {
#       tmpAneu <- "gain"
#     } else if(sum(tmpChr$size[which(tmpChr$seg.mean < -0.2)]) > chrThres){
#       tmpAneu <- "loss"
#     } else{
#       tmpAneu <- "none"
#     }
#     noshadAneploidySub[i,j] <- tmpAneu
#   }
# }


allAscatSegs_filt2 <- allAscatSegs_filt
allAscatSegs_filt2$ploidy <- NA
allAscatSegs_filt2$ploidy <- allSampsCorDf2_bestV2$ploidy[match(allAscatSegs_filt2$sample, allSampsCorDf2_bestV2$string5)]
allAscatSegs_filt2$string <- paste(allAscatSegs_filt2$sample, allAscatSegs_filt2$chr,
                                    allAscatSegs_filt2$startpos, allAscatSegs_filt2$endpos,
                                    allAscatSegs_filt2$ploidy)


allAscatSegs_filt2 <- allAscatSegs_filt2[-which(duplicated(allAscatSegs_filt2$string)), ]
y <- 2
for (y in 2:4) {

  tmpSnpDf <- allAscatSegs_filt2[which(allAscatSegs_filt2$ploidy == y),]

  ### calculating aneuploidy gains/losses
  ascatAneploidy <- matrix(0, ncol = length(unique(tmpSnpDf$chr)), nrow = length(unique(tmpSnpDf$sample)))
  rownames(ascatAneploidy) <- unique(tmpSnpDf$sample)
  colnames(ascatAneploidy) <- paste0("chr", 1:19)

  i <- rownames(ascatAneploidy)[1]
  for (i in rownames(ascatAneploidy)) {
    tmpDf <- tmpSnpDf[grep(i, tmpSnpDf$sample),]
    tmpPloidy <- y
    tmpDf$totalSc_n <- tmpDf$totalSc - tmpPloidy
    j <- 1
    for (j in seq_along(tmpChrList)) {
      tmpChr <- tmpDf[which(tmpDf$chr == tmpChrList[j]), ]
      chrThres <- sum(tmpChr$size) * 0.7
      tmpAneu <- "none"
      if (sum(tmpChr$size[which(tmpChr$totalSc_n > 0)]) > chrThres) {
        tmpAneu <- "gain"
      } else if(sum(tmpChr$size[which(tmpChr$totalSc_n < 0)]) > chrThres){
        tmpAneu <- "loss"
      } else{
        tmpAneu <- "none"
      }
      ascatAneploidy[i,j] <- tmpAneu
    }
  }
  assign(x = paste0("ascatAneploidy", y, "n"), value = ascatAneploidy)
}

rownames(ascatAneploidy2n) <- rownames(noshadAneploidySub)
rownames(ascatAneploidy3n) <- rownames(noshadAneploidySub)
rownames(ascatAneploidy4n) <- rownames(noshadAneploidySub)

library(caret)
# load("/mnt/DATA5/tmp/kev/misc/20230424sampleQcDf.Rdata")
# load("/mnt/DATA5/tmp/kev/misc/20230519sampleQcDf.Rdata")
load("/mnt/DATA5/tmp/kev/misc/20230614sampleQcDf.Rdata")

sampleQcDf
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "14399RT", "14399rt")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "2027_LT-E", "2027lte")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "13085ROT\\(L\\)", "13085lt")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "15774rt", "tmp")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "15774lt", "15774rt")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "tmp", "15774lt")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "2027lte",  "mg4")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "13085lt",  "mg15")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "14399rt",  "mg20")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "14656peritonealmt",  "14656peritnealmt")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "133576rt",  "13576rt")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "14150rt",  "14150lt")
sampleQcDf$sample <- str_replace_all(sampleQcDf$sample, "14154rot",  "14154lt")
goodSamps <- sampleQcDf$sample[which(sampleQcDf$qc == "good")]
badSamps <- sampleQcDf$sample[which(sampleQcDf$qc == "bad")]



### even with differing ploidy states - there are just some bad types of calls, could be b/c of bad regions
### do disjoin to for areas that are bad
### after finalizing this. do best ploidy provided by ASCAT and see comparison to just cnr - might just be the convex algo

tmpGrange1 <-  allAscatSegs_filt2[which(allAscatSegs_filt2$ploidy == 2), c("chr", "startpos", "endpos", "sample")]
allAscat2nGrange <- GRanges(seqnames = tmpGrange1$chr,
                            ranges = IRanges(start = tmpGrange1$startpos, end = tmpGrange1$endpos))
tmpGrange2 <-  allAscatSegs_filt2[which(allAscatSegs_filt2$ploidy == 3), c("chr", "startpos", "endpos", "sample")]
allAscat3nGrange <- GRanges(seqnames = tmpGrange2$chr,
                            ranges = IRanges(start = tmpGrange2$startpos, end = tmpGrange2$endpos))
tmpGrange3 <-  allAscatSegs_filt2[which(allAscatSegs_filt2$ploidy == 4), c("chr", "startpos", "endpos", "sample")]
allAscat4nGrange <- GRanges(seqnames = tmpGrange3$chr,
                            ranges = IRanges(start = tmpGrange3$startpos, end = tmpGrange3$endpos))

tmpGrange4 <- segResFilt3[which(segResFilt3$chrom %in% paste0("chr", 1:19)), c("chrom", "loc.start", "loc.end")]
colnames(tmpGrange4) <- c("chr", "startpos", "endpos")
allNgsGrange <- GRanges(seqnames = tmpGrange4$chr,
                            ranges = IRanges(start =  tmpGrange4$startpos, end = tmpGrange4$endpos))

### find all segments called, then which ones are exclusively found in ngs
tmpGrange5 <- rbind(tmpGrange1, tmpGrange2, tmpGrange3)
allGrange <- GRanges(seqnames = tmpGrange5$chr,
                     ranges = IRanges(start =  tmpGrange5$startpos, end = tmpGrange5$endpos))

library(GenomicRanges)

allNgsDis <- disjoin(GRanges(seqnames = paste0("chr", segResFilt2$chrom), IRanges(start = segResFilt2$loc.start, end = segResFilt2$loc.end)))
disjoinRanges <- disjoin(allGrange)
allNgsDisOverlap <- allNgsDis[unique(subjectHits(findOverlaps(disjoinRanges, allNgsDis))),]

### all there minus x chromosome I forgot to remove

allNgsDisDf <- data.frame(allNgsDisOverlap)
allNgsDisDf <- allNgsDisDf[-which(allNgsDisDf$width < 2),]
allNgsDisReducedGrange <- GRanges(seqnames = allNgsDisDf$seqnames, IRanges(start = allNgsDisDf$start, end = allNgsDisDf$end))
mouseBed2 <- mouseBed[-lqAmpsAllZ$idNum,]

### next thing to do is for 2n, 3n and 4n is to find which areas with copy number change have the most disagreement
### then some type of coverage metric can be calculated i.e number of markers per window ... could correlate with how well aneuploidy
## and large CNAs are called

ngsDisjoin <- NULL
for (i in unique(segResFilt3$ID)) {
  print(i)
  tmpNgs <- segResFilt3[which(segResFilt3$ID == i),]
  tmpRes <- allNgsDisDf[,1:3]
  tmpRes$sample <- i
  queryIdx <- 1:nrow(tmpRes)
  tmpGrange <- GRanges(seqnames = tmpNgs$chrom, IRanges(start = tmpNgs$loc.start, end = tmpNgs$loc.end))
  tmpIdx <- subjectHits(findOverlaps(query = allNgsDisReducedGrange, subject = tmpGrange))
  testQuery <- queryHits(findOverlaps(query = allNgsDisReducedGrange, subject = tmpGrange))
  if (length(which(duplicated(queryHits(findOverlaps(query = allNgsDisReducedGrange, subject = tmpGrange))))) > 0) {
    tmpIdx <-  tmpIdx[-which(duplicated(queryHits(findOverlaps(query = allNgsDisReducedGrange, subject = tmpGrange))))]
    testQuery <- testQuery[-which(duplicated(queryHits(findOverlaps(query = allNgsDisReducedGrange, subject = tmpGrange))))]
  }
  ### error probably happens b/c there is no match for one of the queries
  tmpRes$cn <- NA
  tmpRes$cn[testQuery] <- tmpNgs$seg.mean[tmpIdx]
  tmpRes2 <- tmpRes[, c("sample", "seqnames", "start", "end", "cn")]
  ngsDisjoin  <- rbind(ngsDisjoin , tmpRes2)
}


### the idea of doing this is to figure out why certain areas and
### sizes of CNAs are good/bad concordance

### before I did per chromosome based on status calls, I should do
### something similar but just general correlation of signal .. i.e bad areas will be detected ....
### for areas that are bad I can then use a marker per window metric to gauge how many are necessary ....

ngsDisjoin <- ngsDisjoin[order(ngsDisjoin$sample),]
ngsDisjoin$string <- paste(ngsDisjoin$seqnames, ngsDisjoin$start, ngsDisjoin$end)

### for paper figure about concordance of calls, it'll be a barplot with 5 different measurements per chromosome - in chromosome order
### spearman for genes on the chromosome, recall for aneuploidy losses and gain, and recall for large CNAs losses and gain


library(Repitools)
library(BSgenome.Mmusculus.UCSC.mm10)
mm10_1Mb <- genomeBlocks(genome = BSgenome.Mmusculus.UCSC.mm10, seqnames(BSgenome.Mmusculus.UCSC.mm10),
                         width = 1e6)
mm10_5Mb <- genomeBlocks(genome = BSgenome.Mmusculus.UCSC.mm10, seqnames(BSgenome.Mmusculus.UCSC.mm10),
                         width = 5e6)
mm10_10Mb <- genomeBlocks(genome = BSgenome.Mmusculus.UCSC.mm10, seqnames(BSgenome.Mmusculus.UCSC.mm10),
                         width = 1e7)
mm10_20Mb <- genomeBlocks(genome = BSgenome.Mmusculus.UCSC.mm10, seqnames(BSgenome.Mmusculus.UCSC.mm10),
                          width = 2e7)

mm10_1Mb <- data.frame(mm10_1Mb)
mm10_1Mb <- mm10_1Mb[which(mm10_1Mb$seqnames %in% paste0("chr", 1:19)),]
mm10_1MbGr  <- GRanges(seqnames = mm10_1Mb$seqnames, IRanges(start = mm10_1Mb$start, end = mm10_1Mb$end))

mm10_5Mb <- data.frame(mm10_5Mb)
mm10_5Mb <- mm10_5Mb[which(mm10_5Mb$seqnames %in% paste0("chr", 1:19)),]
mm10_5MbGr  <- GRanges(seqnames = mm10_5Mb$seqnames, IRanges(start = mm10_5Mb$start, end = mm10_5Mb$end))

mm10_10Mb <- data.frame(mm10_10Mb)
mm10_10Mb <- mm10_10Mb[which(mm10_10Mb$seqnames %in% paste0("chr", 1:19)),]
mm10_10MbGr  <- GRanges(seqnames = mm10_10Mb$seqnames, IRanges(start = mm10_10Mb$start, end = mm10_10Mb$end))


mm10_20Mb <- data.frame(mm10_20Mb)
mm10_20Mb <- mm10_20Mb[which(mm10_20Mb$seqnames %in% paste0("chr", 1:19)),]
mm10_20MbGr  <- GRanges(seqnames = mm10_20Mb$seqnames, IRanges(start = mm10_20Mb$start, end = mm10_20Mb$end))

### need to reduce mouse bed into marker file where every gene only counts for 1 marker before I overlap
### With the windows, is there a minimum I need to detect aneuploidy 

mouseBedMarker <- mouseBed2[-which(duplicated(mouseBed2$V8)),]
removeMarkers <- c("Esr1303K", "Esr1537Y", "Ar875H", "Gnaq209Q", "TertUtr146",
                   "Gna11183R", "Pdgfra842D", "KIT", "Fgfr3249S", "Fgfr3375Y",
                   "Pdgfra842D", "Sf3b1625R")

mouseBedMarker <- mouseBedMarker[-which(mouseBedMarker$V8 %in% removeMarkers),]
mouseBedMarker <- mouseBedMarker[which(mouseBedMarker$V1 %in% paste0("chr", 1:19)),]
mouseBedMarkerGrange <- GRanges(seqnames = mouseBedMarker$V1, IRanges(start = mouseBedMarker$V2, end = mouseBedMarker$V3))


mm10_1Mb$markerCount <- 0
for (i in 1:length(mm10_1MbGr)) {
  tmpOverlaps <- subjectHits(findOverlaps(mm10_1MbGr[i], mouseBedMarkerGrange))
  mm10_1Mb$markerCount[i] <- length(tmpOverlaps)
}


mm10_5Mb$markerCount <- 0
for (i in 1:length(mm10_5MbGr)) {
  tmpOverlaps <- subjectHits(findOverlaps(mm10_5MbGr[i], mouseBedMarkerGrange))
  mm10_5Mb$markerCount[i] <- length(tmpOverlaps)
}


mm10_10Mb$markerCount <- 0
for (i in 1:length(mm10_10MbGr)) {
  tmpOverlaps <- subjectHits(findOverlaps(mm10_10MbGr[i], mouseBedMarkerGrange))
  mm10_10Mb$markerCount[i] <- length(tmpOverlaps)
}



mm10_20Mb$markerCount <- 0
for (i in 1:length(mm10_20MbGr)) {
  tmpOverlaps <- subjectHits(findOverlaps(mm10_20MbGr[i], mouseBedMarkerGrange))
  mm10_20Mb$markerCount[i] <- length(tmpOverlaps)
}

### okay I have windows with markers in them 2 things to look at
### (1) look at distribution of correlation for minimal common regions and how it relates to the region they lie in i.e # of markers per window
### this might be harder to interpret b/c 
### (2)
### then test the calling of both aneuploidy and and CNAsand seeing how the calls are by with each filter

### doing simple filtering first for aneuploidy calls - only need to filter ngs
colnames(ngsDisjoin)[c(2,5)] <- c("chrom", "seg.mean")
ngsDisjoinGr <- GRanges(seqnames = ngsDisjoin$chrom, IRanges(start = ngsDisjoin$start, end = ngsDisjoin$end))
ngsDisjoin$size <- ngsDisjoin$end - ngsDisjoin$start

### loop for the above for multiple tables
for (i in c(0:4)) {
  tmp1Mb <- data.frame(mm10_1MbGr[which(mm10_1Mb$markerCount > i)])
  tmp5Mb <- data.frame(mm10_5MbGr[which(mm10_5Mb$markerCount > i)])
  tmp10Mb <- data.frame(mm10_10MbGr[which(mm10_10Mb$markerCount > i)])
  tmp20Mb <- data.frame(mm10_10MbGr[which(mm10_20Mb$markerCount > i)])
  
  write.table(tmp1Mb, paste0("/mnt/DATA5/tmp/kev/misc/20230606mm10_", i, "_1Mb.bed"), sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(tmp5Mb, paste0("/mnt/DATA5/tmp/kev/misc/20230606mm10_", i, "_5Mb.bed"), sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(tmp10Mb, paste0("/mnt/DATA5/tmp/kev/misc/20230606mm10_", i, "_10Mb.bed"), sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(tmp20Mb, paste0("/mnt/DATA5/tmp/kev/misc/20230606mm10_", i, "_20Mb.bed"), sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}



### below is obsolete - made into loop to get more variations
### only possibly important pieces of code are (1) per chromosome barplot
### (2) plot for precision and recall stats 
### before runnign make sure I have the proper bedtool results 

tmpDf2n <- allAscatSegs_filt2[which(allAscatSegs_filt2$ploidy == 2), ]
tmpDf3n <- allAscatSegs_filt2[which(allAscatSegs_filt2$ploidy == 3), ]
tmpDf4n <- allAscatSegs_filt2[which(allAscatSegs_filt2$ploidy == 4), ]

### i is for Mb and j is for number of markers

i <- 1
j <- 1
allCnaStatsDf2 <- NULL
for (i in c(1, 5, 10, 20)) {
  
  for (j in c(0:4)) {
    ### read in tables after they've gone through bedtools
    tmpDisJoin <- read.table(paste0("/mnt/DATA5/tmp/kev/misc/20230606ngsInter_", j, "_", i, "Mb.bed"),
                             sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    
    # colnames(tmpDisJoin) <- colnames(ngsDisjoinBed)
    colnames(tmpDisJoin) <- c("chrom", "start", "end","sample", "seg.mean", "str", "size")
    
    tmpDisJoinAneuploidy <- matrix(0, nrow = length(unique(sampleNames2)), ncol = length(unique(segResFilt3$chrom)))
    rownames(tmpDisJoinAneuploidy) <- unique(sampleNames2)
    colnames(tmpDisJoinAneuploidy) <- unique(segResFilt3$chrom)
    
    tmp_cna <- NULL
    k <- rownames(tmpDisJoinAneuploidy)[1]
    for (k in rownames(tmpDisJoinAneuploidy)) {
      tmpDf <- tmpDisJoin[grep(k, tmpDisJoin$sample),]
      
      for (l in seq_along(tmpChrList)) {
        tmpChr <- tmpDf[which(tmpDf$chrom == tmpChrList[l]), ]
        tmpChr2 <- tmpChr
        chrThres <- sum(tmpChr$size) * 0.7
        tmpAneu <- "none"
        if (sum(tmpChr$size[which(tmpChr$seg.mean > 0.3)]) > chrThres) {
          tmpAneu <- "gain"
          tmpChr2$seg.mean[which(tmpChr2$seg.mean > 0.3)] <- 0
        } else if(sum(tmpChr$size[which(tmpChr$seg.mean < -0.4)]) > chrThres){
          tmpAneu <- "loss"
          tmpChr2$seg.mean[which(tmpChr2$seg.mean < -0.4)] <- 0
        } else{
          tmpAneu <- "none"
        }
        tmpDisJoinAneuploidy[k, l] <- tmpAneu
        tmp_cna <- rbind(tmp_cna, tmpChr2)
      }
    }
    
    ### commented out to compare good samples to all
    
    tmpGood <- tmpDisJoinAneuploidy[which(rownames(tmpDisJoinAneuploidy) %in% goodSamps), 1:19]
    # tmpGood <- tmpDisJoinAneuploidy[, 1:19]
    tmpGoodVec <- as.vector(tmpGood)
    tmpGoodVec <- factor(tmpGoodVec, levels = c("none", "gain", "loss"))
    # 
    goodAscatVec2n <- factor(ascatAneploidy2n[which(rownames(ascatAneploidy2n) %in% goodSamps), ], levels = c("none", "gain", "loss"))
    goodAscatVec3n <- factor(ascatAneploidy3n[which(rownames(ascatAneploidy3n) %in% goodSamps), ], levels = c("none", "gain", "loss"))
    goodAscatVec4n <- factor(ascatAneploidy4n[which(rownames(ascatAneploidy4n) %in% goodSamps), ], levels = c("none", "gain", "loss"))

    # goodAscatVec2n <- factor(ascatAneploidy2n, levels = c("none", "gain", "loss"))
    # goodAscatVec3n <- factor(ascatAneploidy3n, levels = c("none", "gain", "loss"))
    # goodAscatVec4n <- factor(ascatAneploidy4n, levels = c("none", "gain", "loss"))
    # 
    goodConMat2n_pr <- confusionMatrix(tmpGoodVec, reference = goodAscatVec2n, mode = "prec_recall")
    goodConMat3n_pr <- confusionMatrix(tmpGoodVec, reference = goodAscatVec3n, mode = "prec_recall")
    goodConMat4n_pr <- confusionMatrix(tmpGoodVec, reference = goodAscatVec4n, mode = "prec_recall")
    
    ### comment out for all samples
    tmp_cna <-  tmp_cna[which(tmp_cna$sample %in% goodSamps), ]
    
    if (length(which(is.na(tmp_cna$seg.mean))) > 0) {
      tmp_cna <- tmp_cna[-which(is.na(tmp_cna$seg.mean)),]
    }
    
    tmp_cna$seg.mean[which(abs(tmp_cna$seg.mean) < 0.2)] <- 0
    
    tmp_cna2 <- NULL
    k <- unique(tmp_cna$sample)[1]
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
    
    if (length(which(tmp_cna2$seg.mean != 0)) == 0) {
     next() 
    }
    
    
    ### ascat portion - need tmpGrange 1 to change if I use different ploidy
    for (q in 2:4) {
      if(q == 2){
        tmpGrange <- tmpGrange1
        tmpDfPloidy <- tmpDf2n
      } else if(q == 3){
        tmpGrange <- tmpGrange2
        tmpDfPloidy <- tmpDf3n
      } else if(q == 4){
        tmpGrange <- tmpGrange3
        tmpDfPloidy <- tmpDf4n
      }
      tmpAscatCna <- NULL
      k <- unique(tmpGrange$sample)[1]
      for (k in unique(tmpGrange$sample)) {
        tmpN <- tmpDfPloidy[which(tmpDfPloidy$sample == k),]
        tmpRes <- tmp_cna2[which(tmp_cna2$seg.mean != 0) ,1:4]
        
        if (grepl("133576rt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "13576rt"), ]
        } else if (grepl("14150rt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "14150lt"), ]
        } else if (grepl("14154rot", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "14154lt"), ]
        } else if (grepl("133576rt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "13576rt"), ]
        } else if (grepl("14656peritonealmt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "14656peritnealmt"), ] 
        } else{
          tmpRes <- tmpRes[which(tmpRes$sample == allSampsCorDf2_bestV2$sampleStripped2[which(allSampsCorDf2_bestV2$string5 == k)]), ]
        }
        
        if (nrow(tmpRes) == 0) {
          next()
        }
        
        queryIdx <- 1:nrow(tmpRes)
        tmpGr <- GRanges(seqnames = tmpN$chr, IRanges(start = tmpN$startpos, end = tmpN$endpos))
        tmpGr2 <- GRanges(seqnames = tmpRes$chrom, IRanges(start = tmpRes$start, end = tmpRes$end))
        tmpIdx <- subjectHits(findOverlaps(query = tmpGr2, subject = tmpGr))
        testQuery <- queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))
        
        if (length(which(duplicated(queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))))) > 0) {
          tmpIdx <-  tmpIdx[-which(duplicated(queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))))]
          testQuery <- testQuery[-which(duplicated(queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))))]
        }
        ### error probably happens b/c there is no match for one of the queries
        tmpRes$cn <- NA
        tmpRes$cn[testQuery] <- tmpN$totalSc[tmpIdx]
        tmpRes2 <- tmpRes[, c("sample", "chrom", "start", "end", "cn")]
        tmpAscatCna <- rbind(tmpAscatCna, tmpRes2)
      }
      
      tmpAscatCna$absCn <- tmpAscatCna$cn - q
      tmpAscatCna$type <- "none"
      tmpAscatCna$type[which(tmpAscatCna$absCn > 0)] <- "gain"
      tmpAscatCna$type[which(tmpAscatCna$absCn < 0)] <- "loss"
      tmpAscatCna$size <- tmpAscatCna$end - tmpAscatCna$start
      
      ### for only good samples
      tmpAscatCna <- tmpAscatCna[which(tmpAscatCna$sample %in% goodSamps), ]
      
      
      tmp_cna3 <- tmp_cna2[which(tmp_cna2$seg.mean != 0),]
      if (nrow(tmp_cna3) == 0) {
        next()
      }
      tmp_cna3 <- tmp_cna3[order(match(tmp_cna3$sample, unique(tmpAscatCna$sample))),]
      tmp_cna3$type <- "none"
      tmp_cna3$type[which(tmp_cna3$seg.mean >  0.2)] <- "gain"
      tmp_cna3$type[which(tmp_cna3$seg.mean <  -0.2)] <- "loss"
      
      tmpLess1Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size < 1e6)], levels = c("none", "gain", "loss")),
                                    reference = factor(tmpAscatCna$type[which(tmpAscatCna$size < 1e6)], levels = c("none", "gain", "loss")),
                                    mode = "prec_recall")
      
      tmpGreater1Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 1e6)], levels = c("none", "gain", "loss")),
                                       reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 1e6)], levels = c("none", "gain", "loss")),
                                       mode = "prec_recall")
      
      tmpGreater5Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 5e6)], levels = c("none", "gain", "loss")),
                                       reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 5e6)], levels = c("none", "gain", "loss")),
                                       mode = "prec_recall")
      
      tmpGreater10Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 10e6)], levels = c("none", "gain", "loss")),
                                        reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 10e6)], levels = c("none", "gain", "loss")),
                                        mode = "prec_recall")
      
      tmpGreater20Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 20e6)], levels = c("none", "gain", "loss")),
                                        reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 20e6)], levels = c("none", "gain", "loss")),
                                        mode = "prec_recall")
      
      tmpCnaStatsDf <- NULL
      
      goodConMat2n_pr
      
      tmpCnaStatsDf <- rbind(tmpCnaStatsDf,
                             data.frame("type" = rep("aneuploidy", 3), "change"  = rownames(eval(parse(text = paste0("goodConMat", q,"n", "_pr")))$byClass),
                                        eval(parse(text = paste0("goodConMat", q, "n", "_pr")))$byClass),
                             data.frame("type" = rep("<1Mb", 3), "change"  = rownames(tmpLess1Mb$byClass), tmpLess1Mb$byClass),
                             data.frame("type" = rep(">1Mb", 3), "change"  = rownames(tmpGreater1Mb$byClass), tmpGreater1Mb$byClass),
                             data.frame("type" = rep(">5Mb", 3), "change"  = rownames(tmpGreater5Mb$byClass), tmpGreater5Mb$byClass),
                             data.frame("type" = rep(">10Mb", 3), "change"  = rownames(tmpGreater10Mb$byClass), tmpGreater10Mb$byClass),
                             data.frame("type" = rep(">20Mb", 3), "change"  = rownames(tmpGreater20Mb$byClass), tmpGreater20Mb$byClass))
      tmpCnaStatsDf$windowSize <- paste0(i, "Mb")
      tmpCnaStatsDf$ploidy <- paste0(q, "n")
      tmpCnaStatsDf$markerFilt <- j
      allCnaStatsDf2 <- rbind(allCnaStatsDf2, tmpCnaStatsDf)
    }
  }
}

allCnaStatsDf2_goodSamps <- allCnaStatsDf2

### doing it with all variables

allCnaStatsDf3 <- allCnaStatsDf2[-which(allCnaStatsDf2$change ==  "Class: none"), ]
allCnaStatsDf3$type <- factor(allCnaStatsDf3$type, levels = c("aneuploidy", ">20Mb", ">10Mb", ">5Mb", ">1Mb", "<1Mb"))
allCnaStatsDf3$markerFilt <- factor(allCnaStatsDf3$markerFilt)
allCnaStatsDf3$windowSize <- factor(allCnaStatsDf3$windowSize, levels = c("1Mb", "5Mb", "10Mb", "20Mb"))
mixedPlot <- ggplot(allCnaStatsDf3, aes(x = type, y = F1)) + geom_boxplot() + 
  geom_point(aes(shape = markerFilt, color = windowSize), position=position_jitterdodge(jitter.width = 4, dodge.width = 0)) +
  facet_wrap(~ploidy, ncol = 3) + 
  scale_color_manual(values = c("1Mb" = "darkred",
                               "5Mb" = "darkblue",
                               "10Mb" = "darkorange",
                               "20Mb" = "purple")) 

### choosing the best performing ones and then comparing below

allCnaStatsDf4 <- allCnaStatsDf2[-which(allCnaStatsDf2$change ==  "Class: none"), ]
allCnaStatsDf4 <- allCnaStatsDf4[which(allCnaStatsDf4$markerFilt == 3 & allCnaStatsDf4$windowSize == "5Mb"), ]
allCnaStatsDf4$type <- factor(allCnaStatsDf4$type, levels = c("aneuploidy", ">20Mb", ">10Mb", ">5Mb", ">1Mb", "<1Mb"))

a3_5_R <- ggplot(allCnaStatsDf4, aes(x = type, y = Recall)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 3 with 5Mb Window Recall") + ylim(0, 1)

a3_5_P <- ggplot(allCnaStatsDf4, aes(x = type, y = Precision)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 3 with 5Mb Window Precision") + ylim(0, 1)


allCnaStatsDf4 <- allCnaStatsDf2[-which(allCnaStatsDf2$change ==  "Class: none"), ]
allCnaStatsDf4 <- allCnaStatsDf4[which(allCnaStatsDf4$markerFilt == 4 & allCnaStatsDf4$windowSize == "20Mb"), ]
allCnaStatsDf4$type <- factor(allCnaStatsDf4$type, levels = c("aneuploidy", ">20Mb", ">10Mb", ">5Mb", ">1Mb", "<1Mb"))

a4_20_R <- ggplot(allCnaStatsDf4, aes(x = type, y = Recall)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 4 with 20Mb Window Recall") + ylim(0, 1)

a4_20_P <- ggplot(allCnaStatsDf4, aes(x = type, y = Precision)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 4 with 20Mb Window Precision") + ylim(0, 1)


allCnaStatsDf4 <- allCnaStatsDf2[-which(allCnaStatsDf2$change ==  "Class: none"), ]
allCnaStatsDf4 <- allCnaStatsDf4[which(allCnaStatsDf4$markerFilt == 3 & allCnaStatsDf4$windowSize == "10Mb"), ]
allCnaStatsDf4$type <- factor(allCnaStatsDf4$type, levels = c("aneuploidy", ">20Mb", ">10Mb", ">5Mb", ">1Mb", "<1Mb"))

a3_10_R <- ggplot(allCnaStatsDf4, aes(x = type, y = Recall)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 3 with 10Mb Window Recall") + ylim(0, 1)

a3_10_P <- ggplot(allCnaStatsDf4, aes(x = type, y = Precision)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 3 with 10Mb Window Precision") + ylim(0, 1)

grid.arrange(a3_5_R, a3_5_P,
             a3_10_R, a3_10_P,
             a4_20_R, a4_20_P,
             nrow = 3)


### also change windowSize in the df to total genome size the windowed filter covers
### look at cases where gains dont match up  might be samples with low tc - unlikely



allCnaStatsDfAll <- NULL
for (i in c(1, 5, 10, 20)) {
  
  for (j in c(0:4)) {
    ### read in tables after they've gone through bedtools
    tmpDisJoin <- read.table(paste0("/mnt/DATA5/tmp/kev/misc/20230606ngsInter_", j, "_", i, "Mb.bed"),
                             sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    
    colnames(tmpDisJoin) <- colnames(ngsDisjoinBed)
    
    tmpDisJoinAneuploidy <- matrix(0, nrow = length(unique(sampleNames2)), ncol = length(unique(segResFilt3$chrom)))
    rownames(tmpDisJoinAneuploidy) <- unique(sampleNames2)
    colnames(tmpDisJoinAneuploidy) <- unique(segResFilt3$chrom)
    
    tmp_cna <- NULL
    k <- rownames(tmpDisJoinAneuploidy)[1]
    for (k in rownames(tmpDisJoinAneuploidy)) {
      tmpDf <- tmpDisJoin[grep(k, tmpDisJoin$sample),]
      
      for (l in seq_along(tmpChrList)) {
        tmpChr <- tmpDf[which(tmpDf$chrom == tmpChrList[l]), ]
        tmpChr2 <- tmpChr
        chrThres <- sum(tmpChr$size) * 0.7
        tmpAneu <- "none"
        if (sum(tmpChr$size[which(tmpChr$seg.mean > 0.2)]) > chrThres) {
          tmpAneu <- "gain"
          tmpChr2$seg.mean[which(tmpChr2$seg.mean > 0.2)] <- 0
        } else if(sum(tmpChr$size[which(tmpChr$seg.mean < -0.2)]) > chrThres){
          tmpAneu <- "loss"
          tmpChr2$seg.mean[which(tmpChr2$seg.mean < -0.2)] <- 0
        } else{
          tmpAneu <- "none"
        }
        tmpDisJoinAneuploidy[k, l] <- tmpAneu
        tmp_cna <- rbind(tmp_cna, tmpChr2)
      }
    }
    
    ### commented out to compare good samples to all
    
    tmpGood <- tmpDisJoinAneuploidy[which(rownames(tmpDisJoinAneuploidy) %in% badSamps), 1:19]
    # tmpGood <- tmpDisJoinAneuploidy[, 1:19]
    tmpGoodVec <- as.vector(tmpGood)
    tmpGoodVec <- factor(tmpGoodVec, levels = c("none", "gain", "loss"))
    # 
    goodAscatVec2n <- factor(ascatAneploidy2n[which(rownames(ascatAneploidy2n) %in% badSamps), ], levels = c("none", "gain", "loss"))
    goodAscatVec3n <- factor(ascatAneploidy3n[which(rownames(ascatAneploidy3n) %in% badSamps), ], levels = c("none", "gain", "loss"))
    goodAscatVec4n <- factor(ascatAneploidy4n[which(rownames(ascatAneploidy4n) %in% badSamps), ], levels = c("none", "gain", "loss"))

    # 
    goodConMat2n_pr <- confusionMatrix(tmpGoodVec, reference = goodAscatVec2n, mode = "prec_recall")
    goodConMat3n_pr <- confusionMatrix(tmpGoodVec, reference = goodAscatVec3n, mode = "prec_recall")
    goodConMat4n_pr <- confusionMatrix(tmpGoodVec, reference = goodAscatVec4n, mode = "prec_recall")
    
    tmp_cna <-  tmp_cna[which(tmp_cna$sample %in% badSamps), ]
    
    if (length(which(is.na(tmp_cna$seg.mean))) > 0) {
      tmp_cna <- tmp_cna[-which(is.na(tmp_cna$seg.mean)),]
    }
    
    tmp_cna$seg.mean[which(abs(tmp_cna$seg.mean) < 0.2)] <- 0
    
    tmp_cna2 <- NULL
    k <- unique(tmp_cna$sample)[1]
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
    
    if (length(which(tmp_cna2$seg.mean != 0)) == 0) {
      next() 
    }
    
    
    ### ascat portion - need tmpGrange 1 to change if I use different ploidy
    for (q in 2:4) {
      if(q == 2){
        tmpGrange <- tmpGrange1
        tmpDfPloidy <- tmpDf2n
      } else if(q == 3){
        tmpGrange <- tmpGrange2
        tmpDfPloidy <- tmpDf3n
      } else if(q == 4){
        tmpGrange <- tmpGrange3
        tmpDfPloidy <- tmpDf4n
      }
      tmpAscatCna <- NULL
      k <- unique(tmpGrange$sample)[1]
      for (k in unique(tmpGrange$sample)) {
        tmpN <- tmpDfPloidy[which(tmpDfPloidy$sample == k),]
        tmpRes <- tmp_cna2[which(tmp_cna2$seg.mean != 0) ,1:4]
        
        if (grepl("133576rt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "13576rt"), ]
        } else if (grepl("14150rt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "14150lt"), ]
        } else if (grepl("14154rot", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "14154lt"), ]
        } else if (grepl("133576rt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "13576rt"), ]
        } else if (grepl("14656peritonealmt", k)) {
          tmpRes <- tmpRes[which(tmpRes$sample == "14656peritnealmt"), ] 
        } else{
          tmpRes <- tmpRes[which(tmpRes$sample == allSampsCorDf2_bestV2$sampleStripped2[which(allSampsCorDf2_bestV2$string5 == k)]), ]
        }
        
        if (nrow(tmpRes) == 0) {
          next()
        }
        
        queryIdx <- 1:nrow(tmpRes)
        tmpGr <- GRanges(seqnames = tmpN$chr, IRanges(start = tmpN$startpos, end = tmpN$endpos))
        tmpGr2 <- GRanges(seqnames = tmpRes$chrom, IRanges(start = tmpRes$start, end = tmpRes$end))
        tmpIdx <- subjectHits(findOverlaps(query = tmpGr2, subject = tmpGr))
        testQuery <- queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))
        
        if (length(which(duplicated(queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))))) > 0) {
          tmpIdx <-  tmpIdx[-which(duplicated(queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))))]
          testQuery <- testQuery[-which(duplicated(queryHits(findOverlaps(query = tmpGr2, subject = tmpGr))))]
        }
        ### error probably happens b/c there is no match for one of the queries
        tmpRes$cn <- NA
        tmpRes$cn[testQuery] <- tmpN$totalSc[tmpIdx]
        tmpRes2 <- tmpRes[, c("sample", "chrom", "start", "end", "cn")]
        tmpAscatCna <- rbind(tmpAscatCna, tmpRes2)
      }
      
      tmpAscatCna$absCn <- tmpAscatCna$cn - q
      tmpAscatCna$type <- "none"
      tmpAscatCna$type[which(tmpAscatCna$absCn > 0)] <- "gain"
      tmpAscatCna$type[which(tmpAscatCna$absCn < 0)] <- "loss"
      tmpAscatCna$size <- tmpAscatCna$end - tmpAscatCna$start
      
      tmpAscatCna <- tmpAscatCna[which(tmpAscatCna$sample %in% badSamps), ]
      
      tmp_cna3 <- tmp_cna2[which(tmp_cna2$seg.mean != 0),]
      if (nrow(tmp_cna3) == 0) {
        next()
      }
      
      tmp_cna3 <- tmp_cna3[order(match(tmp_cna3$sample, unique(tmpAscatCna$sample))),]
      tmp_cna3$type <- "none"
      tmp_cna3$type[which(tmp_cna3$seg.mean >  0.2)] <- "gain"
      tmp_cna3$type[which(tmp_cna3$seg.mean <  -0.2)] <- "loss"
      
      tmpLess1Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size < 1e6)], levels = c("none", "gain", "loss")),
                                    reference = factor(tmpAscatCna$type[which(tmpAscatCna$size < 1e6)], levels = c("none", "gain", "loss")),
                                    mode = "prec_recall")
      
      tmpGreater1Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 1e6)], levels = c("none", "gain", "loss")),
                                       reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 1e6)], levels = c("none", "gain", "loss")),
                                       mode = "prec_recall")
      
      tmpGreater5Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 5e6)], levels = c("none", "gain", "loss")),
                                       reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 5e6)], levels = c("none", "gain", "loss")),
                                       mode = "prec_recall")
      
      tmpGreater10Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 10e6)], levels = c("none", "gain", "loss")),
                                        reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 10e6)], levels = c("none", "gain", "loss")),
                                        mode = "prec_recall")
      
      tmpGreater20Mb <- confusionMatrix(factor(tmp_cna3$type[which(tmpAscatCna$size > 20e6)], levels = c("none", "gain", "loss")),
                                        reference = factor(tmpAscatCna$type[which(tmpAscatCna$size > 20e6)], levels = c("none", "gain", "loss")),
                                        mode = "prec_recall")
      
      tmpCnaStatsDf <- NULL
      
      goodConMat2n_pr
      
      tmpCnaStatsDf <- rbind(tmpCnaStatsDf,
                             data.frame("type" = rep("aneuploidy", 3), "change"  = rownames(eval(parse(text = paste0("goodConMat", q,"n", "_pr")))$byClass),
                                        eval(parse(text = paste0("goodConMat", q, "n", "_pr")))$byClass),
                             data.frame("type" = rep("<1Mb", 3), "change"  = rownames(tmpLess1Mb$byClass), tmpLess1Mb$byClass),
                             data.frame("type" = rep(">1Mb", 3), "change"  = rownames(tmpGreater1Mb$byClass), tmpGreater1Mb$byClass),
                             data.frame("type" = rep(">5Mb", 3), "change"  = rownames(tmpGreater5Mb$byClass), tmpGreater5Mb$byClass),
                             data.frame("type" = rep(">10Mb", 3), "change"  = rownames(tmpGreater10Mb$byClass), tmpGreater10Mb$byClass),
                             data.frame("type" = rep(">20Mb", 3), "change"  = rownames(tmpGreater20Mb$byClass), tmpGreater20Mb$byClass))
      tmpCnaStatsDf$windowSize <- paste0(i, "Mb")
      tmpCnaStatsDf$ploidy <- paste0(q, "n")
      tmpCnaStatsDf$markerFilt <- j
      allCnaStatsDfAll <- rbind(allCnaStatsDfAll, tmpCnaStatsDf)
    }
  }
}




allCnaStatsDfAll2 <- allCnaStatsDfAll[-which(allCnaStatsDfAll$change ==  "Class: none"), ]
allCnaStatsDfAll2$type <- factor(allCnaStatsDfAll2$type, levels = c("aneuploidy", ">20Mb", ">10Mb", ">5Mb", ">1Mb", "<1Mb"))
allCnaStatsDfAll2$markerFilt <- factor(allCnaStatsDfAll2$markerFilt)
allCnaStatsDfAll2$windowSize <- factor(allCnaStatsDfAll2$windowSize, levels = c("1Mb", "5Mb", "10Mb", "20Mb"))
mixedPlot2 <- ggplot(allCnaStatsDfAll2, aes(x = type, y = F1)) + geom_boxplot() + 
  geom_point(aes(shape = markerFilt, color = windowSize), position=position_jitterdodge(jitter.width = 4, dodge.width = 0)) +
  facet_wrap(~ploidy, ncol = 3) + 
  scale_color_manual(values = c("1Mb" = "darkred",
                                "5Mb" = "darkblue",
                                "10Mb" = "darkorange",
                                "20Mb" = "purple")) 

### choosing the best performing ones and then comparing below

allCnaStatsDfAll3 <- allCnaStatsDfAll[-which(allCnaStatsDfAll$change ==  "Class: none"), ]
allCnaStatsDfAll3 <- allCnaStatsDfAll3[which(allCnaStatsDfAll3$markerFilt == 3 & allCnaStatsDfAll3$windowSize == "5Mb"), ]
allCnaStatsDfAll3$type <- factor(allCnaStatsDfAll3$type, levels = c("aneuploidy", ">20Mb", ">10Mb", ">5Mb", ">1Mb", "<1Mb"))

a3_5_R_all <- ggplot(allCnaStatsDfAll3, aes(x = type, y = Recall)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 3 with 5Mb Window Recall") + ylim(0, 1)

a3_5_P_all <- ggplot(allCnaStatsDfAll3, aes(x = type, y = Precision)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 3 with 5Mb Window Precision") + ylim(0, 1)


allCnaStatsDfAll3 <- allCnaStatsDfAll[-which(allCnaStatsDfAll$change ==  "Class: none"), ]
allCnaStatsDfAll3 <- allCnaStatsDfAll3[which(allCnaStatsDfAll3$markerFilt == 4 & allCnaStatsDfAll3$windowSize == "10Mb"), ]
allCnaStatsDfAll3$type <- factor(allCnaStatsDfAll3$type, levels = c("aneuploidy", ">20Mb", ">10Mb", ">5Mb", ">1Mb", "<1Mb"))

a4_10_R_all <- ggplot(allCnaStatsDfAll3, aes(x = type, y = Recall)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 4 with 10Mb Window Recall") + ylim(0, 1)

a4_10_P_all <- ggplot(allCnaStatsDfAll3, aes(x = type, y = Precision)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 4 with 10Mb Window Precision") + ylim(0, 1)

allCnaStatsDfAll3 <- allCnaStatsDfAll[-which(allCnaStatsDfAll$change ==  "Class: none"), ]
allCnaStatsDfAll3 <- allCnaStatsDfAll3[which(allCnaStatsDfAll3$markerFilt == 2 & allCnaStatsDfAll3$windowSize == "5Mb"), ]
allCnaStatsDfAll3$type <- factor(allCnaStatsDfAll3$type, levels = c("aneuploidy", ">20Mb", ">10Mb", ">5Mb", ">1Mb", "<1Mb"))

a2_5_R_all <- ggplot(allCnaStatsDfAll3, aes(x = type, y = Recall)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 2 with 5Mb Window Recall") + ylim(0, 1)

a2_5_P_all <- ggplot(allCnaStatsDfAll3, aes(x = type, y = Precision)) + geom_boxplot(aes(fill = change)) +
  ggtitle("Marker >= 2 with 5Mb Window Precision") + ylim(0, 1)

grid.arrange(a2_5_R_all, a2_5_P_all,
             a3_5_R_all, a3_5_P_all,
             a4_10_R_all, a4_10_P_all, nrow = 3)






### easy to expect that the larger number of markers in a smaller window i.e more general coverage == better result
### however a lot of the times that leads to over filtering of results
### the aneuploidy does seem pretty stable between the the different filters meaning maybe only a general coverage is needed
### it doesn't change at all which is pretty strange ... but oh well maybe much of the genome isn't filtered out
### need to get plots that are stacked for overall genome coverage for the windows





