# nameStripper <- function(df){
#   require(stringr)
#   df <- str_remove(df, "_MG_X.*")
#   df <- str_remove(str_remove(df, "^X"), "_X.*")
#   df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
#   df  <- str_remove(df, "o")
#   df  <- str_replace_all(df, " ", "")
#   return(df)
# }

#
# sampleMap <- read.table("/mnt/DATA5/tmp/kev/misc/Sample_Map.txt", header = TRUE,
#                         sep = "\t", stringsAsFactors = FALSE)
# sampleMap$stripped <- nameStripper(sampleMap$ID)
# # sampleMap$stripped[21] <- "13085lt"
# sampleMap$stripped[32] <- "14154lt"
# sampleMap$stripped[25] <- "13576rt"
# sampleMap$stripped[31] <- "14150ltt"
#
# tcDf <- read.table("/mnt/DATA5/tmp/kev/misc/20220321hgscTc.txt", sep = "\t",
#                    header = TRUE, stringsAsFactors = FALSE)
#
#
# tcDf2 <- tcDf[which(tcDf$sample %in% sampleMap$stripped),]
#
# write.table(tcDf2, "/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
#


### 20220521: getting all samples with correct Tc

# mouseNormal <- c("EF_D03_MG_X14", "MG_18X50", "MG_21X53", "MG_6X38",
#                  "MG_8X40","MG_11X43", "MG_13X45")
# 
# nameStripper <- function(df){
#   require(stringr)
#   df <- str_remove(df, "_MG_X.*")
#   df <- str_remove(str_remove(df, "^X"), "_X.*")
#   df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
#   df  <- str_remove(df, "o")
#   df  <- str_replace_all(df, " ", "")
#   return(df)
# }
# 
# cnCalls_1 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-138-MG_cho_20210621_354_343/cnMatrix_gene.txt",
#                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
# cnCalls_2 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/cnMatrix_gene.txt",
#                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
# cnCalls_3 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-142-MG_cho_20210701_357_353/cnMatrix_gene.txt",
#                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
# cnCalls_4 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-76-MG_test1_255_185/cnMatrix_gene.txt",
#                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
# cnCalls_5 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/cnMatrix_gene.txt",
#                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
# cnCalls_6 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384/cnMatrix_gene.txt",
#                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
# cnCalls_7 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/cnMatrix_gene.txt",
#                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
# 
# 
# all_cnCalls <- cbind(cnCalls_1[grep("Del", cnCalls_1$Gene),], cnCalls_2[grep("Del", cnCalls_2$Gene),2:ncol(cnCalls_2)],
#                      cnCalls_3[grep("Del", cnCalls_3$Gene),2:ncol(cnCalls_3)], cnCalls_4[grep("Del", cnCalls_4$Gene),2:ncol(cnCalls_4)],
#                      cnCalls_5[grep("Del", cnCalls_5$Gene),2:ncol(cnCalls_5)], cnCalls_6[grep("Del", cnCalls_6$Gene),2:ncol(cnCalls_6)],
#                      cnCalls_7[grep("Del", cnCalls_7$Gene),2:ncol(cnCalls_7)])
# all_cnCalls <- all_cnCalls[,-which(duplicated(colnames(all_cnCalls)))]
# 
# 
# tc <- 1 - apply(all_cnCalls[grep("Trp53Del", all_cnCalls$Gene), 2:ncol(all_cnCalls)], 2, min)
# tc[which(names(tc) %in% mouseNormal)] <- 1
# names(tc) <- nameStripper(names(tc))
# names(tc) <- str_remove(names(tc), "x.*")
# tcDf <- data.frame("sample" = names(tc), "tc" = tc, stringsAsFactors = FALSE)
# tcHgscDfNonZero <- tcDf
# 
# tcDf$tc[which(tcDf$tc < 0.40)] <- .40
# 
# tcApc <- 1 - all_cnCalls[grep("ApcDel", all_cnCalls$Gene), 2:ncol(all_cnCalls)]
# tcApc[which(names(tcApc) %in% mouseNormal)] <- 1
# tcApcZero  <- tcApc
# tcApc[tcApc < 0.40] <- .40
# names(tcApc) <- nameStripper(names(tcApc))
# names(tcApc) <- str_remove(names(tcApc), "x.*")
# tcApcDf <- data.frame("sample" = names(tcApc), "tc" = unlist(tcApc), stringsAsFactors = FALSE)
# 
# names(tcApcZero) <- nameStripper(names(tcApcZero))
# names(tcApcZero) <- str_remove(names(tcApcZero), "x.*")
# tcApcZeroDf <- data.frame("sample" = names(tcApcZero), "tc" = unlist(tcApcZero), stringsAsFactors = FALSE, check.names = TRUE)
# tcApcZeroDf <- tcApcZeroDf[grep("efd", tcApcZeroDf$sample),]
# tcDf$tc[grep("efd", tcApcDf$sample)] <- tcApcDf$tc[grep("efd", tcApcDf$sample)]
# 
# tcHgscDfNonZero <- tcHgscDfNonZero[-grep("efd", tcHgscDfNonZero$sample), ]
# 
# tcEfdDfNonZero <- tcHgscDfNonZero[-grep("efd", tcHgscDfNonZero$sample), ]


# write.table(tcApcZeroDf, "/mnt/DATA5/tmp/kev/misc/20220720allSampsApc.txt", sep = "\t",
#             col.names = TRUE, row.names = FALSE)

# write.table(tcDf, "/mnt/DATA5/tmp/kev/misc/20220718allSampsDf40.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# 
# write.table(tcApcZeroDf, "/mnt/DATA5/tmp/kev/misc/20220718onlyEfdApcTc50.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# 
# write.table(tcHgscDfNonZero, "/mnt/DATA5/tmp/kev/misc/20220718onlyHgscNonZeroTc40.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)

library(data.table)
library(doParallel)
library(foreach)
library(GenomicRanges)
library(gtools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(raster)
library(ggpubr)
library(jointseg)
library(stringr)

load("/mnt/DATA5/tmp/kev/misc/20220224noshadFunctionsLoaded.RData")



extraCols <- c("AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene", "NumGC",
               "Length", "GC", "TotalPool", "Weights")

# genes = read_rds("~/Codes/TPO/tpo-sane1/genes.rds")
gobj <- getGobj("mm10", NULL, FALSE)

# mouseNormal <- c("MG_17X49", "MG_18X50", "MG_23X55", "MG_6X38",
#                  "MG_8X40","MG_11X43", "MG_13X45")

mouseNormal <- c("EF_D03_MG_X14", "MG_18X50", "MG_21X53", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

eros <- c("/mnt/DATA6/mouseData/copynumber/")
listOfDirectories <- c("Auto_user_AUS5-138-MG_cho_20210621_354_343", "Auto_user_AUS5-76-MG_test1_255_185",
                       "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                       "Auto_user_AUS5-120-MG_EFD4_BBN_334_304", "Auto_user_AUS5-156-MG_Fearon_20210809_374_382",
                       "Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384")

listOfDirectories <- paste0(eros, listOfDirectories)


### 20220224: create tc of only releveant cell line derived samples to look at
### 14806, 12167, 12402, 14908 and 3867 (Cho) & 19821, 19396, 19904, 20390 (need efd numbers) - original

### 14806rt, 12167rt, 12167mt, 12402rt, 12402lt, 14908lt 14908rt and 3867lt (Cho) &
### 19821, 19396, efd47-51, efd8-11 (need efd numbers) - sample names

# tmpTable <- fread("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt")
# tmpTable <- tmpTable[which(tmpTable$sample %in% c("14806rt","12167rt", "12167mt",
#                                                   "12402rt", "12402lt", "14908lt",
#                                                   "14908rt", "3867lt"))]
# 
# efNames <- c("efd08", "efd09", "efd10", "efd11", "efd47", "efd48", "efd49", "efd51")
# efTc <- c(.95, .95, .95, .95, .43, .37, .65, .32)
# tmpTable2 <- cbind(efNames, efTc)
# colnames(tmpTable2) <- colnames(tmpTable)
# tmpTable3 <- rbind(tmpTable, tmpTable2)
# tmpTable3$tc <- as.numeric(tmpTable3$tc)
# 
# write.table(tmpTable3, "/mnt/DATA5/tmp/kev/misc/20220224cellLinesSubSet.txt", quote = FALSE,
#             sep = "\t", row.names = FALSE, col.names = TRUE)
# 


tableRes <- NULL

for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  # segs    <- fread("segResults.txt")
  segs    <- fread("pbsSegRes.txt")
  lrs     <- fread("mouseAmplicons.txt")
  
  probloc <- fread("mouseProbeLoc.txt")
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  Pp$sex  <- "female"
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  #segs$seg.mean[segs$seg.mean > 3] <- 0
  # segs$seg.mean[segs$seg.mean < log2(5/4) & segs$seg.mean > log2(3/4)] <- 0
  # segs$seg.mean[abs(segs$seg.mean) < 0.2] <- 0
  # segs$seg.mean[segs$length < 5e3] <-0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(5)
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnp/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  # print(i)
  # print(head(res))
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}

write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/20220420snpPloidyPrediction.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


### create filtered version of pbsSeg

# tmpTc <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
# tmpTc2 <- data.frame("sample" = c("mg4", "14150lt", "mg20", "14109rt", "mg15"),
#                      "tc" = c(0.96, 0.89, 0.84, 0.69, 0.78))
# tmpTc <- rbind(tmpTc, tmpTc2)
# write.table(tmpTc, file = "/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt", sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)



tableRes <- NULL
for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  # segs    <- fread("segResults.txt")
  segs    <- fread("pbsSegRes.txt")
  lrs     <- fread("mouseAmplicons.txt")
  
  probloc <- fread("mouseProbeLoc.txt")
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  Pp$sex  <- "female"
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  segs$seg.mean[which(abs(segs$seg.mean) < 0.2)] <- 0
  # segs$seg.mean[segs$seg.mean < log2(5/4) & segs$seg.mean > log2(3/4)] <- 0
  segs$seg.mean[segs$length < 5e3] <- 0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(5)
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFilt/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  # print(i)
  # print(head(res))
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}


write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPrediction.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


### running the prediction on the different undo.sd 1.25



tableRes <- NULL

for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  # segs    <- fread("segResults.txt")
  segs    <- fread("pbsSegRes_1.25.txt")
  lrs     <- fread("mouseAmplicons.txt")
  
  probloc <- fread("mouseProbeLoc.txt")
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  Pp$sex  <- "female"
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  segs$seg.mean[which(abs(segs$seg.mean) < 0.2)] <- 0
  # segs$seg.mean[segs$seg.mean < log2(5/4) & segs$seg.mean > log2(3/4)] <- 0
  segs$seg.mean[segs$length < 5e3] <- 0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(5)
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFilt1.25/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  # print(i)
  # print(head(res))
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}


write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPrediction1.25.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


### undo 1.50


tableRes <- NULL

for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  # segs    <- fread("segResults.txt")
  segs    <- fread("pbsSegRes_1.5.txt")
  lrs     <- fread("mouseAmplicons.txt")
  
  probloc <- fread("mouseProbeLoc.txt")
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  Pp$sex  <- "female"
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  segs$seg.mean[which(abs(segs$seg.mean) < 0.2)] <- 0
  # segs$seg.mean[segs$seg.mean < log2(5/4) & segs$seg.mean > log2(3/4)] <- 0
  segs$seg.mean[segs$length < 5e3] <- 0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(5)
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFilt1.50/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  # print(i)
  # print(head(res))
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}


write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPrediction1.5.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


### sd 2


tableRes <- NULL

for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  # segs    <- fread("segResults.txt")
  segs    <- fread("pbsSegRes_2.txt")
  lrs     <- fread("mouseAmplicons.txt")
  
  probloc <- fread("mouseProbeLoc.txt")
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  Pp$sex  <- "female"
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  segs$seg.mean[which(abs(segs$seg.mean) < 0.2)] <- 0
  # segs$seg.mean[segs$seg.mean < log2(5/4) & segs$seg.mean > log2(3/4)] <- 0
  segs$seg.mean[segs$length < 5e3] <- 0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(5)
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFilt2/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  # print(i)
  # print(head(res))
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}


write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPrediction2.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



### sd 3




tableRes <- NULL

for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  # segs    <- fread("segResults.txt")
  segs    <- fread("pbsSegRes_3txt")
  lrs     <- fread("mouseAmplicons.txt")
  
  probloc <- fread("mouseProbeLoc.txt")
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  Pp$sex  <- "female"
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  segs$seg.mean[which(abs(segs$seg.mean) < 0.2)] <- 0
  # segs$seg.mean[segs$seg.mean < log2(5/4) & segs$seg.mean > log2(3/4)] <- 0
  segs$seg.mean[segs$length < 5e3] <- 0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(5)
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFilt3/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  # print(i)
  # print(head(res))
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}


write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPrediction3.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



### segmenting with log2 ratio z scores



tableRes <- NULL

for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  # segs    <- fread("segResults.txt")
  segs    <- fread("pbsSegResLogZ.txt")
  lrs     <- fread("mouseAmplicons.txt")
  
  probloc <- fread("mouseProbeLoc.txt")
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  Pp$sex  <- "female"
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  segs$seg.mean[which(abs(segs$seg.mean) < 0.2)] <- 0
  # segs$seg.mean[segs$seg.mean < log2(5/4) & segs$seg.mean > log2(3/4)] <- 0
  segs$seg.mean[segs$length < 5e3] <- 0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(5)
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFiltLogZ/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  # print(i)
  # print(head(res))
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}


write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPredictionLogZ.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# hclust 

source("/home/kevhu/scripts/20220427clusterSeg.R")
tableRes <- NULL
file.remove("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15/absoluteRes.txt")
for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  # segs    <- fread("segResults.txt")
  segs    <- fread("pbsSegResPrune15.txt")
  lrs     <- fread("mouseAmplicons.txt")
  
  probloc <- fread("mouseProbeLoc.txt")
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  # tmpPp <- Pp$tc[which(Pp$sample  == "15774lt")]
  # Pp$tc[which(Pp$sample  == "15774lt")] <- Pp$tc[which(Pp$sample  == "15774rt")]
  # Pp$tc[which(Pp$sample  == "15774rt")] <- tmpPp
  Pp$sex  <- "female"
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  # segs$seg.mean <- log2(segs$seg.mean)
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  segs$seg.mean[which(abs(segs$seg.mean) < 0.2)] <- 0
  #  longest on panel is apc at 
  segs$seg.mean[segs$length < 5e3] <- 0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(5)
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  # print(i)
  # print(head(res))
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}


write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPredictionHclust15.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

### try to apply gaps to lrs data


source("/home/kevhu/scripts/20220427clusterSeg.R")
tableRes <- NULL
file.remove("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15All/absoluteRes.txt")
for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  # segs    <- fread("segResults.txt")
  segs    <- fread("pbsSegResPrune15.txt")
  lrs     <- fread("mouseAmplicons.txt")
  
  probloc <- fread("mouseProbeLoc.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt")
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  Pp$sex  <- "female"
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  segs$seg.mean[which(abs(segs$seg.mean) < 0.2)] <- 0
  #  longest on panel is apc at 
  segs$seg.mean[segs$length < 5e3] <- 0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(5)
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15All/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}


write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPredictionHclust15All.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



### all samps correct tc for fearon samples


# eros <- c("/mnt/DATA6/mouseData/copynumber/")
# listOfDirectories <- c("Auto_user_AUS5-156-MG_Fearon_20210809_374_382", "Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384")
# listOfDirectories <- paste0(eros, listOfDirectories)
# 


source("/home/kevhu/scripts/20220427clusterSeg.R")
tableRes <- NULL
file.remove("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV2/absoluteRes.txt")
for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  # segs    <- fread("segResults.txt")
  segs    <- fread("pbsSegResPrune15.txt")
  lrs     <- fread("mouseAmplicons.txt")
  
  probloc <- fread("mouseProbeLoc.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220521allSampsDf.txt")
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  Pp$sex  <- "female"
  Pp$sex[grep("efd", Pp$sample)] <- "male"
  
  Pp$tc[which(Pp$tc > 0.95)] <- 0.95
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  if(length(which(duplicated(colnames(lrs)))) > 0){
    lrs <- data.frame(lrs, stringsAsFactors = FALSE)
    lrs <- data.table(lrs)
  }
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  
  if (z == "/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384") {
    colnames(lrs)[2:ncol(lrs)] <- str_remove(colnames(lrs)[2:ncol(lrs)], "x.*")
  }
  

  
  if (length(which(duplicated(colnames(lrs)))) > 0) {
    lrs <- data.frame(lrs, stringsAsFactors = FALSE)
    lrs <- lrs[, -grep("\\.1", colnames(lrs))]
    lrs <- data.table(lrs)
  }
  
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  segs$seg.mean[which(abs(segs$seg.mean) < 0.2)] <- 0
  #  longest on panel is apc at 
  segs$seg.mean[segs$length < 5e3] <- 0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(5)
  i <- 26
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV2/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}


write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/ploidyPredictionHclust15All.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


###
###
### making the gene filter a bit smaller ...


load("/mnt/DATA5/tmp/kev/misc/20220628badIdx.RData")
listOfDirectories  <- listOfDirectories[5:7]
source("/home/kevhu/scripts/20220427clusterSeg.R")
tableRes <- NULL
# file.remove("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genes/absoluteRes.txt")
z <- listOfDirectories[1]
for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  
  segs    <- fread("pbsSegResPrune20.txt")
  
  lrs     <- fread("mouseAmplicons.txt")
  lrs <- lrs[-badAmps,]
  
  probloc <- fread("mouseProbeLoc.txt")
  probloc <- probloc[-badAmps,]
  
  probloc <- fread("mouseProbeLoc.txt")
  ### 20220718 new tc df has 40% cutoff instead of 50%
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220521allSampsDf.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220718allSampsDf40.txt")
  Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  Pp$sex  <- "female"
  Pp$sex[grep("efd", Pp$sample)] <- "male"
  
  Pp$tc[which(Pp$tc > 0.95)] <- 0.95
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  if(length(which(duplicated(colnames(lrs)))) > 0){
    lrs <- data.frame(lrs, stringsAsFactors = FALSE)
    lrs <- data.table(lrs)
  }
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  
  if (z == "/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384") {
    colnames(lrs)[2:ncol(lrs)] <- str_remove(colnames(lrs)[2:ncol(lrs)], "x.*")
  }
  
  
  
  if (length(which(duplicated(colnames(lrs)))) > 0) {
    lrs <- data.frame(lrs, stringsAsFactors = FALSE)
    lrs <- lrs[, -grep("\\.1", colnames(lrs))]
    lrs <- data.table(lrs)
  }
  
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  segs$seg.mean[which(abs(segs$seg.mean) < 0.2)] <- 0
  #  longest on panel is apc at 
  # segs$seg.mean[segs$length < 5e3] <- 0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(15)
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genes/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}


write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/ploidyPredictionHclust15AllV3_genes.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



### no cnr filter but qval Filter
load("/mnt/DATA5/tmp/kev/misc/20220628badIdx.RData")
# listOfDirectories <- listOfDirectories[6]
source("/home/kevhu/scripts/20220427clusterSeg.R")
tableRes <- NULL
file.remove("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genes_qval/absoluteRes.txt")
z <- listOfDirectories[1]
for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  # segs    <- fread("segResults.txt")
  
  segs    <- fread("pbsSegResPrune15.txt")
  
  lrs     <- fread("mouseAmplicons.txt")
  lrs <- lrs[-badAmps,]
  
  probloc <- fread("mouseProbeLoc.txt")
  probloc <- probloc[-badAmps,]
  
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220521allSampsDf.txt")
  Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220420snpPanelTcDf.txt")
  Pp$sex  <- "female"
  Pp$sex[grep("efd", Pp$sample)] <- "male"
  
  Pp$tc[which(Pp$tc > 0.95)] <- 0.95
  
  mouseNormal <- nameStripper(mouseNormal)
  # Pp <- Pp[-which(duplicated(Pp$sample)),]
  # Pp <- Pp[-which(Pp$sample %in% mouseNormal),]
  
  if(length(which(duplicated(colnames(lrs)))) > 0){
    lrs <- data.frame(lrs, stringsAsFactors = FALSE)
    lrs <- data.table(lrs)
  }
  
  lrs <- lrs[, -which(colnames(lrs) %in% extraCols), with = FALSE]
  lrs[,2:ncol(lrs)] <- log2(lrs[,2:ncol(lrs)])
  colnames(lrs)[2:ncol(lrs)] <- nameStripper(colnames(lrs)[2:ncol(lrs)])
  
  if (z == "/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384") {
    colnames(lrs)[2:ncol(lrs)] <- str_remove(colnames(lrs)[2:ncol(lrs)], "x.*")
  }
  
  
  
  if (length(which(duplicated(colnames(lrs)))) > 0) {
    lrs <- data.frame(lrs, stringsAsFactors = FALSE)
    lrs <- lrs[, -grep("\\.1", colnames(lrs))]
    lrs <- data.table(lrs)
  }
  
  namesInter <- intersect(colnames(lrs)[2:ncol(lrs)], Pp$sample)
  
  if (isEmpty(namesInter)) {
    next()
  }
  
  
  Pp <- Pp[match(namesInter, Pp$sample), ]
  if (any(is.na(Pp$sample))) {
    Pp <- Pp[-which(is.na(Pp$sample)), ]
  }
  
  lrs <- lrs[, c(1, match(namesInter, colnames(lrs))), with = FALSE]
  segs$ID <- nameStripper(segs$ID)
  segs <- segs[-which(segs$ID %in% mouseNormal),]
  segs$length <- segs$loc.end - segs$loc.start
  segs$seg.mean[which(abs(segs$`q-val2`) > 0.05)] <- 0
  #  longest on panel is apc at 
  # segs$seg.mean[segs$length < 5e3] <- 0
  
  ## melt lr to make format consistent
  lrs <- melt.data.table(lrs, id.vars = "AmpliconId")
  
  ## add prob loc to lrs
  probloc <- probloc[, AmpliconId := sprintf("AMP_%d", AmpliconId)]
  lrs     <- merge.data.table(lrs, probloc, by = "AmpliconId", all.x = TRUE)
  
  
  ## sort locations
  lrs <- lrs[order(nchar(lrs$ChromNum), lrs$ChromNum)]
  probloc <- lrs[order(nchar(probloc$ChromNum), probloc$ChromNum)]
  
  ## convert sex chromosome to X
  segs[, chrom := sprintf("chr%s", chrom)]
  segs[chrom == "chr23", chrom := "chrX"]
  segs[chrom == "chr20", chrom := "chrX"]
  
  probloc[, ChromNum := sprintf("chr%s", ChromNum)]
  probloc[ChromNum == "chr23", ChromNum := "chrX"]
  probloc[ChromNum == "chr20", ChromNum := "chrX"]
  
  lrs[ , ChromNum := sprintf("chr%s", ChromNum)]
  lrs[ChromNum == "chr23", ChromNum := "chrX"]
  lrs[ChromNum == "chr20", ChromNum := "chrX"]
  
  ### getting lrs lr values to match filtered seg
  
  
  lrs2 <- NULL
  for (i in unique(segs$ID)) {
    tmpSeg <- segs[which(segs$ID == i),]
    tmpLrs <- lrs[which(lrs$variable == i),]
    lrsSubject <- GRanges(seqnames = tmpLrs$ChromNum,
                          IRanges(start = tmpLrs$StartPos, end = tmpLrs$EndPos))
    
    for (j in 1:nrow(tmpSeg)) {
      queryGrange <- GRanges(seqnames = tmpSeg$chrom[j],
                             IRanges(start = tmpSeg$loc.start[j], end = tmpSeg$loc.end[j]))
      tmpLrs$value[subjectHits(findOverlaps(queryGrange, lrsSubject))] <- tmpSeg$seg.mean[j]
    }
    lrs2 <- rbind(lrs2, tmpLrs)
  }
  
  lrs <- lrs2
  
  
  
  ## making CNVEX format CNV object and then run the analysis
  samples <- unique(Pp$sample)
  registerDoParallel(5)
  i <- 26
  res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
    print(i)
    lr     <- lrs[variable == samples[i], ]
    
    if (isEmpty(lr$variable)) {
      res <- NULL
      return(res)
    }
    
    seg    <- segs[ID == samples[i]]
    purity <- Pp[sample == samples[i], tc]
    sex    <- Pp$sex[i]
    ## make tiled genome object
    cnv      <- list()
    cnv$var  <- GRanges()
    cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
    seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                       lr.mean = seg$seg.mean, nlr = seg$num.mark)
    seqinfo(seg) <- gobj$seqi
    ## model cnv
    mcnv             <- cnv
    mcnv$tile$arm    <- seqnames(cnv$tile)
    mcnv$tile$target <- TRUE
    mcnv$tile$hq     <- TRUE
    mcnv$tile$gap    <- 0
    mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
    ## normal copy number
    if( sex == "male") {
      mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
    } else {
      mcnv$tile$nC  <- 2
    }
    mcnv$tile$baf       <- NA_real_
    mcnv$tile$baf.depth <- NA_real_
    mcnv$tile$baf.n     <- NA_real_
    
    data <- .opt.data(mcnv, seg, opts)
    opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                    error = function(e) NULL)
    
    if (is.null(opt)) {
      res <- NULL
      return(res)
    }
    
    # this is where the final likelihoods are so after sort, append to df
    opt <- opt[order(-L)]
    opt$Lrank <- 1:nrow(opt)
    
    lapply(1:nrow(opt), function(j) {
      purity <- opt$p0[j]
      ploidy <- opt$P0[j]
      LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genes_qval/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}


write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/ploidyPredictionHclust15AllV3_qval_genes.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




ploidyTableRaw <- read.table("/mnt/DATA5/tmp/kev/misc/20220420snpPloidyPrediction.txt", 
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)


ploidyTableFilt <- read.table("/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPrediction.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)


ploidyTableFilt2 <- read.table("/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPrediction2.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)



ploidyTableFilt3 <- read.table("/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPrediction3.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)

ploidyTableFiltAns <- read.table("/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPredictionLogZAns.txt",
                               sep = "\t", stringsAsFactors = FALSE, header = TRUE)


ploidyExpected <- read.table("/mnt/DATA5/tmp/kev/misc/20220328ploidyChart.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)

ploidyExpected <- ploidyExpected[-which(ploidyExpected$change == "hom loss"),]

ploidyTableFilt <- ploidyTableFilt[-which(ploidyTableFilt$P0 > 5),]

ploidyFilt15 <- read.table("/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPredictionHclust15.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)

ploidyFilt25 <- read.table("/mnt/DATA5/tmp/kev/misc/20220420snpFiltPloidyPredictionHclust25.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)

