nameStripper <- function(df){
  require(stringr)
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
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
  
  # uniqIdx <- unique(match(str_replace_all(segRes$ID, "[[:punct:]]", ""), 
  #                         str_replace_all(str_replace_all(colnames(zscore_gc_oe_ratios), "[[:punct:]]", ""), " ", "")))
  
  uniqIdx <- unique(match(segRes$ID, colnames(zscore_gc_oe_ratios)))
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


calcZscore_log2 <- function(segRes){
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
  
  # uniqIdx <- unique(match(str_replace_all(segRes$ID, "[[:punct:]]", ""), 
  #                         str_replace_all(str_replace_all(colnames(zscore_gc_oe_ratios), "[[:punct:]]", ""), " ", "")))
  
  uniqIdx <- unique(match(segRes$ID, colnames(zscore_gc_oe_ratios)))
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
        
        tmpTest <- log2(tmpTest)
        tmpNormal2 <- log2(tmpNormal2)
        
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

whiteListGrange <- readRDS("/mnt/DATA5/tmp/kev/misc/20221212whiteListGRangesMousePanel.rds")

source("/home/kevhu/scripts/20220427clusterSeg.R")
library("PSCBS")
library(copynumber)
mainDir <- c("/mnt/DATA6/mouseData/copynumber/")
listOfReports <- c("Auto_user_AUS5-120-MG_EFD4_BBN_334_304", "Auto_user_AUS5-138-MG_cho_20210621_354_343",
                   "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                   "Auto_user_AUS5-156-MG_Fearon_20210809_374_382", "Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384",
                   "Reanalysis_AUS5-76-MG_test1_217", "Auto_user_AUS5-76-MG_test1_255_185")


allReportSegs <- NULL
i <- listOfReports[2]
for (i in listOfReports) {
  amplicons <- read.table(paste0(mainDir, i, "/cnAmplicon_matrix.txt"),
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  zscore_gc_oe_ratios <- read.table(paste0(mainDir, i, "/gcCorrectedCounts_matrix.txt"),
                                    sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  bed <- read.table(paste0(mainDir, i, "/bed.txt"),
                    sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  normal <- read.table(paste0(mainDir, i, "/normals.txt"),
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  
  mouseAmplicons <- amplicons
  mouseBedFile <- bed
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
  colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)] <- nameStripper(colnames(mouseAmplicons2)[2:ncol(mouseAmplicons2)])
  allMouseCalls <- colnames(mouseAmplicons2)[2:(ncol(mouseAmplicons2))]
  colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)] <- nameStripper(colnames(zscore_gc_oe_ratios)[2:(ncol(zscore_gc_oe_ratios)-10)])
  
  all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                               ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                                end = zscore_gc_oe_ratios$EndPos))
  mouseNormal <- normal$V1
  mouseNormal <- nameStripper(mouseNormal)
  mouseBed2$avgPos <- apply(mouseBed2[,3:4], 1, function(x) sum(x)/2)
  
  ### edited this out for my newer version of filtering based on whitelisting
  # noiseAmps <- mouseAmplicons2[,which(colnames(mouseAmplicons2) %in% mouseNormal)]
  # if (any(grepl("\\.", colnames(noiseAmps)))) {
  #   noiseAmps <- noiseAmps[, -grep("\\.", colnames(noiseAmps))]
  # }
  # noiseAmps  <- 2^noiseAmps
  # quantile(apply(noiseAmps, 1, sd), seq(0,1,0.01))
  # mouseBed2 <- mouseBed2[-which(apply(noiseAmps, 1, sd) > 0.27),]
  # mouseAmplicons2 <- mouseAmplicons2[-which(apply(noiseAmps, 1, sd) > 0.27),]
  # 
  ###
  
  all_probes_grange_chr <- GRanges(seqnames = Rle(paste0("chr", zscore_gc_oe_ratios$ChromNum)),
                               ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                                end = zscore_gc_oe_ratios$EndPos))
  
  goodAmps <- unique(queryHits(findOverlaps(all_probes_grange_chr, whiteListGrange)))
  mouseBed2 <- mouseBed2[goodAmps, ]
  mouseAmplicons2 <- mouseAmplicons2[goodAmps,]
  rownames(mouseBed2) <- NULL
  rownames(mouseAmplicons2) <- NULL
  
  j <- allMouseCalls[1]
  
  cl <- makeCluster(25)
  registerDoParallel(cl)
  allSeg <- foreach(j = allMouseCalls, .combine = "rbind",
                    .packages = c("PSCBS", "copynumber")) %dopar% {
                      source("/home/kevhu/scripts/20220427clusterSeg.R")
                      tmp <- data.frame("chromosome"  = mouseBed2$ChromNum, "x" = mouseBed2$StartPos,
                                        "y" = mouseAmplicons2[[j]])
                      tmpName <- j
                      colnames(tmp) <- c("chromosome", "x", "y")
                      rownames(tmp) <- mouseBed2$AmpliconId
                      tmpGaps <- findLargeGaps(tmp, minLength = 1.5e+07)
                      tmpKnownSegs <- gapsToSegments(tmpGaps)
                      tmpFit <- segmentByCBS(tmp, undo = 1, knownSegments = tmpKnownSegs,
                                             p.method = c("hybrid"))
                      tmpSegRes <- tmpFit$output
                      tmpSegRes$sampleName <- tmpName
                      tmpSegRes <- tmpSegRes[-which(is.na(tmpSegRes$mean)),]
                      tmpSegRes2 <- clustSeg(tmpSegRes, distVar = 0.2)
                      tmpSegRes2
                    }
  stopCluster(cl)
  colnames(allSeg) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
  tmpZScore <- calcZscore(allSeg)
  segZfilt <- segZscoresFilt(allSeg, tmpZScore)
  segZfilt$ID <- nameStripper(segZfilt$ID)
  segZfilt$ID <- str_remove(segZfilt$ID, "x.*")
  
  write.table(segZfilt, file = paste0(mainDir, i,"/pbsSegResPrune20_whitelist.txt"), sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  allReportSegs <- rbind(allReportSegs, segZfilt)
}

# write.table(allReportSegs, paste0("/mnt/DATA5/tmp/kev/misc/2022041915MbGapAllSegPrune20.txt"),
#             sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


### running convex pipeline
###
###

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




load("/mnt/DATA5/tmp/kev/misc/20220628badIdx.RData")
source("/home/kevhu/scripts/20220427clusterSeg.R")
tableRes <- NULL
# file.remove("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genes/absoluteRes.txt")
z <- listOfDirectories
for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  
  segs    <- fread("pbsSegResPrune20_whitelist.txt")
  
  lrs     <- fread("mouseAmplicons.txt")
  lrs <- lrs[goodAmps,]
  
  probloc <- fread("mouseProbeLoc.txt")
  probloc <- probloc[goodAmps,]
  
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
  registerDoParallel(7)
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
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genes_whitelist/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  stopImplicitCluster()
  
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}

### need to redo snp plots but with filtered out regions ... 

write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/ploidyPredictionHclust15AllV3_whitelist_genes.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


