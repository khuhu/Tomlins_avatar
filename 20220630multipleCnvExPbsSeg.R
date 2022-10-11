# original script at /home/kevhu/scripts/20220420noshadMultipleSnp.R

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


gobj <- getGobj("mm10", NULL, FALSE)


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
file.remove("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genes/absoluteRes.txt")
z <- listOfDirectories[1]
for (z in listOfDirectories[1]) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  
  segs    <- fread("pbsSegResPrune20.txt")
  
  lrs     <- fread("mouseAmplicons.txt")
  lrs <- lrs[-badAmps,]
  
  probloc <- fread("mouseProbeLoc.txt")
  probloc <- probloc[-badAmps,]
  
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


### testing global filter

tableRes <- NULL
z <- listOfDirectories[1]
for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  
  segs    <- fread("pbsSegResPrune20G.txt")
  
  lrs     <- fread("mouseAmplicons.txt")
  lrs <- lrs[-badAmps,]
  
  probloc <- fread("mouseProbeLoc.txt")
  probloc <- probloc[-badAmps,]
  
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
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genesG/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}



write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/ploidyPredictionHclust15AllV3_genesG.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




### testing global and chromosomal filter

tableRes <- NULL
z <- listOfDirectories[1]
for (z in listOfDirectories) {
  print(z)
  opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
  setwd(z)
  
  segs    <- fread("pbsSegResPrune20G_2step.txt")
  
  lrs     <- fread("mouseAmplicons.txt")
  lrs <- lrs[-badAmps,]
  
  probloc <- fread("mouseProbeLoc.txt")
  probloc <- probloc[-badAmps,]
  
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
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV3_genesG_2step/")
      id     <- samples[i]
      plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
    })
    
    tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
    return(tmpOpt)
  }
  
  tableRes <- rbind(tableRes, data.frame(res, stringsAsFactors = FALSE))
}



write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/ploidyPredictionHclust15AllV3_genesG_step.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

