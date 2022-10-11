# load("/mnt/DATA5/tmp/kev/misc/20210831noshadFunctionsLoaded.RData")
# 
# files.sources = list.files("/mnt/DATA5/tmp/kev/programs/cnvex/R/");
# files.sources <- paste0("/mnt/DATA5/tmp/kev/programs/cnvex/R/", files.sources);
# sapply(files.sources, source);
# opts <- readRDS("/mnt/DATA5/tmp/kev/misc/opts.rds");
# library(data.table)
# library(doParallel)
# library(foreach)
# library(GenomicRanges)
# library(gtools)
# library(BSgenome.Mmusculus.UCSC.mm10)
# library(raster)
# library(ggpubr)
# library(jointseg)
# 
# extraCols <- c("AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene", "NumGC",
#                "Length", "GC", "TotalPool", "Weights")
# 
# # genes = read_rds("~/Codes/TPO/tpo-sane1/genes.rds")
# gobj <- getGobj("mm10", NULL, FALSE)

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

mouseNormal <- c("MG_17X49", "MG_18X50", "MG_23X55", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

eros <- c("/mnt/DATA6/mouseData/copynumber/")
listOfDirectories <- c("Auto_user_AUS5-138-MG_cho_20210621_354_343", "Auto_user_AUS5-76-MG_test1_255_185",
                       "Auto_user_AUS5-141-MG_cho_202106_3TS_356_349", "Auto_user_AUS5-142-MG_cho_20210701_357_353",
                       "Auto_user_AUS5-120-MG_EFD4_BBN_334_304")

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
  segs    <- fread("segResults.txt")
  lrs     <- fread("mouseAmplicons.txt")
  
  probloc <- fread("mouseProbeLoc.txt")
  # Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt")
  Pp      <- fread("/mnt/DATA5/tmp/kev/misc/20220224cellLinesSubSet.txt")
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
  segs$seg.mean[abs(segs$seg.mean) < 0.2] <- 0
  segs$seg.mean[segs$length < 5e3] <-0
  
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
      w.fn   <- sprintf("/mnt/DATA5/tmp/kev/testNoshad4/")
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

write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/20220224cellLinePloidyPrediction.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# write.table(tableRes, "/mnt/DATA5/tmp/kev/misc/20210917hgscPloidyPredictionV3.txt",
#             sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

tableRes <- read.table("/mnt/DATA5/tmp/kev/misc/20210917hgscPloidyPredictionV3.txt", sep = "\t",
                       header = TRUE, stringsAsFactors = FALSE)

tableOrig <- read.table("/mnt/DATA5/tmp/kev/misc/20210913hgscPloidyPrediction.txt", sep = "\t",
                        header = TRUE, stringsAsFactors = FALSE)

tableTest2 <- read.table("/mnt/DATA5/tmp/kev/misc/20210916hgscPloidyPrediction5kb.txt", sep = "\t",
                         header = TRUE, stringsAsFactors = FALSE)

ppTblRes <- Pp[which(Pp$sample %in% unique(tableRes$Sample)),]

tmpTable <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20210914ploidyAssesmentsHgovc.xlsx")
tmpTable2 <- tmpTable[which(tmpTable$Sample %in% ppTblRes$sample[which(ppTblRes$tc > 0.50)]),]

a <- ggplot(tmpTable, aes(x = IntegerP)) + geom_density() + xlab("ploidy") + ggtitle("HGSC(GEMM) Plodies n = 160") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0,1))

b <- ggplot(tmpTable2, aes(x = IntegerP)) + geom_density() + xlab("ploidy") + ggtitle("HGSC(GEMM) Plodies n = 136; > 50% purity") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(c(0,1))

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210914hgscGemmPloidyDensity.pdf")
grid.arrange(a,b)
dev.off()


c <- ggplot(tmpTable, aes(x = P0)) + geom_histogram(binwidth = 0.1, fill = "white", color = "black")  + 
  xlab("ploidy") + ggtitle("HGSC(GEMM) Plodies n = 169") +
  theme(plot.title = element_text(hjust = 0.5)) 

d <- ggplot(tmpTable2, aes(x = P0)) + geom_histogram(binwidth = 0.1, fill = "white", color = "black")  +
  xlab("ploidy") + ggtitle("HGSC(GEMM) Plodies n = 136; > 50% purity") +
  theme(plot.title = element_text(hjust = 0.5)) 

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210914hgscGemmPloidyDensityNonInt.pdf")
grid.arrange(c,d)
dev.off()

tableRes <- read.table("/mnt/DATA5/tmp/kev/misc/20210913hgscPloidyPrediction.txt", sep = "\t",
                       stringsAsFactors = FALSE, header = TRUE)

tableRes$Ploidy <- round2(tableRes$P0)
tableRes2 <- tableRes[,c(1:3,13, 16, 17)]

tableOrig$Ploidy <- round2(tableOrig$P0)
tableOrig2 <- tableOrig[,c(1:3,13, 16, 17)]


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210921hgscGemmPloidyDensityNonInt.pdf")
ggplot(mouse_ploidy, aes(x = P0)) + geom_histogram(binwidth = 0.1, fill = "white", color = "black")  + 
  xlab("ploidy") + ggtitle("HGSC(GEMM) Plodies n = 100") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

