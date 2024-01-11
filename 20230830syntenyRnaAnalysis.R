### looking the rna expression patterns for samples in coad chr18 loss vs no change


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
library(TCGAbiolinks)

# BiocManager::install("ExperimentHub", force = TRUE)
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks", force = TRUE)


# query2 <- GDCquery(
#   project = "TCGA-COAD",
#   data.category = "Transcriptome Profiling",
#   data.type = "Gene Expression Quantification",
#   workflow.type = "STAR - Counts"
# )
# 
# GDCdownload(query = query2, directory = "/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/")
# coadTxn <- GDCprepare(query = query2, directory = "/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/",
#                       save = TRUE,
#                       save.filename = "coad.rda")
# 
# save(coadTxn, file = "/mnt/DATA4/test_nextflow/20230830tcgaBiolinksCoadTxn.RData")
# load("/mnt/DATA4/test_nextflow/20230830tcgaBiolinksCoadTxn.RData")
# 
# apcCoadcolnamesDash2[which(apcCoadcolnamesDash2 %in% coadTxn@colData[,3])]
# 
# coadRawCounts <- coadTxn@assays@data@listData$unstranded
# colnames(coadRawCounts) <- coadTxn$sample
# rownames(coadRawCounts) <- coadTxn@rowRanges$gene_name
# 
# coadRawCounts_red <- coadRawCounts[, which(colnames(coadRawCounts) %in% apcCoadcolnamesDash2)]
# coadRawCounts_red <- coadRawCounts_red[,-which(duplicated(colnames(coadRawCounts_red)))]


# save.image(file = "/mnt/DATA4/test_nextflow/20230829allSyntenyEnivorn.RData")
# load(file = "/mnt/DATA4/test_nextflow/20230829allSyntenyEnivorn.RData")


# query2 <- GDCquery(
#   project = "TCGA-READ",
#   data.category = "Transcriptome Profiling",
#   data.type = "Gene Expression Quantification",
#   workflow.type = "STAR - Counts"
# )
# 
# GDCdownload(query = query2, directory = "/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/")
# readTxn <- GDCprepare(query = query2, directory = "/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/",
#                       save = TRUE,
#                       save.filename = "read.rda")
# 
# save(readTxn, file = "/mnt/DATA4/test_nextflow/20230830tcgaBiolinksReadTxn.RData")
# load("/mnt/DATA4/test_nextflow/20230830tcgaBiolinksReadTxn.RData")




### need to get 18 arm loss statuses ... maybe i have this df already
### just arm loss status in general like the mice


### for colorectal get aneuploidy statuses 
coadArmStatus <- broadBySampleGisticCoad
colnames(coadArmStatus)[2:ncol(coadArmStatus)] <- gsub("(^.*?-.{3}?)-.*", "\\1", colnames(coadArmStatus)[2:ncol(coadArmStatus)])
coadArmStatusMat <- coadArmStatus[, 2:ncol(coadArmStatus)]
rownames(coadArmStatusMat) <- coadArmStatus$`Chromosome Arm`

# coad5q <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == "5q"), ] < -0.2)]
# coad5p <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == "5q"), ] < -0.2)]
coad7q <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == "7q"), ] > 0.2)]
coad7p <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == "7p"), ] > 0.2)]
coad13q <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == "13q"), ] > 0.2)]
# coad16q <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == "16q"), ] > 0.2)]
# coad16p <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == "16p"), ] > 0.2)]
# coad18q <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == "18q"), ] < -0.2)]
# coad18p <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == "18p"), ] < -0.2)]
# coad17p <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == "17p"), ] < -0.2)]

tmpcgs_adenoCar1 <- cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = adenoCar_arm_allSynTable_amp$h_chr, IRanges(adenoCar_arm_allSynTable_amp$h_start, end = adenoCar_arm_allSynTable_amp$h_end))
                                                              ,cancerGeneCensusGr)),]
tmpcgs_adenoCar1$type <- "amp"
tmpcgs_adenoCar2 <- cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = adenoCar_arm_allSynTable_del$h_chr, IRanges(adenoCar_arm_allSynTable_del$h_start, end = adenoCar_arm_allSynTable_del$h_end))
                                                              ,cancerGeneCensusGr)),]
tmpcgs_adenoCar2$type <- "del"
cgs_adenoCar <- rbind(tmpcgs_adenoCar1, tmpcgs_adenoCar2)
cgs_adenoCarGr <- GRanges(seqnames = cgs_adenoCar$chr, IRanges(as.numeric(cgs_adenoCar$start), as.numeric(cgs_adenoCar$end)))
coadArmGr <- GRanges(seqnames = apcCoadArms$chrStripped, IRanges(apcCoadArms$start, apcCoadArms$end))

cgs_adenoCar$arm <- broadBySampleGisticCoad$Chromosome.Arm[queryHits(findOverlaps(coadArmGr, cgs_adenoCarGr))]

### how i got about this for the figure. I'll do each comparison, and grab all cancer census genes found

normalSamples <- colnames(combinedRsemRawMat)[grep("11A", colnames(combinedRsemRawMat))]

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

listOfArms <- c("7p", "7q", "13q")

library(foreach)
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)
start <- Sys.time()
allArmGains <- foreach(i=seq_along(listOfArms),.combine = 'comb', .multicombine = TRUE,
                        .init = list(list(), list()), .packages = c("edgeR", "stringr")) %dopar% {
                          
                          tmpCoadArmNames <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == listOfArms[i]), ] > 0.2)]
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
                          lrtTmp <- glmLRT(fitTmp,coef=2)
                          lrtTabTmp <- data.frame("genes" = rownames(lrtTmp), lrtTmp$table)
                          lrtTabTmp$qval <- p.adjust(lrtTabTmp$PValue, method = "BH", n = length(lrtTabTmp$PValue))
                          lrtTabTmp$arm <- listOfArms[i]
                          
                          tmpRawCounts2 <- combinedRsemRawMat[, -which(colnames(combinedRsemRawMat) %in% tmpCoadArmNames)]
                          tmpGroup2 <- rep(1, ncol(tmpRawCounts2))
                          tmpGroup2[-which(colnames(tmpRawCounts2) %in% normalSamples)] <- 2
                          dgeTmp2 <- DGEList(counts=tmpRawCounts2, group= tmpGroup2)
                          keep2 <- filterByExpr(dgeTmp2)
                          dgeTmp2 <- dgeTmp2[keep2,,keep.lib.sizes=FALSE]
                          dgeTmp2 <- calcNormFactors(dgeTmp2)
                          designTmp2 <- model.matrix(~tmpGroup2)
                          dgeTmp2 <- estimateDisp(dgeTmp2, designTmp2)
                          fitTmp2 <- glmFit(dgeTmp2, designTmp2)
                          lrtTmp2 <- glmLRT(fitTmp2,coef=2)
                          lrtTabTmp2 <- data.frame("genes" = rownames(lrtTmp2), lrtTmp2$table)
                          lrtTabTmp2$qval <- p.adjust(lrtTabTmp2$PValue, method = "BH", n = length(lrtTabTmp2$PValue))
                          lrtTabTmp2$arm <- listOfArms[i]
                          
                          return(list(lrtTabTmp, lrtTabTmp2))
                        }


stopCluster(cl)
print( Sys.time() - start )


allArmGainsDf <- do.call(rbind, allArmGains[[1]])
allArmGainsRecipDf <- do.call(rbind, allArmGains[[2]])

allArmGainsDf <- allArmGainsDf[-grep("\\?", allArmGainsDf$genes), ]
allArmGainsRecipDf <- allArmGainsRecipDf[-grep("\\?", allArmGainsRecipDf$genes), ]

allArmGainsDf$genes <- str_remove(allArmGainsDf$genes, "\\|.*")
allArmGainsRecipDf$genes <- str_remove(allArmGainsRecipDf$genes, "\\|.*")

### running for losses

listOfArms <- c("18p", "18q")
cl <- makeCluster(2)
registerDoParallel(cl)
start <- Sys.time()
allArmLosses <- foreach(i=seq_along(listOfArms),.combine = 'comb', .multicombine = TRUE,
                       .init = list(list(), list()), .packages = c("edgeR", "stringr")) %dopar% {
                         
                         tmpCoadArmNames <- colnames(coadArmStatusMat)[which(coadArmStatusMat[which(rownames(coadArmStatusMat) == listOfArms[i]), ] < 0.2)]
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
                         lrtTmp <- glmLRT(fitTmp,coef=2)
                         lrtTabTmp <- data.frame("genes" = rownames(lrtTmp), lrtTmp$table)
                         lrtTabTmp$qval <- p.adjust(lrtTabTmp$PValue, method = "BH", n = length(lrtTabTmp$PValue))
                         lrtTabTmp$arm <- listOfArms[i]
                         
                         tmpRawCounts2 <- combinedRsemRawMat[, -which(colnames(combinedRsemRawMat) %in% tmpCoadArmNames)]
                         tmpGroup2 <- rep(1, ncol(tmpRawCounts2))
                         tmpGroup2[-which(colnames(tmpRawCounts2) %in% normalSamples)] <- 2
                         dgeTmp2 <- DGEList(counts=tmpRawCounts2, group= tmpGroup2)
                         keep2 <- filterByExpr(dgeTmp2)
                         dgeTmp2 <- dgeTmp2[keep2,,keep.lib.sizes=FALSE]
                         dgeTmp2 <- calcNormFactors(dgeTmp2)
                         designTmp2 <- model.matrix(~tmpGroup2)
                         dgeTmp2 <- estimateDisp(dgeTmp2, designTmp2)
                         fitTmp2 <- glmFit(dgeTmp2, designTmp2)
                         lrtTmp2 <- glmLRT(fitTmp2,coef=2)
                         lrtTabTmp2 <- data.frame("genes" = rownames(lrtTmp2), lrtTmp2$table)
                         lrtTabTmp2$qval <- p.adjust(lrtTabTmp2$PValue, method = "BH", n = length(lrtTabTmp2$PValue))
                         lrtTabTmp2$arm <- listOfArms[i]
                         
                         return(list(lrtTabTmp, lrtTabTmp2))
                       }


stopCluster(cl)
print( Sys.time() - start )


allArmLossesDf <- do.call(rbind, allArmLosses[[1]])
allArmLossesRecipDf <- do.call(rbind, allArmLosses[[2]])

allArmLossesDf <- allArmLossesDf[-grep("\\?", allArmLossesDf$genes), ]
allArmLossesRecipDf <- allArmLossesRecipDf[-grep("\\?", allArmLossesRecipDf$genes), ]

allArmLossesDf$genes <- str_remove(allArmLossesDf$genes, "\\|.*")
allArmLossesRecipDf$genes <- str_remove(allArmLossesRecipDf$genes, "\\|.*")


write.table(allArmGainsDf, "/mnt/DATA5/tmp/kev/misc/20230922allArmGainsDf.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(allArmGainsRecipDf, "/mnt/DATA5/tmp/kev/misc/20230922allArmGainsRecipDf.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(allArmLossesDf, "/mnt/DATA5/tmp/kev/misc/20230922allArmLossesDf.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(allArmLossesRecipDf, "/mnt/DATA5/tmp/kev/misc/20230922allArmLossesRecipDf.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

### looking to create graphs of pathway changes for with and without arm changes
### gsea done on brahm, from the script, 20230922syntenyExpressionCoad.R
### take armGisticCoadApc_freq matrix, get matrix of change states: gain , loss, none
### then count occurences and intersection of different events

# coadArmStatusFactor <- broadBySampleGisticMeltCoad2
# coadArmStatusFactor$status <- 0
# coadArmStatusFactor$status[which(coadArmStatusFactor$mean > 0.2)] <- 1
# coadArmStatusFactor$status[which(coadArmStatusFactor$mean < -0.2)] <- -1
# coadArmStatusFactor$arm <- rep(coadArmStatus$`Chromosome Arm`, length(unique(coadArmStatusFactor$sampleID)))
# coadArmStatusFactor$status2 <- paste(coadArmStatusFactor$arm, coadArmStatusFactor$status)


### exlucded regions where it's both gained and loss
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
write.table(coadArmMemebershipDf, "/mnt/DATA5/tmp/kev/misc/20230926coadArmChangeMembership.txt",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# library(UpSetR)
# 
# colnames(coadArmMemebershipMat)[2:9]
# coadArmMemebershipMatUpset <- coadArmMemebershipMat[, c(2, 4:7, 9)]
# colnames(coadArmMemebershipMatUpset) <- c("chr1pLoss", "chr7Gain", "chr13qGain",
#                                           "chr14qLoss", "chr15qLoss", "chr18qLoss")
# 
# upset(data.frame(coadArmMemebershipMatUpset),
#       nsets = 100,
#       mb.ratio = c(0.6, 0.4),
#       number.angles = 0, 
#       text.scale = 1.1, 
#       point.size = 2.8, 
#       line.size = 1
# )
# 
# 
# pdf("/mnt/DATA5/tmp/kev/misc/20231011upsetPlotReduced.pdf", useDingbats = FALSE,
#     width = 14, height = 7)
# upset(data.frame(coadArmMemebershipMatUpset),
#       mb.ratio = c(0.6, 0.4),
#       number.angles = 0, 
#       text.scale = 2.5, 
#       point.size = 3.5, 
#       line.size = 0.5,
#       intersections = list(list("chr18qLoss"),
#                            list("chr18qLoss", "chr1pLoss"), list("chr18qLoss", "chr7Gain"),
#                            list("chr18qLoss", "chr13qGain"), list("chr18qLoss", "chr14qLoss"),
#                            list("chr18qLoss", "chr15qLoss"),
#                            list("chr18qLoss", "chr7Gain", "chr13qGain"),
#                            list("chr18qLoss", "chr7Gain", "chr14qLoss"),
#                            list("chr18qLoss", "chr7Gain", "chr15qLoss"),
#                            list("chr18qLoss", "chr14qLoss", "chr15qLoss"),
#                            list("chr18qLoss",  "chr13qGain", "chr15qLoss"),
#                            list("chr18qLoss", "chr7Gain", "chr13qGain", "chr15qLoss"),
#                            list("chr18qLoss", "chr13qGain", "chr14qLoss", "chr15qLoss"),
#                            list("chr18qLoss", "chr7Gain", "chr13qGain", "chr14qLoss"),
#                            list("chr18qLoss", "chr7Gain", "chr13qGain", "chr14qLoss", "chr15qLoss"),
#                            list("chr18qLoss", "chr1pLoss", "chr7Gain", "chr13qGain", "chr14qLoss", "chr15qLoss"))
# )
# 
# dev.off()


# looking at genes within the syntenic regions, only interested in 7, 15q, 18q




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

with_7_15_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7q1 == 1 &
                                       coadArmMemebershipDf$status15q1 == 1 &
                                       coadArmMemebershipDf$status18q1 == 1)]
no_7_15_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7q1 == 0 &
                                       coadArmMemebershipDf$status15q1 == 0 &
                                       coadArmMemebershipDf$status18q1 == 0)]

with_1_7_15_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status1p1 == 1 &
                                                       coadArmMemebershipDf$status7q1 == 1 &
                                                       coadArmMemebershipDf$status15q1 == 1 &
                                                       coadArmMemebershipDf$status18q1 == 1)]
no_1_7_15_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status1p1  == 0 &
                                                     coadArmMemebershipDf$status7q1 == 0 &
                                                     coadArmMemebershipDf$status15q1 == 0 &
                                                     coadArmMemebershipDf$status18q1 == 0)]


with_7_13_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7q1 == 1 &
                                                       coadArmMemebershipDf$status13q1 == 1 &
                                                       coadArmMemebershipDf$status18q1 == 1)]

no_7_13_18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status7q1 == 0 &
                                                  coadArmMemebershipDf$status13q1 == 0 &
                                                  coadArmMemebershipDf$status18q1 == 0)]


###
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




# listOfSamples2 <- list(with_1, with_7, with_13, with_14, with_15, with_18, with_7_15_18, with_1_7_15_18)
# listOfRecip2 <- list(no_1, no_7, no_13, no_14, no_15, no_18, no_7_15_18, no_1_7_15_18 )
# listOfSamplesNames2 <- c("1", "7", "13", "14", "15","18", "7_15_18", "1_7_15_18")


listOfSamples2 <- list(with_18, with_1_18, with_7_18, with_13_18, with_14_18, with_15_18, with_7_13_18, with_7_18_no_13, with_13_18_no_7, with_18_no_7_13, no_18)
listOfRecip2 <- list(no_18, no_1_18, no_7_18, no_13_18, no_14_18, no_15_18, with_18_no_7_13, with_18_no_7_13, with_18_no_7_13, with_18, with_18)
listOfSamplesNames2 <- c("with_18", "with_1_18", "with_7_18", "with_13_18", "with_14_18", "with_15_18",
                         "with_7_13_18", "with_7_18_no_13", "with_13_18_no_7", "with_18_no_7_13", "no_18")

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


write.table(diffArmLossComboDf2, "/mnt/DATA5/tmp/kev/misc/20231019diffArmLossComboDf2.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(diffArmLossComboRecipDf2, "/mnt/DATA5/tmp/kev/misc/20231019diffArmLossRecipComboDf2.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

### looking at the log fold change of arm change vs no arm change

no18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status18q1 == 0)]

only18 <- rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status1p1 == 0 &
                                                 coadArmMemebershipDf$status18q1 == 1 &
                                                 coadArmMemebershipDf$status13q1 == 0 & 
                                                 coadArmMemebershipDf$status15q1 == 0 &
                                                 coadArmMemebershipDf$status14q1 == 0 &
                                                 coadArmMemebershipDf$status7q1 == 0)]


only18_13 <-  rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status1p1 == 0 &
                                                     coadArmMemebershipDf$status18q1 == 1 &
                                                     coadArmMemebershipDf$status13q1 == 1 & 
                                                     coadArmMemebershipDf$status15q1 == 0 &
                                                     coadArmMemebershipDf$status14q1 == 0 &
                                                     coadArmMemebershipDf$status7q1 == 0)]


only18_7 <-  rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status1p1 == 0 &
                                                    coadArmMemebershipDf$status18q1 == 1 &
                                                    coadArmMemebershipDf$status13q1 == 0 & 
                                                    coadArmMemebershipDf$status15q1 == 0 &
                                                    coadArmMemebershipDf$status14q1 == 0 &
                                                    coadArmMemebershipDf$status7q1 == 1)]


only18_15 <-  rownames(coadArmMemebershipDf)[which(coadArmMemebershipDf$status1p1 == 0 &
                                                    coadArmMemebershipDf$status18q1 == 1 &
                                                    coadArmMemebershipDf$status13q1 == 0 & 
                                                    coadArmMemebershipDf$status15q1 == 1 &
                                                    coadArmMemebershipDf$status7q1 == 0)]





listOfSamples3 <- list(only18, only18_13, only18_7, only18_15)
listOfRecip3<- list(no18, only18, only18, only18)
listOfSamplesNames3 <- c("only18", "only18_13", "only18_7", "only18_15")

library(foreach)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
start <- Sys.time()
diffArmLossCombo3 <- foreach(i=seq_along(listOfSamples3),.combine = 'comb', .multicombine = TRUE,
                             .init = list(list(), list()), .packages = c("edgeR", "stringr")) %dopar% {
                               
                               tmpCoadArmNames <- unlist(listOfSamples3[i])
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
                               lrtTabTmp$arm <- listOfSamplesNames3[i]
                               
                               tmpCoadArmNames2 <- unlist(listOfRecip3[i])
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
                               lrtTabTmp2$arm <- listOfSamplesNames3[i]
                               
                               return(list(lrtTabTmp, lrtTabTmp2))
                             }


stopCluster(cl)
print( Sys.time() - start )


diffArmLossComboDf3 <- do.call(rbind, diffArmLossCombo3[[1]])
diffArmLossComboRecipDf3 <- do.call(rbind, diffArmLossCombo3[[2]])

diffArmLossComboDf3$genes <- str_remove(diffArmLossComboDf3$genes, "\\|.*")
diffArmLossComboRecipDf3$genes <- str_remove(diffArmLossComboRecipDf3$genes, "\\|.*")


write.table(diffArmLossComboDf3, "/mnt/DATA5/tmp/kev/misc/20231023diffArmLossComboDf3.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(diffArmLossComboRecipDf3, "/mnt/DATA5/tmp/kev/misc/20231023diffArmLossRecipComboDf3.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)





listOfSamples3 <- list(only18, only18_13, only18_7, only18_15)
listOfRecip3<- list(no18, only18, only18, only18)
listOfSamplesNames3 <- c("only18", "only18_13", "only18_7", "only18_15")

library(foreach)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
start <- Sys.time()
diffArmLossCombo3 <- foreach(i=seq_along(listOfSamples3),.combine = 'comb', .multicombine = TRUE,
                             .init = list(list(), list()), .packages = c("edgeR", "stringr")) %dopar% {
                               
                               tmpCoadArmNames <- unlist(listOfSamples3[i])
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
                               lrtTabTmp$arm <- listOfSamplesNames3[i]
                               
                               tmpCoadArmNames2 <- unlist(listOfRecip3[i])
                               tmpRawCounts2 <- combinedRsemRawMat[, which(colnames(combinedRsemRawMat) %in% c(normalSamples, tmpCoadArmNames2))]
                               tmpGroup2 <- rep(1, ncol(tmpRawCounts2))
                               tmpGroup2[which(colnames(tmpRawCounts2) %in% tmpCoadArmNames2)] <- 2
                               dgeTmp2 <- DGEList(counts=tmpRawCounts2, group= tmpGroup2)
                               keep2 <- filterByExpr(dgeTmp2)
                               dgeTmp2 <- dgeTmp2[keep2,,keep.lib.sizes=FALSE]
                               dgeTmp2 <- calcNormFactors(dgeTmp)
                               designTmp2 <- model.matrix(~tmpGroup2)
                               dgeTmp2 <- estimateDisp(dgeTmp2, designTmp2)
                               fitTmp2 <- glmFit(dgeTmp2, designTmp2)
                               lrtTmp2 <- glmLRT(fitTmp2,coef="tmpGroup2")
                               lrtTabTmp2 <- data.frame("genes" = rownames(lrtTmp2), lrtTmp2$table)
                               lrtTabTmp2$qval <- p.adjust(lrtTab2Tmp$PValue, method = "BH", n = length(lrtTabTmp2$PValue))
                               lrtTabTmp2$arm <- listOfSamplesNames3[i]
                               
                               return(list(lrtTabTmp, lrtTabTmp2))
                             }


stopCluster(cl)
print( Sys.time() - start )


diffArmLossComboDf3 <- do.call(rbind, diffArmLossCombo3[[1]])
diffArmLossComboRecipDf3 <- do.call(rbind, diffArmLossCombo3[[2]])

diffArmLossComboDf3$genes <- str_remove(diffArmLossComboDf3$genes, "\\|.*")
diffArmLossComboRecipDf3$genes <- str_remove(diffArmLossComboRecipDf3$genes, "\\|.*")


write.table(diffArmLossComboDf3, "/mnt/DATA5/tmp/kev/misc/20231023diffArmLossComboDf3.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(diffArmLossComboRecipDf3, "/mnt/DATA5/tmp/kev/misc/20231023diffArmLossRecipComboDf3.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)






adenoCar_arm_allSynTable_amp2 <- adenoCar_arm_allSynTable_amp[which(adenoCar_arm_allSynTable_amp$sigCheck == "good"), ]
adenoCar_arm_allSynTable_amp2$str <- paste(adenoCar_arm_allSynTable_amp2$h_chr, adenoCar_arm_allSynTable_amp2$m_chr)
i <- unique(adenoCar_arm_allSynTable_amp2$h_chr)
affectedGenesGains <- NULL
for (i in unique(adenoCar_arm_allSynTable_amp2$h_chr)) {
  tmp <- adenoCar_arm_allSynTable_amp2[which(adenoCar_arm_allSynTable_amp2$h_chr == i),]
  tmpGr <- GRanges(seqnames = tmp$h_chr, IRanges(start = tmp$h_start, end = tmp$h_end))
  tmpGene <- hg19AllGeneExonLocs$gene[queryHits(findOverlaps(hg19AllGeneExonLocsGr, tmpGr))]
  affectedGenesGains <- rbind(affectedGenesGains, data.frame("h_chr" = tmp$h_chr[subjectHits(findOverlaps(hg19AllGeneExonLocsGr, tmpGr))],
                                                             "mchr" = tmp$m_chr[subjectHits(findOverlaps(hg19AllGeneExonLocsGr, tmpGr))],
                                                             "gene" = tmpGene,
                                                             "type" = rep("gains", length(tmpGene))))
}

adenoCar_arm_allSynTable_del2 <- adenoCar_arm_allSynTable_del[which(adenoCar_arm_allSynTable_del$sigCheck == "good"), ]
adenoCar_arm_allSynTable_del2$str <- paste(adenoCar_arm_allSynTable_del2$h_chr, adenoCar_arm_allSynTable_del2$m_chr)
affectedGenesLosses <- NULL
for (i in unique(adenoCar_arm_allSynTable_del2$h_chr)) {
  tmp <- adenoCar_arm_allSynTable_del2[which(adenoCar_arm_allSynTable_del2$h_chr == i),]
  tmpGr <- GRanges(seqnames = tmp$h_chr, IRanges(start = tmp$h_start, end = tmp$h_end))
  tmpGene <- hg19AllGeneExonLocs$gene[queryHits(findOverlaps(hg19AllGeneExonLocsGr, tmpGr))]
  affectedGenesLosses <- rbind(affectedGenesLosses, data.frame("h_chr" = tmp$h_chr[subjectHits(findOverlaps(hg19AllGeneExonLocsGr, tmpGr))],
                                                               "mchr" = tmp$m_chr[subjectHits(findOverlaps(hg19AllGeneExonLocsGr, tmpGr))],
                                                               "gene" = tmpGene,
                                                               "type" = rep("losses", length(tmpGene))))
}




affectedGenesDf <- rbind(affectedGenesGains, affectedGenesLosses)
affectedGenesDf <- affectedGenesDf[which(affectedGenesDf$h_chr %in% c(1, 7, 8, 13, 14, 15, 18)), ]

### looking at different pathways
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
msig_hallmark <- msig_hallmark[, c("gs_name", "gene_symbol")]

affectedGenesDf_apoptosis <- affectedGenesDf[which(affectedGenesDf$gene %in% msig_hallmark$gene_symbol[which(msig_hallmark$gs_name == "HALLMARK_APOPTOSIS")]), ]
# affectedGenesDf_apoptosis <- affectedGenesDf[which(affectedGenesDf$gene %in% msig_hallmark$gene_symbol[which(msig_hallmark$gs_name == "HALLMARK_IL6_JAK_STAT3_SIGNALING")]), ]

expression_1 <- diffArmLossComboDf2[which(diffArmLossComboDf2$arm == "1"), ]
expression_1_apop <- expression_1[match(affectedGenesDf_apoptosis$gene, expression_1$genes), ]
expression_1_apop$location <- affectedGenesDf_apoptosis$h_chr

expression_no1 <- diffArmLossComboRecipDf2[which(diffArmLossComboRecipDf2$arm == "1"), ]
expression_no1_apop <- expression_no1[match(affectedGenesDf_apoptosis$gene, expression_no1$genes), ]
expression_no1_apop$location <- affectedGenesDf_apoptosis$h_chr



expression_7 <- diffArmLossComboDf2[which(diffArmLossComboDf2$arm == "7"), ]
expression_7_apop <- expression_7[match(affectedGenesDf_apoptosis$gene, expression_7$genes), ]
expression_7_apop$location <- affectedGenesDf_apoptosis$h_chr

expression_no7 <- diffArmLossComboRecipDf2[which(diffArmLossComboRecipDf2$arm == "7"), ]
expression_no7_apop <- expression_no7[match(affectedGenesDf_apoptosis$gene, expression_no7$genes), ]
expression_no7_apop$location <- affectedGenesDf_apoptosis$h_chr


expression_15 <- diffArmLossComboDf2[which(diffArmLossComboDf2$arm == "15"), ]
expression_15_apop <- expression_15[match(affectedGenesDf_apoptosis$gene, expression_15$genes), ]
expression_15_apop$location <- affectedGenesDf_apoptosis$h_chr

expression_no15 <- diffArmLossComboRecipDf2[which(diffArmLossComboRecipDf2$arm == "15"), ]
expression_no15_apop <- expression_no15[match(affectedGenesDf_apoptosis$gene, expression_no15$genes), ]
expression_no15_apop$location <- affectedGenesDf_apoptosis$h_chr


expression_18 <- diffArmLossComboDf2[which(diffArmLossComboDf2$arm == "18"), ]
expression_18_apop <- expression_18[match(affectedGenesDf_apoptosis$gene, expression_18$genes), ]
expression_18_apop$location <- affectedGenesDf_apoptosis$h_chr

expression_no18 <- diffArmLossComboRecipDf2[which(diffArmLossComboRecipDf2$arm == "18"), ]
expression_no18_apop <- expression_no18[match(affectedGenesDf_apoptosis$gene, expression_no18$genes), ]
expression_no18_apop$location <- affectedGenesDf_apoptosis$h_chr


expression_7_15_18 <- diffArmLossComboDf2[which(diffArmLossComboDf2$arm == "7_15_18"), ]
expression_7_15_18_apop <- expression_7_15_18[match(affectedGenesDf_apoptosis$gene, expression_7_15_18$genes), ]
expression_7_15_18_apop$location <- affectedGenesDf_apoptosis$h_chr

expression_no7_15_18 <- diffArmLossComboRecipDf2[which(diffArmLossComboRecipDf2$arm == "7_15_18"), ]
expression_no7_15_18_apop <- expression_no7_15_18[match(affectedGenesDf_apoptosis$gene, expression_no7_15_18$genes), ]
expression_no7_15_18_apop$location <- affectedGenesDf_apoptosis$h_chr

expression_1_7_15_18 <- diffArmLossComboDf2[which(diffArmLossComboDf2$arm == "1_7_15_18"), ]
expression_1_7_15_18_apop <- expression_1_7_15_18[match(affectedGenesDf_apoptosis$gene, expression_1_7_15_18$genes), ]
expression_1_7_15_18_apop$location <- affectedGenesDf_apoptosis$h_chr

expression_no1_7_15_18 <- diffArmLossComboRecipDf2[which(diffArmLossComboRecipDf2$arm == "1_7_15_18"), ]
expression_no1_7_15_18_apop <- expression_no1_7_15_18[match(affectedGenesDf_apoptosis$gene, expression_no1_7_15_18$genes), ]
expression_no1_7_15_18_apop$location <- affectedGenesDf_apoptosis$h_chr

### just looking at overrall heatmap of genes by arm change

affectedGenesDf2 <- rbind(affectedGenesGains, affectedGenesLosses)
affectedGenesDf2 <- affectedGenesDf2[-which(affectedGenesDf2$h_chr == 8),]
diffArmLossComboDf3$arm2 <- str_remove(diffArmLossComboDf3$arm, "with\\_")

interestedGenes <- NULL
allAffectedGeneDf <- NULL
for (i in seq_along(unique(diffArmLossComboDf3$arm))) {
  tmpDf <- diffArmLossComboDf3[which(diffArmLossComboDf3$arm == unique(diffArmLossComboDf3$arm)[i]), ]
  tmpGeneDf <- affectedGenesDf2[which(affectedGenesDf2$h_chr == unique(diffArmLossComboDf3$arm2)[i]), ]
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
for (i in seq_along(unique(diffArmLossComboDf3$arm))) {
  tmpDf <- diffArmLossComboDf3[which(diffArmLossComboDf3$arm == unique(diffArmLossComboDf3$arm)[i]), ]
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
         color = heatmapColGenes, breaks = colorsBreaksGenes, cellwidth = 10, cellheight = 8)

### for each of these genes, look for overlaps in pathways affected

msigDf <- data.frame("genes" = interestedGenes)
msigDf$sig <- "none"
msigDf$sig <- msig_hallmark$gs_name[match(msigDf$genes, msig_hallmark$gene_symbol)]
msigDf$cgs <- "none"
msigDf$cgs <- ifelse(!is.na(match(msigDf$genes, cancerGeneCensus2$Gene.Symbol)), "yes", "no")
msigDf$cgsAll <- "none"
msigDf$cgsAll <- ifelse(!is.na(match(msigDf$genes, cancerGeneCensusAll2$Gene.Symbol)), "yes", "no")
msigDf$chr <- hg19AllGeneExonLocs$chrom[match(msigDf$genes, hg19AllGeneExonLocs$gene)]

### for all genes graph volcano-like plot; color based on gene location ... i.e only show color if
### (1) gene passes fold-change and significance threshold (2) possibly if it's going in the direction
allAffectedGeneDf$log10Q <- -1 * log10(allAffectedGeneDf$qval)
allAffectedGeneDf$color <- "#000000"
allAffectedGeneDf$color[which(allAffectedGeneDf$arm == "with_1")] <- "#ADD8E6"
allAffectedGeneDf$color[which(allAffectedGeneDf$arm == "with_14")] <- "#008080"
allAffectedGeneDf$color[which(allAffectedGeneDf$arm == "with_15")] <- "#800080"
allAffectedGeneDf$color[which(allAffectedGeneDf$arm == "with_18")] <- "#00008B"
allAffectedGeneDf$color[which(allAffectedGeneDf$arm == "with_7")] <- "#ff8c00"
allAffectedGeneDf$color[which(allAffectedGeneDf$arm == "with_13")] <- "#8b0000"
allAffectedGeneDf$color[which(allAffectedGeneDf$logFC > (0.9 * log2(1/2)) & allAffectedGeneDf$logFC < (0.9 * log2(3/2)))] <- "#000000"
allAffectedGeneDf$color[which(allAffectedGeneDf$qval > 0.05)] <- "#000000"


# save.image("/mnt/DATA4/test_nextflow/20231024aneuploidyRna.RData")
# load("/mnt/DATA4/test_nextflow/20231024aneuploidyRna.RData")
library(ggrepel)
pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20231024volcanoAneuploidy", ".pdf"),
    useDingbats = FALSE, width = 10, height = 5)
ggplot(allAffectedGeneDf) + geom_point(aes(x = logFC, y = log10Q), color = allAffectedGeneDf$color)  +
  geom_text_repel(data = allAffectedGeneDf[-which(allAffectedGeneDf$color == "#000000"), ],
            aes(x = logFC, y = log10Q, label = genes),
            inherit.aes = FALSE, vjust = 0.5,size = 2) + 
  scale_x_continuous(breaks = seq(-4, 4, 1)) + ylab("-log10QValue") + xlab("log2FC") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

### easier comparison with boxplots since it's only a few genes
### taking code from 20230913syntenyAllCancerComparison

# armSamples = with_7_15_18
# noChangeSamples =  no_7_15_18
# df =  rownames(expression_7_15_18_apop)
geneZCalc2 <- function(df, armSamples = NULL, noChangeSamples = NULL, normals = normals){
  if (is.null(noChangeSamples) | is.null(armSamples)) {
    stop()
    print("need both samples (var = armSamples) with changes and without (noChangeSamples)")
  }
  
  if (is.null(normals)) {
    stop()
    print("need normal sample IDs")
  }
  
  geneZDf <- NULL
  
  ### takes in geneDf and then calcualtes them one at a time
  for (i in seq_along(df)) {
    tmpDf <- combinedRsemIllumina[which(combinedRsemIllumina$`Hybridization REF` %in% df[i]),]
    tmpDf[,2:ncol(tmpDf)] <- lapply(tmpDf[,2:ncol(tmpDf)], function(x) log2(as.numeric(unlist(x)) + 1))
    colnames(tmpDf) <- gsub("(^.*?-.{3}?)-.*", "\\1", colnames(tmpDf))
    tmpMean <- mean(as.numeric(tmpDf[1, which(colnames(tmpDf) %in% normals)]))
    tmpSd <- sd(as.numeric(tmpDf[1, which(colnames(tmpDf) %in% normals)]))
    tmpGain <-  (as.numeric(tmpDf[1, which(colnames(tmpDf) %in% armSamples)]) - tmpMean)/tmpSd
    tmpDiploid <- (as.numeric(tmpDf[1, which(colnames(tmpDf) %in% noChangeSamples)]) - tmpMean)/tmpSd
    geneZDf <- rbind(geneZDf,
                     data.frame("type" = rep("diploid", length(tmpDiploid)),
                                "gene" = rep(df[i], length(tmpDiploid)),
                                "value" = tmpDiploid),
                     data.frame("type" = rep("gain", length(tmpGain)),
                                "gene" = rep(df[i], length(tmpGain)),
                                "value" = tmpGain))
  }
  return(geneZDf)
}


listOfGenes <- rownames(expression_7_15_18_apop)
listOfGenes <- substr(listOfGenes, 1, nchar(listOfGenes)-1)
chr7_15_18BoxplotDf <- geneZCalc2(listOfGenes, with_7_15_18, no_7_15_18, normals)
chr7_15_18BoxplotDf$gene2 <- str_remove(chr7_15_18BoxplotDf$gene, "\\|.*")
chr7_15_18BoxplotDf$gene2 <- factor(chr7_15_18BoxplotDf$gene2, levels = unique(chr7_15_18BoxplotDf$gene2))

synChr7_15_18 <- ggplot(chr7_15_18BoxplotDf, aes(x = type, y = value)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(color = "black", size=0.4, alpha=0.9) + facet_wrap(~ gene2, nrow = 1) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab("COADREAD-TP chr13q status") +
  ylab("z-scores of noramlized RSEM values relative to normal samples") + 
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 1)) +
  stat_compare_means(method = "t.test")


pdf("/mnt/DATA5/tmp/kev/misc/20231011boxplotApop7_15_18.pdf", useDingbats = FALSE,
    width = 15, height = 7)
synChr7_15_18
dev.off()


listOfGenes <- rownames(expression_1_7_15_18_apop)
listOfGenes <- substr(listOfGenes, 1, nchar(listOfGenes)-1)
chr1_7_15_18BoxplotDf <- geneZCalc2(listOfGenes, with_1_7_15_18, no_1_7_15_18, normals)
chr1_7_15_18BoxplotDf$gene2 <- str_remove(chr1_7_15_18BoxplotDf$gene, "\\|.*")
chr1_7_15_18BoxplotDf$gene2 <- factor(chr1_7_15_18BoxplotDf$gene2, levels = unique(chr1_7_15_18BoxplotDf$gene2))

synChr1_7_15_18 <- ggplot(chr1_7_15_18BoxplotDf, aes(x = type, y = value)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(color = "black", size=0.4, alpha=0.9) + facet_wrap(~ gene2, nrow = 1) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + xlab("COADREAD-TP chr13q status") +
  ylab("z-scores of noramlized RSEM values relative to normal samples") + 
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 1)) +
  stat_compare_means(method = "t.test")

### loop for individual arm changes




### doing similar but for ov



### last thing to check is for any relationship to the staple mutations after APC
### simple oncoprint at first and layer in the anueploidies
coadreadMutTable <- read.table("/mnt/DATA5/tmp/kev/misc/20231012coadread_all_relevant_mut_sample_matrixV2.txt",
                               sep = "\t", header = TRUE)
# library("ComplexHeatmap")
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

coadreadMutTableMat <- coadreadMutTable2[, 2:ncol(coadreadMutTable2)]
rownames(coadreadMutTableMat) <- coadreadMutTable2$studyID.sampleId
coadreadMutTableMat <- t(coadreadMutTableMat)
coadreadMutTableMat <- coadreadMutTableMat[grep("chr", rownames(coadreadMutTableMat)), ]
pdf(paste0("/mnt/DATA5/tmp/kev/misc/20231013aneuploidyFeaturesHeatmap", ".pdf"), useDingbats = FALSE, width = 15)
pheatmap(coadreadMutTableMat, show_colnames = FALSE, cellheight = 20, cellwidth = 2
         ,color = c("white","darkgreen"))
dev.off()




### function to is do pairwise comparisons for fisher's exact method for mutual exclusion or co-occurrence
### https://stackoverflow.com/questions/54543665/looping-all-pairwise-comparisons-from-a-list-in-r



pairwiseCoadRes <- NULL
compMat <- t(combn(colnames(coadreadMutTable2[2:ncol(coadreadMutTable2)]), m = (2)))
for (i in 1:nrow(compMat)) {
  tmpFisherRes <- fisher.test(coadreadMutTable2[ , compMat[i, 1]], coadreadMutTable2[ , compMat[i, 2]],
                              alternative = "two.sided")
  tmpZ <- qnorm(tmpFisherRes$p.value) * sign(-log2(tmpFisherRes$estimate))
  pairwiseCoadRes <- rbind(pairwiseCoadRes, data.frame("Var1" = compMat[i, 1], "Var2" = compMat[i, 2], "Z" = tmpZ))
}

### same as above but removed extra variables used for interactions i.e KRAS, TP53 -> KRAS_TP53

coadreadMutTable3 <- coadreadMutTable2[, -which(colnames(coadreadMutTable2) %in% c("NRAS", "HRAS"))]
colnames(coadreadMutTable3)[2] <- "KRAS/NRAS/HRAS"
coadreadMutTable3$`KRAS/NRAS/HRAS` <- 0
coadreadMutTable3$`KRAS/NRAS/HRAS`[which(apply(coadreadMutTable2[c("KRAS", "NRAS", "HRAS")],1, sum) > 0)] <- 1


pairwiseCoadRes2 <- NULL
compMat2 <- t(combn(colnames(coadreadMutTable3[2:ncol(coadreadMutTable3)]), m = (2)))
for (i in 1:nrow(compMat2)) {
  tmpFisherRes <- fisher.test(coadreadMutTable3[ , compMat2[i, 1]], coadreadMutTable3[ , compMat2[i, 2]],
                              alternative = "two.sided")
  tmpZ <- qnorm(tmpFisherRes$p.value) * sign(-log2(tmpFisherRes$estimate))
  pairwiseCoadRes2 <- rbind(pairwiseCoadRes2, data.frame("Var1" = compMat2[i, 1], "Var2" = compMat2[i, 2], "Z" = tmpZ))
}


coadreadMutTableMat2 <- coadreadMutTable3[, 2:ncol(coadreadMutTable3)]
rownames(coadreadMutTableMat2) <- coadreadMutTable3$studyID.sampleId
coadreadMutTableMat2 <- t(coadreadMutTableMat2)
pdf(paste0("/mnt/DATA5/tmp/kev/misc/20231013molecularFeaturesHeatmap", ".pdf"), useDingbats = FALSE, width = 15)
pheatmap(coadreadMutTableMat2, show_colnames = FALSE, cellheight = 20, cellwidth = 2
         ,color = c("white","darkgreen"))
dev.off()

### redoing the pathway stuff with the genotype data



coadArmStatusFactor <- coadArmStatus
coadArmStatusFactor$`Chromosome Arm` <- coadArmStatus$`Chromosome Arm`
coadArmStatusFactorRed <- coadArmStatusFactor[which(coadArmStatusFactor$`Chromosome Arm` %in% c("1p", "7p", "7q", "8p","13q",
                                                                                             "14q", "15q","18p", "18q")), ]

coadArmStatusFactorRed[1, 2:ncol(coadArmStatusFactorRed)]  <- ifelse(coadArmStatusFactorRed[1, 2:ncol(coadArmStatusFactorRed)] < -0.2, 1, 0)
coadArmStatusFactorRed[2, 2:ncol(coadArmStatusFactorRed)]  <- ifelse(coadArmStatusFactorRed[2, 2:ncol(coadArmStatusFactorRed)] > 0.2, 1, 0)
coadArmStatusFactorRed[3, 2:ncol(coadArmStatusFactorRed)]  <- ifelse(coadArmStatusFactorRed[3, 2:ncol(coadArmStatusFactorRed)] > 0.2, 1, 0)
coadArmStatusFactorRed[4, 2:ncol(coadArmStatusFactorRed)]  <- ifelse(coadArmStatusFactorRed[4, 2:ncol(coadArmStatusFactorRed)] < -0.2, 1, 0)
coadArmStatusFactorRed[5, 2:ncol(coadArmStatusFactorRed)]  <- ifelse(coadArmStatusFactorRed[5, 2:ncol(coadArmStatusFactorRed)] > 0.2, 1, 0)
coadArmStatusFactorRed[6, 2:ncol(coadArmStatusFactorRed)]  <- ifelse(coadArmStatusFactorRed[6, 2:ncol(coadArmStatusFactorRed)] < -0.2, 1, 0)
coadArmStatusFactorRed[7, 2:ncol(coadArmStatusFactorRed)]  <- ifelse(coadArmStatusFactorRed[7, 2:ncol(coadArmStatusFactorRed)] < -0.2, 1, 0)
coadArmStatusFactorRed[8, 2:ncol(coadArmStatusFactorRed)]  <- ifelse(coadArmStatusFactorRed[8, 2:ncol(coadArmStatusFactorRed)] < -0.2, 1, 0)
coadArmStatusFactorRed[9, 2:ncol(coadArmStatusFactorRed)]  <- ifelse(coadArmStatusFactorRed[9, 2:ncol(coadArmStatusFactorRed)] < -0.2, 1, 0)

status1pRed <- factor(unlist(coadArmStatusFactorRed[1, 2:ncol(coadArmStatusFactorRed)]))
status7pRed <- factor(unlist(coadArmStatusFactorRed[2, 2:ncol(coadArmStatusFactorRed)]))
status7qRed <- factor(unlist(coadArmStatusFactorRed[3, 2:ncol(coadArmStatusFactorRed)]))
status8pRed <- factor(unlist(coadArmStatusFactorRed[4, 2:ncol(coadArmStatusFactorRed)]))
status13qRed <- factor(unlist(coadArmStatusFactorRed[5, 2:ncol(coadArmStatusFactorRed)]))
status14qRed <- factor(unlist(coadArmStatusFactorRed[6, 2:ncol(coadArmStatusFactorRed)]))
status15qRed <- factor(unlist(coadArmStatusFactorRed[7, 2:ncol(coadArmStatusFactorRed)]))
status18pRed <- factor(unlist(coadArmStatusFactorRed[8, 2:ncol(coadArmStatusFactorRed)]))
status18qRed <- factor(unlist(coadArmStatusFactorRed[9, 2:ncol(coadArmStatusFactorRed)]))

# coadreadMutTable4 <- coadreadMutTable3[match(names(status1pRed), coadreadMutTable3$studyID.sampleId),]
# statusRas <- factor(coadreadMutTable4$`KRAS/NRAS/HRAS`)
# statusTp53 <- factor(coadreadMutTable4$TP53)
# statusPik3ca <- factor(coadreadMutTable4$PIK3CA)

coadArmMemebershipMatRed <- model.matrix(~ 0 + status1pRed + status7pRed + status7qRed +status8pRed + status13qRed +
                                        status14qRed + status15qRed + status18pRed + status18qRed)

coadArmMemebershipDf2 <- data.frame(coadArmMemebershipMatRed)
rownames(coadArmMemebershipDf2) <- colnames(coadArmStatusFactorRed)[2:ncol(coadArmStatusFactorRed)]



chr18_no13 <- rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status18pRed1 == 1 &
                                                        coadArmMemebershipDf2$status13qRed1 == 0)]

chr18_no7 <- rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status18pRed1 == 1 &
                                                      coadArmMemebershipDf2$status7pRed1 == 0)]

chr18_no15 <- rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status18pRed1 == 1 &
                                                     coadArmMemebershipDf2$status15qRed1 == 0)]

chr18_no14 <- rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status18pRed1 == 1 &
                                                      coadArmMemebershipDf2$status14qRed1 == 0)]

chr18_no13 <- rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status18pRed1 == 1 &
                                                      coadArmMemebershipDf2$status13qRed1 == 0)]

chr18_no1 <- rownames(coadArmMemebershipDf2)[which(coadArmMemebershipDf2$status18pRed1 == 1 &
                                                      coadArmMemebershipDf2$status1pRed1 == 0)]

listOfSamples <- list(chr18_no13, chr18_no7, chr18_no15, chr18_no14, chr18_no13, chr18_no1)
listOfSamplesNames <- c("chr18_no13", "chr18_no7", "chr18_no15", "chr18_no14", "chr18_no13", "chr18_no1")

library(foreach)
library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)
start <- Sys.time()
diffArmLossCombo <- foreach(i=seq_along(listOfSamples),.combine = 'comb', .multicombine = TRUE,
                            .init = list(list()), .packages = c("edgeR", "stringr")) %dopar% {
                              
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
                              
                              return(list(lrtTabTmp))
                            }


stopCluster(cl)
print( Sys.time() - start )


diffArmLossComboDf <- do.call(rbind, diffArmLossCombo[[1]])
diffArmLossComboDf$genes <- str_remove(diffArmLossComboDf$genes, "\\|.*")
write.table(diffArmLossComboDf, "/mnt/DATA5/tmp/kev/misc/20231013diffArmLossComboDf.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)





# GDCdownload(query = query3, directory = "/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/")
# ovTxn <- GDCprepare(query = query3, directory = "/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/")

# save(ovTxn, file = "/mnt/DATA4/test_nextflow/20230830tcgaBiolinksOvTxn.RData")
# load("/mnt/DATA4/test_nextflow/20230830tcgaBiolinksOvTxn.RData")

ovRawCounts <- ovTxn@assays@data@listData$unstranded
colnames(ovRawCounts) <- ovTxn$sample
rownames(ovRawCounts) <- ovTxn@rowRanges$gene_name

tp53OvcolnamesDash <- str_replace_all(colnames(tp53OvArms), "\\.", "\\-")
tp53OvcolnamesDash2 <- gsub("(^.*?-.{3}?)-.*", "\\1", tp53OvcolnamesDash)

ovRawCounts_red <- ovRawCounts[, which(colnames(ovRawCounts) %in% tp53OvcolnamesDash2)]

ovArmStatus <- broadBySampleGistic
colnames(ovArmStatus)[2:ncol(ovArmStatus)] <- stringr::str_replace_all(colnames(ovArmStatus)[2:ncol(ovArmStatus)], "\\.", "\\-")
colnames(ovArmStatus)[2:ncol(ovArmStatus)] <- gsub("(^.*?-.{3}?)-.*", "\\1", colnames(ovArmStatus)[2:ncol(ovArmStatus)])
ovArmStatusMat <- ovArmStatus[, 2:ncol(ovArmStatus)]
rownames(ovArmStatusMat) <- ovArmStatus$Chromosome.Arm

ovArmVector <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "7q", "8q", "9q", "10p", "10q", "11p", "11q",
                 "12p", "12q", "13p", "14p", "15p", "16p", "17p", "17q", "20p", "20q")
ovArmChangeVector <- c("gain", "gain", "gain", "gain", "gain", "gain",
                    "loss",  "gain", "gain", "loss", "gain", "loss",
                    "loss", "loss", "gain", "gain", "loss", "loss",
                    "loss", "loss",  "loss", "loss", "gain", "gain")


tmpcgs_hgsc_bprn1 <- cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bprn_arm_allSynTable_amp$h_chr, IRanges(hgsc_bprn_arm_allSynTable_amp$h_start, end = hgsc_bprn_arm_allSynTable_amp$h_end))
                                                              ,cancerGeneCensusGr)),]
tmpcgs_hgsc_bprn1$type <- "amp"
tmpcgs_hgsc_bprn2 <- cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_bprn_arm_allSynTable_del$h_chr, IRanges(hgsc_bprn_arm_allSynTable_del$h_start, end = hgsc_bprn_arm_allSynTable_del$h_end))
                                                              ,cancerGeneCensusGr)),]
tmpcgs_hgsc_bprn2$type <- "del"
cgs_hgsc_bprn <- rbind(tmpcgs_hgsc_bprn1, tmpcgs_hgsc_bprn2)
cgs_hgsc_bprnGr <- GRanges(seqnames = cgs_hgsc_bprn$chr, IRanges(as.numeric(cgs_hgsc_bprn$start), as.numeric(cgs_hgsc_bprn$end)))
ovArmGr <- GRanges(seqnames = tp53OvArms$chrStripped, IRanges(tp53OvArms$start, tp53OvArms$end))
cgs_hgsc_bprn$arm <- broadBySampleGistic$Chromosome.Arm[queryHits(findOverlaps(ovArmGr, cgs_hgsc_bprnGr))]


ovLogFc <- NULL
ovQval <- NULL
ovAllLrtTable <- list()
for (i in seq_along(ovArmChangeVector)) {
  tmpGroup <- rep(1, ncol(ovRawCounts_red))
  if (ovArmChangeVector[i] == "gain") {
    tmpArmVar <- colnames(ovArmStatusMat)[which(ovArmStatusMat[which(rownames(ovArmStatusMat) == ovArmVector[i]), ] > 0.2)]
  } else if(ovArmChangeVector[i] == "loss"){
    tmpArmVar <- colnames(ovArmStatusMat)[which(ovArmStatusMat[which(rownames(ovArmStatusMat) == ovArmVector[i]), ] < -0.2)]
  }
  tmpGroup[which(colnames(ovRawCounts_red) %in% tmpArmVar)] <- 2
  dgeTmp <- DGEList(counts=ovRawCounts_red,group= tmpGroup)
  keep <- filterByExpr(dgeTmp)
  dgeTmp <- dgeTmp[keep,,keep.lib.sizes=FALSE]
  dgeTmp <- calcNormFactors(dgeTmp)
  designTmp <- model.matrix(~tmpGroup)
  dgeTmp <- estimateDisp(dgeTmp, designTmp)
  fitTmp <- glmFit(dgeTmp, designTmp)
  lrtTmp <- glmLRT(fitTmp,coef=2)
  lrtTabTmp <- data.frame("genes" = rownames(lrtTmp), lrtTmp$table)
  lrtTabTmp$qval <- p.adjust(lrtTabTmp$PValue, method = "BH", n = length(lrtTabTmp$PValue) * length(ovArmChangeVector))
  lrtTabTmp2 <- lrtTabTmp[which(lrtTabTmp$genes %in% cgs_adenoCar$Gene.Symbol), ]
  
  ovLogFc <- cbind(ovLogFc, lrtTabTmp2$logFC)
  ovQval <- cbind(ovQval, lrtTabTmp2$qval)
  
  ovAllLrtTable[[i]] <- lrtTabTmp
}

colnames(ovLogFc) <- paste0(ovArmVector[1:16], "-", ovArmChangeVector[1:16])
rownames(ovLogFc) <- lrtTabTmp2$genes

colnames(ovQval) <- paste0(ovArmVector[1:16], "-", ovArmChangeVector[1:16])
rownames(ovQval) <- lrtTabTmp2$genes

names(allLrtTable) <- ovArmChangeVector

ovLogFc2 <- ovLogFc
ovLogFc2[which(ovQval > 0.05)] <- 0
ovLogFc2 <- ovLogFc2[-which(apply(ovLogFc2, 1, sum) == 0),]

heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-2,2,4/1000)
pheatmap(ovLogFc2, cluster_cols = FALSE, breaks = colors.breaks, color = heatMapCol)
