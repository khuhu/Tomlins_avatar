### analysis looking at segments that have the entire region gained
### filter out extremely small scnasthat have high gains and losses 
### easier to process b/c less multiple overlaps that I need to average out 
### 24kb is median human gene length, so get rid of single gene events



coadSegFile <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/gdac.broadinstitute.org_COADREAD.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/COADREAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt",
                          sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
coadSegFile$Sample <- gsub("(^.*?-.{3}?)-.*", "\\1",coadSegFile$Sample)
coadSegFile$Sample <- substr(coadSegFile$Sample, 1, nchar(coadSegFile$Sample) -1)
coadSegFile <- coadSegFile[which(coadSegFile$Sample %in% coad_anno2$Sample.ID),]
### looks at distribution of probes - again for things that are possible aneuploidy or large changes, but not quite focal
coadSegFile$lengthKb <- (coadSegFile$End - coadSegFile$Start)/1e4
ggplot(coadSegFile) + geom_point(aes(lengthKb, Segment_Mean))

coadSegFile2 <- coadSegFile[which(coadSegFile$lengthKb > 24), ]

adenoCar_arm_allSynTable_amp_good <- adenoCar_arm_allSynTable_amp[which(adenoCar_arm_allSynTable_amp$sigCheck == "good"), ]
adenoCar_arm_allSynTable_del_good <- adenoCar_arm_allSynTable_del[which(adenoCar_arm_allSynTable_del$sigCheck == "good"), ]

adenoCar_arm_allSynTable_amp_good$length <- (adenoCar_arm_allSynTable_amp_good$h_end - adenoCar_arm_allSynTable_amp_good$h_start)/1e6
adenoCar_arm_allSynTable_del_good$length <- (adenoCar_arm_allSynTable_del_good$h_end - adenoCar_arm_allSynTable_del_good$h_start)/1e6

adenoCar_arm_allSynTable_amp_good$block <- paste0("block", 1:nrow(adenoCar_arm_allSynTable_amp_good))
adenoCar_arm_allSynTable_del_good$block <- paste0("block", 1:nrow(adenoCar_arm_allSynTable_del_good))


### find a way to iterate through each significant syntenic region - either 50 indivudal blocks or I combine
### should combine or too many variables for downstream analysis
### nevermind; too much work to combine, need to do individual



i <- unique(adenoCar_arm_allSynTable_amp_good$h_chr)[1]
hGoodDiff <- NULL
for (i in unique(adenoCar_arm_allSynTable_amp_good$h_chr)) {
  tmpDf <- adenoCar_arm_allSynTable_amp_good[which(adenoCar_arm_allSynTable_amp_good$h_chr == i), ]
  tmpDiff <- NULL
  if (nrow(tmpDf) > 1) {
    for (j in 2:nrow(tmpDf)) {
      tmpDiff <- c(tmpDiff, tmpDf$h_start[j] - tmpDf$h_end[j-1])
    }
    tmpDiff <- c(0, tmpDiff)
    hGoodDiff <- c(hGoodDiff, tmpDiff)
  } else{
    hGoodDiff <- c(hGoodDiff, 0)
  }
}
hGoodDiff <- hGoodDiff/1e6
adenoCar_arm_allSynTable_amp_good$diff <- hGoodDiff

### for each region, which samples have these regions amplified or deleted in the same direction
### do differential expression analysis for those samples against tumor samples without;
### make list of all samples then iterate 

coadSegFile_amp <- coadSegFile[which(coadSegFile$Segment_Mean > 0.2), ]
coadSegFile_ampGr <- GRanges(seqnames = coadSegFile_amp$Chromosome,
                             IRanges(start = coadSegFile_amp$Start,
                                     end = coadSegFile_amp$End))

allPercentDf <- NULL
listAllSampsWithGain <- list()
for (i in 1:nrow(adenoCar_arm_allSynTable_amp_good)) {
  tmpHgGrange <- GRanges(seqnames = adenoCar_arm_allSynTable_amp_good$h_chr[i],
                         IRanges(start = adenoCar_arm_allSynTable_amp_good$h_start[i],
                                 end = adenoCar_arm_allSynTable_amp_good$h_end[i]))
  tmpOverlap <- findOverlaps(tmpHgGrange, coadSegFile_ampGr)
  tmpIntersect <- pintersect(tmpHgGrange[queryHits(tmpOverlap)],
                             coadSegFile_ampGr[subjectHits(tmpOverlap)])
  ### getting percent of overlap with the syntenic contig
  tmpPercent <- width(tmpIntersect) / width(tmpHgGrange[queryHits(tmpOverlap)])
  
  tmpAmp <- coadSegFile_amp[subjectHits(tmpOverlap),]
  tmpAmp$percentOverlap <- tmpPercent
  
  tmpPercentDf <- NULL
  for (j in unique(tmpAmp$Sample)) {
    tmpPercentDf <- rbind(tmpPercentDf, data.frame("sample" = j,
                                                   "percent" = sum(tmpAmp$percentOverlap[which(tmpAmp$Sample == j)])))
  }
  
  tmpKeepSamps <- tmpPercentDf$sample[which(tmpPercentDf$percent > 0.99)]
  listAllSampsWithGain[[i]] <- tmpKeepSamps
  
  
  # quantile(tmpPercentDf$percent, seq(0, 1, 0.01))
  ### from all of the overlaps, give a tiny bit of leeway
  ### only keep things >= 99%
  allPercentDf <- rbind(allPercentDf, tmpPercentDf)
}

quantile(allPercentDf$percent, seq(0, 1, 0.01))
hist(allPercentDf$percent, breaks = 99)
names(listAllSampsWithGain) <- paste0("block_amp_", 1:nrow(adenoCar_arm_allSynTable_amp_good))
listAllSampsWithGain2 <- listAllSampsWithGain[-which(unlist(lapply(listAllSampsWithGain, length)) == 0)]
### iterate over rna comparisons. samps with vs without


combinedRsemRawMat2 <- combinedRsemRawMat[, which(colnames(combinedRsemRawMat) %in% paste0(coad_anno2$Sample.ID, "A"))]

normalSamples <- colnames(combinedRsemRawMat)[grep("11A", colnames(combinedRsemRawMat2))]


library(foreach)
library(doParallel)
cl <- makeCluster(15)
registerDoParallel(cl)
start <- Sys.time()
diffArmAmp <- foreach(i=seq_along(listAllSampsWithGain2),.combine = 'comb', .multicombine = TRUE,
                             .init = list(list()), .packages = c("edgeR", "stringr")) %dopar% {
                               
                               tmpCoadArmNames <- unlist(listAllSampsWithGain2[i])
                               tmpMat <- combinedRsemRawMat2
                               # tmpMat <- tmpMat[,-which(colnames(tmpMat) %in% normalSamples)]
                               colnames(tmpMat) <- substr(colnames(tmpMat), 1, nchar(colnames(tmpMat)) -1)
                               tmpRawCounts <- tmpMat
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
                               lrtTabTmp$block <- names(listAllSampsWithGain2)[i]

                               return(list(lrtTabTmp))
                             }


stopCluster(cl)
print( Sys.time() - start )


diffArmAmpDf <- do.call(rbind, diffArmAmp[[1]])
diffArmAmpDf$genes <- str_remove(diffArmAmpDf$genes, "\\|.*")


# write.table(diffArmAmpDf, "/mnt/DATA5/tmp/kev/misc/20231102diffArmAmp.txt",sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
### for deletions
###
###

coadSegFile_del <- coadSegFile[which(coadSegFile$Segment_Mean < -0.2), ]
coadSegFile_delGr <- GRanges(seqnames = coadSegFile_del$Chromosome,
                             IRanges(start = coadSegFile_del$Start,
                                     end = coadSegFile_del$End))


listAllSampsWithLoss <- list()
for (i in 1:nrow(adenoCar_arm_allSynTable_del_good)) {
  tmpHgGrange <- GRanges(seqnames = adenoCar_arm_allSynTable_del_good$h_chr[i],
                         IRanges(start = adenoCar_arm_allSynTable_del_good$h_start[i],
                                 end = adenoCar_arm_allSynTable_del_good$h_end[i]))
  tmpOverlap <- findOverlaps(tmpHgGrange, coadSegFile_delGr)
  tmpIntersect <- pintersect(tmpHgGrange[queryHits(tmpOverlap)],
                             coadSegFile_delGr[subjectHits(tmpOverlap)])
  ### getting percent of overlap with the syntenic contig
  tmpPercent <- width(tmpIntersect) / width(tmpHgGrange[queryHits(tmpOverlap)])
  
  tmpDel <- coadSegFile_del[subjectHits(tmpOverlap),]
  tmpDel$percentOverlap <- tmpPercent
  
  tmpPercentDf <- NULL
  for (j in unique(tmpDel$Sample)) {
    tmpPercentDf <- rbind(tmpPercentDf, data.frame("sample" = j,
                                                   "percent" = sum(tmpDel$percentOverlap[which(tmpDel$Sample == j)])))
  }
  
  tmpKeepSamps <- tmpPercentDf$sample[which(tmpPercentDf$percent > 0.99)]
  listAllSampsWithLoss[[i]] <- tmpKeepSamps
  
}


names(listAllSampsWithLoss) <- paste0("block_loss_", 1:nrow(adenoCar_arm_allSynTable_del_good))
listAllSampsWithLoss2 <- listAllSampsWithLoss[-which(unlist(lapply(listAllSampsWithLoss, length)) == 0)]
### iterate over rna comparisons. samps with vs without


library(foreach)
library(doParallel)
cl <- makeCluster(15)
registerDoParallel(cl)
start <- Sys.time()
diffArmDel <- foreach(i=seq_along(listAllSampsWithLoss2),.combine = 'comb', .multicombine = TRUE,
                      .init = list(list()), .packages = c("edgeR", "stringr")) %dopar% {
                        
                        tmpCoadArmNames <- unlist(listAllSampsWithLoss2[i])
                        tmpMat <- combinedRsemRawMat2
                        # tmpMat <- tmpMat[,-which(colnames(tmpMat) %in% normalSamples)]
                        colnames(tmpMat) <- substr(colnames(tmpMat), 1, nchar(colnames(tmpMat)) -1)
                        tmpRawCounts <- tmpMat
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
                        lrtTabTmp$block <- names(listAllSampsWithLoss2)[i]
                        
                        return(list(lrtTabTmp))
                      }


stopCluster(cl)
print( Sys.time() - start )


diffArmDelDf <- do.call(rbind, diffArmDel[[1]])
diffArmDelDf$genes <- str_remove(diffArmDelDf$genes, "\\|.*")


# write.table(diffArmDelDf, "/mnt/DATA5/tmp/kev/misc/20231102diffArmDel.txt",sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)


### creating volcano plot for loss of block 3 and gain of 22

hg19GeneLocations <- read.table("/mnt/DATA5/tmp/kev/misc/20231018hg19GeneLocations.txt", sep = "\t",
                                header = TRUE, stringsAsFactors = FALSE)

hg19ExonGr <- GRanges(seqnames = hg19GeneLocations$chr, 
                      IRanges(start = hg19GeneLocations$start, hg19GeneLocations$end))
affectedGenesGains <- NULL
for (i in unique(adenoCar_arm_allSynTable_amp_good$block)) {
  tmp <- adenoCar_arm_allSynTable_amp_good[which(adenoCar_arm_allSynTable_amp_good$block == i),]
  tmpGr <- GRanges(seqnames = tmp$h_chr, IRanges(start = tmp$h_start, end = tmp$h_end))
  tmpGene <- hg19GeneLocations[queryHits(findOverlaps(hg19ExonGr, tmpGr)),]
  affectedGenesGains <- rbind(affectedGenesGains, data.frame(tmpGene, "type" = rep(paste0("amp_", i), nrow(tmpGene))))
}


affectedGenesLosses <- NULL
for (i in unique(adenoCar_arm_allSynTable_del_good$block)) {
  tmp <- adenoCar_arm_allSynTable_del_good[which(adenoCar_arm_allSynTable_del_good$block == i),]
  tmpGr <- GRanges(seqnames = tmp$h_chr, IRanges(start = tmp$h_start, end = tmp$h_end))
  tmpGene <- hg19GeneLocations[queryHits(findOverlaps(hg19ExonGr, tmpGr)),]
  affectedGenesLosses <- rbind(affectedGenesLosses, data.frame(tmpGene, "type" = rep(paste0("del_", i), nrow(tmpGene))))
}



affectedGenesDf <- rbind(affectedGenesGains, affectedGenesLosses)

### 

lrtTab_bl3 <- diffArmDelDf[which(diffArmDelDf$block == "block_loss_29"),]

lrtTab_bl3$log10Q <- -1 * log10(lrtTab_bl3$qval)
lrtTab_bl3$color <- "#000000"
lrtTab_bl3$color[which(lrtTab_bl3$qval < 0.05 & lrtTab_bl3$logFC > (log2(1.5 * 0.7)))] <- "#8b0000"
lrtTab_bl3$color[which(lrtTab_bl3$qval < 0.05 & lrtTab_bl3$logFC < (log2(0.5/0.7)))] <- "#00008B"
lrtTab_bl3$on_block3 <- "no"
lrtTab_bl3$on_block3[which(lrtTab_bl3$genes %in% affectedGenesDf$gene[which(affectedGenesDf$type == "del_block29")])] <- "yes"
lrtTab_bl3$color2 <- lrtTab_bl3$color
lrtTab_bl3$color2[which(lrtTab_bl3$on_block3 == "no")] <- "#000000"
lrtTab_bl3$color2[which(lrtTab_bl3$color == "#8b0000")] <- "#000000"
lrtTab_bl3$alpha <- 1
lrtTab_bl3$alpha[which(lrtTab_bl3$color2 == "#000000")] <- 0.1

library(ggrepel)


# pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20231106volcanoBlockLoss29", ".pdf"),
#     useDingbats = FALSE, width = 10, height = 10)
# ggplot(lrtTab_bl3) + geom_point(aes(x = logFC, y = log10Q), color = lrtTab_bl3$color) + 
#   ylab("-log10QValue") + xlab("log2FC") + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         text = element_text(size=20)) +
#   xlim(c(-3, 3))
# dev.off()

# pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20231106volcanoBlockLoss29Filt", ".pdf"),
#     useDingbats = FALSE, width = 10, height = 10)
ggplot(lrtTab_bl3) + geom_point(aes(x = logFC, y = log10Q),
                                color = lrtTab_bl3$color2, alpha = lrtTab_bl3$alpha)  + 
  geom_text_repel(data = lrtTab_bl3[-which(lrtTab_bl3$color2 == "#000000"), ],
                  aes(x = logFC, y = log10Q, label = genes),
                  inherit.aes = FALSE, vjust = 0.5,
                  max.overlaps = 5, size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  xlim(c(-3, 3)) + ylab("-log10QValue") + xlab("log2FC")
# dev.off()



lrtTab_bl22 <- diffArmAmpDf[which(diffArmAmpDf$block == "block_amp_22"),]

lrtTab_bl22$log10Q <- -1 * log10(lrtTab_bl22$qval)
lrtTab_bl22$color <- "#000000"
lrtTab_bl22$color[which(lrtTab_bl22$qval < 0.05 & lrtTab_bl22$logFC > (log2(1.5 * 0.7)))] <- "#8b0000"
lrtTab_bl22$color[which(lrtTab_bl22$qval < 0.05 & lrtTab_bl22$logFC < (log2(0.5/0.7)))] <- "#00008B"
lrtTab_bl22$on_block3 <- "no"
lrtTab_bl22$on_block3[which(lrtTab_bl22$genes %in% affectedGenesDf$gene[which(affectedGenesDf$type == "amp_block22")])] <- "yes"
lrtTab_bl22$color2 <- lrtTab_bl22$color
lrtTab_bl22$color2[which(lrtTab_bl22$on_block3 == "no")] <- "#000000"
lrtTab_bl22$color2[which(lrtTab_bl22$color == "#00008B")] <- "#000000"
lrtTab_bl22$alpha <- 1
lrtTab_bl22$alpha[which(lrtTab_bl22$color2 == "#000000")] <- 0.1



# pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20231106volcanoBlockAmp22", ".pdf"),
#     useDingbats = FALSE, width = 10, height = 10)
ggplot(lrtTab_bl22) + geom_point(aes(x = logFC, y = log10Q), color = lrtTab_bl22$color) +
  ylab("-log10QValue") + xlab("log2FC") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  xlim(c(-3, 3))
# dev.off()

# pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20231106volcanoBlockLoss29Filt", ".pdf"),
#     useDingbats = FALSE, width = 10, height = 10)
ggplot(lrtTab_bl22) + geom_point(aes(x = logFC, y = log10Q),
                                color = lrtTab_bl22$color2, alpha = lrtTab_bl22$alpha)  + 
  geom_text_repel(data = lrtTab_bl22[-which(lrtTab_bl22$color2 == "#000000"), ],
                  aes(x = logFC, y = log10Q, label = genes),
                  inherit.aes = FALSE, vjust = 0.5,
                  max.overlaps = 5, size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  xlim(c(-3, 3)) + ylab("-log10QValue") + xlab("log2FC")
# dev.off()

### using either cox-proportional for all synteny blocks


library(TCGAbiolinks)
library(survminer)
library(survival)

clinical_coad <- GDCquery_clinic("TCGA-COAD")
clinical_read <- GDCquery_clinic("TCGA-READ")

clinical_coad2 <- clinical_coad[, which(colnames(clinical_coad) %in% c("submitter_id", "vital_status", "days_to_last_follow_up", "days_to_death"))]
clinical_read2 <- clinical_read[, which(colnames(clinical_read) %in% c("submitter_id", "vital_status", "days_to_last_follow_up", "days_to_death"))]

clinical_coadread <- rbind(clinical_coad2, clinical_read2)
clinical_coadread$deceased <- ifelse(clinical_coadread$vital_status == "Alive", FALSE, TRUE)
clinical_coadread$overall_survival <- ifelse(clinical_coadread$vital_status == "Alive",
                                             clinical_coadread$days_to_last_follow_up,
                                             clinical_coadread$days_to_death)
clinical_coadread$submitter_id2 <- paste0(clinical_coadread$submitter_id, "-01")


clinical_coadread <- clinical_coadread[which(paste0(clinical_coadread$submitter_id2, "A") %in% colnames(combinedRsemRawMat2)),]



### another way to look at this is do survival analysis for all blocks - see which are interesting

blockStatusLossDf <- NULL
for (i in seq_along(names(listAllSampsWithLoss2))) {
  tmpNames <- listAllSampsWithLoss2[[i]]
  tmpVector <- ifelse(clinical_coadread$submitter_id2 %in% tmpNames, "Yes", "No")
  blockStatusLossDf <- cbind(blockStatusLossDf, tmpVector)
}
colnames(blockStatusLossDf) <- names(listAllSampsWithLoss2)

blockStatusGainDf <- NULL
for (i in seq_along(names(listAllSampsWithGain2))) {
  tmpNames <- listAllSampsWithGain2[[i]]
  tmpVector <- ifelse(clinical_coadread$submitter_id2 %in% tmpNames, "Yes", "No")
  blockStatusGainDf <- cbind(blockStatusGainDf, tmpVector)
}
colnames(blockStatusGainDf) <- names(listAllSampsWithGain2)

clinical_coadread <- cbind(clinical_coadread, blockStatusGainDf, blockStatusLossDf)

fit <- survfit(Surv(overall_survival, deceased) ~ block_amp_20, data = clinical_coadread)
pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20231109logRankBlockAmp20", ".pdf"),
    useDingbats = FALSE, width = 10, height = 5)
ggsurvplot(fit,
           data = clinical_coadread,
           pval = T,
           risk.table = T)
dev.off()
### iterating for all arms
survivalBlocksDf <- NULL
for (i in 8:ncol(clinical_coadread)) {
  tmpFit <- survdiff(Surv(overall_survival, deceased) ~ clinical_coadread[,i], data = clinical_coadread)
  tmpPVal <- pchisq(tmpFit$chisq, length(tmpFit$n)-1, lower.tail = FALSE)
  survivalBlocksDf <- rbind(survivalBlocksDf, data.frame("var" = colnames(clinical_coadread)[i], "pval" = tmpPVal))
}


clinical_coadread2 <- clinical_coadread[which(clinical_coadread$deceased), ]
survivalBlocksDf2 <- NULL
for (i in 8:ncol(clinical_coadread2)) {
  tmpFit <- survdiff(Surv(overall_survival, deceased) ~ clinical_coadread2[,i], data = clinical_coadread2)
  tmpPVal <- pchisq(tmpFit$chisq, length(tmpFit$n)-1, lower.tail = FALSE)
  survivalBlocksDf2 <- rbind(survivalBlocksDf2, data.frame("var" = colnames(clinical_coadread2)[i], "pval" = tmpPVal))
}

fit <- survfit(Surv(overall_survival, deceased) ~ block_amp_20, data = clinical_coadread2)
ggsurvplot(fit,
           data = clinical_coadread2,
           pval = T,
           risk.table = T)




### notes from scott, move legend for B to left, and then arrow to signify meth od


### for each of these genes, look for overlaps in pathways affected
diffArmDelDf$type <- "losses"
diffArmAmpDf$type <- "gains"
affectedGenesDf$type2 <- affectedGenesDf$type
affectedGenesDf$type2 <- str_replace(affectedGenesDf$type2, "del_block", "block_loss_")
affectedGenesDf$type2 <- str_replace(affectedGenesDf$type2, "amp_block", "block_amp_")

allArmLossGainDf <- rbind(diffArmDelDf, diffArmAmpDf)
interestedGenes <- NULL
for (i in seq_along(unique(allArmLossGainDf$block))) {
  tmpDf <- allArmLossGainDf[which(allArmLossGainDf$block == unique(allArmLossGainDf$block)[i]), ]
  tmpGeneDf <- affectedGenesDf[which(affectedGenesDf$type2 == unique(allArmLossGainDf$block)[i]), ]
  tmpDf <- tmpDf[which(tmpDf$genes %in%  tmpGeneDf$gene), ]
  if (tmpDf$type[1] == "gains") {
    tmpGenes <- tmpDf$genes[which(tmpDf$logFC > (log2(1.5 * 0.9)) & tmpDf$qval < 0.05)]
  } else if(tmpDf$type[1] == "losses"){
    tmpGenes <- tmpDf$genes[which(tmpDf$logFC < (log2(0.5/0.9)) & tmpDf$qval < 0.05)]
  }
  interestedGenes <- c(interestedGenes, tmpGenes)
}


library(msigdbr)
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
msig_hallmark <- msig_hallmark[, c("gs_name", "gene_symbol")]

msigDf <- data.frame("genes" = interestedGenes)
msigDf$sig <- "none"
msigDf$sig <- msig_hallmark$gs_name[match(msigDf$genes, msig_hallmark$gene_symbol)]
msigDf$chr <- hg19GeneLocations$chr[match(msigDf$genes, hg19GeneLocations$gene)]
msigDf$block <- affectedGenesDf$type2[match(msigDf$genes, affectedGenesDf$gene)]



### similar to above, but instead of edgeR use DEseq2

cl <- makeCluster(20)
registerDoParallel(cl)
start <- Sys.time()
diffArmAmpDeSeq <- foreach(i=seq_along(listAllSampsWithGain2),.combine = 'comb', .multicombine = TRUE,
                           .init = list(list()), .packages = c("DESeq2", "stringr")) %dopar% {
                             
                             tmpCoadArmNames <- unlist(listAllSampsWithGain2[i])
                             tmpMat <- combinedRsemRawMat2
                             colnames(tmpMat) <- substr(colnames(tmpMat), 1, nchar(colnames(tmpMat)) -1)
                             tmpRawCounts <- tmpMat
                             tmpGroup <- rep(1, ncol(tmpRawCounts))
                             tmpGroup[which(colnames(tmpRawCounts) %in% tmpCoadArmNames)] <- 2
                             
                             tmpColData <- data.frame("block" = factor(tmpGroup))
                             rownames(tmpColData) <- colnames(tmpRawCounts)
                             
                             dds <-  DESeqDataSetFromMatrix(countData = round(tmpRawCounts),
                                                            colData = tmpColData,
                                                            design = ~block)
                             smallestGroupSize <- 3
                             keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
                             dds <- dds[keep,]
                             
                             dds$block <- relevel(dds$block, ref = "1")
                             
                             dds <- DESeq(dds)
                             # res <- results(dds)
                             res2 <- lfcShrink(dds, coef = "block_2_vs_1",type="apeglm")
                             
                             # resTable <- data.frame(cbind("genes" = rownames(res), 
                             #                   res))
                             resTable2 <- data.frame(cbind("genes" = rownames(res2), 
                                                           res2))
                             # or to shrink log fold changes association with condition:
                             
                             resTable2$block <- names(listAllSampsWithGain2)[i]
                             
                             return(list(resTable2))
                           }


stopCluster(cl)
print( Sys.time() - start )


diffArmAmpDeSeqDf <- do.call(rbind, diffArmAmpDeSeq[[1]])
diffArmAmpDeSeqDf$genes <- str_remove(diffArmAmpDeSeqDf$genes, "\\|.*")

# write.table(diffArmAmpDeSeqDf, "/mnt/DATA5/tmp/kev/misc/20231109diffArmAmpDeSeq.txt",sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)



### same DEseq but for losses


cl <- makeCluster(20)
registerDoParallel(cl)
start <- Sys.time()
diffArmDelDeSeq <- foreach(i=seq_along(listAllSampsWithLoss2),.combine = 'comb', .multicombine = TRUE,
                           .init = list(list()), .packages = c("DESeq2", "stringr")) %dopar% {
                             
                             tmpCoadArmNames <- unlist(listAllSampsWithLoss2[i])
                             tmpMat <- combinedRsemRawMat2
                             colnames(tmpMat) <- substr(colnames(tmpMat), 1, nchar(colnames(tmpMat)) -1)
                             tmpRawCounts <- tmpMat
                             tmpGroup <- rep(1, ncol(tmpRawCounts))
                             tmpGroup[which(colnames(tmpRawCounts) %in% tmpCoadArmNames)] <- 2
                             
                             tmpColData <- data.frame("block" = factor(tmpGroup))
                             rownames(tmpColData) <- colnames(tmpRawCounts)
                             
                             dds <-  DESeqDataSetFromMatrix(countData = round(tmpRawCounts),
                                                            colData = tmpColData,
                                                            design = ~block)
                             smallestGroupSize <- 3
                             keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
                             dds <- dds[keep,]
                             
                             dds$block <- relevel(dds$block, ref = "1")
                             
                             dds <- DESeq(dds)
                             # res <- results(dds)
                             res2 <- lfcShrink(dds, coef = "block_2_vs_1",type="apeglm")
                             
                             # resTable <- data.frame(cbind("genes" = rownames(res), 
                             #                   res))
                             resTable2 <- data.frame(cbind("genes" = rownames(res2), 
                                                           res2))
                             # or to shrink log fold changes association with condition:
                             
                             resTable2$block <- names(listAllSampsWithLoss2)[i]
                             
                             return(list(resTable2))
                           }


stopCluster(cl)
print( Sys.time() - start )


diffArmDelDeSeqDf <- do.call(rbind, diffArmDelDeSeq[[1]])
diffArmDelDeSeqDf$genes <- str_remove(diffArmDelDeSeqDf$genes, "\\|.*")

# write.table(diffArmDelDeSeqDf, "/mnt/DATA5/tmp/kev/misc/20231109diffArmDelDeSeq.txt",sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)







lrtTab_bl20 <- diffArmAmpDeSeqDf[which(diffArmAmpDeSeqDf$block == "block_amp_20"),]

lrtTab_bl20$log10Q <- -1 * log10(lrtTab_bl20$padj)
lrtTab_bl20$color <- "#000000"
lrtTab_bl20$color[which(lrtTab_bl20$padj < 0.05 & lrtTab_bl20$log2FoldChange > (log2(1.5 * 0.9)))] <- "#8b0000"
lrtTab_bl20$color[which(lrtTab_bl20$padj < 0.05 & lrtTab_bl20$log2FoldChange < (log2(0.5/0.9)))] <- "#00008B"
lrtTab_bl20$on_block <- "no"
lrtTab_bl20$on_block[which(lrtTab_bl20$genes %in% affectedGenesDf$gene[which(affectedGenesDf$type == "amp_block20")])] <- "yes"
lrtTab_bl20$color2 <- lrtTab_bl20$color
lrtTab_bl20$color2[which(lrtTab_bl20$on_block == "no")] <- "#000000"
lrtTab_bl20$color2[which(lrtTab_bl20$color == "#00008B")] <- "#000000"
lrtTab_bl20$alpha <- 1
lrtTab_bl20$alpha[which(lrtTab_bl20$color2 == "#000000")] <- 0.1



# pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20231109volcanoBlockAmp20", ".pdf"),
#     useDingbats = FALSE, width = 10, height = 10)
ggplot(lrtTab_bl20) + geom_point(aes(x = log2FoldChange, y = log10Q), color = lrtTab_bl20$color) +
  ylab("-log10QValue") + xlab("log2FC") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  xlim(c(-3, 3))
# dev.off()

pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20231109volcanoBlockAmp20Filt", ".pdf"),
    useDingbats = FALSE, width = 10, height = 10)
ggplot(lrtTab_bl20) + geom_point(aes(x = log2FoldChange, y = log10Q),
                                 color = lrtTab_bl20$color2, alpha = lrtTab_bl20$alpha)  + 
  geom_text_repel(data = lrtTab_bl20[-which(lrtTab_bl20$color2 == "#000000"), ],
                  aes(x = log2FoldChange, y = log10Q, label = genes),
                  inherit.aes = FALSE, vjust = 0.5,
                  max.overlaps = 5, size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  xlim(c(-3, 3)) + ylab("-log10QValue") + xlab("log2FC")
dev.off()





lrtTab_bl30 <- diffArmDelDeSeqDf[which(diffArmDelDeSeqDf$block == "block_loss_30"),]

lrtTab_bl30$log10Q <- -1 * log10(lrtTab_bl30$padj)
lrtTab_bl30$color <- "#000000"
lrtTab_bl30$color[which(lrtTab_bl30$padj < 0.05 & lrtTab_bl30$log2FoldChange > (log2(1.5 * 0.9)))] <- "#8b0000"
lrtTab_bl30$color[which(lrtTab_bl30$padj < 0.05 & lrtTab_bl30$log2FoldChange < (log2(0.5/0.9)))] <- "#00008B"
lrtTab_bl30$on_block <- "no"
lrtTab_bl30$on_block[which(lrtTab_bl30$genes %in% affectedGenesDf$gene[which(affectedGenesDf$type == "del_block30")])] <- "yes"
lrtTab_bl30$color2 <- lrtTab_bl30$color
lrtTab_bl30$color2[which(lrtTab_bl30$on_block == "no")] <- "#000000"
lrtTab_bl30$color2[which(lrtTab_bl30$color == "#8b0000")] <- "#000000"
lrtTab_bl30$alpha <- 1
lrtTab_bl30$alpha[which(lrtTab_bl30$color2 == "#000000")] <- 0.1



# pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20231106volcanoBlockAmp22", ".pdf"),
#     useDingbats = FALSE, width = 10, height = 10)
ggplot(lrtTab_bl30) + geom_point(aes(x = log2FoldChange, y = log10Q), color = lrtTab_bl30$color) +
  ylab("-log10QValue") + xlab("log2FC") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  xlim(c(-3, 3))
# dev.off()

# pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/20231106volcanoBlockLoss30Filt", ".pdf"),
#     useDingbats = FALSE, width = 10, height = 10)
ggplot(lrtTab_bl30) + geom_point(aes(x = log2FoldChange, y = log10Q),
                                 color = lrtTab_bl30$color2, alpha = lrtTab_bl30$alpha)  + 
  geom_text_repel(data = lrtTab_bl30[-which(lrtTab_bl30$color2 == "#000000"), ],
                  aes(x = log2FoldChange, y = log10Q, label = genes),
                  inherit.aes = FALSE, vjust = 0.5,
                  max.overlaps = 5, size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  xlim(c(-3, 3)) + ylab("-log10QValue") + xlab("log2FC")
# dev.off()

