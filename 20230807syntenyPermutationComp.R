### is permutation important i.e do it for just genes based on aneuploidy
### (1) make table for all genes found in  both human and mouse
### (2) use the frequency * amplitude values for values
### (3) do permutation for human and mouse should be 30,000 x 300 and 30,000 by 30,000 x 50 matrix


### test task chunking since it's a lot of small tasks with previous permutation



# start <- Sys.time()
# humanOvRegionMat2 <- humanOvRegionMat[ , 4:ncol(humanOvRegionMat)]
# humanPermAmp <- NULL
# humanPermDel <- NULL
# humanPermAmpStat2 <- NULL
# humanPermDelStat2 <- NULL
# for (i in 1:10000) {
#   tmp <- apply(humanOvRegionMat2, 2, function(x) sample(x, 396))
#   tmpAmp <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
#   tmpDel <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
#   tmpAmpStat2 <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
#   tmpDelStat2 <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
#   
#   humanPermAmp <- cbind(humanPermAmp, tmpAmp)
#   humanPermDel <- cbind(humanPermDel, tmpDel)
#   humanPermAmpStat2 <- cbind(humanPermAmpStat2,  tmpAmpStat2)
#   humanPermDelStat2 <- cbind(humanPermDelStat2, tmpDelStat2)
# }
# print( Sys.time() - start )





comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

### works better in chunks 16.5 minutes vs 36 seconds
### check to see if the results are consistent they should be i.e check the perm tables



hg19biomartTable <- read.table("/home/kevhu/data/20230720hg19KnownCanBioMartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

mm10biomartTable <- read.table("/home/kevhu/data/20210203Mm10KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

geneNameDf <- read.table("/home/kevhu/data/20230315geneNameBiomart.txt", header = TRUE,
                         stringsAsFactors = FALSE, sep = "\t")


allMouseHomologsMm <- unique(geneNameDf$mmusculus_homolog_associated_gene_name)
allMouseHomologsHg <- unique(geneNameDf$external_gene_name[which(geneNameDf$mmusculus_homolog_associated_gene_name %in% allMouseHomologsMm)])

geneNameDf2 <- geneNameDf[-which(geneNameDf$mmusculus_homolog_associated_gene_name == ""), ]
geneNameDf2 <- geneNameDf2[which(geneNameDf2$mmusculus_homolog_associated_gene_name %in% allMouseHomologsMm),]
geneNameDf2 <- geneNameDf2[-which(duplicated(geneNameDf2$external_gene_name)),]

hg19biomartTable2 <- hg19biomartTable[which(hg19biomartTable$external_gene_name %in% allMouseHomologsHg), ]
hg19biomartTable2 <- hg19biomartTable2[-which(duplicated(hg19biomartTable2$external_gene_name)), ]
mm10biomartTable2 <- mm10biomartTable[which(mm10biomartTable$external_gene_name %in% allMouseHomologsMm), ]
mm10biomartTable2 <- mm10biomartTable2[-which(duplicated(mm10biomartTable2$external_gene_name)), ]

homologDf <- data.frame("human_gene" = hg19biomartTable2$external_gene_name, "h_chr" = hg19biomartTable2$chromosome_name,
                        "h_start" = hg19biomartTable2$exon_chrom_start, "h_end" = hg19biomartTable2$exon_chrom_end)
homologDf$mouse_gene <- geneNameDf2$mmusculus_homolog_associated_gene_name[match(homologDf$human_gene, geneNameDf2$external_gene_name)]
homologDf <- homologDf[-which(is.na(homologDf$mouse_gene)), ]
homologDf$m_chr <- mm10biomartTable2$chromosome_name[match(homologDf$mouse_gene, mm10biomartTable2$external_gene_name)]
homologDf$m_start <- mm10biomartTable2$exon_chrom_start[match(homologDf$mouse_gene, mm10biomartTable2$external_gene_name)]
homologDf$m_end <- mm10biomartTable2$exon_chrom_end[match(homologDf$mouse_gene, mm10biomartTable2$external_gene_name)]
homologDf <- homologDf[-which(is.na(homologDf$m_chr)), ]
homologDf <- homologDf[-which(homologDf$h_chr %in% c("X", "Y")),]
homologDf <- homologDf[-which(homologDf$m_chr %in% c("X", "Y")),]

ovAneuGrange <- GRanges(seqnames = broadBySampleGistic2$chrStripped,
                        IRanges(start = broadBySampleGistic2$start, end = broadBySampleGistic2$end))

tmp_allMouseAneuploidy <- allMouseAneuploidy
tmp_allMouseAneuploidy[tmp_allMouseAneuploidy == "none"] <- 0
tmp_allMouseAneuploidy[tmp_allMouseAneuploidy == "gain"] <- 1
tmp_allMouseAneuploidy[tmp_allMouseAneuploidy == "loss"] <- -1
tmp_allMouseAneuploidy <- apply(tmp_allMouseAneuploidy, 2, as.numeric)
rownames(tmp_allMouseAneuploidy) <- rownames(allMouseAneuploidy)

allMouseAneuploidy_chr <- data.frame(tmpChrBoundaries, t(tmp_allMouseAneuploidy), check.names = FALSE)
allMouseAneuploidy_chr_bprn <- allMouseAneuploidy_chr[, c(1:3 ,which(colnames(allMouseAneuploidy_chr) %in% annoTable2$mouse_id[which(annoTable2$geno == "BPRN")]))]

bprnAneuGrange <- GRanges(seqnames = str_remove(allMouseAneuploidy_chr_bprn$chr, "chr"),
                        IRanges(start = allMouseAneuploidy_chr_bprn$start,
                                end = allMouseAneuploidy_chr_bprn$end))

colnames(allMouseAneuploidy_chr) <- str_remove(colnames(allMouseAneuploidy_chr), "x.*")
allMouseAneuploidy_chr_coad <- allMouseAneuploidy_chr[, c(1:3 ,which(colnames(allMouseAneuploidy_chr) %in% annoTableCoad2$V6[which(annoTableCoad2$histology == "Adenocarcinoma")]))]


bprnAneuGrange <- GRanges(seqnames = str_remove(allMouseAneuploidy_chr_bprn$chr, "chr"),
                          IRanges(start = allMouseAneuploidy_chr_bprn$start,
                                  end = allMouseAneuploidy_chr_bprn$end))

mmCoadAneuGrange <- GRanges(seqnames = str_remove(allMouseAneuploidy_chr_coad$chr, "chr"),
                          IRanges(start = allMouseAneuploidy_chr_coad$start,
                                  end = allMouseAneuploidy_chr_coad$end))

### used to filter out the set of genes not found in the aneuploidy ranges

hg19HomologGranges <- GRanges(seqnames = homologDf$h_chr, IRanges(start =  homologDf$h_start, end = homologDf$h_end))
mm10HomologGranges <- GRanges(seqnames = homologDf$m_chr, IRanges(start =  homologDf$m_start, end = homologDf$m_end))

tmpVec <- 1:15425
badMouseIdx <- tmpVec[-which(tmpVec %in% queryHits(findOverlaps(mm10HomologGranges, bprnAneuGrange)))]
badHumanIdx <- tmpVec[-which(tmpVec %in% queryHits(findOverlaps(hg19HomologGranges, ovAneuGrange)))]

homologDf <- homologDf[-c(badMouseIdx, badHumanIdx), ]

hg19HomologGranges <- GRanges(seqnames = homologDf$h_chr, IRanges(start =  homologDf$h_start, end = homologDf$h_end))
mm10HomologGranges <- GRanges(seqnames = homologDf$m_chr, IRanges(start =  homologDf$m_start, end = homologDf$m_end))

### creating gene x sample matrices for cnr

broadBySampleGistic3 <- broadBySampleGistic2[, 4:ncol(broadBySampleGistic2)]
tmpHitsOv <- subjectHits(findOverlaps(hg19HomologGranges, ovAneuGrange))
ovGeneCnDf <- NULL
for (i in 1:ncol(broadBySampleGistic3)) {
  ovGeneCnDf <- cbind(ovGeneCnDf, broadBySampleGistic3[[i]][tmpHitsOv])
}
rownames(ovGeneCnDf) <- homologDf$human_gene
colnames(ovGeneCnDf) <- colnames(broadBySampleGistic3)


coadAneuGrange <- GRanges(seqnames = broadBySampleGisticCoad2$chrStripped,
                          IRanges(start = broadBySampleGisticCoad2$start,
                                  end = broadBySampleGisticCoad2$end))

broadBySampleGisticCoad3 <- broadBySampleGisticCoad2[, 4:ncol(broadBySampleGisticCoad2)]
tmpHitsCoad <- subjectHits(findOverlaps(hg19HomologGranges, coadAneuGrange))
coadGeneCnDf <- NULL
for (i in 1:ncol(broadBySampleGisticCoad3)) {
  coadGeneCnDf <- cbind(coadGeneCnDf, broadBySampleGistic3[[i]][tmpHitsCoad])
}
rownames(coadGeneCnDf) <- homologDf$human_gene
colnames(coadGeneCnDf) <- colnames(broadBySampleGisticCoad3)


allMouseAneuploidy_chr_bprn2 <- allMouseAneuploidy_chr_bprn[, 4:ncol(allMouseAneuploidy_chr_bprn)]
tmpHitsBprn <- subjectHits(findOverlaps(mm10HomologGranges, bprnAneuGrange))
bprnGeneCnDf <- NULL
for (i in 1:ncol(allMouseAneuploidy_chr_bprn2)) {
  bprnGeneCnDf <- cbind(bprnGeneCnDf, allMouseAneuploidy_chr_bprn2[[i]][tmpHitsBprn])
}

rownames(bprnGeneCnDf) <- homologDf$mouse_gene
colnames(bprnGeneCnDf) <- colnames(allMouseAneuploidy_chr_bprn2)

allMouseAneuploidy_chr_coad2 <- allMouseAneuploidy_chr_coad[, 4:ncol(allMouseAneuploidy_chr_coad)]
tmpHitsMmCoad <- subjectHits(findOverlaps(mm10HomologGranges, mmCoadAneuGrange))
mmCoadGeneCnDf <- NULL
for (i in 1:ncol(allMouseAneuploidy_chr_coad2)) {
  mmCoadGeneCnDf <- cbind(mmCoadGeneCnDf, allMouseAneuploidy_chr_coad2[[i]][tmpHitsMmCoad])
}

rownames(mmCoadGeneCnDf) <- homologDf$mouse_gene
colnames(mmCoadGeneCnDf) <- colnames(allMouseAneuploidy_chr_coad2)

### test perm 

cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
hgscPerm <- foreach(i=1:25,
                    .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                      tmpAmpAll <- NULL
                      tmpDelAll <- NULL
                      tmpAmpFreq <- NULL
                      tmpDelFreq <- NULL
                      j <- 0
                      while (j < 400) {
                        tmp <- apply(ovGeneCnDf, 2, function(x) sample(x, nrow(ovGeneCnDf)))
                        ampFreq <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
                        delFreq <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
                        tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
                        tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
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

ovGeneAmp <- do.call(cbind, hgscPerm[[1]])
ovGeneDel <- do.call(cbind, hgscPerm[[2]])
ovGeneAmpFreq <- do.call(cbind, hgscPerm[[3]])
ovGeneDelFreq <- do.call(cbind, hgscPerm[[4]])


cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
bprnPerm <- foreach(i=1:25,
                    .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                      tmpAmpAll <- NULL
                      tmpDelAll <- NULL
                      tmpAmpFreq <- NULL
                      tmpDelFreq <- NULL
                      j <- 0
                      while (j < 400) {
                        tmp <- apply(bprnGeneCnDf, 2, function(x) sample(x, nrow(bprnGeneCnDf)))
                        ampFreq <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
                        delFreq <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
                        tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
                        tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
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

bprnGeneAmp <- do.call(cbind, bprnPerm[[1]])
bprnGeneDel <- do.call(cbind, bprnPerm[[2]])
bprnGeneAmpFreq <- do.call(cbind, bprnPerm[[3]])
bprnGeneDelFreq <- do.call(cbind, bprnPerm[[4]])

cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
coadPerm <- foreach(i=1:25,
                    .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                      tmpAmpAll <- NULL
                      tmpDelAll <- NULL
                      tmpAmpFreq <- NULL
                      tmpDelFreq <- NULL
                      j <- 0
                      while (j < 400) {
                        tmp <- apply(coadGeneCnDf, 2, function(x) sample(x, nrow(coadGeneCnDf)))
                        ampFreq <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
                        delFreq <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
                        tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
                        tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
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

coadGeneAmp <- do.call(cbind, coadPerm[[1]])
coadGeneDel <- do.call(cbind, coadPerm[[2]])
coadGeneAmpFreq <- do.call(cbind, coadPerm[[3]])
coadGeneDelFreq <- do.call(cbind, coadPerm[[4]])


cl <- makeCluster(25)
registerDoParallel(cl)
start <- Sys.time()
mmCoadPerm <- foreach(i=1:25,
                      .combine = 'comb', .multicombine = TRUE, .init = list(list(), list(), list(), list())) %dopar% {
                        tmpAmpAll <- NULL
                        tmpDelAll <- NULL
                        tmpAmpFreq <- NULL
                        tmpDelFreq <- NULL
                        j <- 0
                        while (j < 400) {
                          tmp <- apply(mmCoadGeneCnDf, 2, function(x) sample(x, nrow(mmCoadGeneCnDf)))
                          ampFreq <- apply(tmp, 1, function(x) length(which(x > 0))/length(x))
                          delFreq <- apply(tmp, 1, function(x) length(which(x < 0))/length(x))
                          tmpAmp <- apply(tmp, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
                          tmpDel <- apply(tmp, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
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

mmCoadGeneAmp <- do.call(cbind, mmCoadPerm[[1]])
mmCoadGeneDel <- do.call(cbind, mmCoadPerm[[2]])
mmCoadGeneAmpFreq <- do.call(cbind, mmCoadPerm[[3]])
mmCoadGeneDelFreq <- do.call(cbind, mmCoadPerm[[4]])

# save.image(file = "/mnt/DATA4/test_nextflow/20230809allSyntenyEnivorn.RData")

### below is the permutation pvalue analysis
###


ovGeneAmp <- do.call(cbind, hgscPerm[[1]])
ovGeneDel <- do.call(cbind, hgscPerm[[2]])
ovGeneAmpFreq <- do.call(cbind, hgscPerm[[3]])
ovGeneDelFreq <- do.call(cbind, hgscPerm[[4]])

ampSynTable_hgsc <- data.frame("gene" = rownames(ovGeneCnDf))
delSynTable_hgsc <- data.frame("gene" = rownames(ovGeneCnDf))
ampSynTable_hgsc$h_freq <- apply(ovGeneCnDf, 1, function(x) length(which(x > 0))/length(x))
ampSynTable_hgsc$m_freq <- apply(bprnGeneCnDf, 1, function(x) length(which(x > 0))/length(x))
delSynTable_hgsc$h_freq <- apply(ovGeneCnDf, 1, function(x) length(which(x < 0))/length(x))
delSynTable_hgsc$m_freq <- apply(bprnGeneCnDf, 1, function(x) length(which(x < 0))/length(x))
ampSynTable_hgsc$h_stat <- apply(ovGeneCnDf, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
ampSynTable_hgsc$m_stat <- apply(bprnGeneCnDf, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
delSynTable_hgsc$h_stat <- apply(ovGeneCnDf, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
delSynTable_hgsc$m_stat <- apply(bprnGeneCnDf, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))

which(is.nan(ampSynTable_hgsc$m_freq))

ampSynTable_hgsc$h_freq[which(ampSynTable_hgsc$h_freq < 0.1)] <- 0
ampSynTable_hgsc$m_freq[which(ampSynTable_hgsc$m_freq < 0.1)] <- 0
delSynTable_hgsc$h_freq[which(delSynTable_hgsc$h_freq < 0.1)] <- 0
delSynTable_hgsc$m_freq[which(delSynTable_hgsc$m_freq < 0.1)] <- 0
ampSynTable_hgsc$h_stat[which(ampSynTable_hgsc$h_freq < 0.1)] <- 0
ampSynTable_hgsc$m_stat[which(ampSynTable_hgsc$m_freq < 0.1)] <- 0
delSynTable_hgsc$h_stat[which(delSynTable_hgsc$h_freq < 0.1)] <- 0
delSynTable_hgsc$m_stat[which(delSynTable_hgsc$m_freq < 0.1)] <- 0

humanPermAmp2_hgsc <- ovGeneAmpFreq
mousePermAmp2_hgsc <- bprnGeneAmpFreq
humanPermDel2_hgsc <- ovGeneDelFreq
mousePermDel2_hgsc <- bprnGeneDelFreq
mousePermAmpStat2_hgsc <- bprnGeneAmp
mousePermDelStat2_hgsc <- bprnGeneDel
humanPermAmpStat2v2 <- ovGeneAmp
humanPermDelStat2v2 <- ovGeneDel

humanPermAmp2_hgsc[humanPermAmp2_hgsc < 0.1] <- 0
mousePermAmp2_hgsc[mousePermAmp2_hgsc < 0.1] <- 0
humanPermDel2_hgsc[humanPermDel2_hgsc < 0.1] <- 0
mousePermDel2_hgsc[mousePermDel2_hgsc < 0.1] <- 0
humanPermAmpStat2v2[humanPermAmp2 < 0.1] <- 0
humanPermDelStat2v2[humanPermDel2 < 0.1] <- 0
mousePermAmpStat2_hgsc[mousePermAmp2_hgsc < 0.1] <- 0
mousePermDelStat2_hgsc[mousePermDel2_hgsc < 0.1] <- 0

permStatResTbl_hgsc <- NULL
for (i in 1:15424) {
  tmpAmp <- humanPermAmp2_hgsc[i, ] + mousePermAmp2_hgsc[i, ]
  tmpDel <- humanPermDel2_hgsc[i, ] + mousePermDel2_hgsc[i, ]
  
  tmpAmpStat2 <- humanPermAmpStat2v2[i, ] + mousePermAmpStat2_hgsc[i, ]
  tmpDelStat2 <- humanPermDelStat2v2[i, ] + mousePermDelStat2_hgsc[i, ]
  
  if (ampSynTable_hgsc$h_freq[i] == 0 | ampSynTable_hgsc$m_freq[i] == 0) {
    tmpAmpZ <- 0
    tmpAmpZStat2 <- 0
  } else{
    tmpAmpZ <- ((ampSynTable_hgsc$h_freq[i] + ampSynTable_hgsc$m_freq[i]) - mean(tmpAmp, na.rm = TRUE))/sd(tmpAmp, na.rm = TRUE)
    tmpAmpZStat2 <- ((ampSynTable$h_stat2[i] + ampSynTable$m_stat2[i]) - mean(tmpAmpStat2, na.rm = TRUE))/sd(tmpAmpStat2, na.rm = TRUE)
  }
  
  if (delSynTable_hgsc$h_freq[i] == 0 | delSynTable_hgsc$m_freq[i] == 0) {
    tmpDelZ <- 0
    tmpDelZStat2 <- 0
  } else{
    tmpDelZ <- ((delSynTable_hgsc$h_freq[i] + delSynTable_hgsc$m_freq[i]) -  mean(tmpDel, na.rm = TRUE))/sd(tmpDel, na.rm = TRUE)
    tmpDelZStat2 <- ((delSynTable$h_stat2[i] + delSynTable$m_stat2[i]) - mean(tmpDelStat2, na.rm = TRUE))/sd(tmpDelStat2, na.rm = TRUE)
  }
  
  permStatResTbl_hgsc <- rbind(permStatResTbl_hgsc, data.frame("ampZ" = tmpAmpZ, "delZ" = tmpDelZ, "ampStat2Z" = tmpAmpZStat2, "delStat2Z" = tmpDelZStat2))
}


permStatResTbl_hgsc$ampStat2Z[which(is.na(permStatResTbl_hgsc$ampStat2Z))] <- 0
permStatResTbl_hgsc$delStat2Z[which(is.na(permStatResTbl_hgsc$delStat2Z))] <- 0
permStatResTbl_hgsc$ampP2 <- pnorm(q=permStatResTbl_hgsc$ampStat2Z, lower.tail=FALSE)
permStatResTbl_hgsc$delP2 <- pnorm(q=permStatResTbl_hgsc$delStat2Z, lower.tail=FALSE)

permStatResTbl_hgsc_amp <- permStatResTbl_hgsc
permStatResTbl_hgsc_del <- permStatResTbl_hgsc

permStatResTbl_hgsc_amp$ampPAdj2 <- p.adjust(permStatResTbl_hgsc_amp$ampP2, method = "BH", n = (length(permStatResTbl_hgsc_amp$ampP2) + length(permStatResTbl_hgsc_del$ampP2)))
permStatResTbl_hgsc_del$delPAdj2 <- p.adjust(permStatResTbl_hgsc_del$delP2, method = "BH", n = (length(permStatResTbl_hgsc_amp$ampP2) + length(permStatResTbl_hgsc_del$ampP2)))

tmpSynTable2_hgsc <- data.frame("Gene" = homologDf$human_gene, 
                                    "chr" = homologDf$h_chr, "start" = homologDf$h_start, 
                                    "end" = homologDf$h_end)

tmpSynTable2_hgsc$ampPAdj2 <- permStatResTbl_hgsc_amp$ampPAdj2
tmpSynTable2_hgsc$delPAdj2 <- permStatResTbl_hgsc_del$delPAdj2

hgsc_hgsc_arm_allSynTable_amp <- tmpSynTable2_hgsc_amp[which(permStatResTbl_hgsc_amp$ampPAdj2 < 0.05), ]
hgsc_hgsc_arm_allSynTable_del <- tmpSynTable2_hgsc_del[which(permStatResTbl_hgsc_del$delPAdj2 < 0.05), ]

cancerGeneCensusGr <- GRanges(seqnames = cancerGeneCensus$chr,
                              IRanges(start = as.numeric(cancerGeneCensus$start), end = as.numeric(cancerGeneCensus$end)))
cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_hgsc_arm_allSynTable_amp$h_chr, IRanges(hgsc_hgsc_arm_allSynTable_amp$h_start, end = hgsc_hgsc_arm_allSynTable_amp$h_end))
                                          ,cancerGeneCensusGr)),]
cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = hgsc_hgsc_arm_allSynTable_del$h_chr, IRanges(hgsc_hgsc_arm_allSynTable_del$h_start, end = hgsc_hgsc_arm_allSynTable_del$h_end))
                                          ,cancerGeneCensusGr)),]



### same but for coad

ampSynTable_coad <- data.frame("gene" = rownames(coadGeneCnDf))
delSynTable_coad <- data.frame("gene" = rownames(coadGeneCnDf))
ampSynTable_coad$h_freq <- apply(coadGeneCnDf, 1, function(x) length(which(x > 0))/length(x))
ampSynTable_coad$m_freq <- apply(mmCoadGeneCnDf, 1, function(x) length(which(x > 0))/length(x))
delSynTable_coad$h_freq <- apply(coadGeneCnDf, 1, function(x) length(which(x < 0))/length(x))
delSynTable_coad$m_freq <- apply(mmCoadGeneCnDf, 1, function(x) length(which(x < 0))/length(x))
ampSynTable_coad$h_stat <- apply(coadGeneCnDf, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
ampSynTable_coad$m_stat <- apply(mmCoadGeneCnDf, 1, function(x) mean(x[which(x > 0)]) * length(which(x > 0))/length(x))
delSynTable_coad$h_stat <- apply(coadGeneCnDf, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))
delSynTable_coad$m_stat <- apply(mmCoadGeneCnDf, 1, function(x) abs(mean(x[which(x < 0)])) * length(which(x < 0))/length(x))


ampSynTable_coad$h_freq[which(ampSynTable_coad$h_freq < 0.1)] <- 0
ampSynTable_coad$m_freq[which(ampSynTable_coad$m_freq < 0.1)] <- 0
delSynTable_coad$h_freq[which(delSynTable_coad$h_freq < 0.1)] <- 0
delSynTable_coad$m_freq[which(delSynTable_coad$m_freq < 0.1)] <- 0
ampSynTable_coad$h_stat[which(ampSynTable_coad$h_freq < 0.1)] <- 0
ampSynTable_coad$m_stat[which(ampSynTable_coad$m_freq < 0.1)] <- 0
delSynTable_coad$h_stat[which(delSynTable_coad$h_freq < 0.1)] <- 0
delSynTable_coad$m_stat[which(delSynTable_coad$m_freq < 0.1)] <- 0

humanPermAmp2_coad <- coadGeneAmpFreq
mousePermAmp2_coad <- mmCoadGeneAmpFreq
humanPermDel2_coad <- coadGeneDelFreq
mousePermDel2_coad <- mmCoadGeneDelFreq
mousePermAmpStat2_coad <- mmCoadGeneAmp
mousePermDelStat2_coad <- mmCoadGeneDel
humanPermAmpStat2v2 <- coadGeneAmp
humanPermDelStat2v2 <- coadGeneDel

humanPermAmp2_coad[humanPermAmp2_coad < 0.1] <- 0
mousePermAmp2_coad[mousePermAmp2_coad < 0.1] <- 0
humanPermDel2_coad[humanPermDel2_coad < 0.1] <- 0
mousePermDel2_coad[mousePermDel2_coad < 0.1] <- 0
humanPermAmpStat2v2[humanPermAmp2 < 0.1] <- 0
humanPermDelStat2v2[humanPermDel2 < 0.1] <- 0
mousePermAmpStat2_coad[mousePermAmp2_coad < 0.1] <- 0
mousePermDelStat2_coad[mousePermDel2_coad < 0.1] <- 0

permStatResTbl_coad <- NULL
for (i in 1:15424) {
  tmpAmp <- humanPermAmp2_coad[i, ] + mousePermAmp2_coad[i, ]
  tmpDel <- humanPermDel2_coad[i, ] + mousePermDel2_coad[i, ]
  
  tmpAmpStat2 <- humanPermAmpStat2v2[i, ] + mousePermAmpStat2_coad[i, ]
  tmpDelStat2 <- humanPermDelStat2v2[i, ] + mousePermDelStat2_coad[i, ]
  
  if (ampSynTable_coad$h_freq[i] == 0 | ampSynTable_coad$m_freq[i] == 0) {
    tmpAmpZ <- 0
    tmpAmpZStat2 <- 0
  } else{
    tmpAmpZ <- ((ampSynTable_coad$h_freq[i] + ampSynTable_coad$m_freq[i]) - mean(tmpAmp, na.rm = TRUE))/sd(tmpAmp, na.rm = TRUE)
    tmpAmpZStat2 <- ((ampSynTable$h_stat2[i] + ampSynTable$m_stat2[i]) - mean(tmpAmpStat2, na.rm = TRUE))/sd(tmpAmpStat2, na.rm = TRUE)
  }
  
  if (delSynTable_coad$h_freq[i] == 0 | delSynTable_coad$m_freq[i] == 0) {
    tmpDelZ <- 0
    tmpDelZStat2 <- 0
  } else{
    tmpDelZ <- ((delSynTable_coad$h_freq[i] + delSynTable_coad$m_freq[i]) -  mean(tmpDel, na.rm = TRUE))/sd(tmpDel, na.rm = TRUE)
    tmpDelZStat2 <- ((delSynTable$h_stat2[i] + delSynTable$m_stat2[i]) - mean(tmpDelStat2, na.rm = TRUE))/sd(tmpDelStat2, na.rm = TRUE)
  }
  
  permStatResTbl_coad <- rbind(permStatResTbl_coad, data.frame("ampZ" = tmpAmpZ, "delZ" = tmpDelZ, "ampStat2Z" = tmpAmpZStat2, "delStat2Z" = tmpDelZStat2))
}

permStatResTbl_coad$ampStat2Z[which(is.na(permStatResTbl_coad$ampStat2Z))] <- 0
permStatResTbl_coad$delStat2Z[which(is.na(permStatResTbl_coad$delStat2Z))] <- 0
permStatResTbl_coad$ampP2 <- pnorm(q=permStatResTbl_coad$ampStat2Z, lower.tail=FALSE)
permStatResTbl_coad$delP2 <- pnorm(q=permStatResTbl_coad$delStat2Z, lower.tail=FALSE)

permStatResTbl_coad_amp <- permStatResTbl_coad
permStatResTbl_coad_del <- permStatResTbl_coad

permStatResTbl_coad_amp$ampPAdj2 <- p.adjust(permStatResTbl_coad_amp$ampP2, method = "BH", n = (length(permStatResTbl_coad_amp$ampP2) + length(permStatResTbl_coad_del$ampP2)))
permStatResTbl_coad_del$delPAdj2 <- p.adjust(permStatResTbl_coad_del$delP2, method = "BH", n = (length(permStatResTbl_coad_amp$ampP2) + length(permStatResTbl_coad_del$ampP2)))

tmpSynTable2_coad <- data.frame("Gene" = homologDf$human_gene, 
                                "chr" = homologDf$h_chr, "start" = homologDf$h_start, 
                                "end" = homologDf$h_end)

tmpSynTable2_coad$ampPAdj2 <- permStatResTbl_coad_amp$ampPAdj2
tmpSynTable2_coad$delPAdj2 <- permStatResTbl_coad_del$delPAdj2

coad_coad_arm_allSynTable_amp <- tmpSynTable2_coad_amp[which(permStatResTbl_coad_amp$ampPAdj2 < 0.05), ]
coad_coad_arm_allSynTable_del <- tmpSynTable2_coad_del[which(permStatResTbl_coad_del$delPAdj2 < 0.05), ]

cancerGeneCensusGr <- GRanges(seqnames = cancerGeneCensus$chr,
                              IRanges(start = as.numeric(cancerGeneCensus$start), end = as.numeric(cancerGeneCensus$end)))
cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = coad_coad_arm_allSynTable_amp$h_chr, IRanges(coad_coad_arm_allSynTable_amp$h_start, end = coad_coad_arm_allSynTable_amp$h_end))
                                          ,cancerGeneCensusGr)),]
cancerGeneCensus[subjectHits(findOverlaps(GRanges(seqnames = coad_coad_arm_allSynTable_del$h_chr, IRanges(coad_coad_arm_allSynTable_del$h_start, end = coad_coad_arm_allSynTable_del$h_end))
                                          ,cancerGeneCensusGr)),]




### find pathway overlaps with general hypergeomtric test for signficance
### simple hyper geometric test here: https://montilab.github.io/BS831/articles/docs/HyperEnrichment.html
### be sure to do the permutation test with coad too and replace nan with 0's
### then do the general large CNAs and aneuploidy stats



human_hallmarks <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/MSigDb/h.all.v7.4.symbols.gmt",
                              sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

human_hallmarks_list <- NULL
for (i in 1:nrow(human_hallmarks)) {
  tmpList <- list(human_hallmarks[i, 3:ncol(human_hallmarks)])
  names(tmpList) <- human_hallmarks[i, 1]
  human_hallmarks_list[[i]] <- tmpList
}

### the hyper geometric test is done correctly
### maybe just too much signal, look at the colorectal data next
### could indicate one of two things no pathways altered -> aneuploidues
### target singlar oncogenes or TSGs

bprnHgscHyper <- pathwayHyper(syntenyGains = hgsc_bprn_arm_allSynTable_amp, 
                              syntenyLosses = hgsc_bprn_arm_allSynTable_del,
                              pathwayList = human_hallmarks_list)


coadHgscHyper <- pathwayHyper(syntenyGains = adenoCar_arm_allSynTable_amp, 
                              syntenyLosses = adenoCar_arm_allSynTable_del,
                              pathwayList = human_hallmarks_list)

