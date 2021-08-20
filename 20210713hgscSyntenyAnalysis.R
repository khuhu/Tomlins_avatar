### /home/kevhu/scripts/20210629syntenyDocumentation.R used as reference

source("/home/kevhu/scripts/20210802syntenyFunctions.R")
synteny_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_syntenyDf.txt",sep = "\t", stringsAsFactors = FALSE, header = TRUE)

cyto_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_cytoDf.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = FALSE)
tcga_seg <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/TCGA_mastercalls.abs_segtabs.fixed.txt",sep = "\t", stringsAsFactors = FALSE, header = TRUE)

tcga_ploidy <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",sep = "\t", stringsAsFactors = FALSE, header = TRUE)

segResults_1 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/segResults.txt",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
segResults_2 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/segResults.txt",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
segResults_3 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/segResults.txt",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
segResults_4 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-76-MG_test1_255_185/segResults.txt",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")

cnCalls_1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-138-MG_cho_20210621_354_343/cnMatrix_gene.txt",
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t")
cnCalls_2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/cnMatrix_gene.txt",
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t")
cnCalls_3 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-142-MG_cho_20210701_357_353/cnMatrix_gene.txt",
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t")
cnCalls_4 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-76-MG_test1_255_185/cnMatrix_gene.txt",
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t")

zscore_gc_oe_ratios <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/gcCorrectedCounts_matrix.txt", sep = "\t",
                                  stringsAsFactors = FALSE, header = TRUE)

mouseNormal <- c("MG_17X49", "MG_18X50", "MG_23X55", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

mm10ChromSize <- read.table("/mnt/DATA5/tmp/kev/misc/20210809mm10ChromSizes.txt", sep = "\t", stringsAsFactors = FALSE,
                            header = TRUE)

hg38biomartTable <- read.table("/home/kevhu/data/20201030hg38KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

mm10biomartTable <- read.table("/home/kevhu/data/20201030Mm10KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

geneNameDf <- read.table("/home/kevhu/data/20201021geneNameBiomart.txt", header = TRUE,
                         stringsAsFactors = FALSE, sep = "\t")


human_hallmarks <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/MSigDb/h.all.v7.4.symbols.gmt",
                              sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

mouseBedFile <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.gc.bed", stringsAsFactors = FALSE, sep = "\t",
                           header = FALSE)


hgsc_anno <- read.table("/mnt/DATA5/tmp/kev/misc/PATIENT_DATA_oncoprint_hgsc_Tp53_Brca.tsv", sep = "\t",
           stringsAsFactors = FALSE, header = TRUE)

human_arm <- read.table("/mnt/DATA5/tmp/kev/misc/20210730hg19chromSizeTable.txt", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)
mouse_arm <- read.table("/mnt/DATA5/tmp/kev/misc/20210713mouse_arm_syn.txt", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)

finalTable <- read.table("/mnt/DATA5/tmp/kev/misc/20210628panCancerAllSegReduced_hg19.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")

annoTable <- read.table("/mnt/DATA5/tmp/kev/misc/20210713allHgscMouseAnno_noCrisp.txt", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)

human_hallmarks <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/MSigDb/h.all.v7.4.symbols.gmt",
                              sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

hgsc_anno2 <- hgsc_anno[4:nrow(hgsc_anno),3:ncol(hgsc_anno)]
rownames(hgsc_anno2) <- paste0(hgsc_anno$track_name[4:nrow(hgsc_anno)],
                               "_", hgsc_anno$track_type[4:nrow(hgsc_anno)])

hgsc_anno2 <- data.frame(t(hgsc_anno2), stringsAsFactors = FALSE)

apply(hgsc_anno2, 2, unique)

tp53 <- which(hgsc_anno2$TP53_CNA == "homdel_rec" | hgsc_anno2$TP53_MUTATIONS == "Truncating mutation (putative driver)" |
        hgsc_anno2$TP53_MUTATIONS == "Inframe Mutation (putative driver)" | hgsc_anno2$TP53_MUTATIONS == "Missense Mutation (putative driver)" |
        hgsc_anno2$TP53_MUTATIONS == "splice_rec")

hgsc_anno2_tp53 <- hgsc_anno2[tp53,]
hgsc_anno2_tp53_brca12 <- hgsc_anno2_tp53[which(hgsc_anno2_tp53$BRCA1_CNA == "homdel_rec" | hgsc_anno2_tp53$BRCA2_CNA == "homdel_rec" |
                                                  hgsc_anno2_tp53$BRCA1_MUTATIONS == "Truncating mutation (putative driver)" |
                                                  hgsc_anno2_tp53$BRCA1_MUTATIONS == "Missense Mutation (putative driver)" |
                                                  hgsc_anno2_tp53$BRCA1_MUTATIONS == "splice_rec" |
                                                  hgsc_anno2_tp53$BRCA2_MUTATIONS == "Truncating mutation (putative driver)" |
                                                  hgsc_anno2_tp53$BRCA2_MUTATIONS == "Missense Mutation (putative driver)"),]

#all_segRes <- rbind(segResults_1, segResults_2, segResults_3, segResults_4)

all_cnCalls <- cbind(cnCalls_1[grep("Del", cnCalls_1$Gene),], cnCalls_2[grep("Del", cnCalls_2$Gene),2:ncol(cnCalls_2)],
                     cnCalls_3[grep("Del", cnCalls_3$Gene),2:ncol(cnCalls_3)], cnCalls_4[grep("Del", cnCalls_4$Gene),2:ncol(cnCalls_4)])
all_cnCalls <- all_cnCalls[,-which(duplicated(colnames(all_cnCalls)))]

all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                             ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                              end = zscore_gc_oe_ratios$EndPos))


### like Scott said, low tc values used for correction can overcorrect - crazy high amp
tc <- 1 - apply(all_cnCalls[grep("Del", all_cnCalls$Gene), 2:ncol(all_cnCalls)], 2, min)
tc[which(names(tc) %in% mouseNormal)] <- 1
tc[tc < 0.50] <- .50
names(tc) <- str_remove(names(tc), "_MG_X.*")
names(tc) <- str_remove(str_remove(names(tc), "^X"), "_X.*")
names(tc) <- tolower(str_remove_all(str_remove_all(str_remove(names(tc), "X.*"), "_"), "\\."))
names(tc)  <- str_remove(names(tc), "o")

tcDf <- data.frame("sample" = names(tc), "tc" = tc, stringsAsFactors = FALSE)
# write.table(tcDf, "/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt", sep = "\t",
#               quote = FALSE, row.names = FALSE, col.names = TRUE)


finalTable$length <- finalTable$End - finalTable$Start
finalTable2 <- finalTable[which(finalTable$Chromosome %in% c(1:22)),]
finalTable2 <- finalTable2[order(finalTable2$Chromosome, finalTable2$Start),]
# armRes <- separateSegments_tcgaV3(finalTable2, tcga_ploidy)
# save(armRes, file = "/mnt/DATA5/tmp/kev/misc/20210802tcgaArmRes.Robj")
load("/mnt/DATA5/tmp/kev/misc/20210802tcgaArmRes.Robj")

armTable <- do.call(rbind, armRes[[1]])
cnaTable <- do.call(rbind, armRes[[2]])
countTable <- do.call(rbind, armRes[[3]])


tcga_arm <- armTable
tcga_cna <- cnaTable

hgsc_anno2_tp53_names <- str_replace_all(rownames(hgsc_anno2_tp53), "\\.", "\\-")
hgsc_anno2_tp53_brca12_names <- str_replace_all(rownames(hgsc_anno2_tp53_brca12), "\\.", "\\-")
hgsc_anno2_tp53_names <- hgsc_anno2_tp53_names[-which(hgsc_anno2_tp53_names %in% hgsc_anno2_tp53_brca12_names)]


hgsc_arm <- tcga_arm[grep(paste(hgsc_anno2_tp53_names, collapse = "|"), tcga_arm$Sample),]
hgsc_cna <- tcga_cna[grep(paste(hgsc_anno2_tp53_names, collapse = "|"), tcga_arm$Sample),]

brca_arm <- tcga_arm[grep(paste(hgsc_anno2_tp53_brca12_names, collapse = "|"), tcga_arm$Sample),]
brca_cna <- tcga_cna[grep(paste(hgsc_anno2_tp53_brca12_names, collapse = "|"), tcga_cna$Sample),]



hgsc_arm2 <- data.frame("sampleID" = hgsc_arm$Sample, "chrom" = hgsc_arm$Chromosome,
                        "start.pos" = hgsc_arm$Start, "end.pos" = hgsc_arm$End,
                        "n.probes" = NA, "mean" = hgsc_arm$newTotalCN, stringsAsFactors = FALSE)

hgsc_cna2 <- data.frame("sampleID" = hgsc_cna$Sample, "chrom" = hgsc_cna$Chromosome,
                        "start.pos" = hgsc_cna$Start, "end.pos" = hgsc_cna$End,
                        "n.probes" = NA, "mean" = hgsc_cna$newTotalCN, stringsAsFactors = FALSE)


brca_arm2 <- data.frame("sampleID" = brca_arm$Sample, "chrom" = brca_arm$Chromosome,
                       "start.pos" = brca_arm$Start, "end.pos" = brca_arm$End,
                       "n.probes" = NA, "mean" = brca_arm$newTotalCN, stringsAsFactors = FALSE)
brca_cna2 <- data.frame("sampleID" = brca_cna$Sample, "chrom" = brca_cna$Chromosome,
                       "start.pos" = brca_cna$Start, "end.pos" = brca_cna$End,
                       "n.probes" = NA, "mean" = brca_cna$newTotalCN, stringsAsFactors = FALSE)



hgsc_arm_freq_out <- getFreqData(hgsc_arm2)
ampDels_hgsc_arm <- ampsDels(hgsc_arm_freq_out)
hgsc_arm_amp_bed <- reducingFreqBed(ampDels_hgsc_arm[[1]], ampDels_hgsc_arm[[2]])
hgsc_arm_del_bed <- reducingFreqBed(ampDels_hgsc_arm[[3]], ampDels_hgsc_arm[[4]])

hgsc_cna_freq_out <- getFreqData(hgsc_cna2)
ampDels_hgsc_cna <- ampsDels(hgsc_cna_freq_out)
hgsc_cna_amp_bed <- reducingFreqBed(ampDels_hgsc_cna[[1]], ampDels_hgsc_cna[[2]])
hgsc_cna_del_bed <- reducingFreqBed(ampDels_hgsc_cna[[3]], ampDels_hgsc_cna[[4]])


brca_arm_freq_out <- getFreqData(brca_arm2)
ampDels_brca_arm <- ampsDels(brca_arm_freq_out)
brca_arm_amp_bed <- reducingFreqBed(ampDels_brca_arm[[1]], ampDels_brca_arm[[2]])
brca_arm_del_bed <- reducingFreqBed(ampDels_brca_arm[[3]], ampDels_brca_arm[[4]])

brca_cna_freq_out <- getFreqData(brca_cna2)
ampDels_brca_cna <- ampsDels(brca_cna_freq_out)
brca_cna_amp_bed <- reducingFreqBed(ampDels_brca_cna[[1]], ampDels_brca_cna[[2]])
brca_cna_del_bed <- reducingFreqBed(ampDels_brca_cna[[3]], ampDels_brca_cna[[4]])


### mouse data - need to filt out non-hgsc samples
### make sure sampels used are of good quality too
###

all_segRes <- rbind(segResults_1, segResults_2, segResults_3, segResults_4)
all_segRes <- all_segRes[-which(all_segRes$ID %in% mouseNormal),]
all_segRes$ID <- str_remove(all_segRes$ID, "_MG_X.*")
all_segRes$ID <- str_remove(str_remove(all_segRes$ID, "^X"), "_X.*")
all_segRes$ID <- tolower(str_remove_all(str_remove_all(str_remove(all_segRes$ID, "X.*"), "_"), "\\."))
all_segRes$ID  <- str_remove(all_segRes$ID, "o")
matchingTc <- tcDf$tc[match(tolower(all_segRes$ID), tcDf$sample)]
all_segRes$seg.mean[which(!is.na(matchingTc))] <- all_segRes$seg.mean[which(!is.na(matchingTc))]/matchingTc[which(!is.na(matchingTc))]
all_segRes$chrom <- as.numeric(str_replace(all_segRes$chrom, "23", "20"))

annoTable$mouse_id <- str_remove(str_remove(str_remove_all(tolower(annoTable$mouse_id), "_"),"-"), "o")
annoTable$type <- tolower(annoTable$type)
annoTable$type[142] <- "ehgsc"
annoTable$mouse_id <- str_remove(str_replace_all(annoTable$mouse_id, " ", ""), "x.*")
annoTable2 <- annoTable[which(annoTable$type %in% c("hgsc", "mmmt")),]


unique(all_segRes$ID)[which(unique(all_segRes$ID) %in% annoTable$mouse_id[which(annoTable$geno %in% c("BPRN", "BPN", "BPP"))])]
unique(all_segRes$ID)[which(unique(all_segRes$ID) %in% annoTable$mouse_id[which(annoTable$geno %in% c("PRN"))])]
unique(all_segRes$ID)[-which(unique(all_segRes$ID) %in% annoTable$mouse_id)]

all_segRes_BPRN <- all_segRes[which(all_segRes$ID %in% annoTable2$mouse_id[which(annoTable2$geno %in% c("BPRN"))]),]
all_segRes_BPP <- all_segRes[which(all_segRes$ID %in% annoTable2$mouse_id[which(annoTable2$geno %in% c("BPP"))]),]
all_segRes_BPN <- all_segRes[which(all_segRes$ID %in% annoTable2$mouse_id[which(annoTable2$geno %in% c("BPN"))]),]
all_segRes_PRN <- all_segRes[which(all_segRes$ID %in% annoTable2$mouse_id[which(annoTable2$geno %in% c("PRN"))]),]


# for (i in seq_along(tc)) {
#   tc_names <- names(tc)[i]
#   all_segRes$seg.mean[which(all_segRes$ID == i)] <- all_segRes$seg.mean[which(all_segRes$ID == i)]/tc[i]
# }
                                                

#segZscores_BPRN <- calcZscore(all_segRes_BPRN)
segZfilt_BPRN <- segZscoresFilt_zeroOut(all_segRes_BPRN)
tmpSeg_BPRN <- cbind(segZfilt_BPRN[,1:4], NA, segZfilt_BPRN[,5])
colnames(tmpSeg_BPRN) <- c("sampleID","chrom", "start.pos","end.pos", "n.probes", "mean")
tmpSeg_BPRN$string <- paste0(tmpSeg_BPRN$sampleID, tmpSeg_BPRN$chrom, tmpSeg_BPRN$start.pos, tmpSeg_BPRN$end.pos)
tmpSeg_BPRN <- tmpSeg_BPRN[-which(duplicated(tmpSeg_BPRN$string)),]
tmpSeg_BPRN <- tmpSeg_BPRN[order(tmpSeg_BPRN$chrom, tmpSeg_BPRN$start.pos),]

#segZscores_BPP <- calcZscore(all_segRes_BPP)
segZfilt_BPP <- segZscoresFilt_zeroOut(all_segRes_BPP)
tmpSeg_BPP <- cbind(segZfilt_BPP[,1:4], NA, segZfilt_BPP[,5])
colnames(tmpSeg_BPP) <- c("sampleID","chrom", "start.pos","end.pos", "n.probes", "mean")
tmpSeg_BPP$string <- paste0(tmpSeg_BPP$sampleID, tmpSeg_BPP$chrom, tmpSeg_BPP$start.pos, tmpSeg_BPP$end.pos)
#tmpSeg_BPP <- tmpSeg_BPP[-which(duplicated(tmpSeg_BPP$string)),]
tmpSeg_BPP <- tmpSeg_BPP[order(tmpSeg_BPP$chrom, tmpSeg_BPP$start.pos),]


#segZscores_BPN <- calcZscore(all_segRes_BPN)
segZfilt_BPN <- segZscoresFilt_zeroOut(all_segRes_BPN)
tmpSeg_BPN <- cbind(segZfilt_BPN[,1:4], NA, segZfilt_BPN[,5])
colnames(tmpSeg_BPN) <- c("sampleID","chrom", "start.pos","end.pos", "n.probes", "mean")
tmpSeg_BPN$string <- paste0(tmpSeg_BPN$sampleID, tmpSeg_BPN$chrom, tmpSeg_BPN$start.pos, tmpSeg_BPN$end.pos)
#tmpSeg_BPN <- tmpSeg_BPN[-which(duplicated(tmpSeg_BPN$string)),]
tmpSeg_BPN <- tmpSeg_BPN[order(tmpSeg_BPN$chrom, tmpSeg_BPN$start.pos),]

#segZscores_PRN <- calcZscore(all_segRes_PRN)
segZfilt_PRN <- segZscoresFilt_zeroOut(all_segRes_PRN)
tmpSeg_PRN <- cbind(segZfilt_PRN[,1:4], NA, segZfilt_PRN[,5])
colnames(tmpSeg_PRN) <- c("sampleID","chrom", "start.pos","end.pos", "n.probes", "mean")
tmpSeg_PRN$string <- paste0(tmpSeg_PRN$sampleID, tmpSeg_PRN$chrom, tmpSeg_PRN$start.pos, tmpSeg_PRN$end.pos)
#tmpSeg_PRN <- tmpSeg_PRN[-which(duplicated(tmpSeg_PRN$string)),]
tmpSeg_PRN <- tmpSeg_PRN[order(tmpSeg_PRN$chrom, tmpSeg_PRN$start.pos),]




#armRes <- separateSegments_m(tmpSeg, mouse_arm)
armRes_bprn <- separateSegments_mV4(tmpSeg_BPRN)
bprn_arm <- do.call(rbind, armRes_bprn[[1]])
bprn_cna <- do.call(rbind, armRes_bprn[[2]])
bprn_count <- do.call(rbind, armRes_bprn[[3]])

bprn_arm_freq <- getFreqData(bprn_arm)
bprn_arm_res <- ampsDels(bprn_arm_freq)
bprn_arm_amp_bed <- reducingFreqBed(bprn_arm_res[[1]], bprn_arm_res[[2]])
bprn_arm_del_bed <- reducingFreqBed(bprn_arm_res[[3]], bprn_arm_res[[4]])

bprn_cna_freq <- getFreqData(bprn_cna)
bprn_cna_res <- ampsDels(bprn_cna_freq)
bprn_cna_amp_bed <- reducingFreqBed(bprn_cna_res[[1]], bprn_cna_res[[2]])
bprn_cna_del_bed <- reducingFreqBed(bprn_cna_res[[3]], bprn_cna_res[[4]])



hgsc_bprn_arm_allSynTable <- circosFreq(bprn_arm_amp_bed, bprn_arm_del_bed, hgsc_arm_amp_bed,
                                   hgsc_arm_del_bed, filename = "20210805kathyHgsc_arm_hgsc_bprn")

hgsc_bprn_cna_allSynTable <- circosFreq(bprn_cna_amp_bed, bprn_cna_del_bed, hgsc_cna_amp_bed, hgsc_cna_del_bed,
                                   filename = "20210805kathyHgsc_cna_hgsc_bprn")

brca_bprn_arm_allSynTable <- circosFreq(bprn_arm_amp_bed, bprn_arm_del_bed, brca_arm_amp_bed,
                                   brca_arm_del_bed, filename = "20210805kathyHgsc_arm_brca_bprn")

brca_bprn_cna_allSynTable <- circosFreq(bprn_cna_amp_bed, bprn_cna_del_bed, hgsc_cna_amp_bed, hgsc_cna_del_bed,
                                   filename = "20210805kathyHgsc_cna_brca_bprn")



### look at freqs
gi_count_BPRN <- do.call(rbind, armRes_bprn[[3]])
colnames(gi_count_BPRN) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")
gi_count_BPRN_df <- data.frame("samples" = rownames(gi_count_BPRN))
gi_count_BPRN_df <- cbind(gi_count_BPRN_df, gi_count_BPRN)
gi_count_BPRN_df <- melt(gi_count_BPRN_df)


armRes_bpn <- separateSegments_mV4(tmpSeg_BPN)
gi_count_BPN <- do.call(rbind, armRes_bpn[[3]])
colnames(gi_count_BPN) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")
gi_count_BPN_df <- data.frame("samples" = rownames(gi_count_BPN))
gi_count_BPN_df <- cbind(gi_count_BPN_df, gi_count_BPN)
gi_count_BPN_df <- melt(gi_count_BPN_df)

armRes_bpp <- separateSegments_mV4(tmpSeg_BPP)
gi_count_BPP <- do.call(rbind, armRes_bpp[[3]])
colnames(gi_count_BPP) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")
gi_count_BPP_df <- data.frame("samples" = rownames(gi_count_BPP))
gi_count_BPP_df <- cbind(gi_count_BPP_df, gi_count_BPP)
gi_count_BPP_df <- melt(gi_count_BPP_df)


armRes_prn <- separateSegments_mV4(tmpSeg_PRN)
gi_count_PRN <- do.call(rbind, armRes_prn[[3]])
colnames(gi_count_PRN) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")
gi_count_PRN_df <- data.frame("samples" = rownames(gi_count_PRN))
gi_count_PRN_df <- cbind(gi_count_PRN_df, gi_count_PRN)
gi_count_PRN_df <- melt(gi_count_PRN_df)


gi_count_WGS <- do.call(rbind, armRes_WGS[[3]])
colnames(gi_count_WGS) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")
gi_count_WGS_df <- data.frame("samples" = rownames(gi_count_WGS))
gi_count_WGS_df <- cbind(gi_count_WGS_df, gi_count_WGS)
gi_count_WGS_df <- melt(gi_count_WGS_df)


gi_count_human <- do.call(rbind, armRes[[3]])
colnames(gi_count_human) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")

hgsc_anno2_tp53_names <- str_replace_all(rownames(hgsc_anno2_tp53), "\\.", "\\-")
hgsc_anno2_tp53_brca12_names <- str_replace_all(rownames(hgsc_anno2_tp53_brca12), "\\.", "\\-")
hgsc_anno2_tp53_names <- hgsc_anno2_tp53_names[-which(hgsc_anno2_tp53_names %in% hgsc_anno2_tp53_brca12_names)]

gi_count_hgsc <- gi_count_human[grep(paste(hgsc_anno2_tp53_names, collapse = "|"), rownames(gi_count_human)),]
gi_count_hgsc_df <- data.frame("samples" = rownames(gi_count_hgsc))
gi_count_hgsc_df <- cbind(gi_count_hgsc_df, gi_count_hgsc)
gi_count_hgsc_df <- melt(gi_count_hgsc_df)

gi_count_brca <- gi_count_human[grep(paste(hgsc_anno2_tp53_brca12_names, collapse = "|"), rownames(gi_count_human)),]
gi_count_brca_df <- data.frame("samples" = rownames(gi_count_brca))
gi_count_brca_df <- cbind(gi_count_brca_df, gi_count_brca)
gi_count_brca_df <- melt(gi_count_brca_df)
# make boxplot for both human and mouse samples

a <- ggplot(gi_count_BPRN_df) + geom_boxplot(aes(x = variable, y = value)) + ylim(c(0,30)) + ggtitle("BPRN n = 78") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
b <- ggplot(gi_count_BPP_df) + geom_boxplot(aes(x = variable, y = value)) + ylim(c(0,30)) + ggtitle("BPP n = 11") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
c <- ggplot(gi_count_BPN_df) + geom_boxplot(aes(x = variable, y = value)) + ylim(c(0,30)) + ggtitle("BPN n = 7") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
d <- ggplot(gi_count_PRN_df) + geom_boxplot(aes(x = variable, y = value)) + ylim(c(0,30)) + ggtitle("PRN n = 10") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
e <- ggplot(gi_count_WGS_df) + geom_boxplot(aes(x = variable, y = value)) + ylim(c(0,30)) + ggtitle("WGS n = 10") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f <- ggplot(gi_count_hgsc_df) + geom_boxplot(aes(x = variable, y = value)) + ylim(c(0,30)) + ggtitle("hgsc n = 331") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
g <- ggplot(gi_count_brca_df) + geom_boxplot(aes(x = variable, y = value)) + ylim(c(0,30)) + ggtitle("hgsc(brca) n = 33") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210809countingGI.pdf", useDingbats = FALSE)
grid.arrange(a,b,c,d,e,f,g, ncol =3)
dev.off()


finalTable_hgsc <- finalTable2
finalTable_hgsc <- finalTable_hgsc[grep(paste(hgsc_anno2_tp53_names, collapse = "|"), finalTable_hgsc$Sample),]

finalTable_brca <- finalTable2
finalTable_brca <- finalTable_brca[grep(paste(hgsc_anno2_tp53_brca12_names, collapse = "|"), finalTable_brca$Sample),]

fga_bpn <- fgaCalculator_amp(tmpSeg_BPN)
fga_bprn <- fgaCalculator_amp(tmpSeg_BPRN)
fga_bpp <- fgaCalculator_amp(tmpSeg_BPP)
fga_prn <- fgaCalculator_amp(tmpSeg_PRN)
fga_wgs <- fgaCalculator_amp(allSamps2, mm10ChromSize2)
fga_hgsc <- fgaCalculator_tcga(finalTable_hgsc, tcga_ploidy)
fga_brca <- fgaCalculator_tcga(finalTable_brca, tcga_ploidy)

fga_bpn$type <- "bpn"
fga_bprn$type <- "bprn"
fga_bpp$type <- "bpp"
fga_prn$type <- "prn"
fga_wgs$type <- "wgs"
fga_brca$type <- "brca"
fga_hgsc$type <- "hgsc"

fga_all <- rbind(fga_bpn, fga_bprn, fga_bpp,
                 fga_prn, fga_wgs, fga_hgsc,
                 fga_brca)

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210810_fga.pdf", useDingbats = FALSE)
ggplot(fga_all, aes(x = type, y = X2)) + geom_boxplot() + xlab("genotype") + ylab("Fraction of genome altered")
dev.off()


bpp_arm <- do.call(rbind, armRes_bpp[[1]])
bpp_arm_freq <- getFreqData(bpp_arm)
bpp_arm_res <- ampsDels(bpp_arm_freq)
bpp_arm_amp_bed <- reducingFreqBed(bpp_arm_res[[1]], bpp_arm_res[[2]])
bpp_arm_del_bed <- reducingFreqBed(bpp_arm_res[[3]], bpp_arm_res[[4]])
bpp_cna <- do.call(rbind,armRes_bpp[[2]])
bpp_cna_freq <- getFreqData(bpp_cna)
bpp_cna_res <- ampsDels(bpp_cna_freq)
bpp_cna_amp_bed <- reducingFreqBed(bpp_cna_res[[1]], bpp_cna_res[[2]])
bpp_cna_del_bed <- reducingFreqBed(bpp_cna_res[[3]], bpp_cna_res[[4]])

bpn_arm <- do.call(rbind, armRes_bpn[[1]])
bpn_arm_freq <- getFreqData(bpn_arm)
bpn_arm_res <- ampsDels(bpn_arm_freq)
bpn_arm_amp_bed <- reducingFreqBed(bpn_arm_res[[1]], bpn_arm_res[[2]])
bpn_arm_del_bed <- reducingFreqBed(bpn_arm_res[[3]], bpn_arm_res[[4]])
bpn_cna <- do.call(rbind, armRes_bpn[[2]])
bpn_cna_freq <- getFreqData(bpn_cna)
bpn_cna_res <- ampsDels(bpn_cna_freq)
bpn_cna_amp_bed <- reducingFreqBed(bpn_cna_res[[1]], bpn_cna_res[[2]])
bpn_cna_del_bed <- reducingFreqBed(bpn_cna_res[[3]], bpn_cna_res[[4]])

prn_arm <- do.call(rbind, armRes_prn[[1]])
prn_arm_freq <- getFreqData(prn_arm)
prn_arm_res <- ampsDels(prn_arm_freq)
prn_arm_amp_bed <- reducingFreqBed(prn_arm_res[[1]], prn_arm_res[[2]])
prn_arm_del_bed <- reducingFreqBed(prn_arm_res[[3]], prn_arm_res[[4]])
prn_cna <- do.call(rbind, armRes_prn[[2]])
prn_cna_freq <- getFreqData(prn_cna)
prn_cna_res <- ampsDels(prn_cna_freq)
prn_cna_amp_bed <- reducingFreqBed(prn_cna_res[[1]], prn_cna_res[[2]])
prn_cna_del_bed <- reducingFreqBed(prn_cna_res[[3]], prn_cna_res[[4]])


wgs_arm <- do.call(rbind, armRes_WGS[[1]])
wgs_arm_freq <- getFreqData(wgs_arm)
wgs_arm_res <- ampsDels(wgs_arm_freq)
wgs_arm_amp_bed <- reducingFreqBed(wgs_arm_res[[1]], wgs_arm_res[[2]])
wgs_arm_del_bed <- reducingFreqBed(wgs_arm_res[[3]], wgs_arm_res[[4]])
wgs_cna <- do.call(rbind, armRes_WGS[[2]])
wgs_cna_freq <- getFreqData(wgs_cna)
wgs_cna_res <- ampsDels(wgs_cna_freq)
wgs_cna_amp_bed <- reducingFreqBed(wgs_cna_res[[1]], wgs_cna_res[[2]])
wgs_cna_del_bed <- reducingFreqBed(wgs_cna_res[[3]], wgs_cna_res[[4]])









a <- freqPlot(hgsc_arm_amp_bed, hgsc_arm_del_bed, speciesType = "human", main = "hgsc n = 331")
b <- freqPlot(brca_arm_amp_bed, brca_arm_del_bed, speciesType = "human", main = "hgsc w/ brca n = 33")
c <- freqPlot(bprn_arm_amp_bed, bprn_arm_del_bed, speciesType = "mouse", main = "bprn n = 78")
d <- freqPlot(bpp_arm_amp_bed, bpp_arm_del_bed, speciesType = "mouse", main = "bpp n = 11")
e <- freqPlot(bpn_arm_amp_bed, bpn_arm_del_bed, speciesType = "mouse", main = "bpn n = 7")
f <- freqPlot(prn_arm_amp_bed, prn_arm_del_bed, speciesType = "mouse", main = "prn n = 10")
g <- freqPlot(wgs_arm_amp_bed, wgs_arm_del_bed, speciesType = "mouse", main = "wgs n = 10")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210809freqPlotsAllArms.pdf", height = 9, width = 9, useDingbats = FALSE)
grid.arrange(a,b,c,d,e,f,g, ncol = 1)
dev.off()

a <- freqPlot(hgsc_cna_amp_bed, hgsc_cna_del_bed, speciesType = "human", main = "hgsc n = 331")
b <- freqPlot(brca_cna_amp_bed, brca_cna_del_bed, speciesType = "human", main = "hgsc w/ brca n = 33")
c <- freqPlot(bprn_cna_amp_bed, bprn_cna_del_bed, speciesType = "mouse", main = "bprn n = 78")
d <- freqPlot(bpp_cna_amp_bed, bpp_cna_del_bed, speciesType = "mouse", main = "bpp n = 11")
e <- freqPlot(bpn_cna_amp_bed, bpn_cna_del_bed, speciesType = "mouse", main = "bpn n = 7")
f <- freqPlot(prn_cna_amp_bed, prn_cna_del_bed, speciesType = "mouse", main = "prn n = 10")
g <- freqPlot(wgs_cna_amp_bed, wgs_cna_del_bed, speciesType = "mouse", main = "wgs n = 10")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210809freqPlotsAllCnas.pdf", height = 9, width = 9, useDingbats = FALSE)
grid.arrange(a,b,c,d,e,f,g, ncol = 1)
dev.off()


### pathway analysis
###
###
human_hallmarks_list <- NULL
for (i in 1:nrow(human_hallmarks)) {
  tmpList <- list(human_hallmarks[i, 3:ncol(human_hallmarks)])
  names(tmpList) <- human_hallmarks[i, 1]
  human_hallmarks_list[[i]] <- tmpList
}


hgsc_bprn_arm_geneList <- getGeneList(hgsc_bprn_arm_allSynTable)
hgsc_bprn_cna_geneList <- getGeneList(hgsc_bprn_cna_allSynTable) 
brca_bprn_arm_geneList <- getGeneList(brca_bprn_arm_allSynTable)
brca_bprn_cna_geneList <- getGeneList(brca_bprn_cna_allSynTable)

hgsc_bprn_arm_enrichmentStats <- enrichmentStats(hgsc_bprn_arm_geneList, human_hallmarks_list)
hgsc_bprn_arm_path_res <- hgsc_bprn_arm_enrichmentStats[[1]]
hgsc_bprn_arm_path_df <- hgsc_bprn_arm_enrichmentStats[[2]]
hgsc_bprn_cna_enrichmentStats <- enrichmentStats(hgsc_bprn_cna_geneList, human_hallmarks_list)
hgsc_bprn_cna_path_res <- hgsc_bprn_cna_enrichmentStats[[1]]
hgsc_bprn_cna_path_df <- hgsc_bprn_cna_enrichmentStats[[2]]

enrichmentBarplot(hgsc_bprn_arm_path_res)
enrichmentBarplot(hgsc_bprn_cna_path_res)

brca_bprn_arm_enrichmentStats <- enrichmentStats(brca_bprn_arm_geneList, human_hallmarks_list)
brca_bprn_cna_enrichmentStats <- enrichmentStats(brca_bprn_cna_geneList, human_hallmarks_list)

### only thing to add is the segmentation data from the other mouse model