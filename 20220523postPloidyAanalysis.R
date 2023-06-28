source("/home/kevhu/scripts/20210802syntenyFunctions.R")

library(ggpubr)

nameStripper <- function(df){
  require(stringr)
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
}



synteny_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_syntenyDf.txt",sep = "\t", stringsAsFactors = FALSE, header = TRUE)

allPanCancerSampleAnno <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/merged_sample_quality_annotations.tsv",
                                     sep = "\t", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
allHumanHgscAnno <- allPanCancerSampleAnno[which(allPanCancerSampleAnno$cancer.type == 'OV'),]


humanHrdCalc <- read.table("/mnt/DATA5/tmp/kev/misc/bmc2015HrdRes.txt", sep = "\t",
                           stringsAsFactors = FALSE, header = TRUE)

tcga_ploidy <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",sep = "\t", stringsAsFactors = FALSE, header = TRUE)
hgsc_anno <- read.table("/mnt/DATA5/tmp/kev/misc/PATIENT_DATA_oncoprint_hgsc_Tp53_Brca.tsv", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)

tcga2012crcMethBeta <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/tcga_crc_2012_beta_val.txt", sep = "\t",
                                  stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)


allNoshadSeg <- read.table("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV2/absoluteRes.txt", sep = "\t",
                           stringsAsFactors = FALSE, header = TRUE, fill = TRUE)

allPloidyCalls <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220517allHgscPloidy.xlsx")
allPloidyCalls <- allPloidyCalls[-which(allPloidyCalls$ploidy_int == "NA"), ]

allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "mg4", "2027lte")
allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "mg15", "13085lt")
allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "mg20", "14399rt")



allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "14150lt", "14150rt")
allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "14154lt", "14154rot")
allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "14656peritnealmt", "14656peritonealmt")
allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "13576rt", "133576rt")





allNoshadSeg_filt <- allNoshadSeg[grep(paste0(paste(allPloidyCalls$sample, round(as.numeric(allPloidyCalls$purityV2), 3),
                                                    round(as.numeric(allPloidyCalls$ploidyV2), 2), sep = "_"), collapse = "|"), allNoshadSeg$sample),]
allNoshadSeg_filt$sC <- as.numeric(allNoshadSeg_filt$sC)

### check
paste(allPloidyCalls$sample, round(as.numeric(allPloidyCalls$purityV2), 3),
      round(as.numeric(allPloidyCalls$ploidyV2), 2), sep = "_")[-which(paste(allPloidyCalls$sample, round(as.numeric(allPloidyCalls$purityV2), 3),
                                                                          round(as.numeric(allPloidyCalls$ploidyV2), 2), sep = "_") %in% unique(allNoshadSeg_filt$sample))]
i <- 25
allNoshadSeg_filt$absCn <- NA
for (i in 1:nrow(allPloidyCalls)) {
  tmpString <- paste(allPloidyCalls$sample[i], round(as.numeric(allPloidyCalls$purityV2[i]), 3),
                     round(as.numeric(allPloidyCalls$ploidyV2[i]), 2), sep = "_")
  tmpMatch <- grep(tmpString, allNoshadSeg_filt$sample)
  tmpDf <- allNoshadSeg_filt[tmpMatch,]
  allNoshadSeg_filt$absCn[tmpMatch] <- round(allNoshadSeg_filt$sC[tmpMatch] - as.numeric(allPloidyCalls$ploidy_int[i]))
}

### prepping the absCn data for downstream pipelines
noshadAbsCnFreqDf <- allNoshadSeg_filt[,c("sample", "chr", "start", "end", "K", "absCn")]
noshadAbsCnFreqDf$K <- NA
colnames(noshadAbsCnFreqDf) <- c("sampleID", "chrom", "start.pos", "end.pos", "n.probes", "mean")
noshadAbsCnFreqDf$string <- paste0(noshadAbsCnFreqDf$sampleID, noshadAbsCnFreqDf$chrom,noshadAbsCnFreqDf$start.pos, noshadAbsCnFreqDf$end.pos)
noshadAbsCnFreqDf$chrom <- str_remove(noshadAbsCnFreqDf$chrom, "chr")
noshadAbsCnFreqDf <- noshadAbsCnFreqDf[which(noshadAbsCnFreqDf$chrom %in% c(1:19)),]
noshadAbsCnFreqDf[,c("chrom", "start.pos", "end.pos", "mean")] <- lapply(noshadAbsCnFreqDf[,c("chrom", "start.pos", "end.pos", "mean")], as.numeric)

### from documentation https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
### splitting sample strings to grab normal IDS
### use zscore method from tcga supplement https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3401966/

beta_names <- str_split(colnames(tcga2012crcMethBeta), "\\-")
beta_names_sample <- unlist(lapply(beta_names, '[', 4))
beta_names_sample[grep(paste0(c(10:19), collapse = "|"), beta_names_sample)]

controlBetaCrc_mlh1 <- tcga2012crcMethBeta[grep(paste0(c("MLH1", "EPM2AIP1"), collapse = "|"), tcga2012crcMethBeta$Gene),
                                      c(1:4, grep(paste0(c(10:19), collapse = "|"),
                                                  beta_names_sample))]


allKathyAnno <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220523allHgscNewPanelAnnoV2.xlsx")
allKathyAnno$Sample <- nameStripper(allKathyAnno$Sample)
allKathyAnno$Sample <- str_remove(allKathyAnno$Sample, "\\-")
allKathyAnno$Sample[9] <- "2519lt"

allKathyAnno <- allKathyAnno[,1:3]
annoTable <- read.table("/mnt/DATA5/tmp/kev/misc/20210713allHgscMouseAnno_noCrisp.txt", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)
annoTable <- annoTable[, c(1,3,2)]
annoTable$mouse_id <- str_remove_all(annoTable$mouse_id, "ovc")
annoTable$mouse_id <- nameStripper(annoTable$mouse_id)
annoTable$mouse_id <- str_remove_all(annoTable$mouse_id, "\\-")
colnames(allKathyAnno) <- colnames(annoTable)
annoTable_combined <- rbind(annoTable, allKathyAnno)
annoTable_combined$type2 <- tolower(annoTable_combined$type)
annoTable_combined <- annoTable_combined[which(annoTable_combined$type2 %in% c("hgsc", "mmmt")), ]

# allPloidyCalls <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220517allHgscPloidy.xlsx")
colnames(allPloidyCalls)[5] <- "altPloidy"


hgscNonzeroTc <- read.table("/mnt/DATA5/tmp/kev/misc/20220521onlyHgscNonZeroTc.txt", sep = "\t",
                            header = TRUE, stringsAsFactors = FALSE)
hgscNonzeroTc <- hgscNonzeroTc[-grep("umm", hgscNonzeroTc$sample),]
hgscNonzeroTc <- hgscNonzeroTc[-grep("cet", hgscNonzeroTc$sample),]

otherSamps <- c("mg28", "mg29", "mg30", "mg31", "mg32", "mg33", "mg34", "mg35")
hgscPloidyCalls <- allPloidyCalls[which(allPloidyCalls$sample %in% hgscNonzeroTc$sample),]
hgscPloidyCalls <- hgscPloidyCalls[-which(hgscPloidyCalls$sample %in% otherSamps),]

hgscNonzeroTc <- hgscNonzeroTc[which(hgscNonzeroTc$sample %in% hgscPloidyCalls$sample), ]

dupeIdx <- c(49, 51, 7, 27, 112)
hgscPloidyCalls <- hgscPloidyCalls[-dupeIdx, ]

hgscPloidyCalls$altPloidy <- hgscNonzeroTc$tc[match(hgscPloidyCalls$sample, hgscNonzeroTc$sample)]
hgscPloidyCalls[,2:ncol(hgscPloidyCalls)] <- lapply(hgscPloidyCalls[,2:ncol(hgscPloidyCalls)], function(x) as.numeric(x))
hgscPloidyCalls$geno <- annoTable_combined$geno[match(hgscPloidyCalls$sample, annoTable_combined$mouse_id)]
hgscPloidyCalls <- hgscPloidyCalls[-which(is.na(hgscPloidyCalls$geno)), ]
hgscPloidyCalls <- hgscPloidyCalls[-which(hgscPloidyCalls$geno == "unk"), ]

mouseDensity <- ggplot(hgscPloidyCalls) + geom_density(aes(ploidy)) + xlim(c(1, 5)) + ggtitle("Mouse hgsc ploidy distribution")

hgscPloidyCalls_filt <- hgscPloidyCalls[which(hgscPloidyCalls$geno %in% c("BPRN", "BPN", "BPP")),]
ggplot(hgscPloidyCalls_filt, aes(x = ploidy)) + geom_density(aes(group = geno, fill = geno)) + xlim(c(1, 5)) + 
  ggtitle("Mouse hgsc split by geno")

### human hgsc from 20210916absMouseHgscSynteny

hgsc_anno2 <- hgsc_anno[4:nrow(hgsc_anno),3:ncol(hgsc_anno)]
rownames(hgsc_anno2) <- paste0(hgsc_anno$track_name[4:nrow(hgsc_anno)],
                               "_", hgsc_anno$track_type[4:nrow(hgsc_anno)])

hgsc_anno2 <- data.frame(t(hgsc_anno2), stringsAsFactors = FALSE)
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
rownames(hgsc_anno2_tp53_brca12) <- str_replace_all(rownames(hgsc_anno2_tp53_brca12), pattern = "\\.", replacement = "\\-")

rownames(hgsc_anno2_tp53) <- str_replace_all(rownames(hgsc_anno2_tp53), pattern = "\\.", replacement = "\\-")
hg19_hgsc_ploidy <- tcga_ploidy[grep(paste0(rownames(hgsc_anno2_tp53), collapse = "|"), tcga_ploidy$sample),]
hrd_hg19_tp53 <- humanHrdCalc[which(humanHrdCalc$Tumor %in% rownames(hgsc_anno2_tp53) ),]

hrd_hg19_tp53 <- humanHrdCalc[which(humanHrdCalc$Tumor %in% rownames(hgsc_anno2_tp53) ),]
hrd_hg19_tp53$hrd.score <- hrd_hg19_tp53$NtAI + hrd_hg19_tp53$LST + hrd_hg19_tp53$HRD.LOH
hrd_hg19_tp53$hrd.status <- ifelse(hrd_hg19_tp53$hrd.score > 62, "HRD", "non-HRD")
hrd_hg19_tp53$hrd.status <- factor(hrd_hg19_tp53$hrd.status)

humanDensity <- ggplot(hg19_hgsc_ploidy) + geom_density(aes(ploidy)) + xlim(c(1, 5)) + ggtitle("Human hgsc ploidy distribution")


ggplot(hrd_hg19_tp53, aes(x = Ploidy)) + geom_density(aes(group = hrd.status, fill = hrd.status)) + xlim(c(1, 5)) + 
  ggtitle("Human ploidy split by HRD score > 62")




gridExtra::grid.arrange(humanDensity, mouseDensity, ncol = 1)


quantile(hrd_hg19_tp53$Ploidy, seq(0,1,0.05))
quantile(hg19_hgsc_ploidy$ploidy, seq(0,1,0.05),na.rm = TRUE)
quantile(hgscPloidyCalls_filt$ploidy[which(hgscPloidyCalls_filt$geno == "BPRN")], seq(0,1,0.05),na.rm = TRUE)
quantile(hgscPloidyCalls_filt$ploidy[which(hgscPloidyCalls_filt$geno == "BPP")], seq(0,1,0.05),na.rm = TRUE)
quantile(hgscPloidyCalls_filt$ploidy[which(hgscPloidyCalls_filt$geno == "BPN")], seq(0,1,0.05),na.rm = TRUE)


### think about splitting by LST or just ones with BRCA-like HRD scores b/c things with 62 >= are similar to BRCA mutations
###
###


coad_anno <- read.table("/mnt/DATA5/tmp/kev/misc/COAD_ALL_APC_TP53_KRAS_MLH1_PATIENT_DATA_oncoprint.tsv", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)
coad_anno2 <- data.frame(t(coad_anno[c(4:9, 19:23),]), stringsAsFactors = FALSE)
new_header <- paste0(coad_anno2[1,], "_", coad_anno2[2,])
colnames(coad_anno2) <- new_header
coad_anno2 <- coad_anno2[3:nrow(coad_anno2), ]
coad_anno2[, 7:10] <- lapply(coad_anno2[, 7:10], as.numeric)
rownames(coad_anno2) <- str_replace_all(rownames(coad_anno2), "\\.", "\\-")

listOfMatched <- c("cg02279071", "cg00893636", "cg10990993", "cg18320188")
listOfMatched3 <- c("cg18320188")
zScoreDf <- NULL
pValDf <- NULL
for (i in listOfMatched3) {
  tmpAvg <- mean(as.numeric(controlBetaCrc_mlh1[which(controlBetaCrc_mlh1$Composite == i), 5:ncol(controlBetaCrc_mlh1)]))
  tmpSd <-  sd(as.numeric(controlBetaCrc_mlh1[which(controlBetaCrc_mlh1$Composite == i), 5:ncol(controlBetaCrc_mlh1)]))
  zScoreList <- sapply(coad_anno2[, grep(i, colnames(coad_anno2))], function(x) (as.numeric(x) - tmpAvg)/tmpSd)
  pvalList <- sapply(zScoreList, function(x) 2*pnorm(x, lower.tail=FALSE))
  pAdjList <- sapply(pvalList, function(x) p.adjust(x, method = "BH", n = length(pvalList)))
  
  zScoreDf <- cbind(zScoreDf, zScoreList)
  pValDf <- cbind(pValDf, pAdjList)
}


# colnames(zScoreDf) <- paste0("zscore_", listOfMatched)
# colnames(pValDf) <- paste0("pval_", listOfMatched)

coad_anno2 <- cbind(coad_anno2, zScoreDf, pValDf)
coad_anno2$hyperMethStatus <- "none"
coad_anno2$hyperMethStatus[which(coad_anno2$zScoreList > 0 & coad_anno2$pAdjList < 0.05)] <- "hyperMethyl"

### the general idea of sporadic cancers according to Fearon and Vogelstein
### even the mice show not all 3 are necessary to produce the tumors
### the more mutations, the more penetrant the COAD phenotype is? we don't have the numbers i.e same biological samples for a lot of them
### should add some additional split in the density based on hypermethylated genotype - not exact way to determine hypermethylation but
### we know ~15% of sporadic cancers are hypermethylated ... we can find the cutoff for that in our sample

only_Apc <-  which(coad_anno2$APC_CNA == "homdel_rec" | coad_anno2$APC_MUTATIONS == "Truncating mutation (putative driver)")


sporadic_all_3 <- Reduce(intersect, list(which(coad_anno2$APC_CNA == "homdel_rec" | coad_anno2$APC_MUTATIONS == "Truncating mutation (putative driver)"),
                                    which(coad_anno2$TP53_CNA == "homdel_rec" | coad_anno2$TP53_MUTATIONS == "Inframe Mutation (putative driver)" |
                                            coad_anno2$TP53_MUTATIONS == "Missense Mutation (putative driver)" | coad_anno2$TP53_MUTATIONS == "Truncating mutation (putative driver)"),
                                    which(coad_anno2$KRAS_MUTATIONS == "Inframe Mutation (putative driver)" | coad_anno2$KRAS_MUTATIONS == "Missense Mutation (putative driver)")))

sporadic_apc_tp53 <- Reduce(intersect, list(which(coad_anno2$APC_CNA == "homdel_rec" | coad_anno2$APC_MUTATIONS == "Truncating mutation (putative driver)"),
                                   which(coad_anno2$TP53_CNA == "homdel_rec" | coad_anno2$TP53_MUTATIONS == "Inframe Mutation (putative driver)" |
                                           coad_anno2$TP53_MUTATIONS == "Missense Mutation (putative driver)" | coad_anno2$TP53_MUTATIONS == "Truncating mutation (putative driver)")))

sporadic_apc_kras <- Reduce(intersect, list(which(coad_anno2$APC_CNA == "homdel_rec" | coad_anno2$APC_MUTATIONS == "Truncating mutation (putative driver)"),
                                   which(coad_anno2$KRAS_MUTATIONS == "Inframe Mutation (putative driver)" | coad_anno2$KRAS_MUTATIONS == "Missense Mutation (putative driver)")))

sporadic_apc_tp53 <- sporadic_apc_tp53[-which(sporadic_apc_tp53 %in% sporadic_all_3)]
sporadic_apc_kras <- sporadic_apc_kras[-which(sporadic_apc_kras %in% sporadic_all_3)]

coad_only_apc <- coad_anno2[only_Apc, ]
coad_sporadic_all <- coad_anno2[sporadic_all_3,]
coad_sporadic_tp53 <- coad_anno2[sporadic_apc_tp53,]
coad_sporadic_kras <- coad_anno2[sporadic_apc_kras,]



coad_ploidy <- tcga_ploidy[grep(paste0(rownames(coad_anno2), collapse = "|"),
                                tcga_ploidy$array),]
coad_ploidy$Mlh1Utr5 <- "none"
coad_ploidy$Mlh1Utr5[grep(paste0(rownames(coad_anno2)[which(coad_anno2$hyperMethStatus == "hyperMethyl")], collapse = "|"),
                          coad_ploidy$array)] <- "hyperMethyl"

coad_spor_only_apc <- coad_ploidy[grep(paste0(rownames(coad_only_apc), collapse = "|"), coad_ploidy$sample),]
coad_spor_3_ploidy <- coad_ploidy[grep(paste0(rownames(coad_sporadic_all), collapse = "|"), coad_ploidy$sample),]
coad_spor_tp53_ploidy <- coad_ploidy[grep(paste0(rownames(coad_sporadic_tp53), collapse = "|"), coad_ploidy$sample),]
coad_spor_kras_ploidy <- coad_ploidy[grep(paste0(rownames(coad_sporadic_kras), collapse = "|"), coad_ploidy$sample),]



coad_3 <- ggplot(coad_spor_3_ploidy) + geom_density(aes(ploidy)) + xlim(c(1, 5)) + ggtitle("TCGA COAD Apc, P53, Kras ploidy distribution mut. exclusive")
coad_tp53 <- ggplot(coad_spor_tp53_ploidy) + geom_density(aes(ploidy)) + xlim(c(1, 5)) + ggtitle("TCGA COAD Apc, P53 ploidy distribution mut. exclusive")
coad_kras <- ggplot(coad_spor_kras_ploidy) + geom_density(aes(ploidy)) + xlim(c(1, 5)) + ggtitle("TCGA COAD Apc,Kras ploidy distribution mut. exclusive")

gridExtra::grid.arrange(coad_3, coad_tp53,  coad_kras, ncol = 1)

### for the mouse data ...  separate by histology
### then by genotype  - 2 separate

allFearonAnno <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220531allFearonAnno.xlsx")
allFearonAnno$Sample <- nameStripper(allFearonAnno$`Sequencing name`)
allFearonAnno$Sample <- str_remove(allFearonAnno$Sample, "\\-")
allFearonAnno$Sample[1:9] <- paste0("efd0", 1:9)

coadNonzeroTc <- read.table("/mnt/DATA5/tmp/kev/misc/20220521onlyEfdApcTc.txt", sep = "\t",
                            header = TRUE, stringsAsFactors = FALSE)
coadNonzeroTc <- coadNonzeroTc[which(coadNonzeroTc$tc > 0.2), ]
coadPloidyCalls <- allPloidyCalls[which(allPloidyCalls$sample %in% coadNonzeroTc$sample),]

coadPloidyCalls[,2:4] <- lapply(coadPloidyCalls[,2:4], function(x) as.numeric(x))
coadPloidyCalls$pathology <- allFearonAnno$Pathology[match(coadPloidyCalls$sample, allFearonAnno$Sample)]
coadPloidyCalls$geno <- allFearonAnno$Genotype[match(coadPloidyCalls$sample, allFearonAnno$Sample)]

coadPloidyCalls$ploidy_int <- as.numeric(coadPloidyCalls$ploidy_int)

coadPloidyCalls <- coadPloidyCalls[which(coadPloidyCalls$geno %in% "Adenocarcinoma"),]
mouseDensityCoad <- ggplot(coadPloidyCalls) + geom_density(aes(ploidy_int)) + xlim(c(1, 5)) + ggtitle("Mouse coad ploidy distribution")

hgscPloidyCalls_filt <- hgscPloidyCalls[which(hgscPloidyCalls$geno %in% c("BPRN", "BPN", "BPP")),]





### frequency results for humans. mice might be different and easier - i.e use method from snpRes comp to get anueploidy, discount for cna ... etc
### added options to change minlength
### for these process them in previous script 20210916absMouseHgscSynteny and then load them here



### processing all aneuploidy + cna stuff then

### all tcga first

finalTable <- read.table("/mnt/DATA5/tmp/kev/misc/20210628panCancerAllSegReduced_hg19.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")

finalTable$length <- finalTable$End - finalTable$Start
finalTable2 <- finalTable[which(finalTable$Chromosome %in% c(1:22)),]
finalTable2 <- finalTable2[order(finalTable2$Chromosome, finalTable2$Start),]
armRes <- separateSegments_tcgaV3(finalTable2, tcga_ploidy, minLength = 20000000)

armTable <- do.call(rbind, armRes[[1]])
cnaTable <- do.call(rbind, armRes[[2]])
countTable <- do.call(rbind, armRes[[3]])

tcga_arm <- armTable
tcga_cna <- cnaTable


### human hgsc
###
###

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


### making hgsc hrd split
hrd_hg19_tp53$Tumor[which(hrd_hg19_tp53$hrd.status == "non-HRD")]

hgsc_hrd_arm <- tcga_arm[grep(paste(hrd_hg19_tp53$Tumor[which(hrd_hg19_tp53$hrd.status == "HRD")],
                                    collapse = "|"), tcga_arm$Sample),]
hgsc_hrd_cna <- tcga_cna[grep(paste(hrd_hg19_tp53$Tumor[which(hrd_hg19_tp53$hrd.status == "HRD")],
                                    collapse = "|"), tcga_arm$Sample),]

hgsc_nonhrd_arm <- tcga_arm[grep(paste(hrd_hg19_tp53$Tumor[which(hrd_hg19_tp53$hrd.status == "non-HRD")],
                                       collapse = "|"), tcga_arm$Sample),]
hgsc_nonhrd_cna <- tcga_cna[grep(paste(hrd_hg19_tp53$Tumor[which(hrd_hg19_tp53$hrd.status == "non-HRD")],
                                       collapse = "|"), tcga_arm$Sample),]



hrd_arm2 <- data.frame("sampleID" = hgsc_hrd_arm$Sample, "chrom" = hgsc_hrd_arm$Chromosome,
                        "start.pos" = hgsc_hrd_arm$Start, "end.pos" = hgsc_hrd_arm$End,
                        "n.probes" = NA, "mean" = hgsc_hrd_arm$newTotalCN, stringsAsFactors = FALSE)
hrd_cna2 <- data.frame("sampleID" = hgsc_hrd_cna$Sample, "chrom" = hgsc_hrd_cna$Chromosome,
                        "start.pos" = hgsc_hrd_cna$Start, "end.pos" = hgsc_hrd_cna$End,
                        "n.probes" = NA, "mean" = hgsc_hrd_cna$newTotalCN, stringsAsFactors = FALSE)

nonhrd_arm2 <- data.frame("sampleID" = hgsc_nonhrd_arm$Sample, "chrom" = hgsc_nonhrd_arm$Chromosome,
                        "start.pos" = hgsc_nonhrd_arm$Start, "end.pos" =hgsc_nonhrd_arm$End,
                        "n.probes" = NA, "mean" = hgsc_nonhrd_arm$newTotalCN, stringsAsFactors = FALSE)
nonhrd_cna2 <- data.frame("sampleID" = hgsc_nonhrd_cna$Sample, "chrom" = hgsc_nonhrd_cna$Chromosome,
                        "start.pos" = hgsc_nonhrd_cna$Start, "end.pos" = hgsc_nonhrd_cna$End,
                        "n.probes" = NA, "mean" = hgsc_nonhrd_cna$newTotalCN, stringsAsFactors = FALSE)



hrd_arm_freq_out <- getFreqData(hrd_arm2)
ampDels_hrd_arm <- ampsDels(hrd_arm_freq_out)
hrd_arm_amp_bed <- reducingFreqBed(ampDels_hrd_arm[[1]], ampDels_hrd_arm[[2]])
hrd_arm_del_bed <- reducingFreqBed(ampDels_hrd_arm[[3]], ampDels_hrd_arm[[4]])

hrd_cna_freq_out <- getFreqData(hrd_cna2)
ampDels_hrd_cna <- ampsDels(hrd_cna_freq_out)
hrd_cna_amp_bed <- reducingFreqBed(ampDels_hrd_cna[[1]], ampDels_hrd_cna[[2]])
hrd_cna_del_bed <- reducingFreqBed(ampDels_hrd_cna[[3]], ampDels_hrd_cna[[4]])

nonhrd_arm_freq_out <- getFreqData(nonhrd_arm2)
ampDels_nonhrd_arm <- ampsDels(nonhrd_arm_freq_out)
nonhrd_arm_amp_bed <- reducingFreqBed(ampDels_nonhrd_arm[[1]], ampDels_nonhrd_arm[[2]])
nonhrd_arm_del_bed <- reducingFreqBed(ampDels_nonhrd_arm[[3]], ampDels_nonhrd_arm[[4]])

nonhrd_cna_freq_out <- getFreqData(nonhrd_cna2)
ampDels_nonhrd_cna <- ampsDels(nonhrd_cna_freq_out)
nonhrd_cna_amp_bed <- reducingFreqBed(ampDels_nonhrd_cna[[1]], ampDels_nonhrd_cna[[2]])
nonhrd_cna_del_bed <- reducingFreqBed(ampDels_nonhrd_cna[[3]], ampDels_nonhrd_cna[[4]])



### human coad
###
###

coad_sporadic_all_arm <- tcga_arm[grep(paste(rownames(coad_sporadic_all), collapse = "|"), tcga_arm$Sample),]
coad_sporadic_all_cna <- tcga_cna[grep(paste(rownames(coad_sporadic_all), collapse = "|"), tcga_arm$Sample),]

coad_sporadic_tp53_arm <- tcga_arm[grep(paste(rownames(coad_sporadic_tp53), collapse = "|"), tcga_arm$Sample),]
coad_sporadic_tp53_cna <- tcga_cna[grep(paste(rownames(coad_sporadic_tp53), collapse = "|"), tcga_cna$Sample),]

coad_sporadic_kras_arm <- tcga_arm[grep(paste(rownames(coad_sporadic_kras), collapse = "|"), tcga_arm$Sample),]
coad_sporadic_kras_cna <- tcga_cna[grep(paste(rownames(coad_sporadic_kras), collapse = "|"), tcga_cna$Sample),]



coad_sporadic_all_arm2 <- data.frame("sampleID" = coad_sporadic_all_arm$Sample, "chrom" = coad_sporadic_all_arm$Chromosome,
                        "start.pos" = coad_sporadic_all_arm$Start, "end.pos" = coad_sporadic_all_arm$End,
                        "n.probes" = NA, "mean" = coad_sporadic_all_arm$newTotalCN, stringsAsFactors = FALSE)

coad_sporadic_all_cna2 <- data.frame("sampleID" = coad_sporadic_all_cna$Sample, "chrom" =coad_sporadic_all_cna$Chromosome,
                        "start.pos" = coad_sporadic_all_cna$Start, "end.pos" = coad_sporadic_all_cna$End,
                        "n.probes" = NA, "mean" = coad_sporadic_all_cna$newTotalCN, stringsAsFactors = FALSE)

coad_sporadic_tp53_arm2 <- data.frame("sampleID" = coad_sporadic_tp53_arm$Sample, "chrom" = coad_sporadic_tp53_arm$Chromosome,
                        "start.pos" =coad_sporadic_tp53_arm$Start, "end.pos" = coad_sporadic_tp53_arm$End,
                        "n.probes" = NA, "mean" = coad_sporadic_tp53_arm$newTotalCN, stringsAsFactors = FALSE)
coad_sporadic_tp53_cna2 <- data.frame("sampleID" = coad_sporadic_tp53_cna$Sample, "chrom" = coad_sporadic_tp53_cna$Chromosome,
                        "start.pos" = coad_sporadic_tp53_cna$Start, "end.pos" = coad_sporadic_tp53_cna$End,
                        "n.probes" = NA, "mean" = coad_sporadic_tp53_cna$newTotalCN, stringsAsFactors = FALSE)

coad_sporadic_kras_arm2 <- data.frame("sampleID" = coad_sporadic_kras_arm$Sample, "chrom" = coad_sporadic_kras_arm$Chromosome,
                                      "start.pos" = coad_sporadic_kras_arm$Start, "end.pos" = coad_sporadic_kras_arm$End,
                                      "n.probes" = NA, "mean" = coad_sporadic_kras_arm$newTotalCN, stringsAsFactors = FALSE)
coad_sporadic_kras_cna2 <- data.frame("sampleID" = coad_sporadic_kras_cna$Sample, "chrom" = coad_sporadic_kras_cna$Chromosome,
                                      "start.pos" = coad_sporadic_kras_cna$Start, "end.pos" = coad_sporadic_kras_cna$End,
                                      "n.probes" = NA, "mean" = coad_sporadic_kras_cna$newTotalCN, stringsAsFactors = FALSE)



coad_sporadic_all_freq_out <- getFreqData(coad_sporadic_all_arm2)
ampDels_coad_sporadic_all <- ampsDels(coad_sporadic_all_freq_out)
coad_sporadic_all_amp_bed <- reducingFreqBed(ampDels_coad_sporadic_all[[1]], ampDels_coad_sporadic_all[[2]])
coad_sporadic_all_del_bed <- reducingFreqBed(ampDels_coad_sporadic_all[[3]], ampDels_coad_sporadic_all[[4]])

coad_sporadic_all_cna_freq_out <- getFreqData(coad_sporadic_all_cna2)
ampDels_coad_sporadic_all_cna <- ampsDels(coad_sporadic_all_cna_freq_out)
coad_sporadic_all_cna_amp_bed <- reducingFreqBed(ampDels_coad_sporadic_all_cna[[1]], ampDels_coad_sporadic_all_cna[[2]])
coad_sporadic_all_cna_del_bed <- reducingFreqBed(ampDels_coad_sporadic_all_cna[[3]], ampDels_coad_sporadic_all_cna[[4]])


coad_sporadic_tp53_freq_out <- getFreqData(coad_sporadic_tp53_arm2)
ampDels_coad_sporadic_tp53 <- ampsDels(coad_sporadic_tp53_freq_out)
coad_sporadic_tp53_amp_bed <- reducingFreqBed(ampDels_coad_sporadic_tp53[[1]], ampDels_coad_sporadic_tp53[[2]])
coad_sporadic_tp53_del_bed <- reducingFreqBed(ampDels_coad_sporadic_tp53[[3]], ampDels_coad_sporadic_tp53[[4]])

coad_sporadic_tp53_cna_freq_out <- getFreqData(coad_sporadic_tp53_cna2)
ampDels_coad_sporadic_tp53_cna <- ampsDels(coad_sporadic_tp53_cna_freq_out)
coad_sporadic_tp53_cna_amp_bed <- reducingFreqBed(ampDels_coad_sporadic_tp53_cna[[1]], ampDels_coad_sporadic_tp53_cna[[2]])
coad_sporadic_tp53_cna_del_bed <- reducingFreqBed(ampDels_coad_sporadic_tp53_cna[[3]], ampDels_coad_sporadic_tp53_cna[[4]])


coad_sporadic_kras_freq_out <- getFreqData(coad_sporadic_kras_arm2)
ampDels_coad_sporadic_kras <- ampsDels(coad_sporadic_kras_freq_out)
coad_sporadic_kras_amp_bed <- reducingFreqBed(ampDels_coad_sporadic_kras[[1]], ampDels_coad_sporadic_kras[[2]])
coad_sporadic_kras_del_bed <- reducingFreqBed(ampDels_coad_sporadic_kras[[3]], ampDels_coad_sporadic_kras[[4]])

coad_sporadic_kras_cna_freq_out <- getFreqData(coad_sporadic_kras_cna2)
ampDels_coad_sporadic_kras_cna <- ampsDels(coad_sporadic_kras_cna_freq_out)
coad_sporadic_kras_cna_amp_bed <- reducingFreqBed(ampDels_coad_sporadic_kras_cna[[1]], ampDels_coad_sporadic_kras_cna[[2]])
coad_sporadic_kras_cna_del_bed <- reducingFreqBed(ampDels_coad_sporadic_kras_cna[[3]], ampDels_coad_sporadic_kras_cna[[4]])



a <- freqPlot(coad_sporadic_all_amp_bed, coad_sporadic_all_del_bed, speciesType = "human", main = "coad aneploidy Apc Tp53 Kras n = 100")
b <- freqPlot(coad_sporadic_tp53_amp_bed, coad_sporadic_tp53_del_bed, speciesType = "human", main = "coad aneuploidy Apc Tp53 n = 146")
c <- freqPlot(coad_sporadic_kras_amp_bed, coad_sporadic_kras_del_bed, speciesType = "human", main = "coad aneuploidy Apc Kras n = 80")
d <- freqPlot(coad_sporadic_all_cna_amp_bed , coad_sporadic_all_cna_del_bed , speciesType = "human", main = "coad cna Apc Tp53 Kras n = 100")
e <- freqPlot(coad_sporadic_tp53_cna_amp_bed, coad_sporadic_tp53_cna_del_bed, speciesType = "human", main = "coad cna Apc Tp53 n = 146")
f <- freqPlot(coad_sporadic_kras_cna_amp_bed, coad_sporadic_kras_cna_del_bed, speciesType = "human", main = "coad cna Apc Kras n = 80")

gridExtra::grid.arrange(a, b, c, d, e, f, ncol = 1)


### getting mouse hgsc version of this

absCn_BPRN <- noshadAbsCnFreqDf[grep(paste0(hgscPloidyCalls$sample[which(hgscPloidyCalls$geno == "BPRN")], collapse = "|"),
                                     noshadAbsCnFreqDf$sampleID),]

absCn_BPP <- noshadAbsCnFreqDf[grep(paste0(hgscPloidyCalls$sample[which(hgscPloidyCalls$geno == "BPP")], collapse = "|"),
                                     noshadAbsCnFreqDf$sampleID),]

absCn_BPN <- noshadAbsCnFreqDf[grep(paste0(hgscPloidyCalls$sample[which(hgscPloidyCalls$geno == "BPN")], collapse = "|"),
                                    noshadAbsCnFreqDf$sampleID),]

absCn_PRN <- noshadAbsCnFreqDf[grep(paste0(hgscPloidyCalls$sample[which(hgscPloidyCalls$geno == "PRN")], collapse = "|"),
                                    noshadAbsCnFreqDf$sampleID),]


armRes_bprn <- separateSegments_mV4(absCn_BPRN, minLength = 20000000, percentage = 0.60)
bprn_arm <- do.call(rbind, armRes_bprn[[1]])
bprn_cna <- do.call(rbind, armRes_bprn[[2]])
bprn_count <- do.call(rbind, armRes_bprn[[3]])
colnames(bprn_count) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")
bprn_count_df <- data.frame("samples" = rownames(bprn_count))
bprn_count_df <- cbind(bprn_count_df, bprn_count)
bprn_count_df <- melt(bprn_count)



bprn_arm_freq <- getFreqData(bprn_arm)
bprn_arm_res <- ampsDels(bprn_arm_freq)
bprn_arm_amp_bed <- reducingFreqBed(bprn_arm_res[[1]], bprn_arm_res[[2]])
bprn_arm_del_bed <- reducingFreqBed(bprn_arm_res[[3]], bprn_arm_res[[4]])

bprn_cna_freq <- getFreqData(bprn_cna)
bprn_cna_res <- ampsDels(bprn_cna_freq)
bprn_cna_amp_bed <- reducingFreqBed(bprn_cna_res[[1]], bprn_cna_res[[2]])
bprn_cna_del_bed <- reducingFreqBed(bprn_cna_res[[3]], bprn_cna_res[[4]])


armRes_bpp <- separateSegments_mV4(absCn_BPP, minLength = 20000000, percentage = 0.60)
bpp_arm <- do.call(rbind, armRes_bpp[[1]])
bpp_cna <- do.call(rbind, armRes_bpp[[2]])
bpp_count <- do.call(rbind, armRes_bpp[[3]])
colnames(bpp_count) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")
bpp_count_df <- data.frame("samples" = rownames(bpp_count))
bpp_count_df <- cbind(bpp_count_df, bpp_count)
bpp_count_df <- melt(bpp_count)


bpp_arm_freq <- getFreqData(bpp_arm)
bpp_arm_res <- ampsDels(bpp_arm_freq)
bpp_arm_amp_bed <- reducingFreqBed(bpp_arm_res[[1]], bpp_arm_res[[2]])
bpp_arm_del_bed <- reducingFreqBed(bpp_arm_res[[3]], bpp_arm_res[[4]])

bpp_cna_freq <- getFreqData(bpp_cna)
bpp_cna_res <- ampsDels(bpp_cna_freq)
bpp_cna_amp_bed <- reducingFreqBed(bpp_cna_res[[1]], bpp_cna_res[[2]])
bpp_cna_del_bed <- reducingFreqBed(bpp_cna_res[[3]], bpp_cna_res[[4]])

armRes_bpn <- separateSegments_mV4(absCn_BPN, minLength = 20000000, percentage = 0.60)
bpn_arm <- do.call(rbind, armRes_bpn[[1]])
bpn_cna <- do.call(rbind, armRes_bpn[[2]])
bpn_count <- do.call(rbind, armRes_bpn[[3]])
colnames(bpn_count) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")
bpn_count_df <- data.frame("samples" = rownames(bpn_count))
bpn_count_df <- cbind(bpn_count_df, bpn_count)
bpn_count_df <- melt(bpn_count)

bpn_arm_freq <- getFreqData(bpn_arm)
bpn_arm_res <- ampsDels(bpn_arm_freq)
bpn_arm_amp_bed <- reducingFreqBed(bpn_arm_res[[1]], bpn_arm_res[[2]])
bpn_arm_del_bed <- reducingFreqBed(bpn_arm_res[[3]], bpn_arm_res[[4]])

bpn_cna_freq <- getFreqData(bpn_cna)
bpn_cna_res <- ampsDels(bpn_cna_freq)
bpn_cna_amp_bed <- reducingFreqBed(bpn_cna_res[[1]], bpn_cna_res[[2]])
bpn_cna_del_bed <- reducingFreqBed(bpn_cna_res[[3]], bpn_cna_res[[4]])


armRes_prn <- separateSegments_mV4(absCn_PRN, minLength = 20000000, percentage = 0.60)
prn_arm <- do.call(rbind, armRes_prn[[1]])
prn_cna <- do.call(rbind, armRes_prn[[2]])
prn_count <- do.call(rbind, armRes_prn[[3]])
colnames(prn_count) <- c("aneuGain", "aneuLoss", "cnaGain", "cnaLoss")
prn_count_df <- data.frame("samples" = rownames(prn_count))
prn_count_df <- cbind(prn_count_df, prn_count)
prn_count_df <- melt(prn_count)

prn_arm_freq <- getFreqData(prn_arm)
prn_arm_res <- ampsDels(prn_arm_freq)
prn_arm_amp_bed <- reducingFreqBed(prn_arm_res[[1]], prn_arm_res[[2]])
prn_arm_del_bed <- reducingFreqBed(prn_arm_res[[3]], prn_arm_res[[4]])

prn_cna_freq <- getFreqData(prn_cna)
prn_cna_res <- ampsDels(prn_cna_freq)
prn_cna_amp_bed <- reducingFreqBed(prn_cna_res[[1]], prn_cna_res[[2]])
prn_cna_del_bed <- reducingFreqBed(prn_cna_res[[3]], prn_cna_res[[4]])


### for coad samples

coadPloidyCalls <- allPloidyCalls[which(allPloidyCalls$sample %in% coadNonzeroTc$sample),]
hgscNonzeroTc <- hgscNonzeroTc[which(hgscNonzeroTc$sample %in% hgscPloidyCalls$sample), ]

dupeIdx <- c(49, 51, 7, 27, 112)
hgscPloidyCalls <- hgscPloidyCalls[-dupeIdx, ]

hgscPloidyCalls$altPloidy <- hgscNonzeroTc$tc[match(hgscPloidyCalls$sample, hgscNonzeroTc$sample)]
hgscPloidyCalls[,2:ncol(hgscPloidyCalls)] <- lapply(hgscPloidyCalls[,2:ncol(hgscPloidyCalls)], function(x) as.numeric(x))
hgscPloidyCalls$geno <- annoTable_combined$geno[match(hgscPloidyCalls$sample, annoTable_combined$mouse_id)]
hgscPloidyCalls <- hgscPloidyCalls[-which(is.na(hgscPloidyCalls$geno)), ]
hgscPloidyCalls <- hgscPloidyCalls[-which(hgscPloidyCalls$geno == "unk"), ]
hgscPloidyCalls$sample[which(hgscPloidyCalls$geno == "BPRN")]

absCn_BPRN <- noshadAbsCnFreqDf[grep(paste0(hgscPloidyCalls$sample[which(hgscPloidyCalls$geno == "BPRN")], collapse = "|"),
                                     noshadAbsCnFreqDf$sampleID),]

absCn_BPP <- noshadAbsCnFreqDf[grep(paste0(hgscPloidyCalls$sample[which(hgscPloidyCalls$geno == "BPP")], collapse = "|"),
                                    noshadAbsCnFreqDf$sampleID),]

absCn_BPN <- noshadAbsCnFreqDf[grep(paste0(hgscPloidyCalls$sample[which(hgscPloidyCalls$geno == "BPN")], collapse = "|"),
                                    noshadAbsCnFreqDf$sampleID),]

absCn_PRN <- noshadAbsCnFreqDf[grep(paste0(hgscPloidyCalls$sample[which(hgscPloidyCalls$geno == "PRN")], collapse = "|"),
                                    noshadAbsCnFreqDf$sampleID),]




###
###
### graphs for cross species comparisons


hgsc_aneu <- freqPlot(bprn_arm_amp_bed, hgsc_arm_del_bed, speciesType = "human", main = "Aneuploidy: hgsc n = 331")
hgsc_brca_aneu <- freqPlot(brca_arm_amp_bed, brca_arm_del_bed, speciesType = "human", main = "hgsc w/ brca n = 33")
hgsc_cna <- freqPlot(hgsc_cna_amp_bed, hgsc_cna_del_bed, speciesType = "human", main = "CNA > 20Mb: hgsc n = 331")
hgsc_brca_cna <- freqPlot(brca_cna_amp_bed, brca_cna_del_bed, speciesType = "human", main = "hgsc w/ brca n = 33")

hrd_aneu <- freqPlot(bprn_arm_amp_bed, hrd_arm_del_bed, speciesType = "human", main = "Aneuploidy: hrd n = 152")
hrd_cna <- freqPlot(hrd_cna_amp_bed, hrd_cna_del_bed, speciesType = "human", main = "CNA > 20Mb: hrd n = 152")
hgsc_nonhrd_aneu <- freqPlot(nonhrd_arm_amp_bed, nonhrd_arm_del_bed, speciesType = "human", main = "hgsc w/ nonhrd n = 165")
hgsc_nonhrd_cna <- freqPlot(nonhrd_cna_amp_bed, nonhrd_cna_del_bed, speciesType = "human", main = "hgsc w/ nonhrd n = 165")


bprn_aneu <- freqPlot(bprn_arm_amp_bed,bprn_arm_del_bed, speciesType = "mouse", main = "hgsc BPRN n = 83")
bpp_aneu <- freqPlot(bpp_arm_amp_bed, bpp_arm_del_bed, speciesType = "mouse", main = "hgsc BPP n = 12")
bpn_aneu <- freqPlot(bpn_arm_amp_bed, bpn_arm_del_bed, speciesType = "mouse", main = "hgsc BPN n = 12")
prn_aneu <- freqPlot(prn_arm_amp_bed, prn_arm_del_bed, speciesType = "mouse", main = "hgsc PRN n = 10")

bprn_cna <- freqPlot(bprn_cna_amp_bed,bprn_cna_del_bed, speciesType = "mouse", main = "BPRN n = 83")
bpp_cna <- freqPlot(bpp_cna_amp_bed, bpp_cna_del_bed, speciesType = "mouse", main = "BPP n = 12")
bpn_cna <- freqPlot(bpn_cna_amp_bed, bpn_cna_del_bed, speciesType = "mouse", main = "BPN brca n = 12")
prn_cna <- freqPlot(prn_cna_amp_bed, prn_cna_del_bed, speciesType = "mouse", main = "PRN brca n = 10")

dev.off()
pdf(file = "/mnt/DATA5/tmp/kev/misc/20220608hgscFreqCompAneu.pdf")
gridExtra::grid.arrange(hgsc_aneu, hgsc_brca_aneu, bprn_aneu, bpp_aneu, bpn_aneu, prn_aneu, ncol = 1)
dev.off


dev.off()
pdf(file = "/mnt/DATA5/tmp/kev/misc/20220608hgscFreqCompCNA.pdf")
gridExtra::grid.arrange(hgsc_cna, hgsc_brca_cna, bprn_cna, bpp_cna, bpn_cna, prn_cna, ncol = 1)
dev.off()


dev.off()
pdf(file = "/mnt/DATA5/tmp/kev/misc/20220608hgscHrdFreqCompAneu.pdf")
gridExtra::grid.arrange(hrd_aneu, hgsc_nonhrd_aneu, bprn_aneu, bpp_aneu, bpn_aneu, prn_aneu, ncol = 1)
dev.off


dev.off()
pdf(file = "/mnt/DATA5/tmp/kev/misc/20220608hgscHrdFreqCompCNA.pdf")
gridExtra::grid.arrange(hrd_cna, hgsc_nonhrd_cna, bprn_cna, bpp_cna, bpn_cna, prn_cna, ncol = 1)
dev.off()

### general feature count per genotype doesn't do much. fga more valuable stat or maybe useful when comparing within a sample

# gi_bprn <- ggplot(bprn_count_df) + geom_boxplot(aes(x = Var1, y = value)) + ylim(c(0,20)) + ggtitle("BPRN n = 83") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# gi_bpp <- ggplot(bpp_count_df) + geom_boxplot(aes(x = Var1, y = value)) + ylim(c(0,20)) + ggtitle("BPP n = 12") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# gi_bpn <- ggplot(bpn_count_df) + geom_boxplot(aes(x = Var1, y = value)) + ylim(c(0,20)) + ggtitle("BPN n = 12") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# gi_prn <- ggplot(prn_count_df) + geom_boxplot(aes(x = Var1, y = value)) + ylim(c(0,20)) + ggtitle("PRN n = 10") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# gridExtra::grid.arrange(gi_bprn, gi_bpp, gi_bpn, gi_prn)

### fga

finalTable_hgsc <- finalTable2
finalTable_hgsc <- finalTable_hgsc[grep(paste(hgsc_anno2_tp53_names, collapse = "|"),
                                        finalTable_hgsc$Sample),]

finalTable_brca <- finalTable2
finalTable_brca <- finalTable_brca[grep(paste(hgsc_anno2_tp53_brca12_names, collapse = "|"),
                                        finalTable_brca$Sample),]

finalTable_hrd <- finalTable2
finalTable_hrd <- finalTable_hrd[grep(paste(hrd_hg19_tp53$Tumor[which(hrd_hg19_tp53$hrd.status == "HRD")],
                                            collapse = "|"), finalTable_hrd$Sample),]

finalTable_nonhrd <- finalTable2
finalTable_nonhrd <- finalTable_nonhrd[grep(paste(hrd_hg19_tp53$Tumor[which(hrd_hg19_tp53$hrd.status == "non-HRD")],
                                                  collapse = "|"), finalTable_nonhrd$Sample),]


fga_bpn <- fgaCalculator_amp(absCn_BPN)
fga_bprn <- fgaCalculator_amp(absCn_BPRN)
fga_bpp <- fgaCalculator_amp(absCn_BPP)
fga_prn <- fgaCalculator_amp(absCn_PRN)
# fga_wgs <- fgaCalculator_amp(allSamps2, mm10ChromSize2)
fga_hgsc <- fgaCalculator_tcga(finalTable_hgsc, tcga_ploidy)
fga_brca <- fgaCalculator_tcga(finalTable_brca, tcga_ploidy)
fga_hrd <- fgaCalculator_tcga(finalTable_hrd, tcga_ploidy)
fga_nonhrd <- fgaCalculator_tcga(finalTable_nonhrd, tcga_ploidy)



fga_bpn$type <- "bpn"
fga_bprn$type <- "bprn"
fga_bpp$type <- "bpp"
fga_prn$type <- "prn"
# fga_wgs$type <- "wgs"
fga_brca$type <- "brca"
fga_hgsc$type <- "hgsc"
fga_hrd$type <- "hrd"
fga_nonhrd$type <- "nonhrd"

fga_all <- rbind(fga_bpn, fga_bprn, fga_bpp,
                 fga_prn, fga_hrd,fga_nonhrd)

colnames(fga_all)[1:2] <- c("sample", "fga")
fga_all$type <- factor(fga_all$type)

library(tidyr)
library(rstatix)

stat_test <- fga_all %>% t_test(fga ~ type)
stat_test <- stat_test %>% add_xy_position(x = "type")

ggboxplot(fga_all, x = "type", y = "fga") +
  stat_pvalue_manual(stat_test, label = "p.adj.signif", tip.length = 0.01)


### calculating things for coad so then I can just try and fit story from there


