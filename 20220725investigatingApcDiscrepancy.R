source("/home/kevhu/scripts/20210802syntenyFunctions.R")

annotationList <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20211207yingSecondSet.xlsx")
annotationList$sampleStripped <- str_remove(nameStripper(annotationList$`#`), "\\-")
annotationList <- annotationList[-which(annotationList$sampleStripped == "efd36"),]
annotationList$sampleStripped[which(annotationList$sampleStripped %in% paste0("efd", c(1,2,4:7)))] <- paste0("efd0", c(1,2,4:7))


histologyList <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/AKP&CDX2Braftumor.xlsx", skip = 1)
histologyList$sampleStripped <- str_remove(nameStripper(histologyList$Sample), "\\-")
histologyList <- histologyList[1:61,]


fearonVar1 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-76-MG_test1_255_185_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonVar2 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonVar3 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-156-MG_Fearon_20210809_374_382_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonVar4 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)


allAnno <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220725sampleYingBasicAnno.xlsx")
allAnno$stripped <- str_remove(nameStripper(allAnno$efid), "\\-")
allAnno$stripped[1:9] <- paste0("efd0", 1:9)

allmuts <- rbind(fearonVar1, fearonVar2, fearonVar3, fearonVar4)
allmuts$stripped <- str_remove(str_remove(nameStripper(allmuts$Sample), "x.*"), "\\-")

allmuts_filt <- allmuts[grep("efd", allmuts$stripped),]

allmuts_filt$genotype <- allAnno$genotype[match(allmuts_filt$stripped, allAnno$stripped)]

combinedVars <- allmuts_filt


#### wrote file for Aaron
# write.table(allmuts_filt, "/mnt/DATA5/tmp/kev/misc/20220725fearonAllMutAnno.txt", sep = "\t", 
#             col.names = TRUE, row.names = FALSE, quote = FALSE)


### getting allele Afs

newAnnoTable <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220725fearonGemmsAfDelFreq.xlsx")
newAnnoTable$stripped <- str_remove(str_remove(nameStripper(newAnnoTable$`#`), "x.*"), "\\-")
newAnnoTable$stripped[1:9] <- paste0("efd0", 1:9)

mgpVars <- c("hom", "het", "unknown")
combinedVars <- combinedVars[which(combinedVars$mm10_mpgpv6_Indels == ""),]

fdpFilt <- which(combinedVars$FDP > 100)
faoFilt <- which(combinedVars$FAO > 10)
freqFilt <- which(combinedVars$AF > 0.10)
hrunFilt <- which(combinedVars$HRUN < 4)
strandRatio <- intersect(which(combinedVars$FSAF/combinedVars$FSAR > 0.2),
                         which(combinedVars$FSAF/combinedVars$FSAR < 5))
qualFilt <- which(combinedVars$QUAL >= 100)

goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, qualFilt, hrunFilt))
combinedVars_goodsamps <- combinedVars[goodSamps,]
combinedVars_goodsamps_exon <- combinedVars_goodsamps[which(combinedVars_goodsamps$Func.refGene == "exonic"), ]
combinedVars_goodsamps_exon  <- combinedVars_goodsamps_exon[grep("EF", combinedVars_goodsamps_exon$Sample),]

combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[which(combinedVars_goodsamps_exon$name == "."),]                                         
combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[-grep("nonframeshift",
                                                                 combinedVars_goodsamps_exon$ExonicFunc.refGene),]  
combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[-grep("^synonymous SNV",
                                                                 combinedVars_goodsamps_exon$ExonicFunc.refGene),]  

combinedVars_goodsamps_exon$Sample <- str_remove(str_remove(nameStripper(combinedVars_goodsamps_exon$Sample), "\\-"), "x.*")

goodSamps_anno <- annotationList$sampleStripped[-which(annotationList$Qual == "Bad")]

combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[which(combinedVars_goodsamps_exon$Sample %in% goodSamps_anno),]


mutOfInterest <- c("Kras:NM_021284:exon2:c.G35A:p.G12D", "Trp53:NM_001127233:exon8:c.G809A:p.R270H,Trp53:NM_011640:exon8:c.G809A:p.R270H")
mutDf_kras <- combinedVars[grep(mutOfInterest[1], combinedVars$AAChange.refGene),]
mutDf_p53 <- combinedVars[grep(mutOfInterest[2], combinedVars$AAChange.refGene),]

combinedVars_somatic <- combinedVars_goodsamps_exon[-grep(paste(mutOfInterest, collapse = "|"), combinedVars_goodsamps_exon$AAChange.refGene),]

newAnnoTable$KRASG12D <- signif(mutDf_kras$AF, digits = 2)[match(newAnnoTable$stripped, mutDf_kras$stripped)]
newAnnoTable$TRP53R270H <- signif(mutDf_p53$AF, digits = 2)[match(newAnnoTable$stripped, mutDf_p53$stripped)]



### getting dels

mouseNormal <- c("EF_D03_MG_X14", "MG_18X50", "MG_21X53", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

cnCalls_1 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-138-MG_cho_20210621_354_343/cnMatrix_gene.txt",
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
cnCalls_2 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/cnMatrix_gene.txt",
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
cnCalls_3 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-142-MG_cho_20210701_357_353/cnMatrix_gene.txt",
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
cnCalls_4 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-76-MG_test1_255_185/cnMatrix_gene.txt",
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
cnCalls_5 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/cnMatrix_gene.txt",
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
cnCalls_6 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384/cnMatrix_gene.txt",
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
cnCalls_7 <- read.table("/mnt/DATA5/tmp/kev/mouseDels/mousePanelDelCalls/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/cnMatrix_gene.txt",
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)

all_Apc <- cbind(cnCalls_1[grep("Apc$", cnCalls_1$Gene),], cnCalls_2[grep("Apc$", cnCalls_2$Gene),2:ncol(cnCalls_2)],
                     cnCalls_3[grep("Apc$", cnCalls_3$Gene),2:ncol(cnCalls_3)], cnCalls_4[grep("Apc$", cnCalls_4$Gene),2:ncol(cnCalls_4)],
                     cnCalls_5[grep("Apc$", cnCalls_5$Gene),2:ncol(cnCalls_5)], cnCalls_6[grep("Apc$", cnCalls_6$Gene),2:ncol(cnCalls_6)],
                     cnCalls_7[grep("Apc$", cnCalls_7$Gene),2:ncol(cnCalls_7)])


all_cnCalls <- cbind(cnCalls_1[grep("Del", cnCalls_1$Gene),], cnCalls_2[grep("Del", cnCalls_2$Gene),2:ncol(cnCalls_2)],
                     cnCalls_3[grep("Del", cnCalls_3$Gene),2:ncol(cnCalls_3)], cnCalls_4[grep("Del", cnCalls_4$Gene),2:ncol(cnCalls_4)],
                     cnCalls_5[grep("Del", cnCalls_5$Gene),2:ncol(cnCalls_5)], cnCalls_6[grep("Del", cnCalls_6$Gene),2:ncol(cnCalls_6)],
                     cnCalls_7[grep("Del", cnCalls_7$Gene),2:ncol(cnCalls_7)])
all_cnCalls <- all_cnCalls[,-which(duplicated(colnames(all_cnCalls)))]

trp53Del <- unlist(log2(all_cnCalls[grep("Trp53Del", all_cnCalls$Gene), 2:ncol(all_cnCalls)]))
names(trp53Del) <- str_remove(nameStripper(names(trp53Del)), "x.*")

apcDel <- unlist(log2(all_cnCalls[grep("ApcDel", all_cnCalls$Gene), 2:ncol(all_cnCalls)]))
names(apcDel) <- str_remove(nameStripper(names(apcDel)), "x.*")


apcNonDel <- unlist(log2(all_Apc[, 2:ncol(all_Apc)]))
names(apcNonDel) <- str_remove(nameStripper(names(apcNonDel)), "x.*")


newAnnoTable$trp53Del <- trp53Del[match(newAnnoTable$stripped, names(trp53Del))]
newAnnoTable$apcDel <- apcDel[match(newAnnoTable$stripped, names(apcDel))]
newAnnoTable$apcNonDel <- apcNonDel[match(newAnnoTable$stripped, names(apcNonDel))]


tc <- 1 - apply(all_cnCalls[grep("Trp53Del", all_cnCalls$Gene), 2:ncol(all_cnCalls)], 2, min)
tcApc <- 1 - all_cnCalls[grep("ApcDel", all_cnCalls$Gene), 2:ncol(all_cnCalls)]


### calculating expected cnr of 1 copy loss given the tumor content from homozygous deletion
### for labeling mechanisms of loss setting minimum TC to 30% - thought process behind this is
### below 30% tumor content, the expected 1 copy loss reaches levels of CNR where it can be noise


newAnnoTable$predictedTc <- 1 - 2^newAnnoTable$apcDel

newAnnoTable$expectedChr18Loss <- log2((0.5 * ((1 - 2^newAnnoTable$apcDel) * 100) - 100) / -100) # based on (one loss cnr (0.5))/100 = 1 - x

# newAnnoTable$expectedChr18Loss<- log2(0.5) * (1 - 2^newAnnoTable$apcDel)
# newAnnoTable$expectedChr18Loss2 <- log2((0.5 * ((1 - 2^newAnnoTable$apcDel) * 100) - 100) / -100)

newAnnoTable$chr18Status <- "copy neutral LOH"
### 0.2 cutoff made by graphing all samples, a lot of small change but majority between 0 to 0.2
newAnnoTable$chr18Status[which(abs(newAnnoTable$expectedChr18Loss - newAnnoTable$apcNonDel) < 0.2 & newAnnoTable$apcNonDel < -0.2)] <- "copy loss LOH"
newAnnoTable$chr18Status[which(newAnnoTable$predictedTc < 0.3)] <- "TC too low"
newAnnoTable$chr18Status[which(is.na(newAnnoTable$KRASG12D))] <- "Bad Sample"
newAnnoTable$chr18Status[grep("Wild|wild", newAnnoTable$Genotype)] <- "Normal"


newAnnoTable$predictedTc[grep("Wild|wild", newAnnoTable$Genotype)] <- 0
newAnnoTable$predictedTc[which(newAnnoTable$predictedTc < 0)] <- 0


newAnnoTable_noR270H <- newAnnoTable[-grep("R270H|Wild|wild", newAnnoTable$Genotype),]
newAnnoTable_noR270H <- newAnnoTable_noR270H[grep("p53", newAnnoTable_noR270H$Genotype),]


newAnnoTable_cutoff <- newAnnoTable[-grep("Braf|Wild|wild", newAnnoTable$Genotype),]
hist(abs(newAnnoTable_cutoff$expectedChr18Loss - newAnnoTable_cutoff$apcNonDel))

cor(newAnnoTable_noR270H$trp53Del[-which(is.na(newAnnoTable_noR270H$trp53Del))],
    newAnnoTable_noR270H$apcDel[-which(is.na(newAnnoTable_noR270H$trp53Del))])

write.table(newAnnoTable, "/mnt/DATA5/tmp/kev/misc/20220726FearonAnnoAfDelFreq.txt", sep = "\t",
            row.names = FALSE,  col.names = TRUE, quote = FALSE)


# 
# newAnnoTableV2 <- newAnnoTable
# newAnnoTableV2$KRASG12D_expected <- NA
# newAnnoTableV2$KRASG12D_observed <- NA
# newAnnoTableV2$BRAFV600_expected <- NA
# newAnnoTableV2$BRAFV600_observed <- NA
# newAnnoTableV2$TRP53R270H_expected <- NA
# newAnnoTableV2$TRP53R270H_observed <- NA
# newAnnoTableV2$TRP53_flox_expected <- NA
# newAnnoTableV2$TRP53_flox_observed <- NA
# newAnnoTableV2$Apc_flox_expected <- NA
# newAnnoTableV2$Apc_flox_observed <- NA
# 
# 
# newAnnoTableV2$KRASG12D_expected <- ifelse(grepl("G12D", newAnnoTableV2$Genotype), "+", "-")
# newAnnoTableV2$KRASG12D_observed <- ifelse(newAnnoTableV2$KRASG12D > 0.1, "+", "-")
# newAnnoTableV2$BRAFV600_expected <- ifelse(grepl("Braf", newAnnoTableV2$Genotype), "+", "-")
# newAnnoTableV2$BRAFV600_observed <- ifelse(newAnnoTableV2$BRAFV600_expected == "+", "+", "-") # change this to the vaf later, not enough time for last minute change
# newAnnoTableV2$TRP53R270H_expected <- ifelse(grepl("R270H", newAnnoTableV2$Genotype), "+", "-")
# newAnnoTableV2$TRP53R270H_observed <- ifelse(newAnnoTableV2$TRP53R270H > 0.1, "+", "-")
# newAnnoTableV2$TRP53_flox_expected <- ifelse(grepl("p53", newAnnoTableV2$Genotype), "+", "-")
# newAnnoTableV2$TRP53_flox_observed <- ifelse(newAnnoTableV2$trp53Del < -0.2, "+", "-")
# newAnnoTableV2$Apc_flox_expected <- ifelse(grepl("Apc", newAnnoTableV2$Genotype), "+", "-")
# newAnnoTableV2$Apc_flox_observed <- ifelse(newAnnoTableV2$apcDel < -0.2, "+", "-")
# 
# 
# newAnnoTableV3 <- newAnnoTableV2[, 15:24]
# 
# 
# res <- cbind(newAnnoTableV2$`Sample name`, newAnnoTableV2$Genotype, newAnnoTableV3)
# colnames(res)[1:2] <- c("sample", "genotype")
# 
# write.table(res, "/mnt/DATA5/tmp/kev/misc/20220726FearonAnnoAfDelFreqSimp.txt", sep = "\t",
#             row.names = FALSE,  col.names = TRUE, quote = FALSE)

