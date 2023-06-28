### idea is to use any and all snps because if GEMM is composed of more than one GEMM,
### homs can become het

statsTab1 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-138-MG_cho_20210621_354_343/plugin_out/coverageAnalysis_out.668/Auto_MG_cho_20210621_eros_343.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

statsTab2 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/plugin_out/coverageAnalysis_out.681/Auto_MG_cho_202106_3TS_eros_349.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

statsTab3 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-142-MG_cho_20210701_357_353/plugin_out/coverageAnalysis_out.693/Auto_MG_cho_20210701_eros_353.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)


statsTab4 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/plugin_out/coverageAnalysis_out.577/Auto_MG_EFD4_BBN_eros_304.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)


statsTab5 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-76-MG_test1_255_185/plugin_out/coverageAnalysis_out.303/Auto_MG_test1_eros_185.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

statsTab6 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/plugin_out/coverageAnalysis_out.757/Auto_MG_Fearon_20210809_eros_382.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
statsTab7 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/plugin_out/coverageAnalysis_out.577/Auto_MG_EFD4_BBN_eros_304.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)


allStats <- rbind(statsTab1, statsTab2, statsTab3, statsTab4, statsTab5, statsTab6, statsTab7)
allStats$On.Target <- as.numeric(str_remove(allStats$On.Target, "%"))
allStats$Uniformity <- as.numeric(str_remove(allStats$Uniformity, "%"))
badSamps <- allStats$Sample.Name[which(allStats$On.Target < 80 | allStats$Mean.Depth < 100 | allStats$Uniformity < 80)]
badSamps



varTab1 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-138-MG_cho_20210621_354_343_anno.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
varTab2 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349_anno.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
varTab3 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-142-MG_cho_20210701_357_353_anno.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
varTab4 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-76-MG_test1_255_185_anno.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
varTab5 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_anno.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
varTab6 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-156-MG_Fearon_20210809_374_382_anno.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
varTab7 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384_anno.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

allVars <- rbind(varTab1, varTab2, varTab3,
                 varTab4, varTab4, varTab5,
                 varTab6, varTab7)


fdpFilt <- which(allVars$FDP > 20)
faoFilt <- which(allVars$FAO > 5)
freqFilt <- which(allVars$AF > 0.10)
hrunFilt <- which(allVars$HRUN < 4)
strandRatio <- intersect(which(allVars$FSAF/allVars$FSAR > 0.2),
                         which(allVars$FSAF/allVars$FSAR < 5))
qualFilt <- which(allVars$QUAL >= 100)

goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, qualFilt))
allVars_goodsamps <- allVars[goodSamps,]

allVars_goodsamps$tmpString <- paste0(allVars_goodsamps$Chr, allVars_goodsamps$Start,
                                      allVars_goodsamps$End)

stringNames <- names(table(allVars_goodsamps$tmpString))
stringNames <- stringNames[which(table(allVars_goodsamps$tmpString) > length(unique(allVars_goodsamps$Sample)) * .10)]
allVars_goodsamps <- allVars_goodsamps[which(allVars_goodsamps$tmpString %in% stringNames),]


### 1009 hom SNPs + 5 het snps - ref may whether het or hom may change based on mouse
### TC independent for reference i.e indepedent on reference het or hom status
### 
tmpVars <- allVars_goodsamps[which(allVars_goodsamps$mm10_mpgpv6_Indels == "hom"),]
unique(tmpVars$tmpString)
tmpVars2 <- allVars_goodsamps[which(allVars_goodsamps$mm10_mpgpv6_Indels == "het"),]
unique(tmpVars2$tmpString)

tmpVars3 <- allVars_goodsamps[which(allVars_goodsamps$mm10_mpgpv6_Indels == "het" | 
                                      allVars_goodsamps$mm10_mpgpv6_Indels == "hom"),]

### mouse SNPs are annotated SNPs from MGP, so I'll use those - i.e verified

tmpVars4 <- tmpVars3[which(tmpVars3$Ref == "-"),]
table(tmpVars4$tmpString)
tmpVars5 <- tmpVars3[which(tmpVars3$Alt == "-"),]
table(tmpVars5$tmpString)

### 986 useable SNPs annotated on MGP
tmpVars6 <- tmpVars3[-which(tmpVars3$Ref == "-" |
                              tmpVars3$Alt == "-"),]


### 2021004: I can just use the annotation SNP file I made for annovar 
annovarTable <- read.table("/mnt/DATA6/mouseData/mm10_mgp.v6.combinedMouseFilt.txt", sep = "\t",
                           header = FALSE, stringsAsFactors = FALSE)

snpFilt <- c("A", "C", "T", "G")
annovarTable2 <- annovarTable[-which(annovarTable$V4 == "-" |
                              annovarTable$V5 == "-"),]
annovarTable2 <- annovarTable2[-which(annovarTable2$V7 < 100),]
annovarTable2 <- annovarTable2[which(annovarTable2$V4 %in% snpFilt),]
annovarTable2 <- annovarTable2[which(annovarTable2$V5 %in% snpFilt),]

tmpVcf <- vcfR::read.vcfR("/mnt/DATA6/mouseData/tvcOut/TSVC_variants.vcf")
tmpVcf2 <- data.frame(tmpVcf@fix, stringsAsFactors = FALSE)
tmpVcf2 <- tmpVcf2[which(tmpVcf2$FILTER == "PASS"),]

tmpVcf3 <- NULL
strSplitRes <- strsplit(tmpVcf2$INFO, ";")
for (i in 1:length(strSplitRes)) {
  tmpVector <- unlist(strSplitRes[[i]])[c(1:5,8,9)]
  tmpVector <- str_remove(tmpVector, ".*=")
  tmpVector2 <- unlist(c(tmpVcf2[i,1:7], tmpVector))
  tmpVcf3 <- rbind(tmpVcf3, tmpVector2)
}
tmpVcf3 <- data.frame(tmpVcf3, stringsAsFactors = FALSE)
rownames(tmpVcf3) <- NULL
colnames(tmpVcf3) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                       "FILT", "AF", "AO", "DP", "FAO", "FDP",
                       "FSAF", "FSAR")
tmpVcf3[,c(6, 8:14)] <- lapply(tmpVcf3[,c(6, 8:14)], as.numeric)
