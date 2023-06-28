source("/home/kevhu/scripts/20210802syntenyFunctions.R")
sum_1 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-138-MG_cho_20210621_354_343/plugin_out/coverageAnalysis_out.668/Auto_MG_cho_20210621_eros_343.bc_summary.xls",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
sum_2 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/plugin_out/coverageAnalysis_out.681/Auto_MG_cho_202106_3TS_eros_349.bc_summary.xls",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
sum_3 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-142-MG_cho_20210701_357_353/plugin_out/coverageAnalysis_out.693/Auto_MG_cho_20210701_eros_353.bc_summary.xls",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
sum_4 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-76-MG_test1_255_185/plugin_out/coverageAnalysis_out.303/Auto_MG_test1_eros_185.bc_summary.xls",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")
sum_5 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/plugin_out/coverageAnalysis_out.577/Auto_MG_EFD4_BBN_eros_304.bc_summary.xls",
                    header = TRUE, stringsAsFactors = FALSE, sep = "\t")
sum_6 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/plugin_out/coverageAnalysis_out.757/Auto_MG_Fearon_20210809_eros_382.bc_summary.xls",
                    header = TRUE, stringsAsFactors = FALSE, sep = "\t")
sum_7 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384/plugin_out/coverageAnalysis_out.761/Auto_MG_Fearon_20210809_2_eros_384.bc_summary.xls",
                    header = TRUE, stringsAsFactors = FALSE, sep = "\t")

all_sum <- rbind(sum_1, sum_2, sum_3, sum_4, sum_5, sum_6, sum_7)


# tcDf  <- read.table("/mnt/DATA5/tmp/kev/misc/20211101hgscTcDf_median.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


mouseNormal <- c("EF_D03_MG_X14", "MG_18X50", "MG_21X53", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

mouseNormal <- nameStripper(mouseNormal)

# amplicon_1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-138-MG_cho_20210621_354_343_median/gcCorrectedCounts_matrix.txt",
#                            header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
# amplicon_2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349_median/gcCorrectedCounts_matrix.txt",
#                            header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
# amplicon_3 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-142-MG_cho_20210701_357_353_median/gcCorrectedCounts_matrix.txt",
#                            header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
# amplicon_4 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-76-MG_test1_255_185_median/gcCorrectedCounts_matrix.txt",
#                            header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)

amplicon_1 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/mouseAmplicons_tc.txt",
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
amplicon_2 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/mouseAmplicons_tc.txt",
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
amplicon_3 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/mouseAmplicons_tc.txt",
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
amplicon_4 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-76-MG_test1_255_185/mouseAmplicons_tc.txt",
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
amplicon_5 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/mouseAmplicons_tc.txt",
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
amplicon_6 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384/mouseAmplicons_tc.txt",
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
amplicon_7 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/mouseAmplicons_tc.txt",
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)




fieldNames <- c("AmpliconId", "AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene")
filtNames <- c(fieldNames, "NumGC", "Length", "GC", "MinIndex", "MaxIndex", "NumProbes", "Label", "GeneNum", "Color",
               "Weights", "TotalPool")

paraDf <- amplicon_1[, fieldNames]
amplicon_1 <- amplicon_1[,-which(colnames(amplicon_1) %in% filtNames)]
amplicon_2 <- amplicon_2[,-which(colnames(amplicon_2) %in% filtNames)]
amplicon_3 <- amplicon_3[,-which(colnames(amplicon_3) %in% filtNames)]
amplicon_4 <- amplicon_4[,-which(colnames(amplicon_4) %in% filtNames)]
amplicon_5 <- amplicon_5[,-which(colnames(amplicon_5) %in% filtNames)]
amplicon_6 <- amplicon_6[,-which(colnames(amplicon_6) %in% filtNames)]
amplicon_7 <- amplicon_7[,-which(colnames(amplicon_7) %in% filtNames)]


all_amp <- cbind(paraDf, amplicon_1, amplicon_2, amplicon_3, amplicon_4, amplicon_5, amplicon_6, amplicon_7)
all_amp_normal <- all_amp[,c(1:6, which(colnames(all_amp) %in% mouseNormal))]
all_amp_normal <- all_amp_normal[,1:13]

all_amp <- all_amp[,-which(colnames(all_amp) %in% mouseNormal)]
colnames(all_amp)[7:ncol(all_amp)] <- nameStripper(colnames(all_amp)[7:ncol(all_amp)])
colnames(all_amp_normal)[7:ncol(all_amp_normal)] <- nameStripper(colnames(all_amp_normal)[7:ncol(all_amp_normal)])


all_amp[,7:ncol(all_amp)] <- log2(all_amp[,7:ncol(all_amp)])
all_amp$ChromNum <- as.numeric(str_replace(all_amp$ChromNum, "23", "20"))


write.table(all_amp, "/mnt/DATA5/tmp/kev/misc/20211206all_mouse_amps_tc.txt", sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)

# write.table(all_amp, "/mnt/DATA5/tmp/kev/misc/20211101hgsc_amps.txt", sep = "\t", quote = FALSE, col.names = TRUE,
#             row.names = FALSE)

# write.table(all_amp, "/mnt/DATA5/tmp/kev/misc/20211101hgsc_amps_median.txt", sep = "\t", quote = FALSE, col.names = TRUE,
#             row.names = FALSE)

all_amp_normal[,7:ncol(all_amp_normal)] <- log2(all_amp_normal[,7:ncol(all_amp_normal)])
all_amp_normal$ChromNum <- as.numeric(str_replace(all_amp_normal$ChromNum, "23", "20"))
write.table(all_amp_normal, "/mnt/DATA5/tmp/kev/misc/20211101normals_amps.txt", sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)




### testing for strata binning and norm

strataBinning <- read.table("/mnt/DATA6/mouseData/20211028strataPythonBinning.txt", sep = "\t", stringsAsFactors = FALSE,
                            header = TRUE)

poolBeds <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-142-MG_cho_20210701_357_353/plugin_out/coverageAnalysis_out.693/local_beds/IAD202670_167_Designed.gc.bed",
                       skip = 1, header = FALSE, sep = "\t")


strataBinning$ChromNum <- str_replace(strataBinning$ChromNum, "23", "X")
strataBinning$string <- paste0(strataBinning$ChromNum, strataBinning$StartPos, strataBinning$EndPos)

poolBeds$string <- paste0(str_remove(poolBeds$V1, "chr"), poolBeds$V2, poolBeds$V3)
which(strataBinning$string != poolBeds$string)
poolBeds$pool <- str_remove(poolBeds$V5, ".*Pool=")



for (i in unique(strataBinning$gcbin)) {
  print(paste(i, min(strataBinning$GC[which(strataBinning$gcbin == i)]),
                              max(strataBinning$GC[which(strataBinning$gcbin == i)])))
}

library(pracma)
ampGc142 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/amplicon.GCinput.txt", sep = "\t",
                       stringsAsFactors = FALSE, header = TRUE)
linBreaksGc <- linspace(min(ampGc142$GC), max(ampGc142$GC), 5)
linBreaksLen <- linspace(min(ampGc142$Length), max(ampGc142$Length), 5)


numpy_digitize <- function(vec, breaks){
  res <- rep(0, length(vec))
  for (i in seq_along(breaks)) {
    if (i < length(breaks)) {
      res[which(vec >= breaks[i] & vec < breaks[i + 1])] <- i
      } else {
        res[which(vec == breaks[i])] <- i
      }
    }
  return(res)
}


ampGc142$gcbin <- numpy_digitize(ampGc142$GC, linBreaksGc)
ampGc142$lenbin <- numpy_digitize(ampGc142$Length, linBreaksLen)


length(which(paste0(ampGc142$lenbin, ampGc142$gcbin) == paste0(strataBinning$lenbin, strataBinning$gcbin)))


### comparing old method and new



tmp <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/combinedCalls.txt",
                  sep = "\t", stringsAsFactors = FALSE, header = TRUE)
corList <- NULL
for (i in unique(tmp$Sample)) {
  tmp2 <- tmp[which(tmp$Sample == i),]
  tmpCor <- cor(tmp2$CopyNumberRatio[match(geneEst$Gene, tmp2$Gene)], geneEst[[i]])
  corList <- c(corList, tmpCor)
}








