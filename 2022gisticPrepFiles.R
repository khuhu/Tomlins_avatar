#inflating the markers by double - because the start and end need to be a marker otherwise gistic won't match the
#start and ends of the segments

bedFile <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.gc.noChr.bed", sep = "\t", 
                      header = FALSE)

markerFile_start <- bedFile[,c(1, 2)]
markerFile_start <- cbind(paste0(bedFile$V4, "_s"), markerFile_start)
markerFile_end <- bedFile[,c(1, 3)]
markerFile_end <- cbind(paste0(bedFile$V4, "_e"), markerFile_end)
colnames(markerFile_start) <- c("V1", "V2", "V3")
colnames(markerFile_end) <- c("V1", "V2", "V3")

markerFile <- rbind(markerFile_start, markerFile_end)
#markerFile$V1 <- str_replace(markerFile$V1, "X", "20")
markerFile <- markerFile[order(markerFile$V2, markerFile$V3),]
markerFile <- markerFile[-which(markerFile$V2 == "X"),]


#write.table(markerFile, "/mnt/DATA5/tmp/kev/programs/GISTIC2/testfiles/20220103IAD202670_167.markerFiled.txt", sep = "\t",
#            col.names = FALSE, row.names = FALSE, quote = FALSE)


segColnames <- c("ID", "chrom", "loc.start", "loc.end",
                 "num.mark", "seg.mean")
mouseNormal <- c("EF_D03_MG_X14", "MG_18X50", "MG_21X53", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

testSeg141 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/segResults.txt",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
testSeg141 <- testSeg141[, segColnames]


testSeg142 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/segResults.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
testSeg142 <- testSeg142[, segColnames]


testSeg76 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-76-MG_test1_255_185/segResults.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
testSeg76 <- testSeg76[, segColnames]

testSeg138 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/segResults.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
testSeg138 <- testSeg138[, segColnames]

testSeg156 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/segResults.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
testSeg156 <- testSeg156[, segColnames]

testSeg157 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384/segResults.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
testSeg157 <- testSeg157[, segColnames]


testSegFile <- rbind(testSeg76, testSeg138, testSeg141, testSeg142, testSeg156, testSeg157)

# testSegFile$chrom <- str_replace(testSegFile$chrom, "23", "X")

testSegFile <- testSegFile[-which(testSegFile$chrom == "20"), ]
testSegFile <- testSegFile[-which(testSegFile$chrom == "23"), ]
testSegFile <- testSegFile[-which(testSegFile$ID %in% mouseNormal), ]

### the sample names from 20220104segmentalHeatmaps script
testSegFile$ID <- str_remove(nameStripper(testSegFile$ID), "x.*")

testSegFile_hgsc <- testSegFile[which(testSegFile$ID %in% hgscSampleNames),]
testSegFile_crc <- testSegFile[which(testSegFile$ID %in% crcSampleNames),]
testSegFile_adenoma <- testSegFile[which(testSegFile$ID %in% adenomaSampleNames),]


# crc_dup <- paste0(testSegFile_crc$ID, testSegFile_crc$chrom, testSegFile_crc$loc.start,
#                   testSegFile_crc$loc.end, testSegFile_crc$num.mark, testSegFile_crc$seg.mean)
# testSegFile_crc <- testSegFile_crc[-which(duplicated(crc_dup)),]

# this is just a quick fix .... doesnt work for multiple samples with overlaps --- can do this at thre report level 
# if need be
rownames(testSegFile_hgsc) <- NULL
testSegFile_hgsc <- testSegFile_hgsc[-c(1773:1810),]

write.table(testSegFile_crc, "/mnt/DATA5/tmp/kev/programs/GISTIC2/testfiles/20220103allMouseCrcSeg.txt", sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(testSegFile_hgsc, "/mnt/DATA5/tmp/kev/programs/GISTIC2/testfiles/20220103allMouseHgscSeg.txt", sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(testSegFile_adenoma, "/mnt/DATA5/tmp/kev/programs/GISTIC2/testfiles/20220103allMouseAdenomaSeg.txt", sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

