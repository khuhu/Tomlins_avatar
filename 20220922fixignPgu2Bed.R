### fixing bed file for pangu2

library(GenomicRanges)

pangu <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD203665_173_Designed.gc.bed", sep = "\t",
                    stringsAsFactors = FALSE, header = FALSE)
pangu2 <- read.table("/home/kevhu/data/bedFiles/WG_IAD127899.20170720.designed.forscript.GC.bed",
                     sep = "\t", stringsAsFactors = FALSE, header = FALSE)

pangu_reg <- pangu[grep("REGION", pangu$V8),]

panguregGR <- GRanges(seqnames = pangu_reg$V1, IRanges(start = pangu_reg$V2, end = pangu_reg$V3))
pangu2GR <- GRanges(seqnames = pangu2$V1, IRanges(start = pangu2$V2, end = pangu2$V3))


pgu2Idx <- subjectHits(findOverlaps(panguregGR, pangu2GR))
pgu2Idx <- pgu2Idx[-which(duplicated(queryHits(findOverlaps(panguregGR, pangu2GR))))]


### weird how it's off by one index? maybe one of the new regions doesn't have a match
### so every region minus 778 is not in the previous bed file

tmpIdx <- data.frame(findOverlaps(panguregGR, pangu2GR))
(1:1072)[-which(1:1072 %in% tmpIdx$queryHits)]

pangu_reg$string <- paste0(pangu_reg$V1, pangu_reg$V2, pangu_reg$V3)

pangu$string <- paste0(pangu$V1, pangu$V2, pangu$V3)
pangu$V9 <- pangu$V8
pangu$V9[which(pangu$string %in% pangu_reg$string[c(1:777,779:1072)])] <- pangu2$V8[pgu2Idx]
pangu$V9[which(pangu$V9 == "REGION_2411")]  <- "BRCA2"

pangu$V9[grep("VHL", pangu$V9)] <- "VHL"

write.table(pangu[,c(1:7, 10)], "/home/kevhu/data/bedFiles/IAD203665_173_Designed_fixed_KH.gc.bed", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(pangu[,c(1:7, 10)], "/mnt/DATA6/mouseData/bedFiles/IAD203665_173_Designed_fixed_KH.gc.bed", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)


mouseBed <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed", sep = "\t",
                       header = FALSE)
mouseBed2 <- mouseBed[-grep("SNP", mouseBed$V8), ]
