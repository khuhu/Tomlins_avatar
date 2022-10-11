#!/usr/bin/env Rscript

library(stringr)

### quick way to find directories that are DNA; i.e if RNA will not have cnAmplicon_matrix.txt

tmpBed <- read.table("/mnt/DATA6/mouseData/mouseBedsIdx.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
dnaDirs <- system("find /mnt/DATA6/mouseData/copynumber/ -type f -name '*cnAmplicon_matrix.txt'", intern = TRUE)
dnaDirs <- str_remove(dnaDirs, "\\/cnAmplicon_matrix.txt")
dnaDirs <- str_remove(dnaDirs, "\\/mnt\\/DATA6\\/mouseData\\/copynumber\\/")

tmpBed2 <- tmpBed[which(tmpBed$reports %in% dnaDirs),]
write.table(tmpBed2, "/mnt/DATA6/mouseData/20211203mouseBeds_tc.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)
