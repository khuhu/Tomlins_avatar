### gc bed

pgu2_bed <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-114-PGU2_NL20001_324_287/plugin_out/coverageAnalysis_out.529/local_beds/IAD203665_173_Designed.gc.bed",
                       sep = "\t", stringsAsFactors = FALSE, skip = 1)

otherBed <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed", sep = "\t",
                       stringsAsFactors = FALSE)


v6 <- pgu2_bed$V3 - pgu2_bed$V2
v5 <- round(v6 * (pgu2_bed$V6/100))
v7 <- signif(pgu2_bed$V6/100, digits = 2)
v8 <- str_remove(pgu2_bed$V5, "SUBMITTED_REGION\\=")
v8[grep("REGION", v8)] <- str_remove(v8[grep("REGION", v8)], ";Pool*.")
v8[-grep("REGION", v8)] <- str_remove(v8[-grep("REGION", v8)], "\\_.*")

finalBed <- data.frame("V1" = pgu2_bed$V1, "V2" = pgu2_bed$V2, "V3" = pgu2_bed$V3,
                       "V4" = pgu2_bed$V4, "V5" = v5, "V6" = v6,
                       "V7" = v7, "V8" = v8, stringsAsFactors = FALSE)

write.table(finalBed, "/home/kevhu/data/bedFiles/IAD203665_173_Designed.gc.bed",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

normalFiles <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-114-PGU2_NL20001_324_287/plugin_out/coverageAnalysis_out.529/Auto_PGU2_NL20001_eros_287.bc_summary.xls",
                          sep = "\t", stringsAsFactors = FALSE, header = TRUE)

normalFiles_samps <- normalFiles$Sample.Name[grep("NL", normalFiles$Sample.Name)]
normalFiles_barcodes <- normalFiles$Barcode.ID[grep("NL", normalFiles$Sample.Name)]

writeLines(normalFiles_samps, "/mnt/DATA6/mouseData/normals/20220914pgu2.IDlist.txt")


normalDir <- "/mnt/DATA3/eros_tmp/Auto_user_AUS5-114-PGU2_NL20001_324_287/plugin_out/coverageAnalysis_out.529/"
setwd(normalDir)

normalFilesLoc <- system("find '/mnt/DATA3/eros_tmp/Auto_user_AUS5-114-PGU2_NL20001_324_287/plugin_out/coverageAnalysis_out.529/' -type f -name '*.amplicon.cov.xls'",
       intern = TRUE)

normalFilesLoc2 <- normalFilesLoc[grep(paste0(normalFiles_barcodes, collapse = "|"), normalFilesLoc)]

normalFilesFinal <- data.frame("V1" = normalFiles_samps, "V2" = normalFiles_barcodes, "V3" = normalFilesLoc2,
                               stringsAsFactors = FALSE)

write.table(normalFilesFinal, "/mnt/DATA6/mouseData/normals/20220914pgu2.txt", sep = "\t",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

