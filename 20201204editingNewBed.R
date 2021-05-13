tmpBed <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-76-MG_test1_255_185/plugin_out/coverageAnalysis_out.303/local_beds/IAD202670_167_Designed.gc.bed",
                     sep = "\t", header = FALSE, skip = 1)
tmpBed$V7 <- tmpBed$V3 - tmpBed$V2
tmpBed$V8 <- tmpBed$V6/tmpBed$V7
geneColumn <- str_remove(tmpBed$V5, pattern = "\\;.*")
geneColumn <- str_remove(geneColumn, pattern = "^.*\\=")
geneColumn2 <- geneColumn

geneColumn <- str_remove(geneColumn, pattern = "exon.*")
geneColumn <- str_remove(geneColumn, pattern = "Exon.*")
geneColumn <- str_remove(geneColumn, pattern = "^.*\\,")

tmpNum <- 1:length((which(geneColumn == "SNP_Genotyping")))
newSnpLabels <- paste0(geneColumn[which(geneColumn == "SNP_Genotyping")], "_",tmpNum)
geneColumn[which(geneColumn == "SNP_Genotyping")] <- newSnpLabels

tmpBed$V9 <- geneColumn
newGcBed <- tmpBed[,c(1:4,6:9)]

altBed <- cbind(newGcBed[,1:7], geneColumn2)

write.table(newGcBed, "/home/kevhu/data/bedFiles/IAD202296_167_Designed.gc.bed",
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

write.table(altBed, "/home/kevhu/data/bedFiles/IAD202296_167_Designed.fullName.gc.bed",
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
