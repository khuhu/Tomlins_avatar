library(readxl)
library(stringr)

ocaplus_bed <- read.table("/mnt/DATA6/mouseData/bedFiles/20210726_OCAP_gc_noheader.bed", sep = "\t",
                            stringsAsFactors = FALSE, header = FALSE)

pg1 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/OCApWBC_lob_all.amplicon.cov.xlsx", sheet = 1)
pg2 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/OCApWBC_lob_all.amplicon.cov.xlsx", sheet = 2)
pg3 <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/OCApWBC_lob_all.amplicon.cov.xlsx", sheet = 3)


fp1 <- read.table("/mnt/DATA5/tmp/kev/20220211Andi/fullPanel/combinedCalls_Hayes_132_346.txt", sep = "\t",
                  header = TRUE, stringsAsFactors = FALSE)

fp2 <- read.table("/mnt/DATA5/tmp/kev/20220211Andi/fullPanel/combinedCalls_Hayes_good_3_270.txt", sep = "\t",
                  header = TRUE, stringsAsFactors = FALSE)

fp1_matrix <- read.table("/mnt/DATA5/tmp/kev/20220211Andi/fullPanel/cnMatrix_gene_Hayes_132_346.txt", sep = "\t",
                  header = TRUE, stringsAsFactors = FALSE)

fp2_matrix <- read.table("/mnt/DATA5/tmp/kev/20220211Andi/fullPanel/cnGeneMatrix_Hayes_good_3_270.txt", sep = "\t",
                  header = TRUE, stringsAsFactors = FALSE)

fp1_edited <- fp1[grep("ctD", fp1$Sample),]
fp2_edited <- fp2[grep("SRSQ|HRZ", fp2$Sample),]

fp1_matrix_edited <- fp1_matrix[,c(1, grep("ctD", colnames(fp1_matrix)))]
fp2_matrix_edited <- fp2_matrix[,c(1, grep("SRSQ|HRZ", colnames(fp2_matrix)))]


write.table(fp1_edited, "/mnt/DATA5/tmp/kev/20220211Andi/fullPanel/combinedCalls_Hayes_132_346_edited.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(fp2_edited, "/mnt/DATA5/tmp/kev/20220211Andi/fullPanel/combinedCalls_Hayes_good_3_270_edited.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(fp1_matrix_edited, "/mnt/DATA5/tmp/kev/20220211Andi/fullPanel/cnMatrix_Hayes_132_346_edited.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(fp2_matrix_edited, "/mnt/DATA5/tmp/kev/20220211Andi/fullPanel/cnMatrix_Hayes_good_3_270_edited.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

ocaplus_bed2 <- ocaplus_bed[which(ocaplus_bed$V4 %in% pg2$region_id),]
ocaplus_bed3 <- ocaplus_bed[-which(ocaplus_bed$V4 %in% pg3$region_id),]


write.table(ocaplus_bed2, "/mnt/DATA6/mouseData/bedFiles/20210726_OCAP_gc_Andi2_noheader.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ocaplus_bed3, "/mnt/DATA6/mouseData/bedFiles/20210726_OCAP_gc_Andi3_noheader.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

mouseBedIdx <- read.table("/mnt/DATA6/mouseData/mouseBedsIdx.txt", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)

reportList <- c("Auto_user_AUS5-89-HayesLobOCAplusAC01_294_230", "Auto_user_AUS5-102-HayesLobOCAplusAC01_2_310_261", 
                "HayesLobOCAplusAC01_3_good_270", "Auto_user_AUS5-132-HayesLobRun1_346_329", "Auto_user_AUS5-133-HayesLobRun2_347_331",
                "Auto_user_AUS5-136-HayesLobRun3_348_339", "Auto_user_AUS5-137-HayesLobRun4_349_341", "Auto_user_AUS5-159-Hayes_LobCTC_Run5_373_388")

mouseBedIdx$cnBed[grep(paste(reportList, collapse = "|"), mouseBedIdx$directories)] <- "20210726_OCAP_gc_Andi2_noheader.bed"
write.table(mouseBedIdx, "/mnt/DATA6/mouseData/mouseBedsIdx.txt", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)


### copying over all the pdfs/pngs + getting all the cnMatrices + combinedCalls and combining them
### two times one for each type of bed file

reportList <- c("Auto_user_AUS5-89-HayesLobOCAplusAC01_294_230", "Auto_user_AUS5-102-HayesLobOCAplusAC01_2_310_261", 
                "HayesLobOCAplusAC01_3_good_270", "Auto_user_AUS5-132-HayesLobRun1_346_329", "Auto_user_AUS5-133-HayesLobRun2_347_331",
                "Auto_user_AUS5-136-HayesLobRun3_348_339", "Auto_user_AUS5-137-HayesLobRun4_349_341", "Auto_user_AUS5-159-Hayes_LobCTC_Run5_373_388")


normalIds <- as.vector(read.table("/mnt/DATA6/mouseData/normals/20220211ocap_WBC.IDlist.txt", 
                        stringsAsFactors = FALSE, header = FALSE)$V1)

dir <- "/mnt/DATA6/mouseData/copynumber/"
allCombinedCalls <- NULL
allCnMatrix <- NULL
for (i in reportList) {
  tmpDir <- paste0(dir, i)
  setwd(tmpDir)
  tmpRes <- system(paste0("find ", tmpDir ," -type f -name '*combinedCalls.txt*'"),
                   intern = TRUE)
  tmpCombinedCall <- read.table(tmpRes, sep = "\t", stringsAsFactors = FALSE,
                                header = TRUE)
  tmpRes2 <- system(paste0("find ", tmpDir ," -type f -name '*cnMatrix_gene.txt*'"),
                   intern = TRUE)
  tmpGeneMatrix <- read.table(tmpRes2, sep = "\t", stringsAsFactors = FALSE,
                              header = TRUE)
  if (i == reportList[1]) {
    allCnMatrix <- tmpGeneMatrix
    allCombinedCalls <- tmpCombinedCall
    
  } else{
    tmpGeneMatrix <- tmpGeneMatrix[,-which(colnames(tmpGeneMatrix) %in% normalIds)]
    allCnMatrix <- cbind(allCnMatrix, tmpGeneMatrix[,2:ncol(tmpGeneMatrix)])
    
    tmpCombinedCall <- tmpCombinedCall[-which(tmpCombinedCall$Sample %in% normalIds),]
    allCombinedCalls <- rbind(allCombinedCalls, tmpCombinedCall)
  }
  
  system("cp *pdf* /mnt/DATA5/tmp/kev/20220211Andi/panel2/")
}

write.table(allCombinedCalls, "/mnt/DATA5/tmp/kev/20220211Andi/panel2/combinedCalls_andi2.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(allCnMatrix, "/mnt/DATA5/tmp/kev/20220211Andi/panel2/cnMatrix_andi2.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


### second set of bed


dir <- "/mnt/DATA6/mouseData/copynumber/"
allCombinedCalls <- NULL
allCnMatrix <- NULL
for (i in reportList) {
  tmpDir <- paste0(dir, i)
  setwd(tmpDir)
  tmpRes <- system(paste0("find ", tmpDir ," -type f -name '*combinedCalls.txt*'"),
                   intern = TRUE)
  tmpCombinedCall <- read.table(tmpRes, sep = "\t", stringsAsFactors = FALSE,
                                header = TRUE)
  tmpRes2 <- system(paste0("find ", tmpDir ," -type f -name '*cnMatrix_gene.txt*'"),
                    intern = TRUE)
  tmpGeneMatrix <- read.table(tmpRes2, sep = "\t", stringsAsFactors = FALSE,
                              header = TRUE)
  if (i == reportList[1]) {
    allCnMatrix <- tmpGeneMatrix
    allCombinedCalls <- tmpCombinedCall
    
  } else{
    tmpGeneMatrix <- tmpGeneMatrix[,-which(colnames(tmpGeneMatrix) %in% normalIds)]
    allCnMatrix <- cbind(allCnMatrix, tmpGeneMatrix[,2:ncol(tmpGeneMatrix)])
    
    tmpCombinedCall <- tmpCombinedCall[-which(tmpCombinedCall$Sample %in% normalIds),]
    allCombinedCalls <- rbind(allCombinedCalls, tmpCombinedCall)
  }
  
  system("cp *pdf* /mnt/DATA5/tmp/kev/20220211Andi/panel3/")
}

write.table(allCombinedCalls, "/mnt/DATA5/tmp/kev/20220211Andi/panel3/combinedCalls_andi3.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(allCnMatrix, "/mnt/DATA5/tmp/kev/20220211Andi/panel3/cnMatrix_andi3.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

