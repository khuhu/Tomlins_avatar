
OCPv2.poolIdx <- read.table("/home/kevhu/data/OCP2.20131028.384WellPlateDataSheet.csv", skip = 1, header = TRUE, sep = ",", stringsAsFactors = FALSE)
OCPv2.bed <- read.table("/home/hovelson/lists/bed/OCP2.20131028.designed.noTrack.GC.bed")
masterCovList <- read.table("/mnt/DATA/project_data/lifetech_seq_run/manifest_docs/cna.samp_bc_ampFile.index.20140918.txt", sep = "\t")
andisList <- read.table("/home/kevhu/data/andiListPoolSamps.csv",sep = ",", header = FALSE, skip = 1)


OCPv2.poolIdx$X384Well_PlateID <- gsub(pattern = "_[A-Z]", replacement = "", OCPv2.poolIdx$X384Well_PlateID)
andisList <- andisList[,1:2]
masterCovList$BarcodeNums <- gsub(pattern = "[A-z]", replacement = "",masterCovList$V2)

for(i in seq_along(andisList$V1)){
  if(nchar(andisList$V2[i]) == 1){
    andisList$V2[i] <- paste("00",andisList$V2[i], sep = "")
  }
  if(nchar(andisList$V2[i]) == 2){
    andisList$V2[i] <- paste("0",andisList$V2[i], sep = "")
  }
}

idx <- NULL
for(i in seq_along(andisList$V1)){
  correctMat <- which(grepl(andisList$V1[i], masterCovList$V1) & grepl(andisList$V2[i], masterCovList$V2))
  print(correctMat)
  if(length(correctMat) != 1){
    print("ERROR")
  }
  idx <- c(idx, correctMat)
}

masterCovList.subset <- masterCovList[idx,]

###alternative way 

choice <- c("PR74","PR66","PR69","PR64","LU12","LU23","LU30","LU35")

secondSubset <- NULL
for(i in seq_along(choice)){
  matches <- which(grepl(choice[i],masterCovList$V1))
  secondSubset <- c(secondSubset, matches)
}

masterCovList.subset2 <- masterCovList[secondSubset,]


write.table(masterCovList.subset2, "/home/kevhu/data/OCPv2AnalysisCovListSub.tsv", row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")

###something to note is that the list of primers is larger than the number of amplicons on the panel


counter <- 0
for(i in seq_along(OCPv2.bed$V1)){
  matchIdx <- NULL
  matchIdx <- which(OCPv2.poolIdx$Amplicon_ID == OCPv2.bed$V4[i])
  OCPv2.bed$Pool[i] <- OCPv2.poolIdx$X384Well_PlateID[matchIdx]
  if(isNumeric(matchIdx)){
    counter <- counter + 1
  }
}


write.table(OCPv2.bed, "/home/kevhu/data/OCP2.20131028.designed.noTrack.Pools.GC.bed",sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)

OCPv2.bed.2 <- OCPv2.bed
OCPv2.bed.2$V8 <- paste(OCPv2.bed.2$V8,OCPv2.bed.2$Pool, sep = "")
OCPv2.bed.2 <- OCPv2.bed.2[,1:8]

write.table(OCPv2.bed.2, "/home/kevhu/data/OCP2.20131028.designed.noTrack.Pools.GC2.bed",sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
