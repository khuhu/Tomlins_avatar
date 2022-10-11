#fullTable <- read.table("/mnt/DATA4/kevhu/fullRNACountTable.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
library(xlsx)
fullTable <- read.table("/mnt/DATA4/kevhu/fullRNACountTable2.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

andiList <- read.csv("/mnt/DATA4/kevhu/20180105AndisList.csv",stringsAsFactors = FALSE)


###all this is temporary, just for Andi's samples
#bed.4  <- read.table("/mnt/DATA4/kevhu/WG00196.4_05122017_Designed.noTrack.bed" ,sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#bed.4 <- bed.4[order(bed.4$V4),]
#bed.5 <- read.table("/mnt/DATA4/kevhu/WG00196.5_05122017_Designed.noTrack.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#bed.5 <- bed.5[order(bed.5$V4),]
#listOfAndisBed <- c("bed.5","bed.4")


#test <- NULL
#for(i in 1:nrow(fullTable)){
#  test <- c(test,strsplit(fullTable$ReportID[i],"/", fixed = TRUE)[[1]][5])
#}
#fullTable$ReportID <- test
#write.table(fullTable, "/mnt/DATA4/kevhu/fullRNACountTable2.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


#S <- c("UT80","VCaP_1B_PanGu")
#B <- c("IonXpress_001","IonXpress_096")
#Bed <- c("WG_IAD127899.20170720.designed","WG_IAD127899.20170720.designed")
#Report <- c("","")

key <- paste0(fullTable$SampleName,fullTable$Barcode, fullTable$ReportID)

keyVals <- rep(0, length(key))
hash <- hashmap(keys = key, values = keyVals)
lookup <- paste0(andiList$SampleName, andiList$Barcode, andiList$ReportID)

for(i in seq_along(lookup)){
  if(hash$has_key(lookup[i]) == TRUE){
    counter <- hash$find(lookup[i])
    hash$insert(lookup[i],c(counter + 1))
  }
}


listOfNames<- names(which(c(hash$data()) > 0))
exportTable <- fullTable[which(key %in% listOfNames),]


listOfBedsToSplit <- c(unique(exportTable$Bed))
listOfSplitTables <- NULL


for(i in seq_along(listOfBedsToSplit)){
  a <- NULL
  a <- exportTable[which(exportTable$Bed == listOfBedsToSplit[i]),]
  assign(paste0("RnaCountData", listOfBedsToSplit[i]), a)
  listOfSplitTables <- c(listOfSplitTables, paste0("RnaCountData", listOfBedsToSplit[i]))
}

###time to transform the pulled data into the correct final form

listOfFinalExcelSheets <- NULL

for(i in seq_along(listOfSplitTables)){
  sampCount <- unique(eval(as.name(listOfSplitTables[i]))$SampleName)
  tmpTable <- NULL
  mappedReadsList <- NULL
  for(j in seq_along(sampCount)){
    tmpTable2 <- NULL
    tmpTable2 <- eval(as.name(listOfSplitTables[i]))[which(eval(as.name(listOfSplitTables[i]))$SampleName == sampCount[j]),]
    tmpTable2 <- tmpTable2[order(tmpTable2$AmpliconID),]
    mappedReadsList <- c(mappedReadsList, unique(tmpTable2$NumberOfMappedReads))
    if(j == 1){
      tmpTable <- cbind(tmpTable, tmpTable2[,c("AmpliconID")])
      tmpTable <- cbind(tmpTable, tmpTable2[,c("Gene")])
      tmpTable <- data.frame(tmpTable, stringsAsFactors = FALSE)
      tmpTable <- cbind(tmpTable, tmpTable2[,ncol(tmpTable2)])
      colnames(tmpTable)[1] <- c("AmpliconID")
      colnames(tmpTable)[2] <- c("Gene")
      colnames(tmpTable)[ncol(tmpTable)] <- sampCount[j]
    }
    else{
      tmpTable <- cbind(tmpTable, tmpTable2[,ncol(tmpTable2)])
      colnames(tmpTable)[ncol(tmpTable)] <- sampCount[j]
    }
  }
  ###one - liner for temporary column binding for Andi
  #print(i)
  #colDummy <- c(eval(as.name(listOfAndisBed[i]))$V1)
  #tmpTable <- cbind(colDummy, tmpTable)
  mappedReadsList <- c(NA,NA,NA, mappedReadsList)
  tmpTable <- rbind(mappedReadsList, tmpTable)
  #tmpTable[,4:ncol(tmpTable)] <- as.numeric(tmpTable[,4:ncol(tmpTable)])
  assign(paste0("transformed",listOfBedsToSplit[i]), tmpTable)
  listOfFinalExcelSheets <- c(listOfFinalExcelSheets, paste0("transformed",listOfBedsToSplit[i]))
}




for(i in seq_along(listOfFinalExcelSheets)){
  write.xlsx(x = eval(as.name(listOfFinalExcelSheets[i])), file = "/mnt/DATA4/kevhu/testAndi2.xlsx",sheetName = listOfBedsToSplit[i], row.names = FALSE,append = TRUE)
}




