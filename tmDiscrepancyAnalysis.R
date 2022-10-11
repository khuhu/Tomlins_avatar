tableOfData <- read.table("/home/kevhu/data/allDataTest.txt",
                          header = TRUE,sep = '\t', quote = "", stringsAsFactors = FALSE)

tableOfData <- read.table("/home/kevhu/data/20170810mySQLdata.txt",
                          header = TRUE,sep = '\t', quote = "", stringsAsFactors = FALSE)


a <- tableOfData
tableOfData <- a
unique(tableOfData$BED)
###remapping the names
namesRemap = read.table("/home/kevhu/data/Bed_name_map.July2017.csv", stringsAsFactors = FALSE,sep = ",",header = TRUE)
namesRemap = as.data.frame(namesRemap)
namesRemap$ORIG_BED[which(is.na(namesRemap$ORIG_BED))] = "None"
namesRemap$FINAL_BED[which(is.na(namesRemap$FINAL_BED))] = "None"
namesRemap2 <- namesRemap[19:20,]
namesRemap <- namesRemap[-c(19,20),]
extra <- c("OCP3_20140506_designed_noTrack", "OCP3_20140506_designed")
namesRemap <- rbind(namesRemap, extra)
#col1Remap <- namesRemap[,1]
#col1Remap <- sapply(col1Remap, function(x) paste("\\<",x,"\\>"))
#namesRemap$ORIG_BED <- col1Remap

for(i in seq_along(namesRemap$ORIG_BED)){
  tableOfData$BED <- gsub(pattern = namesRemap$ORIG_BED[i], replacement = namesRemap$FINAL_BED[i], x = tableOfData$BED, fixed = TRUE)
}

OCP_2015 <- which(nchar(tableOfData$BED) == 8)
OCP_201506 <- which(nchar(tableOfData$BED) == 10)

tableOfData$BED[OCP_2015] <- "OCP_20150630_designed"
tableOfData$BED[OCP_201506] <- "OCP_20150630_designed"
unique(tableOfData$BED)
tableOfData$BED[which(tableOfData$BED == "none")] <- "None"
#tableOfData <- tableOfData[-which(tableOfData$BED == "NA"),]
tableOfData$BED[which(is.na(tableOfData$BED))] <- "None"
unique(tableOfData$BED) 

subSample <- tableOfData[,c("SAMPLE","BARCODE","BED","REPORTNUM","REPORT")]

###trying to what is and isnt there in terms of my none names

subSample$SAMPLE <- gsub("[[:punct:]]", "", subSample$SAMPLE)
subSample$SAMPLE <- sapply(subSample$SAMPLE, toupper)
subSample$SAMPLE <- gsub("DNA","", subSample$SAMPLE)
subSample <- subSample[-which(subSample$BED == "prTissue_WG00196_02092016_Designed"),]
#consoliatedDf$Sample <- gsub("[[:punct:]]", "", consoliatedDf$Sample)

#testNames2 <- unique(do.call(paste0,subSample[,c("BARCODE","BED","REPORTNUM")]))
#testNames1 <- unique(do.call(paste0,subSample[,c("SAMPLE","BARCODE","BED","REPORTNUM")]))

numbers_only <- function(x) !grepl("\\D", x)

###this for loop isn't optimal I could probably just use apply and it would be exponentially faster
for(i in which(numbers_only(subSample$BARCODE))){
  if(nchar(subSample$BARCODE[i]) == 1){
    subSample$BARCODE[i] <- paste0("IonXpress_00",subSample$BARCODE[i])
  }
  if(nchar(subSample$BARCODE[i]) == 2){
    subSample$BARCODE[i] <- paste0("IonXpress_0",subSample$BARCODE[i])
  }
}


write.table(subSample, file = "/home/kevhu/data/subSampleData.tsv", row.names = FALSE, sep = '\t')

testNames2 <- unique(apply(X = subSample[,c("BARCODE","BED","REPORTNUM")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))
testNames1 <- unique(apply(X = subSample[,c("SAMPLE","BARCODE","BED","REPORTNUM")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))
testNames3 <- unique(apply(X = subSample[,c("SAMPLE","BED")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))
testNames4 <- unique(apply(X = subSample[,c("SAMPLE","BARCODE","BED")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))
testNames5 <- unique(apply(X = subSample[,c("SAMPLE","BARCODE","BED","REPORT")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))
testNames6 <- unique(apply(X = subSample[,c("SAMPLE","BARCODE","REPORT")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))


testNames2 <- gsub(" ","", testNames2)
testNames1 <- gsub(" ","", testNames1)
testNames3 <- gsub(" ","", testNames3)
testNames4 <- gsub(" ","", testNames4)
testNames5 <- gsub(" ","", testNames5)
testNames6 <- gsub(" ","", testNames6)

subSample.cc <- subSample
subSample.cc$tmp <- NULL
subSample.cc$tmp2 <- NULL

subSample.cc$tmp <- apply(X = subSample[,c("SAMPLE","BARCODE","BED","REPORT")],MARGIN = 1, FUN = function(x) paste(x, collapse = ","))
subSample.cc <- subSample.cc[!duplicated(subSample.cc$tmp),]
subSample.cc$tmp2 <- apply(X = subSample.cc[,c("SAMPLE","BARCODE","REPORT")],MARGIN = 1, FUN = function(x) paste(x, collapse = ","))



#listMatch <- (apply(X = consoliatedDf[,c("Sample","Barcode","Bed","Reportnum")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))
#listMatch2 <- (apply(X = consoliatedDf[,c("Barcode","Bed","Reportnum")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))










#####below is old stuff I don't run anymore
#####
#####



length(which(listMatch2 %in% testNames2))

2901/8222  #testNames2

length(which(testNames1 %in% listMatch))

1994/4592


length(which(listMatch %in% testNames1))
2697/8222

#This shows that there is still a large proportion fo the data not accounted for .....

notMatched <- testNames1[which(!testNames1 %in% listMatch)]

notMatchedf <- NULL
for(i in seq_along(notMatched)){
  a <- unlist(strsplit(notMatched[i], split = ","))
  notMatchedf <- rbind(notMatchedf,a)
}

notMatchedf <- as_data_frame(notMatchedf)
unique(notMatchedf$V2)
notMatchedf$tmp <- apply(notMatchedf[,2:4], MARGIN = 1, FUN = function(x) paste(x,collapse = ","))

consoliatedDf$tmp <- apply(consoliatedDf[,c(2:4)], MARGIN = 1, FUN = function(x) paste(x,collapse = ","))
subSample$tmp <- apply(subSample[,c(2:4)], MARGIN = 1, FUN = function(x) paste(x,collapse = ","))

#~223 with different names

length(which(consoliatedDf$tmp %in% notMatchedf$tmp))
length(which(subSample$tmp %in% notMatchedf$tmp))

duumy <- subSample[which(subSample$tmp %in% notMatchedf$tmp),1:4]
length(which(duumy$SAMPLE == "HorizonDxUM"))
duumy <- duumy[-which(duumy$SAMPLE == "HorizonDxUM"),]

dummy2 <- consoliatedDf[which(consoliatedDf$tmp %in% notMatchedf$tmp),1:4]

#first step is to fix the barcodes in the tableOfData
###thing to note is i need to add IonXpress_00 or IonXpress_0 depending on length of number

copyOfTestTable <- tableOfData
copyOfTestTable <- copyOfTestTable[-grep("IonXpress", copyOfTestTable$BARCODE),]
copyOfTestTable <- copyOfTestTable[-grep("HaloPlex", copyOfTestTable$BARCODE),]
copyOfTestTable <- copyOfTestTable[-grep("TargetSeq", copyOfTestTable$BARCODE),]

copyOfTestTable$BARCODE <- sapply(copyOfTestTable$BARCODE, function(x){
  if(nchar(x) == 2){
    paste0("IonXpress_0",x)
  } else{
    paste0("IonXpress_00",x)
  }
})

rownames(copyOfTestTable)

rownames(tableOfData)

tableOfData <- tableOfData[-which(rownames(tableOfData) %in% rownames(copyOfTestTable)),]
tableOfData <- rbind(tableOfData, copyOfTestTable)
