###Finding duplicates sample ID's on BEDs

library(hashmap)

tableOfData <- read.table("/home/kevhu/data/allDataTest.txt",
                     header = TRUE,sep = '\t', quote = "", stringsAsFactors = FALSE)
a <- tableOfData
tableOfData <- a
unique(tableOfData$BED)
###remapping the names
namesRemap = read.table("/home/kevhu/data/Bed_name_map.July2017.csv", stringsAsFactors = FALSE,sep = ",",header = TRUE)
namesRemap = as.data.frame(namesRemap)
namesRemap$ORIG_BED[which(is.na(namesRemap$ORIG_BED))] = "NA"
namesRemap$FINAL_BED[which(is.na(namesRemap$FINAL_BED))] = "NA"
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
tableOfData <- tableOfData[-which(tableOfData$BED == "none"),]
tableOfData <- tableOfData[-which(tableOfData$BED == "NA"),]
tableOfData <- tableOfData[-which(is.na(tableOfData$BED)),]
unique(tableOfData$BED) 

subSample <- tableOfData[,c("SAMPLE","BARCODE","BED","REPORTNUM")]
#testNames <- do.call(paste0, subSample[,c("SAMPLE","BARCODE","BED")])
#testNames2 <- do.call(paste0, subSample[,c("SAMPLE","BARCODE","BED","REPORTNUM")])
testNames <- apply(subSample[,c("SAMPLE","BARCODE","BED")],1,paste, collapse=",")
testNames2 <- apply(subSample[,c("SAMPLE","BARCODE","BED","REPORTNUM")],1,paste, collapse=",")

length(unique(testNames))
length(unique(testNames2))


testNames.uniq <- unique(testNames)
testNames2.uniq <- unique(testNames2)
tableOfMatches <- NULL
for(i in seq_along(testNames.uniq)){
  tableOfMatches[i] <- list(grep(testNames.uniq[i], testNames2.uniq))
}

multipleMatches <- tableOfMatches[which(lengths(tableOfMatches) > 1)]

tableOfMatches2 <- NULL
for(i in seq_along(multipleMatches)){
  a <- NULL
  for(j in seq_along(multipleMatches[[i]])){
    a <- c(a, testNames2.uniq[multipleMatches[[i]][j]])
  }
  tableOfMatches2[i] <- list(a)
}


save(tableOfMatches2, file = "/home/kevhu/data/listOfDuplicatedSampleRuns.Robj")

##If I were to just subset the rerun guys - I'll just pick the first one .. how many would I be left with
##Let me test thist 

#neat trick below in order to just get first element of each index of list
pickedFirst <- lapply(tableOfMatches2,'[[',1)

pickedFirst <- as.data.frame(pickedFirst, stringsAsFactors = FALSE)
pickedFirst <- t(pickedFirst)
row.names(pickedFirst) <- NULL

tableOfData$tmp <- apply(subSample[,c("SAMPLE","BARCODE","BED","REPORTNUM")],1,paste, collapse=",")

length(which(tableOfData$tmp %in% pickedFirst))
length(which(tableOfData$tmp %in% pickedFirst))/ nrow(tableOfData)

##remember this is number of those matched from the larger uniqe data set so 1767/4668 is non-uniqe
collapsedTableOfMatches <- unlist(which(lengths(tableOfMatches) == 1))
length(unique(collapsedTableOfMatches))
nonDups <- testNames2.uniq[(collapsedTableOfMatches)]
nonDups <- c(nonDups, pickedFirst)
length(which(tableOfData$tmp %in% nonDups))
length(which(tableOfData$tmp %in% nonDups))/ nrow(tableOfData)

###so sample estimates above say we lose ~15-20% of variant calls
###slight problem is i dont get the same number of uniq names with just sample bed and barcode ~26 off


