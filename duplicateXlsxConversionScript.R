library(plyr)
library(WriteXLS)

load("/home/kevhu/data/listOfDuplicatedSampleRuns.Robj")
for(i in seq_along(tableOfMatches2)){
  for(j in 1:length(tableOfMatches2[[i]])){
    tableOfMatches2[[i]][j] <- gsub(",",":", tableOfMatches2[[i]][j])
  }
}



df <- combine(tableOfMatches2)
finalTable <- NULL
for(i in seq_along(tableOfMatches2)){
  a <- NULL
  b <- unlist(tableOfMatches2[i])
  for(j in seq_along(b)){
   a <- cbind(a, b[j])
   a <- data.frame(a, stringsAsFactors = FALSE)
  }
  finalTable <- rbind.fill(finalTable, a)
}
finalTable <- data.frame(lapply(finalTable, as.character), stringsAsFactors=FALSE)
finalTable[is.na(finalTable)] <- "NA"
  

WriteXLS(finalTable, ExcelFileName = "/home/kevhu/data/tableOfDuplicates.xls",
            row.names = FALSE, col.names = FALSE)





