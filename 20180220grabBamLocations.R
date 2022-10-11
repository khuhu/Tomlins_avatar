listOfDirectories <- c("/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-173-Cho_mouse_20171130_362_439/",
                       "/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-174-Cho_mouse_20171213_363_441/")

finalTable <- NULL
finalTable <- c("SampleName","path")
for(j in seq_along(listOfDirectories)){
  directory <- listOfDirectories[j]
  setwd(directory)
  path <- system('find ./ -type f -name "IonXpress*" -maxdepth 1 | grep ".bam" | grep -v ".bai"', intern = TRUE)
  path <- path[order(path)]
  path <- gsub("./","", path)
  path <- paste0(directory, path)
  sampeFile <- system('find . -type f -name "*.bc_summary.xls"', intern = TRUE)
  tableOfInfo <- read.table(sampeFile, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  tableOfInfo <- tableOfInfo[order(tableOfInfo$Barcode.ID),]
  tableOfInfo$Sample.Name <- gsub(" ", "",  tableOfInfo$Sample.Name, fixed = TRUE)
  tmp <- cbind(tableOfInfo$Sample.Name, path)
  finalTable <- rbind(finalTable, tmp)
  setwd("~/")
}

rownames(finalTable) <- NULL
colnames(finalTable) <- finalTable[1,]
finalTable <- finalTable[-1,]
finalTable <- data.frame(finalTable, stringsAsFactors = FALSE)




write.table(x = finalTable, file = "/home/kevhu/data/20180220mouseBams.txt",
            sep = ',',row.names = FALSE, quote = FALSE, col.names = TRUE)


