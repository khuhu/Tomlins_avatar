
library(jsonlite)


### ccp  "/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-79-20170922_PoolA_CCP_AD_HN_ML_151_227/plugin_out/coverageAnalysis_out.378/"

listOfDirectories <- c("/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-173-Cho_mouse_20171130_362_439/plugin_out/coverageAnalysis_out.767/",
                       "/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-174-Cho_mouse_20171213_363_441/plugin_out/coverageAnalysis_out.769/")

finalTable <- NULL
finalTable <- c("SampleName","Barocde","path")
for(j in seq_along(listOfDirectories)){
  directory <- listOfDirectories[j]
  setwd(directory)
  barcodes <- system('find . -type d -name "IonXpress*"', intern = TRUE)
  for(i in seq_along(barcodes)){
    barcodes[i] <- gsub("./","",barcodes[i])
  }
  jsonFile <- fromJSON("./results.json")
  
  
  
  ###making one large index file 
  barcodes <- sort(barcodes)
  for(i in seq_along(barcodes)){
    if(length(jsonFile$barcodes[[i]]$`Uniformity of amplicon coverage`) == 0){
      print(i)
      next()
    }
    sampleName <- jsonFile$barcodes[[i]]$`Sample Name`
    print(sampleName)
    path <- paste(directory, "/",barcodes[i], "/",sep = "")
    print(path)
    setwd(path)
    ampliCovFile <- system('find . -name "*.amplicon.cov.xls" | grep IonXpress',intern = TRUE)
    ampliCovFile <- sub("./","",ampliCovFile)
    ampliCovFile <- paste(path, ampliCovFile, sep = "")
    combined <- c(sampleName, barcodes[i], ampliCovFile)
    finalTable <- rbind(finalTable, combined)
    setwd("~/")
  }
}

rownames(finalTable) <- NULL
colnames(finalTable) <- NULL
finalTable <- finalTable[-1,]
finalTable <- data.frame(finalTable, stringsAsFactors = FALSE)




write.table(x = finalTable, file = "/home/kevhu/data/20180213newMouseIdx.txt",
            sep = '\t',row.names = FALSE, quote = FALSE, col.names = FALSE)


