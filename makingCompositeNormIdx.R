library(jsonlite)


###too lazy to write code to loop through all the directors, but the ones I will use to make a composte are listed below
#/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-94-SHlotanPgu-4_169_261/plugin_out/coverageAnalysis_out.434
#/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-117-KI_DNA_panGU_193_323/plugin_out/coverageAnalysis_out.517
#/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-86-UT_DNA_PanGU_161_241/plugin_out/coverageAnalysis_out.401


#directory <- "/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-94-SHlotanPgu-4_169_261/plugin_out/coverageAnalysis_out.434"
#directory <- "/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-117-KI_DNA_panGU_193_323/plugin_out/coverageAnalysis_out.517"
directory <- "/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-86-UT_DNA_PanGU_161_241/plugin_out/coverageAnalysis_out.401"


setwd(directory)

barcodes <- system('find . -type d -name "IonXpress*"', intern = TRUE)
for(i in seq_along(barcodes)){
  barcodes[i] <- gsub("./","",barcodes[i])
}


jsonFile <- fromJSON("./results.json")



###making one large index file 
###if redoing make sure to wipe the table once
finalTable <- NULL
finalTable <- c("SampleName","Barocde","path")

barcodes <- sort(barcodes)

for(i in seq_along(barcodes)){
  if(length(jsonFile$barcodes[[i]]$`Uniformity of amplicon coverage`) == 0){
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
}

rownames(finalTable) <- NULL
colnames(finalTable) <- NULL
finalTable <- finalTable[-1,]
finalTable <- data.frame(finalTable, stringsAsFactors = FALSE)

write.table(x = finalTable, file = "/home/kevhu/data/normals/PanGU/ampliconCovIdx.PanGU_normals.combined.n77.txt",
            sep = '\t',row.names = FALSE, quote = FALSE, col.names = FALSE)

ID.list <- finalTable[,1]

write.table(x = ID.list, file = "/home/kevhu/data/normals/PanGU/ampliconCovIdx.PanGU_normals.combined.n77.IDlist.txt",
            sep = '\t',row.names = FALSE, quote = FALSE, col.names = FALSE)
