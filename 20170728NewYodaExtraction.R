setwd(dir = "/home/kevhu/")

###Setting up connection to a database

library(jsonlite)
require("RPostgreSQL")

pw2 <- {
  "SAT_1840"
}

drv <- dbDriver("PostgreSQL")

con <- dbConnect(drv, dbname="umst_old",
                 host="localhost", port= 5432,
                 user = "albert", password = pw2)



dbExistsTable(con, c("yoda", "rundb_pluginresult"))
dbExistsTable(con, c("yoda", "rundb_experimentanalysissettings"))

myTable <- dbReadTable(con, c("yoda","rundb_pluginresult"))
myTable.subset <- myTable[,c("result_id","state","store")]
myTable.subset <- myTable.subset[which(myTable.subset$state == "Completed"),]
myTable.subset <- myTable.subset[-which(myTable.subset$store == ""),]

myTable2 <- dbReadTable(con, c("yoda","rundb_experimentanalysissettings"))
myTable2.subset <- myTable2[,c("id","targetRegionBedFile","barcodedSamples")]

tableOfFun <- NULL
tableOfFun <- c("sample", "barcode","bed", "Reportnum","Uniformity","numberMappedReads","percentMapped","avgCov")
for(i in seq_along(myTable.subset$store)){
  #print(i)
  g <- fromJSON(myTable.subset$store[i])
  for(j in seq_along(g$barcodes)){
    sampleName <- g$barcodes[[j]]$`Sample Name`
    bed <- g$barcodes[[j]]$`Targeted Regions`
    uniformity <- g$barcodes[[j]]$`Uniformity of base coverage`
    nameOfbar <- names(g$barcodes[j])
    percentBases <- g$barcodes[[j]]$`Percent base reads on target`
    numberReads <- g$barcodes[[j]]$`Number of mapped reads`
    avgCov <- g$barcodes[[j]]$`Average base coverage depth`
    reportNum <- myTable.subset$result_id[i]
    if(length(uniformity) == 0){
      next()
    }
    if(length(sampleName) == 0){
      #print(i)
      test <- tryCatch(fromJSON(myTable2.subset$barcodedSamples[grep(reportNum, myTable2.subset$id)]), error = function(x) return(NULL))
      #print(test)
      if(!is.null(test)){
        #print(i)
        sampleName <- names(test[grep(nameOfbar, test)])
        #print(sampleName)
      }
      if(length(sampleName) == 0){
        sampleName <- "None"
      }
    }
    if(length(bed) == 0){
      bed <- g$barcodes[[j]]$`Target Regions`
      if(length(bed) == 0){
        test2 <- myTable2.subset$targetRegionBedFile[grep(reportNum, myTable2.subset$id)]
        print(test2)
        print(reportNum)
        if(!is.null(test2)){
          #print(i)
          ssTest <- unlist(strsplit(test2, split = "/"))
          bed <- ssTest[length(ssTest)]
          print(bed)
        }
        if(length(bed) == 0){
          bed <- "None"
        }
      }
      if(length(numberReads) == 0){
        numberReads <- "None"
      }
      if(length(percentBases) == 0){
        percentBases <- "None"
      }
      if(length(avgCov) == 0){
        avgCov <- "None"
      }
    }
    combined <- c(sampleName, nameOfbar, bed, reportNum, uniformity, numberReads, percentBases, avgCov)
    #print(combined)
    tableOfFun <- rbind(tableOfFun, combined)
  }
}


colnames(tableOfFun) <- c("Sample","Barcode","Bed", "Reportnum","Uniformity","numberMappedReads","percentMapped","avgCov")
rownames(tableOfFun) <- NULL
tableOfFun <- tableOfFun[-1,]
tableOfFun <- as.data.frame(tableOfFun, stringsAsFactors = FALSE)

tableOfFun.filtered <- tableOfFun
tableOfFun.filtered$Sample[is.na(tableOfFun$Sample)] <- "None"

length(which(tableOfFun.filtered$Sample == "None"))
184/3276
length(which(tableOfFun.filtered$Bed == "None"))
69/3276


write.table(x = tableOfFun.filtered, file = "/home/kevhu/data/20170728NewYoda.tsv", sep = "\t", row.names = FALSE)
