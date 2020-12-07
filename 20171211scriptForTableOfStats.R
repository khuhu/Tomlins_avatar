library(jsonlite)

### listOfRuns will only be used to access the primary key (result) - found alternative way to do it
#listOfRuns <- system(command = "ls --ignore='*.gz*' /mnt/DATA5/Gandalf_data_keep/", intern = TRUE)

pathToZues <- "/mnt/DATA3/Zeus_data_dump/"
pathToYoda <- "/mnt/DATA3/Yoda_data_dump/"
listofSequencers <- c(pathToYoda, pathToZues)

#setwd(dir = pathToZues)

###used to find results.json to parse everything else
###note this way of creating a connection through the scrip is dangerous b/c it gives the password of the server: future iterations might try to give sudo permissions to the script through the masterscript

listOfJSONs <- NULL
listOfStartPlugins <-NULL

for(i in seq_along(listofSequencers)){
  setwd(listofSequencers[i])
  list1 <- system('sudo -kS find . -name *results.json*  | grep "coverageAnalysis_out"',
                  input="sat1840", intern = TRUE)
  list2 <- system('sudo -kS find . -name *startplugin.json*  | grep "coverageAnalysis_out"',
                  input="sat1840", intern = TRUE)
  
  list1 <- sub("./",listofSequencers[i], list1)
  list2 <- sub("./",listofSequencers[i], list2)
  listOfJSONs <- c(listOfJSONs, list1)
  listOfStartPlugins <- c(listOfStartPlugins, list2)
}

#listOfJSONs <- system('sudo -kS find . -name *results.json*  | grep "coverageAnalysis_out"',
#                      input="sat1840", intern = TRUE)
#listOfStartPlugins <- system('sudo -kS find . -name *startplugin.json*  | grep "coverageAnalysis_out"',
#                             input="sat1840", intern = TRUE)

listOfRNAJsons <- NULL
tableOfFun <- NULL
tableOfFun <- c("sample", "barcode","bed", "report","Uniformity","listJSONkey","numberMappedReads","percentMapped","avgCov")
for(i in seq_along(listOfJSONs)){
  #makeshift tryCatch function for nonaccessible directories - can change later
  readableDir <- tryCatch(fromJSON(listOfJSONs[i]), error = function(x) return(NULL))
  if(is.null(readableDir)){
    next()
  }
  g <- fromJSON(listOfJSONs[i])
  for(j in seq_along(g$barcodes)){
    sampleName <- g$barcodes[[j]]$`Sample Name`
    bed <- g$barcodes[[j]]$`Targeted Regions`
    uniformity <- g$barcodes[[j]]$`Uniformity of base coverage`
    percentBases <- g$barcodes[[j]]$`Percent base reads on target`
    numberReads <- g$barcodes[[j]]$`Number of mapped reads`
    avgCov <- g$barcodes[[j]]$`Average base coverage depth`
    nameOfbar <- names(g$barcodes[j])
    report <- listOfJSONs[i]
    dummy <- unlist(strsplit(listOfJSONs[i],"/"))
    address <- sub("results.json","", listOfJSONs[i])
    address2 <- system(paste("find",address, "-name *bc_summary*"),intern = TRUE)
    address2 <- address2[-which(grepl("scraper", address2))]
    listJSONkey <- i
    if(length(uniformity) == 0){
    listOfRNAJsons <- c(listOfRNAJsons, listOfJSONs[i])
    next()
    }
    if(length(sampleName) == 0 || sampleName == "None"){
      print(i)
      dummyFindSampleVar <- fromJSON(listOfStartPlugins[grep(dummy[2], listOfStartPlugins)][1])
      test <- tryCatch(fromJSON(dummyFindSampleVar$plan$barcodedSample), error = function(x) return(NULL))
      if(!is.null(test)){
        dummyFindSampleVar2 <- fromJSON(dummyFindSampleVar$plan$barcodedSample)
        z <- names(dummyFindSampleVar2[grep(nameOfbar, dummyFindSampleVar2)])
        if(is.list(z)){
          sampleName <- grep(nameOfbar, dummyFindSampleVar2)
        }
        if(length(z) > 0){
          print(z)
          sampleName <- unlist(z)
        }
      }
      if(length(sampleName) == 0 || sampleName == "None"){
        testingtesting <- unlist(strsplit(system(paste("less",address2,"| grep",nameOfbar),
                                                 intern = TRUE),"\t"))
        if(is.character(testingtesting)){
          sampleName <- testingtesting[2]
        }
      }
      if(length(sampleName) == 0 || grepl("IonXpress",sampleName) == TRUE){
        sampleName <- "None"
      }
    }
    if(length(bed) == 0){
      bed <- g$barcodes[[j]]$`Target Regions`
      #if(length(bed) == 0){
      #varJSON <- fromJSON(listOfStartPlugins[grep(dummy[2], listOfStartPlugins)])
      #bed <- varJSON$runinfo$plugin$userInput$sample
      if(length(bed) == 0){
        bed <- "None"
        #}
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
    combined <- c(sampleName, nameOfbar, bed, report, uniformity, listJSONkey, numberReads, percentBases, avgCov)
    #print(combined)
    tableOfFun <- rbind(tableOfFun, combined)
  }
}



tableOfFun <- as.data.frame(tableOfFun, stringsAsFactors = FALSE)
colnames(tableOfFun) <- tableOfFun[1,]
rownames(tableOfFun) <- NULL
tableOfFun <- tableOfFun[-1,]


###this way makes it so it doesn't have to recompute the entire RNA table

#previousRNATab <- read.table("/mnt/DATA4/kevhu/apps/fetchApp/data/fullRNACountTable3.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

previousRNATab <- read.table("/mnt/DATA4/kevhu/apps/fetchApp/data/fullRNACountTable3.txt", sep = "\t", stringsAsFactors = FALSE,
                              header = TRUE, colClasses = c(rep("NULL",10), NA))


listOfRNAJsons <- (unique(listOfRNAJsons))
#listOfPreviousRNAPaths <- previousRNATab$path
listOfPreviousRNAPaths <- previousRNATab
listOfPreviousRNAPaths <- unique(unlist(listOfPreviousRNAPaths))


alreadyProcessed <- which(listOfRNAJsons %in% listOfPreviousRNAPaths)

listOfRNAJsons <- listOfRNAJsons[-c(alreadyProcessed)]

###RNA part
###later add an error file for files not being read


listOfRNAJsons <- (unique(listOfRNAJsons))
tableOfRNA <- NULL
tableOfRNA.nonAmp <- NULL
tableOfRNA <- c("sampleName", "nameOfbar", "bed", "report", "contig_id","region_id","gene","numberReads",
                "percentBases", "tot_e2e")
tableOfRNA.nonAmp <-c("sampleName", "nameOfbar", "bed", "report", "numberReads",
                      "percentBases", "tot_e2e")
  
listofErrors <- NULL

for(i in seq_along(listOfRNAJsons)){
  print(i)
  #setwd(pathToZues)
  #makeshift tryCatch function for nonaccessible directories - can change later
  readableDir <- tryCatch(fromJSON(listOfRNAJsons[i]), error = function(x) return(NULL))
  if(is.null(readableDir)){
    next()
  }
  g <- fromJSON(listOfRNAJsons[i])
  for(j in seq_along(g$barcodes)){
    errorCheck <- g$barcodes[[j]]$`Error`
    if(length(errorCheck) > 0){
      listofErrors <- c(listofErrors, listOfRNAJsons[i]) 
      next()
    }
    sampleName <- g$barcodes[[j]]$`Sample Name`
    bed <- g$barcodes[[j]]$`Targeted Regions`
    percentBases <- g$barcodes[[j]]$`Percent base reads on target`
    numberReads <- g$barcodes[[j]]$`Number of mapped reads`
    nameOfbar <- names(g$barcodes[j])
    report <- listOfRNAJsons[i]
    dummy <- unlist(strsplit(listOfRNAJsons[i],"/"))
    #dummy2 <- paste0(pathToZues,dummy[2],"/primary.key")
    address <- sub("results.json","", listOfJSONs[i])
    address2 <- system(paste("find",address, "-name *bc_summary*"),intern = TRUE)
    address2 <- address2[-which(grepl("scraper", address2))]
    #reportNum <- system(paste("less",dummy2), intern = TRUE)
    listJSONkey <- i
    if(length(sampleName) == 0 || sampleName == "None"){
      #print(i)
      dummyFindSampleVar <- fromJSON(listOfStartPlugins[grep(dummy[2], listOfStartPlugins)][1])
      test <- tryCatch(fromJSON(dummyFindSampleVar$plan$barcodedSample), error = function(x) return(NULL))
      if(!is.null(test)){
        dummyFindSampleVar2 <- fromJSON(dummyFindSampleVar$plan$barcodedSample)
        z <- names(dummyFindSampleVar2[grep(nameOfbar, dummyFindSampleVar2)])
        if(is.list(z)){
          sampleName <- grep(nameOfbar, dummyFindSampleVar2)
        }
        if(length(z) > 0){
          print(z)
          sampleName <- unlist(z)
        }
      }
      if(length(sampleName) == 0 || sampleName == "None"){
        testingtesting <- unlist(strsplit(system(paste("less",address2,"| grep",nameOfbar),
                                                 intern = TRUE),"\t"))
        if(is.character(testingtesting)){
          sampleName <- testingtesting[2]
        }
      }
      if(length(sampleName) == 0 || grepl("IonXpress",sampleName) == TRUE){
        sampleName <- "None"
      }
    }
    if(length(bed) == 0){
      bed <- g$barcodes[[j]]$`Target Regions`
      if(length(bed) == 0){
        bed <- "None"
      }
    }
    if(length(numberReads) == 0){
      numberReads <- "None"
    }
    if(length(percentBases) == 0){
      percentBases <- g$barcodes[[j]]$`Percent reads on target`
      if(length(percentBases) == 0){
        percentBases <- "None"
      }
    }
    
    ###part to obtain the e2e read information
    dirForCovFiles <- sub("results.json","", listOfRNAJsons[i])
    dirForCovFiles <- paste(dirForCovFiles, nameOfbar,"/",sep = "")
    setwd(dirForCovFiles)
    cov.path <- system('find . -name "*amplicon.cov.xls"', intern = TRUE)
    cov.path <- gsub("\\./", dirForCovFiles, cov.path)
    if(length(cov.path) > 1){
      cov.path2 <- cov.path[-c(grep("scraper", cov.path))]
    }
    else{
      cov.path2 <- cov.path
    }
    #dummyTab <- read.table(cov.path, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    dummyTab <- read.table(cov.path2, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    #rowCount <- nrow(dummyTab)
    if(sum(grepl("attributes", colnames(dummyTab))) == 0){
      next()
    }
    
    dummyTab$tot_e2e <- apply(dummyTab[,c("fwd_e2e", "rev_e2e")],1, sum)
    dummyTab$sampleName <- sampleName
    dummyTab$bed <- bed
    dummyTab$report <- report
    dummyTab$numberReads <- numberReads
    dummyTab$percentBases <- percentBases
    dummyTab$nameOfbar <- nameOfbar
    dummyTab$gene <- gsub("\\;.*","", dummyTab$attributes)
    dummyTab$gene <- gsub(".*=","", dummyTab$gene)
    #setwd(pathToZues)
    #print(j)
    
    
    combined <- dummyTab[,c("sampleName", "nameOfbar", "bed", "report", "contig_id","region_id","gene","numberReads",
                            "percentBases", "tot_e2e")]
    tot_e2e <- dummyTab$tot_e2e
    combined.reduced <- c(sampleName, nameOfbar, bed, report, numberReads,
                          percentBases, tot_e2e)
    
    
    tableOfRNA <- rbind(tableOfRNA, combined)
    tableOfRNA.nonAmp <- rbind(tableOfRNA.nonAmp, combined.reduced)
  }
}

###stopped here during processing
tableOfRNA <- data.frame(tableOfRNA, stringsAsFactors = FALSE)
colnames(tableOfRNA) <- tableOfRNA[1,]
tableOfRNA <- tableOfRNA[-1,]
rownames(tableOfRNA) <- NULL

finalRNATab <- tableOfRNA




unique(tableOfRNA$bed)

dnaBedNames <- NULL
dnaBedNames <- c(dnaBedNames, grep("OCP", unique(tableOfRNA$bed)))
dnaBedNames <- c(dnaBedNames, grep("CCP", unique(tableOfRNA$bed)))
dnaBedNames <- unique(tableOfRNA$bed)[dnaBedNames]


listOfDnaBeds <- NULL
for(i in seq_along(dnaBedNames)){
  listOfDnaBeds <- c(listOfDnaBeds,which(tableOfRNA$bed  == dnaBedNames[i]))
}

finalRNATab <- finalRNATab[-listOfDnaBeds,]

colnames(finalRNATab) <- c("SampleName","Barcode","Bed","ReportID","contig_id","AmpliconID","Gene","NumberOfMappedReads",
                           "PercentOfBasesMapped","Total.e2e")

reportName <- NULL
for(i in 1:nrow(finalRNATab)){
  reportName <- c(reportName,strsplit(finalRNATab$ReportID[i],"/", fixed = TRUE)[[1]][5])
}

###switch path of the report to another ID and then readding just the reportname
finalRNATab$path <- finalRNATab$ReportID
finalRNATab$ReportID <- reportName



write.table(finalRNATab, "/mnt/DATA4/kevhu/apps/fetchApp/data/fullRNACountTable3.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)

rdsFile <- read.table("/mnt/DATA4/kevhu/apps/fetchApp/data/fullRNACountTable3.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

saveRDS(rdsFile, "/mnt/DATA4/kevhu/apps/fetchApp/data/data.RDS")
###Below is depreciated version of splitting RNA files and uploading them onto a sql database

###splitting large RNA file into individual tables by bed
#listOfBedsToSplit <- c(unique(finalRNATab$Bed))
#listOfSplitTables <- NULL

#for(i in seq_along(listOfBedsToSplit)){
#  a <- NULL
#  a <- finalRNATab[which(finalRNATab$Bed == listOfBedsToSplit[i]),]
#  assign(paste0("RnaAmpTab_", listOfBedsToSplit[i]), a)
#  listOfSplitTables <- c(listOfSplitTables, paste0("RnaAmpTab_", listOfBedsToSplit[i]))
#}



###pushing all tables 
#library(DBI)
#library(RMySQL)
#m<-dbDriver("MySQL");
#con<-dbConnect(m,user='kevin',password='sat1840',host="localhost",dbname='AnnoDB')

#for(i in seq_along(listOfSplitTables)){
#  
#  dbWriteTable(con, value = eval(as.name(listOfSplitTables[i])), name = listOfSplitTables.names[i], overwrite=TRUE)
#}
#
#dbDisconnect(con)

