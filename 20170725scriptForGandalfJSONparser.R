###This is an Rscript to pull unique sampleIDs along with the stat - uniformity of base cov. - from old Gandalf runs
###Sample name consists of sample name, barcode, bed, and reportnum - all these except reportnum can be found in the 
###following generic path: /mnt/DATA5/Gandalf_data_keep/runName(variable)/plugin_out/coverageAnalysis_out/results.json
###following is path for reportnum: /mnt/DATA5/Gandalf_data_keep/primary.key


###So the problem that I've run into is that not every coverageAnalysis_out dir. is labeled the same
###Also two of the reports don't have a results.json for Cov_analysis 1. 163_re_168 &

library(jsonlite)
### listOfRuns will only be used to access the primary key (result) - found alternative way to do it
#listOfRuns <- system(command = "ls --ignore='*.gz*' /mnt/DATA5/Gandalf_data_keep/", intern = TRUE)

pathToGandalf <- "/mnt/DATA5/Gandalf_data_keep/"
setwd(dir = pathToGandalf)

###used to find results.json to parse everything else
listOfJSONs <- system('find . -name *results.json*  | grep "coverageAnalysis_out"', intern = TRUE)
listOfStartPlugins <- system('find . -name *startplugin.json*  | grep "coverageAnalysis_out"', intern = TRUE)
###temporary fix
blacklist <- c("AC_PR49_51_53_57_good_113","Auto_user_SN2-101-BL211_chp_ocp_129_200",
  "barcode_141")
blacklist.index <- NULL
for(i in seq_along(blacklist)){
  a <- grep(blacklist[i], listOfJSONs)
  blacklist.index <- c(blacklist.index, a)
}
listOfJSONs <- listOfJSONs[-blacklist.index]

tableOfFun <- NULL
tableOfFun <- c("sample", "barcode","bed", "Report","Reportnum","Uniformity","listJSONkey","numberMappedReads","percentMapped","avgCov")
for(i in seq_along(listOfJSONs)){
  #print(i)
  g <- fromJSON(listOfJSONs[i])
  for(j in seq_along(g$barcodes)){
    sampleName <- g$barcodes[[j]]$`Sample Name`
    bed <- g$barcodes[[j]]$`Targeted Regions`
    #print(i)
    #print(bed)
    uniformity <- g$barcodes[[j]]$`Uniformity of base coverage`
    percentBases <- g$barcodes[[j]]$`Percent base reads on target`
    numberReads <- g$barcodes[[j]]$`Number of mapped reads`
    avgCov <- g$barcodes[[j]]$`Average base coverage depth`
    report <- listOfJSONs[i]
    nameOfbar <- names(g$barcodes[j])
    dummy <- unlist(strsplit(listOfJSONs[i],"/"))
    dummy2 <- paste0(pathToGandalf,dummy[2],"/primary.key")
    address <- paste0(dummy[1],"/",dummy[2],"/",dummy[3],"/",dummy[4],"/","*summary*")
    address2 <- system(paste("ls",address),intern = TRUE)
    reportNum <- system(paste("less",dummy2), intern = TRUE)
    listJSONkey <- i
    if(length(uniformity) == 0){
      next()
    }
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
          #print(z)
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
      if(length(sampleName) == 0 || grepl("IonXpress",sampleName) == TRUE || grepl("HaloPlex", sampleName)){
        sampleName <- "None"
      }
    }
    if(length(bed) == 0){
      print(i)
      print(bed)
      bed <- g$barcodes[[j]]$`Target Regions`
      #if(length(bed) == 0){
      #varJSON <- fromJSON(listOfStartPlugins[grep(dummy[2], listOfStartPlugins)])
      #bed <- varJSON$runinfo$plugin$userInput$sample
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
    combined <- c(sampleName, nameOfbar, bed, report,reportNum, uniformity, listJSONkey, numberReads, percentBases, avgCov)
    #print(combined)
    tableOfFun <- rbind(tableOfFun, combined)
  }
}

tableOfFun <- as.data.frame(tableOfFun, stringsAsFactors = FALSE)
colnames(tableOfFun) <- c("Sample","Barcode","Bed", "Report","Reportnum","listOfJSONkey","Uniformity","numberMappedReads","percentMapped","avgCov")
rownames(tableOfFun) <- NULL
tableOfFun <- tableOfFun[-1,]
#haloPlexData <- tableOfFun[grep("HaloPlex",tableOfFun$Barcode),]

###If there is no Uniformity it is useless in our case. So I will subset out ones without it

tableOfFun.filtered <- tableOfFun
tableOfFun.filtered$Sample[is.na(tableOfFun$Sample)] <- "None"
tableOfFun.filtered$tmp <- paste0(tableOfFun.filtered$Sample,tableOfFun.filtered$Barcode,tableOfFun.filtered$Bed,
                                  tableOfFun.filtered$Reportnum)

length(unique(tableOfFun.filtered$tmp))




###If there is no Uniformity it is useless in our case. So I will subset out ones without it

length(which(tableOfFun.filtered$Bed == "None"))
21/1293 #have no bed .. so I did a pretty swell job haha....
length(which(tableOfFun.filtered$Sample == "None")) 
454/1293 # no sample names, update: tried many fixes, but the way I implemented is probably the most complete
# not all names are best labeled. some are just numbers .... Note I did not use expMeta.dat b/c it is inaccurate and
# only displays the first of many sampels in the set - originally used it, but found it to be wrong afterwards
#update 08/03/2017 - dropped from 460 -> ~ 449 (was 60, but ~300 were just barcodes) -> weird truncation with haloplex
tableOfFun.filtered$tmp <- NULL
write.table(x = tableOfFun.filtered, file = "/home/kevhu/data/20170727GandalfTable.tsv",
            sep = '\t',row.names = FALSE)













###to fix sample names. going to pull it from expMeta.dat file - testing extraction
testingtesting <- unlist(strsplit(x = system('grep "Sample" /mnt/DATA5/Gandalf_data_keep/SN2-67_reanalyse_hotspot_156/expMeta.dat',
                                   intern = TRUE), split = " "))
testingtesting[3]


###testing the extraction before I loop it through every file in Zeus
g <- fromJSON(listOfJSONs[100])
length(g$barcodes)
names(g$barcodes)
testTable <- NULL
for(i in seq_along(g$barcodes)){
  sampleName <- g$barcodes[[i]]$`Sample Name`
  bed <- g$barcodes[[i]]$`Targeted Regions`
  uniformity <- g$barcodes[[i]]$`Uniformity of base coverage`
  nameOfbar <- names(g$barcodes[i])
  dummy <- unlist(strsplit(listOfJSONs[i],"/"))
  dummy2 <- paste0(pathToGandalf,dummy[2],"/primary.key")
  reportNum <- system(paste("less",dummy2), intern = TRUE)
  combined <- c(sampleName, nameOfbar, bed, reportNum, uniformity)
  print(combined)
  testTable <- rbind(testTable, combined)
}

colnames(testTable) <- c("Sample","Barcode","Bed", "Reportnum","Uniformity")

