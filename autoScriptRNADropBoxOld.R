#!/usr/bin/env Rscript

.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))

source("/mnt/DATA4/kevhu/scripts/fastReadFile.R")
source("/mnt/DATA4/kevhu/scripts/parseJsons.R")
library(jsonlite);
library(rsconnect);
library(parallel);

### the first part of script figures out which tables/files are DNA, based on statistics that only DNA sequencing have
### some portions of the DNA portion are not optimized and are just artefacts from previous script
### these are just paths to our data dumps and shiny app is where I keep my giant table of data
### part of this script can be ignored since it deals with automatically updating the table on shinyapp the lab uses

pathToZues <- "/mnt/DATA3/Zeus_data_dump/";
pathToYoda <- "/mnt/DATA3/Yoda_data_dump/";
pathToShinyApp <- "/mnt/DATA4/kevhu/apps/fetchApp/";
listofSequencers <- c(pathToYoda, pathToZues);


### collecting all paths to coverage files

listOfJSONs <- NULL;
listOfStartPlugins <-NULL;

for(i in seq_along(listofSequencers)){
  setwd(listofSequencers[i]);
  list1 <- system('find . -name *results.json*  | grep "coverageAnalysis_out"', intern = TRUE);
  list2 <- system('find . -name *startplugin.json*  | grep "coverageAnalysis_out"', intern = TRUE);
  list1 <- sub("./",listofSequencers[i], list1);
  list2 <- sub("./",listofSequencers[i], list2);
  listOfJSONs <- c(listOfJSONs, list1);
  listOfStartPlugins <- c(listOfStartPlugins, list2);
}

print("Done collecting all covs")


### making a table of information from each JSON file to search against

listOfRNAJsons <- NULL;
tableOfFun <- NULL;
tableOfFun <- c("sample", "barcode","bed", "report","Uniformity","listJSONkey","numberMappedReads","percentMapped","avgCov");
for(i in seq_along(listOfJSONs)){
  #makeshift tryCatch function for nonaccessible directories - can change later
  readableDir <- tryCatch(fromJSON(listOfJSONs[i]), error = function(x) return(NULL));
  if(is.null(readableDir)){
    next();
  }
  g <- fromJSON(listOfJSONs[i]);
  for(j in seq_along(g$barcodes)){
    sampleName <- g$barcodes[[j]]$`Sample Name`;
    bed <- g$barcodes[[j]]$`Targeted Regions`;
    uniformity <- g$barcodes[[j]]$`Uniformity of base coverage`;
    percentBases <- g$barcodes[[j]]$`Percent base reads on target`;
    numberReads <- g$barcodes[[j]]$`Number of mapped reads`;
    avgCov <- g$barcodes[[j]]$`Average base coverage depth`;
    nameOfbar <- names(g$barcodes[j]);
    report <- listOfJSONs[i];
    dummy <- unlist(strsplit(listOfJSONs[i],"/"));
    address <- sub("results.json","", listOfJSONs[i]);
    address2 <- system(paste("find",address, "-name *bc_summary*"),intern = TRUE);
    address2 <- address2[-which(grepl("scraper", address2))];
    listJSONkey <- i;
    if(length(uniformity) == 0){
      listOfRNAJsons <- c(listOfRNAJsons, listOfJSONs[i]);
      next();
    }
    if(length(sampleName) == 0 || sampleName == "None"){
      dummyFindSampleVar <- fromJSON(listOfStartPlugins[grep(dummy[2], listOfStartPlugins)][1]);
      test <- tryCatch(fromJSON(dummyFindSampleVar$plan$barcodedSample), error = function(x) return(NULL));
      if(!is.null(test)){
        dummyFindSampleVar2 <- fromJSON(dummyFindSampleVar$plan$barcodedSample);
        z <- names(dummyFindSampleVar2[grep(nameOfbar, dummyFindSampleVar2)]);
        if(is.list(z)){
          sampleName <- grep(nameOfbar, dummyFindSampleVar2);
        }
        if(length(z) > 0){
          print(z);
          sampleName <- unlist(z);
        }
      }
      if(length(sampleName) == 0 || sampleName == "None"){
        testingtesting <- unlist(strsplit(system(paste("less",address2,"| grep",nameOfbar),
                                                 intern = TRUE),"\t"));
        if(is.character(testingtesting)){
          sampleName <- testingtesting[2];
        }
      }
      if(length(sampleName) == 0 || grepl("IonXpress",sampleName) == TRUE){
        sampleName <- "None";
      }
    }
    if(length(bed) == 0){
      bed <- g$barcodes[[j]]$`Target Regions`;
      if(length(bed) == 0){
        bed <- "None";
        #}
      }
    }
    if(length(numberReads) == 0){
      numberReads <- "None";
    }
    if(length(percentBases) == 0){
      percentBases <- "None";
    }
    if(length(avgCov) == 0){
      avgCov <- "None";
    }
    combined <- c(sampleName, nameOfbar, bed, report, uniformity, listJSONkey, numberReads, percentBases, avgCov);
    tableOfFun <- rbind(tableOfFun, combined);
  }
}




tableOfFun <- as.data.frame(tableOfFun, stringsAsFactors = FALSE);
colnames(tableOfFun) <- tableOfFun[1,];
rownames(tableOfFun) <- NULL;
tableOfFun <- tableOfFun[-1,];


### loading a previously made table to see if any entry needs to be added. below should create an empty dataframe if nothign new is found
### there is stop function somwhere below

#previousRNATab <- read.table("/mnt/DATA4/kevhu/apps/fetchApp/data/fullRNACountTable3.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
previousRNATab <- faster.readfile("/mnt/DATA4/kevhu/apps/fetchApp/data/fullRNACountTable3.txt", 0)


### don't get rid of path names or it'll mess up - uncomment this to check for paths afterwords

listOfRNAJsons2 <- (unique(listOfRNAJsons));
listOfPreviousRNAPaths <- previousRNATab$path;
listOfPreviousRNAPaths <- unique(unlist(listOfPreviousRNAPaths));
alreadyProcessed <- which(listOfRNAJsons2 %in% listOfPreviousRNAPaths);
listOfRNAJsons2 <- listOfRNAJsons2[-c(alreadyProcessed)];

###RNA part
###later add an error file for files not being read



#listOfRNAJsons <- (unique(listOfRNAJsons))
tableOfRNA <- NULL;
tableOfRNA.nonAmp <- NULL;
tableOfRNA <- c("sampleName", "nameOfbar", "bed", "report", "contig_id","region_id","gene","numberReads",
                "percentBases", "tot_e2e");
tableOfRNA.nonAmp <-c("sampleName", "nameOfbar", "bed", "report", "numberReads",
                      "percentBases", "tot_e2e");

listofErrors <- NULL;

for(i in seq_along(listOfRNAJsons2)){
  #makeshift tryCatch function for nonaccessible directories - can change later
  readableDir <- tryCatch(fromJSON(listOfRNAJsons2[i]), error = function(x) return(NULL));
  if(is.null(readableDir)){
    next()
  }
  g <- fromJSON(listOfRNAJsons2[i])
  for(j in seq_along(g$barcodes)){
    errorCheck <- g$barcodes[[j]]$`Error`
    if(length(errorCheck) > 0){
      listofErrors <- c(listofErrors, listOfRNAJsons2[i]) 
      next()
    }
    sampleName <- g$barcodes[[j]]$`Sample Name`
    bed <- g$barcodes[[j]]$`Targeted Regions`
    percentBases <- g$barcodes[[j]]$`Percent base reads on target`
    numberReads <- g$barcodes[[j]]$`Number of mapped reads`
    nameOfbar <- names(g$barcodes[j])
    report <- listOfRNAJsons2[i]
    dummy <- unlist(strsplit(listOfRNAJsons2[i],"/"))
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
    dirForCovFiles <- sub("results.json","", listOfRNAJsons2[i])
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
    dummyTab <- read.table(cov.path2, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
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

    combined <- dummyTab[,c("sampleName", "nameOfbar", "bed", "report", "contig_id","region_id","gene","numberReads",
                            "percentBases", "tot_e2e")]
    tot_e2e <- dummyTab$tot_e2e
    combined.reduced <- c(sampleName, nameOfbar, bed, report, numberReads,
                          percentBases, tot_e2e)
    
    
    tableOfRNA <- rbind(tableOfRNA, combined)
    tableOfRNA.nonAmp <- rbind(tableOfRNA.nonAmp, combined.reduced)
  }
}


###20191001: to get full RNA table produced in steps above and environment load file below
#load("/mnt/DATA4/kevhu/apps/fetchApp/data/20191001postFullRnaTab.RData")

tableOfRNA <- data.frame(tableOfRNA, stringsAsFactors = FALSE)
colnames(tableOfRNA) <- tableOfRNA[1,]
tableOfRNA <- tableOfRNA[-1,]
rownames(tableOfRNA) <- NULL

finalRNATab <- tableOfRNA
unique(tableOfRNA$bed)

#dnaBedNames <- NULL
#dnaBedNames <- c(dnaBedNames, grep("OCP", unique(tableOfRNA$bed)))
#dnaBedNames <- c(dnaBedNames, grep("CCP", unique(tableOfRNA$bed)))
#dnaBedNames <- unique(tableOfRNA$bed)[dnaBedNames]


#listOfDnaBeds <- NULL
#for(i in seq_along(dnaBedNames)){
#  listOfDnaBeds <- c(listOfDnaBeds,which(tableOfRNA$bed  == dnaBedNames[i]))
#}

#finalRNATab <- finalRNATab[-listOfDnaBeds,]

if(nrow(finalRNATab) == 0){
  print("Done. Nothing to update.")
  stop()
}

colnames(finalRNATab) <- c("SampleName","Barcode","Bed","ReportID","contig_id","AmpliconID","Gene","NumberOfMappedReads",
                           "PercentOfBasesMapped","Total.e2e")



reportName <- NULL
for(i in 1:nrow(finalRNATab)){
  reportName <- c(reportName,strsplit(finalRNATab$ReportID[i],"/", fixed = TRUE)[[1]][5])
}


###example of faster stringsplit using parallel processes below:
#library(parallel)
#detectCores()
#reportIds <- finalRNATab$ReportID
#strSplitFunc <- function(x) {
#  splitRes <- strsplit(x,"/", fixed = TRUE)[[1]][5]
#  return(splitRes)
#}
#results <- mclapply(reportIds, strSplitFunc,mc.cores = 25)




###switch path of the report to another ID and then readding just the reportname
finalRNATab$path <- finalRNATab$ReportID
finalRNATab$ReportID <- reportName


### 20190729: got rid of the path column, should save some space for shinyapps.io for the meantime .... problem will come up again later
### if problems do persist in the future .... possible fixes are actual paid server, getting rid of barcode and bed columns - need to change search app algorithm
### - blacklisting and get rid of Lorena's whole transcriptome data - to do this portion, I would like to create two separate file indexes -
### one comprehensive with all and one just for targeted data to be used by the app



write.table(finalRNATab, "/mnt/DATA4/kevhu/apps/fetchApp/data/fullRNACountTable3.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)

rdsFile <- faster.readfile("/mnt/DATA4/kevhu/apps/fetchApp/data/fullRNACountTable3.txt", 0)
rdsFile <- rdsFile[,c(1:8,10)]

saveRDS(rdsFile, "/mnt/DATA4/kevhu/apps/fetchApp/data/data.RDS")


setwd(pathToShinyApp);
secret <- unlist(read.table("/mnt/DATA4/kevhu/apps/fetchApp/secret", stringsAsFactors = FALSE))
token <- unlist(read.table("/mnt/DATA4/kevhu/apps/fetchApp/token", stringsAsFactors = FALSE))
setAccountInfo("kevhu", token = token, secret = secret)
deployApp(account = "kevhu", appName = "fetchApp",appFiles = c("app.R", "data/PR_tissue_v1.panel_annotation.xlsx",".httr-oauth","data/data.RDS","data/droptoken.rds"),
        upload = TRUE, server = "shinyapps.io", launch.browser = FALSE, forceUpdate = TRUE) 

restartApp(appName = "fetchApp")

