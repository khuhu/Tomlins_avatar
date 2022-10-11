#!/usr/bin/env Rscript
.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))

library(jsonlite);
library(stringr);

#pathToZues <- "/mnt/DATA3/Zeus_data_dump/";
pathToYoda <- "/mnt/DATA3/Yoda_data_dump/";
pathToEros <- "/mnt/DATA3/eros_tmp/";
pathToCetus <- "/mnt/DATA3/cetus_data_dump/";

listofSequencers <- c(pathToYoda, pathToEros, pathToCetus);


listOfBeds <- NULL

for(i in seq_along(listofSequencers)){
  setwd(listofSequencers[i]);
  list3 <- system('find . -name *.gc.bed*  | grep "local_beds"',
                  intern = TRUE, show.output.on.console = FALSE);
  list3 <- sub("./",listofSequencers[i], list3);
  
  list3_1 <- str_remove_all(list3, "plugin_out.*")
  
  list3_2 <- cbind(list3_1, list3)
  listOfBeds <- rbind(listOfBeds, list3_2)
}

print("getting all dirs - done")

listOfBeds <- as.data.frame(listOfBeds, stringsAsFactors = FALSE)
colnames(listOfBeds) <- c("directories", "bed")
### can intserct the list of directories that use different mouse bed files
mouseIndices <- grep("IAD124056_167_Designed.gc.bed", listOfBeds$bed)
mouseIndices2 <- grep("IAD202670_167_Designed.gc.bed", listOfBeds$bed)
mouseIndices3 <- c(mouseIndices, mouseIndices2)

mouseBeds <- listOfBeds[mouseIndices3,]
mouseBeds$reports <- sapply(mouseBeds$directories, FUN = function(x) unlist(str_split(x, "/"))[5])
mouseBeds$snakemakeInput <- str_remove(mouseBeds$bed, "local_beds.*")

summaryFileList <- NULL
variantDirList <- NULL
for (i in mouseBeds$directories) {
  setwd(i)
  summaryFile <- system('find . -name *bc_summary.xls | grep -v "scraper" | grep -v "dummyCov"', intern = TRUE)
  summaryFile <- sub("./", i, summaryFile)
  
  summaryFileList <- c(summaryFileList, summaryFile)
  
  variantDir <- system('find . -type d -name variantCaller*', intern = TRUE)
  variantDir <- sub("./", i, variantDir)
  variantDir <- paste0(variantDir, "/")
  variantDirList <- c(variantDirList, variantDir)
}

print("getting variant list - done")

# 20210516: KH quick fix for reports not listed in mouse beds (?)
#summaryFileList <- summaryFileList[grep(paste(mouseBeds$directories, collapse = "|"), summaryFileList)]

print(summaryFileList)
mouseBeds$summaryFile <- summaryFileList
mouseBeds$idxFile <- paste0(mouseBeds$reports, "/", "idx.txt")
mouseBeds$variantDir <- variantDirList


### create table for each report I can iterate for file copying
### just use if file exist function for each entry of the should be vcf file

cpTable <- NULL
for (i in seq_along(mouseBeds$summaryFile)) {
  tmpLnTable <- NULL
  tmpTable <- read.table(mouseBeds$summaryFile[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  setwd(mouseBeds$variantDir[i])
  ln1 <- system('find . -name "TSVC_variants.vcf.gz"', intern = TRUE)
  ln1 <- ln1[order(ln1)]
  ln2 <- str_replace(ln1, "./", mouseBeds$variantDir[i])
  ln3 <- NULL
  for (j in seq_along(ln1)) {
    tmpLn1 <- str_remove(ln1[j], "/IonXpress.*")
    tmpLn1 <- paste0(tmpLn1, "/",tmpTable$Sample.Name[j],".vcf.gz")
    #tmpLn1 <- paste0("/home/kevhu/scripts/newMousePanelPipeline/vcfs/", str_replace(tmpLn1, "\\.", mouseBeds$reports[i]))
    tmpLn1 <- paste0("./", str_replace(tmpLn1, "\\.", mouseBeds$reports[i]))
    print(tmpLn1)
    ln3 <- c(ln3, tmpLn1)
  }
  tmpLnTable <- data.frame("from" = ln2, "to" = ln3, stringsAsFactors = FALSE)
  cpTable <- rbind(cpTable, tmpLnTable)
}

cpTable$V3 <- paste0(cpTable$from, ".tbi")
cpTable$V4 <- paste0(cpTable$to, ".tbi")
colnames(cpTable) <- c("fromVcf", "toVcf", "fromTbi", "toTbi")

print("finished list of input for copy - done")

### from list above I iteratre and make a if exist function
### setdir of newMousePipeline - separate dir for vcfs?

setwd("/home/kevhu/scripts/newMousePanelPipeline/vcfs/")

for (i in 1:nrow(cpTable)) {
  print(i)
  if (!file.exists(cpTable$toVcf[i])) {
    #system(sprintf("ln -s %s %s", cpTable$fromVcf[i], cpTable$toVcf[i]))
    system(sprintf("cp %s %s", cpTable$fromVcf[i], cpTable$toVcf[i]))
  }
  if (!file.exists(cpTable$toTbi[i])) {
    #system(sprintf("ln -s %s %s", cpTable$fromVcf[i], cpTable$toVcf[i]))
    system(sprintf("cp %s %s", cpTable$fromTbi[i], cpTable$toTbi[i]))
  }
}

setwd("/home/kevhu/scripts/newMousePanelPipeline/vcfs/")
for (i in 1:nrow(mouseBeds)) {
  setwd(mouseBeds$reports[i])
  print(getwd())
  if (!file.exists("bedfile.txt")) {
    bedName <- sub(x = mouseBeds$bed[i], pattern = ".*local_beds/", replacement = "")
    writeLines(con = "bedfile.txt", text = bedName)
  }
  setwd("/home/kevhu/scripts/newMousePanelPipeline/vcfs/")
}

partFullPath <- getwd()
vcfList <- system('find -maxdepth 2 -name *vcf.gz | grep -v "norm"', intern = TRUE)
vcfList <- sub(vcfList, pattern = "\\.", replacement = partFullPath)
snakeFileOut <- sub(x = vcfList, pattern = "\\.vcf\\.gz", replacement = "")
tableForSnakemake <- data.frame("filename" = snakeFileOut, stringsAsFactors = FALSE)

reportName <- str_remove(tableForSnakemake$filename, "/home/kevhu/scripts/newMousePanelPipeline/vcfs/")
reportName <- str_remove(reportName, "/.*")

tableForSnakemake$report <- reportName
### maybe append instead of rewriting later ..
write.table(tableForSnakemake, "/home/kevhu/scripts/newMousePanelPipeline/mouseVcfTable.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

