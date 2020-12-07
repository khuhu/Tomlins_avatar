#!/usr/bin/env Rscript
.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))

### need to find a way to create an index of only mouse samples

library(jsonlite);
library(stringr);

#pathToZues <- "/mnt/DATA3/Zeus_data_dump/";
pathToYoda <- "/mnt/DATA3/Yoda_data_dump/";
pathToEros <- "/mnt/DATA3/eros_tmp/";
pathToCetus <- "/mnt/DATA3/cetus_data_dump/"

listofSequencers <- c(pathToYoda, pathToEros, pathToCetus);

#listOfJSONs <- NULL;
#listOfStartPlugins <-NULL;
listOfBeds <- NULL

for(i in seq_along(listofSequencers)){
  setwd(listofSequencers[i]);
  list3 <- system('find . -name *.gc.bed*  | grep "local_beds"', intern = TRUE);
  list3 <- sub("./",listofSequencers[i], list3);
  
  list3_1 <- str_remove_all(list3, "plugin_out.*")
  
  list3_2 <- cbind(list3_1, list3)
  listOfBeds <- rbind(listOfBeds, list3_2)
}

listOfBeds <- as.data.frame(listOfBeds, stringsAsFactors = FALSE)
colnames(listOfBeds) <- c("directories", "bed")
### can intserct the list of directories that use different mouse bed files
mouseIndices <- grep("IAD124056_167_Designed.gc.bed", listOfBeds$bed)
mouseIndices2 <- grep("IAD202296_167_Designed.gc.bed", listOfBeds$bed)
mouseIndices3 <- c(mouseIndices, mouseIndices2)

mouseBeds <- listOfBeds[mouseIndices3,]
mouseBeds$reports <- sapply(mouseBeds$directories, FUN = function(x) unlist(str_split(x, "/"))[5])
mouseBeds$snakemakeInput <- str_remove(mouseBeds$bed, "local_beds.*")

### from the listOfBed files, just create a separate dataframe for each bed, then do below.
### can later make below a simple function


summaryFileList <- NULL
for (i in mouseBeds$directories) {
  setwd(i)
  summaryFile <- system('find . -name *bc_summary.xls | grep -v "scraper"', intern = TRUE)
  summaryFile <- sub("./", i, summaryFile)
  summaryFileList <- c(summaryFileList, summaryFile)
}

mouseBeds$summaryFile <- summaryFileList
mouseBeds$idxFile <- paste0(mouseBeds$reports, "/", "idx.txt")


mouseBeds$cnBed <- str_remove(mouseBeds$bed, ".*local_beds/")
mouseBeds$cnBed <- str_replace(mouseBeds$cnBed, "IAD124056_167_Designed.gc.bed", "IAD124056_167_Designed.del.nopool.gc.bed")
mouseBeds$cnBed <- str_replace(mouseBeds$cnBed, "IAD202296_Designed.gc.bed", "placeholder.txt")

mouseBeds$normals <- mouseBeds$cnBed
mouseBeds$normals <- str_replace(mouseBeds$normals, "IAD124056_167_Designed.del.nopool.gc.bed",
                                 "/home/kevhu/data/normals/mousePanel/20191105n9.txt")
mouseBeds$normals <- str_replace(mouseBeds$normals, "IAD202296_Designed.gc.bed",
                                 "/home/kevhu/data/normals/mousePanel/placeholder.txt")

mouseBeds$normalsId <- mouseBeds$cnBed
mouseBeds$normalsId <- str_replace(mouseBeds$normalsId, "IAD124056_167_Designed.del.nopool.gc.bed",
                                   "/home/kevhu/data/normals/mousePanel/20191105n9.IDlist.txt")
mouseBeds$normalsId <- str_replace(mouseBeds$normalsId, "IAD202296_Designed.gc.bed",
                                   "/home/kevhu/data/normals/mousePanel/placeholder.txt")


#fullPathRes <- NULL
#for (i in seq_along(mouseBeds$summaryFile)) {
#  tmpFile <- read.table(mouseBeds$summaryFile[i], sep = "\t", stringsAsFactors = FALSE, header = TRUE)
#  tmpRes <- paste0(mouseBeds$reports, "/", tmpFile$`Sample.Name`, ".pdf")
#  dateTmp <- date()
#  dateTmp <- str_replace_all(dateTmp, pattern = " ", replacement = "_")
#  tmpTable <- cbind(rep(dateTmp, length(tmpRes)),
#                    rep(mouseBeds$bed, length(tmpRes)),
#                    tmpRes)
#  fullPathRes <- rbind(fullPathRes, tmpTable)
#}

#colnames(fullPathRes) <- c("Date", "Bed_name", "fullPath")



### making the necessary directories

for (i in seq_along(mouseBeds$reports)) {
  setwd("/home/kevhu/scripts/newMousePanelPipeline/reports/")
  if(!dir.exists(mouseBeds$reports[i])){
    system(sprintf("mkdir -p %s", mouseBeds$reports[i]))
    setwd(mouseBeds$reports[i])
    tmpDir <- getwd()
    system(sprintf("echo %s > tmp.txt", mouseBeds$snakemakeInput[i]))
    system(sprintf("cp %s %s", paste0("/home/kevhu/data/bedFiles/",mouseBeds$cnBed[i]), "./bed.txt"))
    system(sprintf("cp %s %s", mouseBeds$normalsId[i], "./idList.txt"))
    system(sprintf("cp %s %s", mouseBeds$normals[i], "./normals.txt"))
  }
}


### for vcf sections

for (i in seq_along(mouseBeds$reports)) {
  setwd("/home/kevhu/scripts/newMousePanelPipeline/vcfs/")
  if(!dir.exists(mouseBeds$reports[i])){
    system(sprintf("mkdir -p %s", mouseBeds$reports[i]))
  }
}



### i should make this into append rather than rewrite each time - only way for time to 
### be properly recorded too
write.table(mouseBeds, "/home/kevhu/scripts/newMousePanelPipeline/mouseBedsIdx.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#write.table(fullPathRes, "/home/kevhu/scripts/newMousePanelPipeline/sampleResult.txt"
#            ,sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
