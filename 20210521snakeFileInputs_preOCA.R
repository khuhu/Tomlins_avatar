#!/usr/bin/env Rscript

library(jsonlite);
library(stringr);

#pathToZues <- "/mnt/data_dump/Zeus_data_dump/";
pathToYoda <- "/mnt/DATA3/Yoda_data_dump/";
pathToEros <- "/mnt/DATA3/eros_tmp/";
pathToCetus <- "/mnt/DATA3/cetus_data_dump/"

listofSequencers <- c(pathToYoda, pathToEros, pathToCetus);

listOfBeds <- NULL

for(i in seq_along(listofSequencers)){
  setwd(listofSequencers[i]);
  list3 <- system('find . -name *.gc.bed*  | grep "local_beds" | grep -v "dummyCov"', intern = TRUE);
  list3 <- sub("./",listofSequencers[i], list3);
  
  list3_1 <- str_remove_all(list3, "plugin_out.*")
  
  list3_2 <- cbind(list3_1, list3)
  listOfBeds <- rbind(listOfBeds, list3_2)
}

listOfBeds <- as.data.frame(listOfBeds, stringsAsFactors = FALSE)
colnames(listOfBeds) <- c("directories", "bed")
### can intserct the list of directories that use different mouse bed files
mouseIndices <- grep("IAD124056_167_Designed.gc.bed", listOfBeds$bed)
mouseIndices2 <- grep("IAD202670_167_Designed.gc.bed", listOfBeds$bed)
mouseIndices3 <- c(mouseIndices, mouseIndices2)

mouseBeds <- listOfBeds[mouseIndices3,]
mouseBeds$reports <- sapply(mouseBeds$directories, FUN = function(x) unlist(str_split(x, "/"))[5])
mouseBeds$snakemakeInput <- str_remove(mouseBeds$bed, "local_beds.*")

### from the listOfBed files, just create a separate dataframe for each bed, then do below.
### can later make below a simple function


summaryFileList <- NULL
for (i in mouseBeds$directories) {
  setwd(i)
  summaryFile <- system('find . -name *bc_summary.xls | grep -v "scraper" | grep -v "dummyCov"', intern = TRUE)
  summaryFile <- sub("./", i, summaryFile)
  summaryFileList <- c(summaryFileList, summaryFile)
}

mouseBeds$summaryFile <- summaryFileList
mouseBeds$idxFile <- paste0(mouseBeds$reports, "/", "idx.txt")


mouseBeds$cnBed <- str_remove(mouseBeds$bed, ".*local_beds/")
mouseBeds$cnBed <- str_replace(mouseBeds$cnBed, "IAD124056_167_Designed.gc.bed", "IAD124056_167_Designed.del.nopool.gc.bed")
#mouseBeds$cnBed <- str_replace(mouseBeds$cnBed, "IAD202296_Designed.gc.bed", "placeholder.txt")

mouseBeds$normals <- mouseBeds$cnBed
mouseBeds$normals <- str_replace(mouseBeds$normals, "IAD124056_167_Designed.del.nopool.gc.bed",
                                 "/mnt/DATA6/mouseData/normals/20191105n9.txt")
mouseBeds$normals <- str_replace(mouseBeds$normals, "IAD202670_167_Designed.gc.bed",
                                 "/mnt/DATA6/mouseData/normals/20210518.n7.mouseGeneral.txt")

mouseBeds$normalsId <- mouseBeds$cnBed
mouseBeds$normalsId <- str_replace(mouseBeds$normalsId, "IAD124056_167_Designed.del.nopool.gc.bed",
                                   "/mnt/DATA6/mouseData/normals/20191105n9.IDlist.txt")
mouseBeds$normalsId <- str_replace(mouseBeds$normalsId, "IAD202670_167_Designed.gc.bed",
                                   "/mnt/DATA6/mouseData/normals/20210518.n7.mouseGeneral.IDlist.txt")


# change if I use a shortened version for annotation
#mouseBeds$symlinkFrom <- "empty"
#for (i in 1:nrow(mouseBeds)) {
#  if (mouseBeds$cnBed[i] == "IAD202670_167_Designed.gc.bed") {
#    mouseBeds$symlinkFrom[i] <- "/home/kevhu/programs/annovar/mousedb/IAD202670_167_Designed_mm10_mpgpv6_Indels.avinput"
#  }
#  if (mouseBeds$cnBed[i] == "IAD124056_167_Designed.del.nopool.gc.bed") {
#    mouseBeds$symlinkFrom[i] <- "/home/kevhu/programs/annovar/mousedb/IAD124056_167_Designed_mm10_mpgpv6_Indels.avinput"
#  }
#}




### making the necessary directories

for (i in seq_along(mouseBeds$reports)) {
  setwd("/mnt/DATA6/mouseData/copynumber/")
  if(!dir.exists(mouseBeds$reports[i])){
    system(sprintf("mkdir -p %s", mouseBeds$reports[i]))
    setwd(mouseBeds$reports[i])
    tmpDir <- getwd()
    system(sprintf("echo %s > tmp.txt", mouseBeds$snakemakeInput[i]))
    system(sprintf("cp %s %s", paste0("/mnt/DATA6/mouseData/bedFiles/",mouseBeds$cnBed[i]), "./bed.txt"))
    system(sprintf("cp %s %s", mouseBeds$normalsId[i], "./idList.txt"))
    system(sprintf("cp %s %s", mouseBeds$normals[i], "./normals.txt"))
  }
}


### for vcf sections

for (i in seq_along(mouseBeds$reports)) {
  setwd("/mnt/DATA6/mouseData/vcfs/")
  if(!dir.exists(mouseBeds$reports[i])){
    system(sprintf("mkdir -p %s", mouseBeds$reports[i]))
  }
}



### i should make this into append rather than rewrite each time - only way for time to 
### be properly recorded too
write.table(mouseBeds, "/mnt/DATA6/mouseData/mouseBedsIdx.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

