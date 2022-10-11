#!/usr/bin/env Rscript

library(jsonlite);
library(stringr);

#pathToZues <- "/mnt/data_dump/Zeus_data_dump/";
#pathToYoda <- "/mnt/DATA3/Yoda_data_dump/";
pathToEros <- "/mnt/DATA3/eros_tmp/";
#pathToCetus <- "/mnt/DATA3/cetus_data_dump/"

#listofSequencers <- c(pathToYoda, pathToEros, pathToCetus);
listofSequencers <- c(pathToEros);

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
humanIndice <- grep("OCAPlus.20191203.designed.gc.bed", listOfBeds$bed)
humanIndice2 <- grep("4477685_CCP_designed.gc.bed", listOfBeds$bed)
humanIndice2 <- humanIndice2[-which(humanIndice2 == 80)]
humanIndice3 <- grep("IAD203665_173_Designed.gc.bed", listOfBeds$bed)
mouseIndices3 <- c(mouseIndices, mouseIndices2, humanIndice, humanIndice2, humanIndice3)

mouseBeds <- listOfBeds[mouseIndices3,]
mouseBeds$reports <- sapply(mouseBeds$directories, FUN = function(x) unlist(str_split(x, "/"))[5])
mouseBeds$snakemakeInput <- str_remove(mouseBeds$bed, "local_beds.*")

### from the listOfBed files, just create a separate dataframe for each bed, then do below.
### can later make below a simple function


summaryFileList <- NULL
i <- "/mnt/DATA3/eros_tmp/AUS5-164-A549TS_409/"
for (i in mouseBeds$directories) {
  setwd(i)
  summaryFile <- system('find . -name *bc_summary.xls | grep -v "scraper" | grep -v "dummyCov"', intern = TRUE)
  summaryFile <- sub("./", i, summaryFile)
  if (length(summaryFile) > 1) {
    summaryFile <- summaryFile[which(nchar(summaryFile) == max(nchar(summaryFile)))]
    if (length(summaryFile) > 1) {
      summaryFile <- summaryFile[order(summaryFile, decreasing = TRUE)][1]
    }
  }
  summaryFileList <- c(summaryFileList, summaryFile)
}

mouseBeds$summaryFile <- summaryFileList
mouseBeds$idxFile <- paste0(mouseBeds$reports, "/", "idx.txt")

### stopped here. use this to edit the normal files for human etc

mouseBeds$cnBed <- str_remove(mouseBeds$bed, ".*local_beds/")
mouseBeds$cnBed <- str_replace(mouseBeds$cnBed, "IAD124056_167_Designed.gc.bed", "IAD124056_167_Designed.del.nopool.gc.bed")
mouseBeds$cnBed <- str_replace(mouseBeds$cnBed, "OCAPlus.20191203.designed.gc.bed", "20210726_OCAP_gc_noheader.bed")
mouseBeds$cnBed <- str_replace(mouseBeds$cnBed, "4477685_CCP_designed.gc.bed", "CCP.noTrack.GC.bed")

### 20220211: andis ctc samples - change second line and all following sections for specific
### bed/normals/normalIdx

reportList <- c("Auto_user_AUS5-89-HayesLobOCAplusAC01_294_230", "Auto_user_AUS5-102-HayesLobOCAplusAC01_2_310_261", 
                "HayesLobOCAplusAC01_3_good_270", "Auto_user_AUS5-132-HayesLobRun1_346_329", "Auto_user_AUS5-133-HayesLobRun2_347_331",
                "Auto_user_AUS5-136-HayesLobRun3_348_339", "Auto_user_AUS5-137-HayesLobRun4_349_341", "Auto_user_AUS5-159-Hayes_LobCTC_Run5_373_388",
                "Auto_user_AUS5-231-Hayes_LobCTC_Run6_484_546")

# mouseBeds$cnBed[grep(paste(reportList, collapse = "|"), mouseBeds$directories)] <- "20210726_OCAP_gc_Andi2_noheader.bed"
mouseBeds$cnBed[grep(paste(reportList, collapse = "|"), mouseBeds$directories)] <- "20210726_OCAP_gc_Andi3_noheader.bed"



### Old mouse, New Mouse, OCAP, CCP
mouseBeds$normals <- mouseBeds$cnBed
mouseBeds$normals <- str_replace(mouseBeds$normals, "IAD124056_167_Designed.del.nopool.gc.bed",
                                 "/mnt/DATA6/mouseData/normals/20191105n9.txt")
mouseBeds$normals <- str_replace(mouseBeds$normals, "IAD202670_167_Designed.gc.bed",
                                 "/mnt/DATA6/mouseData/normals/20210518.n7.mouseGeneral.txt")

mouseBeds$normals <- str_replace(mouseBeds$normals, "20210726_OCAP_gc_noheader.bed",
                                 "/mnt/DATA6/mouseData/normals/20210726_OCAP_normals.idx.txt")

mouseBeds$normals <- str_replace(mouseBeds$normals, "CCP.noTrack.GC.bed",
                                 "/mnt/DATA6/mouseData/normals/20211110.n10.CCP.txt")

# mouseBeds$normals <- str_replace(mouseBeds$normals, "20210726_OCAP_gc_Andi2_noheader.bed",
#                                  "/mnt/DATA6/mouseData/normals/20220211ocap_WBC.txt")

mouseBeds$normals <- str_replace(mouseBeds$normals, "20210726_OCAP_gc_Andi3_noheader.bed",
                                 "/mnt/DATA6/mouseData/normals/20220211ocap_WBC.txt")

mouseBeds$normals <- str_replace(mouseBeds$normals, "IAD203665_173_Designed.gc.bed",
                                 "/mnt/DATA6/mouseData/normals/20220914pgu2.txt")


### id list

mouseBeds$normalsId <- mouseBeds$cnBed
mouseBeds$normalsId <- str_replace(mouseBeds$normalsId, "IAD124056_167_Designed.del.nopool.gc.bed",
                                   "/mnt/DATA6/mouseData/normals/20191105n9.IDlist.txt")
mouseBeds$normalsId <- str_replace(mouseBeds$normalsId, "IAD202670_167_Designed.gc.bed",
                                   "/mnt/DATA6/mouseData/normals/20210518.n7.mouseGeneral.IDlist.txt")

mouseBeds$normalsId <- str_replace(mouseBeds$normalsId, "20210726_OCAP_gc_noheader.bed",
                                   "/mnt/DATA6/mouseData/normals/20210726_OCAP_normals.IDlist.txt")

mouseBeds$normalsId <- str_replace(mouseBeds$normalsId, "CCP.noTrack.GC.bed",
                                   "/mnt/DATA6/mouseData/normals/20211110.n10.CCP.IDlist.txt")

# mouseBeds$normalsId <- str_replace(mouseBeds$normalsId, "20210726_OCAP_gc_Andi2_noheader.bed",
#                                    "/mnt/DATA6/mouseData/normals/20220211ocap_WBC.IDlist.txt")

mouseBeds$normalsId <- str_replace(mouseBeds$normalsId, "20210726_OCAP_gc_Andi3_noheader.bed",
                                   "/mnt/DATA6/mouseData/normals/20220211ocap_WBC.IDlist.txt")

mouseBeds$normalsId <- str_replace(mouseBeds$normalsId, "IAD203665_173_Designed.gc.bed",
                                   "/mnt/DATA6/mouseData/normals/20220914pgu2.IDlist.txt")




mouseBeds$segParam <- mouseBeds$cnBed
mouseBeds$segParam <- str_replace(mouseBeds$segParam, "IAD124056_167_Designed.del.nopool.gc.bed", "2")
mouseBeds$segParam<- str_replace(mouseBeds$segParam, "IAD202670_167_Designed.gc.bed", "2")
mouseBeds$segParam <- str_replace(mouseBeds$segParam, "20210726_OCAP_gc_noheader.bed", "5")
mouseBeds$segParam <- str_replace(mouseBeds$segParam, "CCP.noTrack.GC.bed", "5")
# mouseBeds$segParam <- str_replace(mouseBeds$segParam, "20210726_OCAP_gc_Andi2_noheader.bed", "2")
mouseBeds$segParam <- str_replace(mouseBeds$segParam, "20210726_OCAP_gc_Andi3_noheader.bed", "2")
mouseBeds$segParam <- str_replace(mouseBeds$segParam, "IAD203665_173_Designed.gc.bed", "3")


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
  setwd("/mnt/DATA6/mouseData/tvcOut/")
  if(!dir.exists(mouseBeds$reports[i])){
    system(sprintf("mkdir -p %s", mouseBeds$reports[i]))
    setwd(mouseBeds$reports[i])
    system(sprintf("echo %s > directory.txt", mouseBeds$directories[i]))
  }
}


mouseBeds$reports[grep("PGU2", mouseBeds$reports)]
i <- 112

mouseBeds$cnBed[which(mouseBeds$cnBed == "IAD203665_173_Designed.gc.bed")] <- "IAD203665_173_Designed_fixed_KH.gc.bed"

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
    system(sprintf("echo %s > param.txt", mouseBeds$segParam[i]))
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

