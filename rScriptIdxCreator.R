#!/usr/bin/env Rscript
.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))
library(jsonlite)
library(optparse)
#normals <- read.table("/home/kevhu/data/normals/OCP1/normals.n2.20140811.txt", stringsAsFactors = FALSE)

currentDir <- getwd()
finalTable <- NULL;

option_list = list(
  make_option(c("-b", "--blacklist"), type="character", default=NULL, 
              help="blacklist of barcodes, no header", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="list of directories. in the form of /mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-149-AC03_AZ_CTC_1_229_406/plugin_out/coverageAnalysis_out.640", metavar="character"),
  make_option(c("-n", "--normals"), type="character", default=NULL, 
              help="path to list of normal file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="idx.txt", 
              help="output file name - default is idx.txt", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

if (is.null(opt$normals)){
  print_help(opt_parser)
  print("no normal file stated", call.=FALSE)
}

listOfDirectories <- unname(c(unlist(read.table(file = opt$input, stringsAsFactors = FALSE))))
#normals <- read.table(opt$normals, stringsAsFactors = FALSE, sep = "\t")
normals <- tryCatch(read.table(opt$normals, stringsAsFactors = FALSE, sep = "\t"),
                    error = function(x) return(NULL))
print(normals)
blacklist <- tryCatch(read.table(opt$blacklist, header = FALSE, stringsAsFactors = FALSE),
                      error = function(x) return(NULL))
#listOfDirectories <- c("/mnt/DATA3/Zeus_data_dump/ID7_reanalysis2_387/plugin_out/coverageAnalysis_out.604/")


finalTable <- c("SampleName","Barocde","path")
for(j in seq_along(listOfDirectories)){
  directory <- listOfDirectories[j]
  print(directory)
  setwd(directory)
  barcodes <- system('find . -type d -name "IonXpress*"', intern = TRUE)
  for(i in seq_along(barcodes)){
    barcodes[i] <- gsub("./","",barcodes[i])
  }
  jsonFile <- fromJSON("./results.json")
  
  
  
  ###making one large index file 
  barcodes <- sort(barcodes)
  for(i in seq_along(barcodes)){
    if(barcodes[i] %in% blacklist$V1){
      next()
    }
    
    ### added matchd barcode based on results.json misindexing what doesn't show on bc_summary
    ### example in /mnt/DATA3/Yoda_data_dump/ID_OCP_SATProton-139_355/plugin_out/coverageAnalysis_out.615
    
    matchBarcode <- which(names(jsonFile$barcodes[]) %in% barcodes[i])
    
    if(length(jsonFile$barcodes[[matchBarcode]]$`Uniformity of amplicon coverage`) == 0 && length(jsonFile$barcodes[[matchBarcode]]$`Uniformity of base coverage per target`) == 0){
      print(i)
      next()
    }

    sampleName <- jsonFile$barcodes[[matchBarcode]]$`Sample Name`
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
print(finalTable)
finalTable <- finalTable[-1,]
finalTable <- data.frame(finalTable, stringsAsFactors = FALSE)

if(!is.null(normals)){
  colnames(normals) <- colnames(finalTable)
  finalTable <- rbind(finalTable, normals)
}



finalTable[,1] <- apply(finalTable, 2,FUN = function(x) {sub("-","_", x)})[,1]

setwd(currentDir)
write.table(x = finalTable, file = opt$out,
            sep = '\t',row.names = FALSE, quote = FALSE, col.names = FALSE)


print("Done")
