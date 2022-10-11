#!/usr/bin/env Rscript
.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))
library(hashmap)
library(optparse)
library(stringr)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="list of directories", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name - default is out.txt", metavar="character"),
  make_option(c("-b", "--barcode"), type="logical", default=FALSE, 
              help="TRUE for IDs to be barcodes ", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

if (is.null(opt$out)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (out file).n", call.=FALSE)
}

setwd(dir = "/mnt/DATA4/kevhu/apps/fetchApp")




fullTable <- readRDS("data/data.RDS")
fullTable$Total.e2e <- as.numeric(fullTable$Total.e2e)

data <- read.table(file=opt$input, stringsAsFactors = FALSE, header = TRUE)



SampleName <- tryCatch(data$SampleName, error = function(x) return(NULL))
Barcode <- tryCatch(data$Barcode, error = function(x) return(NULL))
Bed <- tryCatch(data$Bed, error = function(x) return(NULL))
Report <- tryCatch(data$ReportID, error = function(x) return(NULL))
tableSamp <- fullTable$SampleName
tableBar <- fullTable$Barcode
tableBed <- fullTable$Bed
tableReport <- fullTable$ReportID
      
##improvement is to make function for below nulls
if(is.null(SampleName)){
  tableSamp <- NULL
}
if(is.null(Barcode)){
  tableBar <- NULL
}
if(is.null(Bed)){
  tableBed <- NULL
}
if(is.null(Report)){
  tableReport <- NULL
}
      
key <- paste0(tableSamp,tableBar, tableBed, tableReport)
keyVals <- rep(0, length(key))
hash <- hashmap(keys = key, values = keyVals)
lookup <- paste0(SampleName, Barcode, Bed, Report)

for(i in seq_along(lookup)){
if(hash$has_key(lookup[i]) == TRUE){
  counter <- hash$find(lookup[i])
  hash$insert(lookup[i],c(counter + 1))
  }
}
      
listOfNames<- names(which(c(hash$data()) > 0))
print("Look up done")
#print(fullTable[which(key %in% listOfNames),])


exportTable <- fullTable[which(key %in% listOfNames),]
exportTable$SampleName <- str_replace_all(string = exportTable$SampleName, pattern = " ", "_")

#print(nrow(exportTable)/length(unique(exportTable$AmpliconID)))

listOfBedsToSplit <- c(unique(exportTable$Bed))
listOfSplitTables <- NULL
      
      
for(i in seq_along(listOfBedsToSplit)){
  a <- NULL
  a <- exportTable[which(exportTable$Bed == listOfBedsToSplit[i]),]
  a$NameAndBar <- paste(a$SampleName, a$Barcode)
  assign(paste0("RnaCountData", listOfBedsToSplit[i]), a)
  listOfSplitTables <- c(listOfSplitTables, paste0("RnaCountData", listOfBedsToSplit[i]))
}

print(paste(listOfSplitTables,"are split into files by bed"))



listOfFinalExcelSheets <- NULL
for(i in seq_along(listOfSplitTables)){
  #print(paste0(i , " for making table"))
  NameAndBar <- paste(eval(as.name(listOfSplitTables[i]))$SampleName, eval(as.name(listOfSplitTables[i]))$Barcode)
  NameAndBar <- unique(NameAndBar)
  sampCount <- lapply(strsplit(NameAndBar, " "), `[[`, 1)
  sampBar <- lapply(strsplit(NameAndBar, " "), `[[`, 2)
  mappedReadsList <- unique(eval(as.name(listOfSplitTables[i]))$NumberOfMappedReads)
  #print(sampCount)
  print(sampBar)
  tmpTable <- NULL
  for(j in seq_along(sampCount)){
    #print(j)
    #print(dim(tmpTable))
    tmpTable2 <- NULL
    tmpTable2 <- eval(as.name(listOfSplitTables[i]))[which(eval(as.name(listOfSplitTables[i]))$NameAndBar == NameAndBar[j]),]
    tmpTable2 <- data.frame(tmpTable2[order(tmpTable2$AmpliconID),], stringsAsFactors = FALSE)
    if(j == 1){
      #print(paste(j, dim(tmpTable2)))
      #print("inital j works")
      tmpTable3 <- NULL
      tmpTable3 <- cbind(tmpTable3, tmpTable2[,c("AmpliconID")])
      tmpTable3 <- cbind(tmpTable3, tmpTable2[,c("Gene")])
      tmpTable3 <- cbind(tmpTable3, tmpTable2[,c("contig_id")])
      tmpTable3 <- data.frame(tmpTable3, stringsAsFactors = FALSE)
      tmpTable <- cbind(tmpTable, tmpTable2[,c("Total.e2e")])
      tmpTable <- data.frame(tmpTable, stringsAsFactors = FALSE)
      colnames(tmpTable3)[1] <- c("AmpliconID")
      colnames(tmpTable3)[2] <- c("Gene")
      colnames(tmpTable3)[3] <- c("Contig_ID")
      if(opt$barcode | sampCount[j] == "None"){
        colnames(tmpTable)[ncol(tmpTable)] <- sampBar[j]
      }
      else{
        colnames(tmpTable)[ncol(tmpTable)] <- sampCount[j]
      }
    }
    else{
      print(paste(j, dim(tmpTable2)))
      tmpTable <- cbind(tmpTable, tmpTable2[,c("Total.e2e")])
      if(opt$barcode | sampCount[j] == "None"){
        colnames(tmpTable)[ncol(tmpTable)] <- sampBar[j]
      }
      else{
        colnames(tmpTable)[ncol(tmpTable)] <- sampCount[j]
      }
    }
  }
  
  tmpTable <- cbind(tmpTable3, tmpTable)
  ###need to add one more NA when I rerun the table to include the contig ID's
  mappedReadsList <- c(NA,NA,NA, mappedReadsList)
  tmpTable <- rbind(mappedReadsList, tmpTable)
  tmpTable[,4:ncol(tmpTable)] <- sapply(tmpTable[,4:ncol(tmpTable)], as.numeric)
  assign(paste0("transformed",listOfBedsToSplit[i]), tmpTable)
  listOfFinalExcelSheets <- c(listOfFinalExcelSheets, paste0("transformed",listOfBedsToSplit[i]))
}

print(dim(eval(as.name(listOfFinalExcelSheets[i]))))

for(i in seq_along(listOfFinalExcelSheets)){
  write.table(x = eval(as.name(listOfFinalExcelSheets[i])), file = paste0("/mnt/DATA4/kevhu/linuxRNApull/",paste0(opt$out,i),".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  print(paste("Wrote txt table",i))
}



