#!/usr/bin/env Rscript
.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))
library(optparse)
library(stringr)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="list of directories, shoudl be full path with header 'Directories'", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default="NULL", 
              help="path to output dir ", metavar="character"),
  make_option(c("-b", "--blacklist"), type="character", default="NULL", 
              help="path to output dir ", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied directory", call.=FALSE)
}


### for the nature of snakemake - should make symbolic links to bams that aren't already there. if an old symbolic link is overwritten, the sample
### may be ran again - should make a table of names at the end of the script and read from the beginning


bamDir <- opt$dir
inputTab <- read.table(opt$input, header = TRUE, stringsAsFactors = FALSE)
directories  <- inputTab$Directories
blacklist <- read.table(opt$blacklist, header = FALSE, stringsAsFactors = FALSE)

snakeMakeFile <- NULL
snakeMakeFile <- "SampleName"
for(i in seq_along(directories)){
  print(i)
  setwd("/")
  setwd(directories[i])
  listOfBams <- system('find . *.bam -maxdepth 1 | grep IonXpress | grep -v "./"', intern = TRUE)
  listOfBais <- system('find . *.bam.bai -maxdepth 1 | grep IonXpress | grep -v "./"', intern = TRUE)
  #listOfVcfs <- listOfVcfs[order(listOfVcfs)]
  locOfNames <- system('find . -name *bc_summary.xls*  | grep -v "link"',
                       intern = TRUE)
  tableOfNames <- read.table(file = locOfNames, stringsAsFactors = FALSE, header = TRUE,sep = "\t")
  tableOfNames.justNames <- NULL
  for(j in seq_along(tableOfNames$Sample.Name)){
    tableOfNames.justNames <- c(tableOfNames.justNames, tableOfNames$Sample.Name)
  }
  tableOfNames.justNames <- sapply(tableOfNames.justNames, FUN = function(x) str_replace_all(x," ",""))
  for(k in seq_along(listOfBams)){
    if(tableOfNames.justNames[k] %in% blacklist$V1){
      next()
    }
    snakeMakeFile <- rbind(snakeMakeFile, tableOfNames.justNames[k])
    system(paste0("ln -s  ", paste0(directories[i], sub("./","",listOfBams[k])), " ", bamDir,tableOfNames.justNames[k],".bam"))
    system(paste0("ln -s  ", paste0(directories[i], sub("./","",listOfBais[k])), " ", bamDir,tableOfNames.justNames[k],".bam.bai"))
  }
  #setwd("/")
}
write.table(snakeMakeFile, "/mnt/DATA4/kevhu/urineRNA/urineSampNames.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")


