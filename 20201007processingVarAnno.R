#!/usr/bin/env Rscript
.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))
library(optparse)
library(stringr)


### actually this script should read in all annoFiles and then output one large table per report

getSeqFields <- function(x){
  tmpFields <- NULL
  for(i in seq_along(x)){
    tmpVar <- unlist(str_split(x[i], "\\:"))
    fieldOneRow <- c(tmpVar[2], tmpVar[4], tmpVar[8], tmpVar[9], tmpVar[15], tmpVar[14])
    tmpFields <- rbind(tmpFields, fieldOneRow)
  }
  return(tmpFields)
}


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="directory where annotations of a specific report lie", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output filename prefix- default is idx.txt", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output prefix).n", call.=FALSE)
}

setwd(paste0(opt$input))
bedFile <- read.table("bedfile.txt" ,stringsAsFactors = FALSE)
bedName <- sub(x = bedFile$V1, pattern = "\\.gc\\.bed", replacement = "")

listOfAnnotationFiles <- system("ls *mm10_multianno.txt", intern = TRUE)
header <- c("Sample","Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene",
            "ExonicFunc.refGene", "AAChange.refGene", "mm10_mpgpv6_Indels",
            "GQ", "FDP", "FAO", "AF", "FSAF","FSAR","Bed.Name")
fullAnnoTable <- NULL
for(i in listOfAnnotationFiles){
  sampleName <- sub(x = i, pattern = "*\\.mm10_multianno\\.txt", replacement = "")
  filename <- paste0(opt$input,i)
  ### fails on empty files, so that's why tryCatch used
  tmpVcf <- tryCatch(read.table(filename, sep = "\t", skip = 1),
                     error=function(x) return(NULL))
  if(is.null(tmpVcf)){
    #print(i)
    next()
  }
  tmpAnno2 <- tmpVcf[,c(1:11, 21)]
  fourCol <- getSeqFields(tmpAnno2$V21)
  rownames(fourCol) <- NULL
  tmpAnno3 <- cbind(rep(sampleName, nrow(fourCol)) ,tmpAnno2[,1:11],
                    fourCol, rep(bedName, nrow(fourCol)))
  colnames(tmpAnno3) <- header
  fullAnnoTable <- rbind(fullAnnoTable, tmpAnno3)
}

setwd("/home/kevhu/scripts/newMousePanelPipeline/reportAnno/")
write.table(fullAnnoTable, paste0(opt$output, "_anno.txt"), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
### notes left behind
### 2:GQ, 4:FDP, 8:FAO, 9:AF
### the scirpt only writes a single reports annotations; just have another write the combined script

