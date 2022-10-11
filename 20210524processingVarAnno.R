library(optparse)
library(stringr)

### add a list here so that I can use bcftoolsc merge for genotype pipeline 
### this script should read in all annoFiles and then output one large table per report

getSeqFields <- function(x,y){
  tmpFields <- NULL
  for(i in seq_along(x)){
    tmpVar <- unlist(str_split(x[i], "\\:"))
    tmpVar2 <- unlist(str_split(y[i], ";"))
    tmpVar2 <- tmpVar2[grep("QD|HRUN", tmpVar2)]
    otherFields <- as.numeric(str_remove(tmpVar2, ".*="))
    
    fieldOneRow <- c(tmpVar[2], tmpVar[4], tmpVar[8], tmpVar[9], tmpVar[15], tmpVar[14],
                     otherFields[1], otherFields[2])
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
### setwd("/mnt/DATA6/mouseData/vcfs/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349")


bedFile <- read.table("bedfile.txt" ,stringsAsFactors = FALSE)
bedName <- sub(x = bedFile$V1, pattern = "\\.gc\\.bed", replacement = "")

listOfAnnotationFiles <- system("ls *mm10_multianno.txt", intern = TRUE)
header <- c("Sample","Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene",
            "ExonicFunc.refGene", "AAChange.refGene", "mm10_mpgpv6_Indels","name","QUAL", 
            "GQ", "FDP", "FAO", "AF", "FSAF","FSAR", "HRUN", "QD", "Bed.Name")
fullAnnoTable <- NULL
for(i in listOfAnnotationFiles){
  sampleName <- sub(x = i, pattern = "*\\.mm10_multianno\\.txt", replacement = "")
  filename <- paste0(opt$input,"/",i)
  
  ### test
  # filename <- paste0("/mnt/DATA6/mouseData/vcfs/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349","/",i)
  ### fails on empty files, so that's why tryCatch used
  tmpVcf <- tryCatch(read.table(filename, sep = "\t", skip = 1),
                     error=function(x) return(NULL))
  if(is.null(tmpVcf)){
    print(i)
    next()
  }
  tmpAnno2 <- tmpVcf[,c(1:11,14, 17, 19, 21)]
  fourCol <- getSeqFields(tmpAnno2$V21, tmpVcf$V19)
  rownames(fourCol) <- NULL
  tmpAnno3 <- cbind(rep(sampleName, nrow(fourCol)) ,tmpAnno2[,1:13],
                    fourCol, rep(bedName, nrow(fourCol)))
  colnames(tmpAnno3) <- header
  fullAnnoTable <- rbind(fullAnnoTable, tmpAnno3)
}

vcfList <- system("ls -d $PWD/*norm.vcf.gz* | grep -v 'tbi'", intern = TRUE)
writeLines(vcfList, "listOfVcfs.txt")

setwd("/mnt/DATA6/mouseData/reportAnno/")
write.table(fullAnnoTable, paste0(opt$output, "_anno.txt"), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
### notes left behind
### 2:GQ, 4:FDP, 8:FAO, 9:AF
### the scirpt only writes a single reports annotations; just have another write the combined script

