bedFile <- read.table("/mnt/DATA4/kevhu/urineRNA/prTissue_WG00196_02092016_Designed.bed",skip = 1, stringsAsFactors = FALSE)


### Code below is used to make ampliseq Txn bed into proper 3 column format for tvc
### USELESS - not the file type I need .... sigh =/

finalTable <- NULL
for(i in seq_along(bedFile$V6)){
  ampName <- bedFile$V4[i]
  name <- bedFile$V1[i]
  aVector <- strsplit(bedFile$V6[i], split = ";")
  if(aVector[[1]][1] == "TYPE=GeneExpression"){
    chrom <- aVector[[1]][6]
    chrom <- str_replace(chrom, ".*=", "")
    dummyStart <- aVector[[1]][7]
    dummyStart <- str_replace(dummyStart, ".*=","")
    dummyStart <- unlist(strsplit(dummyStart, ","))
    dummyEnd <- aVector[[1]][8]
    dummyEnd <- str_replace(dummyEnd, ".*=","")
    dummyEnd <- unlist(strsplit(dummyEnd, ","))
  }
  if(aVector[[1]][1] == "TYPE=Fusion"){
    chrom <- c(aVector[[1]][11],aVector[[1]][14])
    chrom <- str_replace_all(chrom, ".*=", "")
    dummyStart <- c(aVector[[1]][12],aVector[[1]][15])
    dummyStart <- str_replace(dummyStart, ".*=","")
    dummyStart <- unlist(strsplit(dummyStart, ","))
    
    dummyFP <- c(aVector[[1]][12])
    dummyFP <- str_replace(dummyFP, ".*=","")
    dummyFP <- unlist(strsplit(dummyFP,","))
    dummyTP <- c(aVector[[1]][15])
    dummyTP <- str_replace(dummyTP, ".*=","")
    dummyTP <- unlist(strsplit(dummyTP,","))
    
    dummyEnd <- c(aVector[[1]][13],aVector[[1]][16])
    dummyEnd <- str_replace(dummyEnd, ".*=","")
    dummyEnd <- unlist(strsplit(dummyEnd, ","))
    
    if(length(dummyTP) > 1){
      nLength <- length(dummyTP) - 1
      chrom <- c(chrom[1],rep(chrom[2], nLength + 1))
    }
    if(length(dummyFP) > 1){
      nLength <- length(dummyFP) - 1
      chrom <- c(rep(chrom[1], nLength + 1),chrom[2])
    }
  }
  for(j in seq_along(dummyStart)){
      if(length(chrom) > 1){
        dummyVector <- c(chrom[j], dummyStart[j], dummyEnd[j], ampName, name, name) 
        finalTable <- rbind(finalTable, dummyVector)
      }
    else{
      dummyVector <- c(chrom, dummyStart[j], dummyEnd[j], ampName, name, name) 
      finalTable <- rbind(finalTable, dummyVector)
    }
  }
}

finalTable <- finalTable[order(finalTable[,1],finalTable[,2]),]

finalTable <- finalTable[-which(finalTable[,1] %in% c("chr6_ssto_hap7","chrHG104_HG975_PATCH","chrHSCHR6_MHC_MANN")),]

write.table(finalTable,"/mnt/DATA4/kevhu/urineRNA/prTissue_WG00196_02092016_Designed_6col_tvc.bed",sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)
