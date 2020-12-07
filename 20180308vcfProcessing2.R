library(VariantAnnotation)
listOfFiles <- system("ls /mnt/DATA4/kevhu/choLab/vcfs/outLeftAligned/*normed.filtered.vcf.gz | grep -v 'vcf.gz.tbi' | xargs -n1 basename", intern=TRUE)
setwd("/mnt/DATA4/kevhu/choLab/vcfs/outLeftAligned/")

fullAlleleTable <- NULL

for(k in seq_along(listOfFiles)){
  vcf <- readVcf(listOfFiles[k], "mm10")
  alleFreqTab <- c("SampName","Chr","Pos", "Ref", "Alt")
  SampName <- listOfFiles[k]
  SampName <- gsub(".normed.filtered.vcf.gz", "", SampName)
  if(nrow(vcf@info) == 0){
    next()
  }
  for(i in 1:nrow(vcf@fixed)){
    OPOS <- unlist(vcf@info["OPOS"][[1]][i])
    Pos <- as.character(vcf@rowRanges[i,1]@ranges@start)
    Chr <-  as.character(vcf@rowRanges[i,1]@seqnames)
    ALT <- toString(unlist(vcf@fixed["ALT"][[1]][i]))
    REF <- toString(vcf@fixed["REF"][[1]][[i]])
    #OALT <-  unlist(vcf@info["OALT"][[1]][i])
    #OREF <- unlist(vcf@info["OREF"][[1]][i])
    rowOfDat <- list(SampName, Chr,Pos, REF, ALT)
    alleFreqTab <- rbind(alleFreqTab, rowOfDat)
  }
  rownames(alleFreqTab) <- NULL
  colnames(alleFreqTab) <- c("SampName","Chr","Pos", "Ref", "Alt")
  alleFreqTab <- alleFreqTab[-1,]
  alleFreqTab <- data.frame(alleFreqTab, stringsAsFactors = FALSE)
  fullAlleleTable <- rbind(fullAlleleTable, alleFreqTab)
  print(k)
}


allMouseVars <- fullAlleleTable
allMouseVars <- data.frame(apply(allMouseVars, 2, unlist), stringsAsFactors = FALSE)
allMouseVars$tmp <- gsub("chr","",allMouseVars$Chr)
allMouseVars$tmp2 <- paste(allMouseVars$Chr, allMouseVars$Pos, allMouseVars$Ref, allMouseVars$Alt,sep = "/")
allMouseVars <- allMouseVars[order(as.numeric(allMouseVars$tmp),allMouseVars$Pos),]
#allMouseVars$tmp <- NULL


allTable <- read.table("~/programs/annovar/testingtesting.vcf.gz_filtAnno.exonic_variant_function", sep = "\t")
normTable <-  read.table("~/programs/annovar/testingtestingNormals.vcf.gz_filtAnno.exonic_variant_function", sep = "\t")

allTable$tmp <- paste(allTable$V9, allTable$V10,allTable$V12,allTable$V13,sep = "/")
normTable$tmp <- paste(normTable$V9, normTable$V10,normTable$V12,normTable$V13,sep = "/")

allTable <- allTable[-which(allTable$tmp %in% normTable$tmp),]


freqTable <- NULL
namesMatch <- NULL
for(i in seq_along(allTable$tmp)){
  numMatches <- length(which(allMouseVars$tmp2 == allTable$tmp[i]))
  namesMatch <- rbind(namesMatch, paste(allMouseVars$SampName[which(allMouseVars$tmp2 == allTable$tmp[i])], collapse = " "))
  freqTable <- cbind(freqTable, numMatches)
  colnames(freqTable)[i] <- allTable$tmp[i]
}

numRows <- nrow(allTable)

idx <- NULL
for(i in 18:ncol(allTable)){
  if(sum(grepl("\\./",allTable[,i])) == numRows){
    idx <- c(idx,i)
  }
}

finalTable <- allTable[,-idx]

for(i in 1:ncol(finalTable)){
  finalTable[,i] <- sub("\\./.*", "", finalTable[,i])
}

finalTable$tmp <- NULL
idx2 <- c(1,4:8,14:17)
finalTable <- finalTable[,-idx2]

for(i in 8:ncol(finalTable)){
  a <- grepl(":",finalTable[,i])
  b <- finalTable[,i][a]
  d <- which(a)
  for(j in seq_along(b)){
    c <- unlist(strsplit(b[j], ":"))[c(8,4)]
    finalTable[d[j],i] <- paste0(":",c[1],"/",c[2])
  }
}


numCols <- ncol(finalTable)

concatCol <- NULL
for(i in 1:nrow(finalTable)){
  dummyCol <- NULL
  for(j in 8:numCols){
    dummyCol <- c(dummyCol, finalTable[i,j])
    #print(dummyCol)
  }
  concatCol <- rbind(concatCol,paste(dummyCol,collapse = ""))
  #print(concatCol)
}

finalTable2 <- NULL
finalTable2 <- finalTable[,1:7]
finalTable2 <- cbind(finalTable2, concatCol,t(freqTable)[,1])

finalTable2 <- cbind(finalTable2, namesMatch)
#finalTable <- cbind(finalTable, t(freqTable)[,1])
colnames(finalTable2) <- c("mutation","specifics","chromosome","start","type","Ref","Alt","Counts","Recurrence", "Samples")


write.table(finalTable2,"/mnt/DATA4/kevhu/choLab/20180330mouseVariationList.txt", col.names = TRUE,sep = "\t", row.names = FALSE, quote = FALSE)
