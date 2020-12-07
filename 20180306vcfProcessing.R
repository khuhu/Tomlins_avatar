#can make this into a file later




#### 3 steps
###1 Create list of normal mutations that pass the thresholds - only need to do once, then create a table for it

###2 Create list of all calls w/ or w/o normals, doesn't matter as it gets filtered out later
###2.1 Can also count how often the variant is found throughout the samples

###3 filter out and then annotate



library(VariantAnnotation)

normals <- c("2405n.vcf.gz", "2519n.vcf.gz", "2796n.vcf.gz", "3867n.vcf.gz", "13604n.vcf.gz",
              "14104n.vcf.gz", "14154n.vcf.gz", "14433n.vcf.gz", "14085t.vcf.gz",
              "14104t.vcf.gz", "14154t.vcf.gz", "14286t.vcf.gz", "14399t.vcf.gz")
setwd("/mnt/DATA4/kevhu/choLab/vcfs/outLeftAligned/")

###for loop part under to get list for all normal vcfs

#vcf <- readVcf("/mnt/DATA4/kevhu/choLab/vcfs/outLeftAligned/2405n.vcf.gz", "mm10")

fullAlleleTable <- NULL

for(k in seq_along(normals)){
  path <- paste0("/mnt/DATA4/kevhu/choLab/vcfs/outLeftAligned/", normals[k])
  vcf <- readVcf(path, "mm10")
  alleFreqTab <- c("Chr","Pos", "Ref", "Alt","OREF", "OALT", "FAO","FDP","freq","FSAF.1","FSAR.1","ratio")
  for(i in 1:nrow(vcf@fixed)){
    #OPOS <- unlist(vcf@info["OPOS"][[1]][i])
    Pos <- as.character(vcf@rowRanges[i,1]@ranges@start)
    Chr <-  as.character(vcf@rowRanges[i,1]@seqnames)
    for(j in seq_along(OPOS)){
      FAO <- unlist(vcf@info["FAO"][[1]][i])
      FAO.1 <- FAO[j]
      FSAF <- unlist(vcf@info["FSAF"][[1]][i])
      FSAF.1 <- FSAF[j]
      FSAR <- unlist(vcf@info["FSAR"][[1]][i])
      FSAR.1 <- FSAR[j]
      ALT <- toString(unlist(vcf@fixed["ALT"][[1]][i])[[j]])
      REF <- toString(vcf@fixed["REF"][[1]][[i]])
      OALT <-  unlist(vcf@info["OALT"][[1]][i])
      OALT.1 <- OALT[j]
      OREF <- unlist(vcf@info["OREF"][[1]][i])
      OREF.1 <- OREF[j]
      #print(FAO.1)
      ratio <- FSAF.1/FSAR.1
      FDP <- unlist(vcf@info["FDP"][[1]][i])
      #print(FDP.1)
      freq <- FAO.1/FDP
      #print(freq)
      rowOfDat <- list(Chr,Pos, REF, ALT, OREF.1, OALT.1,FAO.1, FDP, freq, FSAF.1, FSAR.1,ratio)
      #print(rowOfDat)
      alleFreqTab <- rbind(alleFreqTab, rowOfDat)
    }
  }
  rownames(alleFreqTab) <- NULL
  colnames(alleFreqTab) <- c("Chr","Pos", "Ref", "Alt","OREF", "OALT", "FAO","FDP","freq","FSAF.1","FSAR.1","ratio")
  alleFreqTab <- alleFreqTab[-1,]
  alleFreqTab <- data.frame(alleFreqTab, stringsAsFactors = FALSE)
  alleFreqTab$ratio[which(alleFreqTab$ratio == "NaN")] <- 0
  ##taken from the UC paper 
  filter <- which(alleFreqTab$FAO < 6 | alleFreqTab$FDP < 20 | alleFreqTab$ratio > 5 | alleFreqTab$ratio < 0.2 | alleFreqTab$freq < 0.1)
  alleFreqTab <- alleFreqTab[-filter,]
  #print(alleFreqTab)
  fullAlleleTable <- rbind(fullAlleleTable, alleFreqTab)
  print(k)
}

fullAlleleTable <- data.frame(apply(fullAlleleTable, 2, unlist), stringsAsFactors = FALSE)
fullAlleleTable2 <- fullAlleleTable
fullAlleleTable2$tmp <- paste(fullAlleleTable2$Chr,fullAlleleTable2$Pos,fullAlleleTable2$Ref,fullAlleleTable2$Alt,sep = "/")
fullAlleleTable2 <- fullAlleleTable2[-c(which(duplicated(fullAlleleTable2$tmp))),]
fullAlleleTable2$tmp <- NULL
normalsList <- fullAlleleTable2


save(normalsList, file = "/mnt/DATA4/kevhu/choLab/vcfs/20180308mouseNormalsVars.Robj")
load("/mnt/DATA4/kevhu/choLab/vcfs/20180308mouseNormalsVars.Robj")

### Portion to grab and QC variants for rest of the mice samples

listOfFiles <- system("ls /mnt/DATA4/kevhu/choLab/vcfs/outLeftAligned/* | grep -v 'vcf.gz.tbi' | xargs -n1 basename", intern=TRUE)
setwd("/mnt/DATA4/kevhu/choLab/vcfs/outLeftAligned/")

fullAlleleTable <- NULL

for(k in seq_along(listOfFiles)){
  vcf <- readVcf(listOfFiles[k], "mm10")
  alleFreqTab <- c("SampName","Chr","Pos", "Ref", "Alt","OREF", "OALT", "FAO","FDP","freq","FSAF.1","FSAR.1","ratio")
  SampName <- listOfFiles[k]
  SampName <- gsub(".vcf.gz", "", SampName)
  for(i in 1:nrow(vcf@fixed)){
    #OPOS <- unlist(vcf@info["OPOS"][[1]][i])
    Pos <- as.character(vcf@rowRanges[i,1]@ranges@start)
    Chr <-  as.character(vcf@rowRanges[i,1]@seqnames)
    for(j in seq_along(OPOS)){
      FAO <- unlist(vcf@info["FAO"][[1]][i])
      FAO.1 <- FAO[j]
      FSAF <- unlist(vcf@info["FSAF"][[1]][i])
      FSAF.1 <- FSAF[j]
      FSAR <- unlist(vcf@info["FSAR"][[1]][i])
      FSAR.1 <- FSAR[j]
      ALT <- toString(unlist(vcf@fixed["ALT"][[1]][i])[[j]])
      REF <- toString(vcf@fixed["REF"][[1]][[i]])
      OALT <-  unlist(vcf@info["OALT"][[1]][i])
      OALT.1 <- OALT[j]
      OREF <- unlist(vcf@info["OREF"][[1]][i])
      OREF.1 <- OREF[j]
      #print(FAO.1)
      ratio <- FSAF.1/FSAR.1
      FDP <- unlist(vcf@info["FDP"][[1]][i])
      #print(FDP.1)
      freq <- FAO.1/FDP
      #print(freq)
      rowOfDat <- list(SampName, Chr,Pos, REF, ALT, OREF.1, OALT.1,FAO.1, FDP, freq, FSAF.1, FSAR.1,ratio)
      #print(rowOfDat)
      alleFreqTab <- rbind(alleFreqTab, rowOfDat)
    }
  }
  rownames(alleFreqTab) <- NULL
  colnames(alleFreqTab) <- c("SampName","Chr","Pos", "Ref", "Alt","OREF", "OALT", "FAO","FDP","freq","FSAF.1","FSAR.1","ratio")
  alleFreqTab <- alleFreqTab[-1,]
  alleFreqTab <- data.frame(alleFreqTab, stringsAsFactors = FALSE)
  alleFreqTab$ratio[which(alleFreqTab$ratio == "NaN")] <- 0
  ##taken from the UC paper 
  filter <- which(alleFreqTab$FAO < 6 | alleFreqTab$FDP < 20 | alleFreqTab$ratio > 5 | alleFreqTab$ratio < 0.2 | alleFreqTab$freq < 0.1)
  alleFreqTab <- alleFreqTab[-filter,]
  #print(alleFreqTab)
  fullAlleleTable <- rbind(fullAlleleTable, alleFreqTab)
  print(k)
}


allMouseVars <- fullAlleleTable
allMouseVars <- data.frame(apply(allMouseVars, 2, unlist), stringsAsFactors = FALSE)
allMouseVars$tmp <- paste(allMouseVars$Chr,allMouseVars$Pos,allMouseVars$Ref,allMouseVars$Alt,sep = "/")
allMouseVars <- allMouseVars[-c(which(duplicated(allMouseVars$tmp))),]
allMouseVars$tmp2 <- gsub("chr","",allMouseVars$Chr)
allMouseVars <- allMouseVars[order(as.numeric(allMouseVars$tmp2),allMouseVars$Pos),]
#allMouseVars$tmp <- NULL
allMouseVars$tmp3 <- paste(allMouseVars$tmp2, allMouseVars$Pos, allMouseVars$OREF, allMouseVars$OALT,sep = "/")

fullAlleleTable$tmp <- paste(fullAlleleTable$Chr,fullAlleleTable$Pos, fullAlleleTable$Ref,fullAlleleTable$Alt,sep = "/")

freqTable <- NULL
for(i in seq_along(allMouseVars$tmp)){
  numMatches <- length(which(fullAlleleTable$tmp == allMouseVars$tmp[i]))
  freqTable <- cbind(freqTable, numMatches)
  colnames(freqTable)[i] <- allMouseVars$tmp3[i]
}


####



aTable <- read.table("~/programs/annovar/allMouseVarsCombined.hotspotLeftAligned.vcf_filtAnno.exonic_variant_function", sep = "\t")
bTable <-  read.table("~/programs/annovar/nonHotspot.merged.vcf_filtAnno.variant_function", sep = "\t")


a.1 <- aTable[,2:8]
colnames(a.1) <- c("mutation type","transcript change", "chromosome","start","end", "ref","alt")
a.1$chromosome <- gsub("chr","", a.1$chromosome)
a.1$chromosome <- as.numeric(a.1$chromosome)
a.1 <- a.1[order(a.1$chromosome,a.1$start,a.1$`mutation type`),]

reduced <- a.1[-c(which(a.1$`mutation type` == "synonymous SNV")),]
reduced$tmp <- paste(reduced$chromosome,reduced$start,reduced$ref,reduced$alt,sep = "/")

###need to use OREF and OALT b/c these are right aligned, and annovar for some reason right aligns the changes
###NEED TO FIGURE OUT WHY THERE ARE SO LITTLE MATCHING POSITONS ... are they just all lowly covered positions and get filtered ?



fullAlleleTable2$tmp2 <- gsub("chr","",fullAlleleTable2$Chr)
fullAlleleTable2$tmp <- paste(fullAlleleTable2$tmp2, fullAlleleTable2$Pos,fullAlleleTable2$OREF, fullAlleleTable2$OALT,sep = "/")

reduced <- reduced[-which(reduced$tmp %in% fullAlleleTable2$tmp),]
which(colnames(freqTable) %in% reduced$tmp)


write.table(a.1,"/mnt/DATA4/kevhu/choLab/mouseComprehensiveVariationList.csv", col.names = TRUE,sep = ",", row.names = FALSE, quote = FALSE)
write.table(reduced,"/mnt/DATA4/kevhu/choLab/mouseReducedVariationList.csv", col.names = TRUE,sep = ",", row.names = FALSE, quote = FALSE)


a.2 <- a.1[,3:5]
for(i in 1:nrow(a.2)){
  a.2$chromosome[i] <- paste0("chr",a.2$chromosome[i]) 
}

reduced.2 <- reduced[,3:5]
for(i in 1:nrow(reduced.2)){
  reduced.2$chromosome[i] <- paste0("chr",reduced.2$chromosome[i])
}


fullAlleleTable2$tmp <- paste(fullAlleleTable2$Chr, fullAlleleTable2$Pos, fullAlleleTable2$OREF, fullAlleleTable2$OALT,sep = "/")
reduced$tmp <- paste(paste0("chr",reduced$chromosome),reduced$start, reduced$ref, reduced$alt, sep = "/")

reduced[-which(reduced$tmp %in% fullAlleleTable2$tmp),]
fullAlleleTable2[-which(fullAlleleTable2$tmp %in% reduced$tmp),]



write.table(a.2,"/mnt/DATA4/kevhu/choLab/mouseComprehensiveVariationList.bed", col.names = FALSE,sep = "\t", row.names = FALSE, quote = FALSE)
write.table(reduced.2,"/mnt/DATA4/kevhu/choLab/mouseReducedVariationList.bed", col.names = FALSE,sep = "\t", row.names = FALSE, quote = FALSE)
