#!/usr/bin/env Rscript

library(hashmap)
library(reshape2)

sampleVarCalls = read.table("/mnt/DATA4/kevhu/allDataTest.txt",
                            header = TRUE,sep = '\t', quote = "", stringsAsFactors = FALSE)


###remapping the names
namesRemap = read.table("/home/kevhu/data/Bed_name_map.July2017.csv", stringsAsFactors = FALSE,sep = ",",header = TRUE)
namesRemap = as.data.frame(namesRemap)
namesRemap$ORIG_BED[which(is.na(namesRemap$ORIG_BED))] = "NA"
namesRemap$FINAL_BED[which(is.na(namesRemap$FINAL_BED))] = "NA"


for(i in seq_along(namesRemap$ORIG_BED)){
  sampleVarCalls$BED <- gsub(pattern = namesRemap$ORIG_BED[i], replacement = namesRemap$FINAL_BED[i], x = sampleVarCalls$BED, fixed = TRUE)
}


###filter the NA's prior to counting
sampleVarCalls <- sampleVarCalls[-which(sampleVarCalls$BED == "NA"),]


###creating hashmaps
KeysTotPos <- unique(do.call(paste0,sampleVarCalls[c("Chr","Start_position","End_position")]))
KeysTotVar <- unique(do.call(paste0,sampleVarCalls[c("Chr","Start_position","End_position","REF_1","ALT_1")]))
KeysBedPos <- unique(do.call(paste0,sampleVarCalls[c("Chr","Start_position","End_position","BED")]))
KeysBedVar <- unique(do.call(paste0,sampleVarCalls[c("Chr","Start_position","End_position","REF_1","ALT_1","BED")]))


Values0 <- rep(0, nrow(sampleVarCalls)) 
totPosHash <- hashmap(KeysTotPos, Values0)
totVarHash <- hashmap(KeysTotVar, Values0)
bedPosHash <- hashmap(KeysBedPos, Values0)
bedVarHash <- hashmap(KeysBedVar, Values0)
#test <- apply(X = sampleVarCalls[,c("SAMPLE","BARCODE","BED","REPORTNUM","RunName")], MARGIN = 1, FUN = paste, collapse = "")
test <- apply(X = sampleVarCalls[,c("SAMPLE","BARCODE","BED","RunName")], MARGIN = 1, FUN = paste, collapse = "")
tot.ns <- length(unique(test))

###way to get counts for bed files
uniqBedNames <- unique(namesRemap$FINAL_BED)
uniqBedNames <- uniqBedNames[-which(uniqBedNames == "NA")]
tableBedCounts <- NULL
for(i in seq_along(uniqBedNames)){
  tableBedCounts <- c(tableBedCounts,length(grep(uniqBedNames[i], unique(test), fixed = TRUE)))
}
bedTable <- data.frame(uniqBedNames, tableBedCounts)
colnames(bedTable)[1] <- c("BED")


###using hashmap to get counts

lookupPos <- do.call(paste0,sampleVarCalls[c("Chr","Start_position","End_position")])
lookupVar <- do.call(paste0,sampleVarCalls[c("Chr","Start_position","End_position","REF_1","ALT_1")])
lookupBedPos <- do.call(paste0,sampleVarCalls[c("Chr","Start_position","End_position","BED")])
lookupBedVar <- do.call(paste0,sampleVarCalls[c("Chr","Start_position","End_position","REF_1","ALT_1","BED")])

for(i in seq_along(lookupPos)){
  if(totPosHash$has_key(lookupPos[i]) == TRUE){
    a <- totPosHash$find(lookupPos[i])
    totPosHash$insert(lookupPos[i], c(a + 1))
  }
}

for(i in seq_along(lookupVar)){
  if(totVarHash$has_key(lookupVar[i])){
    b <- totVarHash$find(lookupVar[i])
    totVarHash$insert(lookupVar[i], c(b + 1))
  }
}

for(i in seq_along(lookupBedPos)){
  if(bedPosHash$has_key(lookupBedPos[i]) == TRUE){
    c <- bedPosHash$find(lookupBedPos[i])
    bedPosHash$insert(lookupBedPos[i], c(c + 1))
  }
}

for(i in seq_along(lookupBedVar)){
  if(bedVarHash$has_key(lookupBedVar[i]) == TRUE){
    e <- bedVarHash$find(lookupBedVar[i])
    bedVarHash$insert(lookupBedVar[i], c(e + 1))
  }
}


###producing the final table
#test <- c(totPosHash$data())
finalTablePos <- data.frame(totPosHash$data())
finalTablePos$stringNamePos <- rownames(finalTablePos)
rownames(finalTablePos) <- NULL
colnames(finalTablePos)[1] <- c("new.tot.obs")
finalTablePos$new.tot.ns <- tot.ns
finalTablePos$new.tot.pos.pct <- round(finalTablePos$new.tot.obs/finalTablePos$new.tot.ns, digits = 4)
finalTableVar <- data.frame(totVarHash$data())
finalTableVar$stringNameVar <- rownames(finalTableVar)
rownames(finalTableVar) <- NULL
colnames(finalTableVar)[1] <- c("tot.var")
finalTableVar$new.tot.ns <- tot.ns
finalTableVar$new.tot.var.pct <- round(finalTableVar$tot.var/finalTableVar$new.tot.ns, digits = 4)
finalTableVar[,c("tot.var")] <- NULL
finalTableVar[,c("new.tot.ns")] <- NULL
finalTableBedPos <- data.frame(bedPosHash$data())
finalTableBedPos$stringNameBedPos <- rownames(finalTableBedPos) 
rownames(finalTableBedPos) <- NULL
colnames(finalTableBedPos)[1] <- c("new.bed.tot.obs")

finalTableBedVar <- data.frame(bedVarHash$data())
finalTableBedVar$stringNameBedVar <- rownames(finalTableBedVar)
rownames(finalTableBedVar) <- NULL
colnames(finalTableBedVar)[1] <- c("new.bed.var.obs")

sampleVarCalls$stringNamePos <- lookupPos
sampleVarCalls$stringNameVar <- lookupVar
sampleVarCalls$stringNameBedPos <- lookupBedPos
sampleVarCalls$stringNameBedVar <- lookupBedVar


###making final table


sampleVarCalls <- merge(sampleVarCalls, finalTablePos, by = c("stringNamePos"))
sampleVarCalls <- merge(sampleVarCalls, finalTableVar, by = c("stringNameVar"))
sampleVarCalls <- merge(sampleVarCalls, finalTableBedPos, by = c("stringNameBedPos"))
sampleVarCalls <- merge(sampleVarCalls, finalTableBedVar, by = c("stringNameBedVar"))
sampleVarCalls <- merge(sampleVarCalls, bedTable, by = c("BED"))

sampleVarCalls$new.bed.tot.pos.pct <- round(sampleVarCalls$new.bed.tot.obs/sampleVarCalls$tableBedCounts, digits = 4) 
sampleVarCalls$new.bed.tot.var.pct <- round(sampleVarCalls$new.bed.var.obs/sampleVarCalls$tableBedCounts, digits = 4)

sampleVarCalls <- sampleVarCalls[,-c(2:5)]


write.table(x = sampleVarCalls, file = "/home/kevhu/data/20170716_recalcAllData.tsv", sep = '\t', row.names = FALSE)


quit()
