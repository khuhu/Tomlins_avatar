oldZues <- read.table("/home/kevhu/data/20170728ZuesOld.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
newZues <- read.table("/home/kevhu/data/20170728ZuesNew.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
gandalf <- read.table("/home/kevhu/data/20170727GandalfTable.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
oldYoda <- read.table("/home/kevhu/data/20170728OldYoda.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
newYoda <- read.table("/home/kevhu/data/20170728NewYoda.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)

###doing some quick cleanup of data so all columns have matching formats

###I mislabled the columsn: listOFJSONkey is suppose to be uniformity and vice versa
which(colnames(gandalf) == "Uniformity")
gandalf <- gandalf[,-6]
colnames(gandalf)[5] <- c("Uniformity")

gandalf$Reportnum <- gsub("ResultsPK = ", "", gandalf$Reportnum)


consoliatedDf <- rbind(gandalf, newYoda, oldYoda, newZues, oldZues)

###formatting the df so numerics are numerics and etcs


consoliatedDf$Uniformity <- gsub("%","", consoliatedDf$Uniformity)
consoliatedDf$percentMapped <- gsub("%","", consoliatedDf$percentMapped)
colnames(consoliatedDf[,4:8])
consoliatedDf[,4:8] <- as.numeric(unlist(consoliatedDf[,4:8]))
###seems percent mapped is only portion with NAs
consoliatedDf[which(consoliatedDf == "None"),]


###remapping the names - need to create new name remapping table as the other one is not perfect
#namesRemap = read.table("/home/kevhu/data/Bed_name_map.July2017.csv", stringsAsFactors = FALSE,sep = ",",header = TRUE)
#namesRemap = as.data.frame(namesRemap)
#namesRemap$ORIG_BED[which(is.na(namesRemap$ORIG_BED))] = "None"
#namesRemap$FINAL_BED[which(is.na(namesRemap$FINAL_BED))] = "None"
#
#for(i in seq_along(namesRemap$ORIG_BED)){
#  consoliatedDf$Bed <- gsub(pattern = namesRemap$ORIG_BED[i], replacement = namesRemap$FINAL_BED[i], x = consoliatedDf$Bed, fixed = TRUE)
#}

listOFBeds <- unique(consoliatedDf$Bed)
listOFreplacements <- c("02867-1358092141_Covered","IAD34847_Designed","OCP3_20140506_designed","OCP_20130724_designed",
                        "OCP3_20140506_designed","02867-1358092141_Covered","IAD46903_31_Designed","CHP2_designed_20120806",
                        "IAD41516_47_Designed","CCP","IAD79597_173_Designed","CCP","KH_IAD41516_47",
                        "SAT_Prostate_2015_1_Covered","CHP2_designed_20120806","None","prTissue_WG00196_02092016_Designed",
                        "OCP_20150630_designed",
                        "TargetSeq-Exome-50Mb-hg19_revA","IAD79597_173_Designed","OCP_20130724_designed","CHP2_designed_20120806",
                        "CCP","OCP_20150630_designed","SAT_prostate_2014_1_Covered","OCP3_20140506_designed",
                        "TargetSeq-Exome-50Mb-hg19_revA","OCP_20130724_designed","TargetSeq-Exome-50Mb-hg19_revA","IAD56687_4_Designed","IAD101642_4_Designed",
                        "prUrine_WG00196.2_02102016_Designed","OCP3_20140506_designed","OCP3_20140506_designed","OCP_20130724_designed",
                        "CCP","hg19_AmpliSeq_Transcriptome_21K_v1")

#listOFBeds <- sapply(listOFBeds, function(x) paste(".*",x,".*"))

namesRemapConsol <- data.frame(listOFBeds, listOFreplacements)


for(i in seq_along(namesRemapConsol$listOFBeds)){
  consoliatedDf$Bed <- gsub(pattern = namesRemapConsol$listOFBeds[i], replacement = namesRemapConsol$listOFreplacements[i], x = consoliatedDf$Bed)
}

unique(consoliatedDf$Bed)
###seing if I can use new consolodated tables to prioritize variants

listMatch <- do.call(paste0, consoliatedDf[,c("Sample","Barcode","Bed","Reportnum")])
listMatch2 <- do.call(paste0, consoliatedDf[,c("Barcode","Bed","Reportnum")])
load("/home/kevhu/data/listOfDuplicatedSampleRuns.Robj")

tableOfMatches2.1 <- gsub(",","",unlist(tableOfMatches2))
tableOfMatches2.1 <- gsub(" ","", tableOfMatches2.1,fixed = TRUE)

#better way to fix problem
tableOfMatches2 <- lapply(tableOfMatches2,function(x) gsub(" ","",x))
tableOfMatches2 <- lapply(tableOfMatches2,function(x) gsub(",","",x))

#only 618 of the 1789 are mathced unfortunately ... and there are 803 uniqe pairs ..... so not even enough data for 1 to 1
#after updated list 931 of the 1789 are found
length(which(tableOfMatches2.1 %in% listMatch))


###finding which ones aren't in the list
tableOfMatches2.1[which(!tableOfMatches2.1 %in% listMatch)]
notMatched <- tableOfMatches2.1[which(!tableOfMatches2.1 %in% listMatch)]


#out of every unique variant: 2190/4663 are found from the list I produced
#after updated list 3831/4663 are found
length(which(listMatch %in% testNames2))

#checking wihtout sample names what i get: still not much better as I only get a bit more 2327/4595
#update 4158/4595 are found 
length(which(listMatch2 %in% testNames))

testList <- NULL
for(i in seq_along(tableOfMatches2)){
  a <- NULL
  for(j in seq_along(tableOfMatches2[[i]])){
    a <- grep(tableOfMatches2[[i]][j],listMatch)
  }
  testList[i] <- list(a)
}

length(which(lengths(testList) == 0))
###now this is where the problem stems from only 419 of the 803 duplicates have matches... fuuu

###creating a cross-tabulation table

library(hashmap)
key3Vars <- unique(do.call(paste0, consoliatedDf[,c("Sample","Barcode","Bed")]))
key2Vars <- unique(do.call(paste0, consoliatedDf[,c("Sample","Bed")]))

Values3 <- rep(0, length(key3Vars)) 
Values2 <- rep(0, length(key2Vars)) 

hashtable3 <- hashmap(keys = key3Vars, values = Values3)
hashtable2 <- hashmap(keys = key2Vars, values = Values2)

lookup3 <- do.call(paste0, consoliatedDf[,c("Sample","Barcode","Bed")])
lookup2 <- do.call(paste0, consoliatedDf[,c("Sample","Bed")])

for(i in seq_along(lookup3)){
  if(hashtable3$has_key(lookup3[i]) == TRUE){
    a <- hashtable3$find(lookup3[i])
    hashtable3$insert(lookup3[i], c(a + 1))
  }
}

for(i in seq_along(lookup2)){
  if(hashtable2$has_key(lookup2[i]) == TRUE){
    a <- hashtable2$find(lookup2[i])
    hashtable2$insert(lookup2[i], c(a + 1))
  }
}

final3 <- data.frame(hashtable3$data())
final3$stringName <- rownames(final3)
rownames(final3) <- NULL
colnames(final3)[1] <- c("counts")

final2   <- data.frame(hashtable2$data())
final2$stringName <- rownames(final2)
rownames(final2) <- NULL
colnames(final2)[1] <- c("counts")


