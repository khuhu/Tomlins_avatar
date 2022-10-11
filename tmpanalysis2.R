gandalf <- read.table("/home/kevhu/data/20170727GandalfTable.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
zues <- read.table("/home/kevhu/data/20170801NonIonDbZues.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
yoda <- read.table("/home/kevhu/data/20170801NonIonDbYoda.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)

LU <- read.csv("/home//kevhu/data/LU.csv")
PR <- read.csv("/home//kevhu/data/PR.csv")
MO <- read.csv("/home//kevhu/data/MO1.csv")

LU <- LU[,c("Sample","Mapped.Reads","Mapped.reads.on.target","Base.coverage..uniformity.","Base.coverage..mean.")]
PR <- PR[,c("Sample","Mapped.Reads","Mapped.reads.on.target","Base.coverage..uniformity.","Base.coverage..mean.")]
MO <- MO[,c("Sample.ID","Mapped.Reads","Mapped.reads.on.target","Base.coverage..uniformity.","Base.coverage..mean.")]
colnames(MO)[1] <- c("Sample")

barcodes1 <- read.csv("/home/kevhu/data/compendiaVCFnames1.csv")
barcodes2 <- read.csv("/home/kevhu/data/compendiaVCFnames2.csv")

barcodes1 <- barcodes1[,c("Sample","IonXpress.Barcode")]
barcodes2 <- barcodes2[,c("Sample","IonXpress.Barcode","PrimerPool")]

#barcodes1$tmp <- paste0(barcodes1$Sample, barcodes1$IonXpress.Barcode)
#barcodes2$tmp <- paste0(barcodes2$Sample, barcodes2$IonXpress.Barcode)

#everyhing in barcodes 1 is in barcodes 
#length(which(barcodes2$tmp %in% barcodes1$tmp))

#barcodes1$tmp <- NULL
#barcodes2$tmp <- NULL

library(reshape2)
test <- NULL
test <- rbind(LU,MO,PR)
test$Sample <- gsub("[[:punct:]]", "", test$Sample)
test$Sample <- sapply(test$Sample, toupper)
testMerge <- merge(test, barcodes2, by = "Sample")

colnames(testMerge) <- c("Sample","numberMappedReads","percentMapped","Uniformity","avgCov",
                         "Barcode","Bed")


numbers_only <- function(x) !grepl("\\D", x)

for(i in which(numbers_only(testMerge$Barcode))){
  if(nchar(testMerge$Barcode[i]) == 1){
    testMerge$Barcode[i] <- paste0("IonXpress_00",testMerge$Barcode[i])
  }
  if(nchar(testMerge$Barcode[i]) == 2){
    testMerge$Barcode[i] <- paste0("IonXpress_0",testMerge$Barcode[i])
  }
}


testMerge$Bed <- as.character(testMerge$Bed)
testMerge$Bed[which(testMerge$Bed == "OCP3")] <- "OCP3_20140506_designed"
testMerge$Bed[which(testMerge$Bed == "OCP2")] <- "OCP3_20140506_designed"


upgradedDups <- system('sudo -kS ls /mnt/DATA3/zeus_temp_2017/',
                       input="sat1840", intern = TRUE)

newExclusion <- NULL
for(i in seq_along(upgradedDups)){
  newExclusion <- c(newExclusion, which(subSample$REPORT == upgradedDups[i]))
}

subSample <- subSample[-newExclusion,]

subSample$tmp <- paste0(subSample$SAMPLE,subSample$BARCODE,subSample$BED)
subSample$tmp2 <- paste0(subSample$SAMPLE,subSample$BARCODE)
testMerge$tmp <- paste0(testMerge$Sample, testMerge$Barcode, testMerge$Bed)
testMerge$tmp2 <- paste0(testMerge$Sample, testMerge$Barcode)
testMerge$Reportnum <- NULL
testMerge$Report <- NULL
for(i in seq_along(testMerge$tmp)){
  dummyVar <- unique(subSample[which(subSample$tmp == testMerge$tmp[i]),c("REPORT")])
  if(length(dummyVar) > 0){
    testMerge$Report[i] <- unique(subSample[which(subSample$tmp == testMerge$tmp[i]),c("REPORT")])
  }
  if(length(dummyVar) == 0){
    dummyVar <- unique(subSample[which(subSample$tmp2 == testMerge$tmp2[i]),c("REPORT")])
    print(dummyVar)
    if(length(dummyVar) > 0){
    testMerge$Report[i] <- unique(subSample[which(subSample$tmp2 == testMerge$tmp2[i]),c("REPORT")])
    }
    if(length(dummyVar) == 0){
    testMerge$Report[i] <- "not found"
    print(i)
    }
  }
}

length(which(testMerge$Reportnum == "not found"))
length(which(testMerge$Report == "not found"))

###checking to see which doubles have a match and which don't


for(i in (which(duplicated(testMerge$Sample)))){
  prevIndex <- i - 1
  bool1 <- testMerge$Report[i] == "not found"
  bool2 <- testMerge$Report[prevIndex] == "not found"
  if(bool1 != bool2){
    print(paste("Index",prevIndex,"and",i,"have a uniqe reportnum"))
  }
  if(bool1 == bool2){
    print(paste("Index",prevIndex,"and",i,"both have or both missing"))
  }
}


#since from above code, i found out i was only missing 2 sets I can just set against database
#LU10 bc 18
testMerge$Reportnum[which(testMerge$Sample == "LU10" & testMerge$Barcode == "18")] <- 57
testMerge$Reportnum[which(testMerge$Sample == "LU13" & testMerge$Barcode == "71")] <- 51

testMerge$Report[which(testMerge$Sample == "LU101" & testMerge$Barcode == "IonXpress_023")] <- "OCPv3"
testMerge$Report[which(testMerge$Sample == "LU92" & testMerge$Barcode == "IonXpress_022")] <- "OCPv3"

testMerge <- testMerge[-c(6,105),]

#these are non-dups that are still missing values - reasonI removed them is probably b/c they had no beds maybe
testMerge$Reportnum[which(testMerge$Sample == "MO55" & testMerge$Barcode == "54")] <- 55
testMerge$Reportnum[which(testMerge$Sample == "MO57" & testMerge$Barcode == "42")] <- 57

#testMerge$Reportnum[which(testMerge$Sample == "MO87" & testMerge$Barcode == "8")] <- 57
#dropped b/c not found in database -> probably b/c it has bad cov and unifomrity too

#testMerge$Reportnum[which(testMerge$Sample == "PR118" & testMerge$Barcode == "28")] <- 57
#dropped b/c it has no values from table for anything



length(which(testMerge$Reportnum == "not found"))

##remove the tmp ones -> most of them are there with a double 
#testMerge <- testMerge[-which(testMerge$Reportnum == "not found"),]
testMerge <- testMerge[-which(testMerge$Report == "not found"),]
testMerge <- data.frame(testMerge,stringsAsFactors = FALSE)
testMerge$tmp <- NULL
testMerge$tmp2 <- NULL
write.table(testMerge, file = "/home/kevhu/data/compendiaVCFfinalTable.tsv", sep = '\t', row.names = FALSE)



unique(testMerge$Bed)

orderOfNonYoda <- colnames(zues)
yoda <- yoda[orderOfNonYoda]
newdf <- rbind(gandalf,zues,yoda)
newdf$Uniformity <- NULL
newdf$Reportnum <- NULL
colnames(newdf)[5] <- c("Uniformity")
#testing the report extraction
unlist(strsplit(newdf$Report[1],"/"))[2]
for(i in seq_along(newdf$Report)){
  newdf$Report[i] <- unlist(strsplit(newdf$Report[i],"/"))[2]
}

orderOfColnames <- colnames(newdf)
testMerge <- testMerge[orderOfColnames]

newdf <- rbind(newdf, testMerge)
newdf$Uniformity <- gsub("%","", newdf$Uniformity)
newdf$percentMapped <- gsub("%","", newdf$percentMapped)
#newdf$Reportnum <- gsub("ResultsPK = ","", newdf$Reportnum)
colnames(newdf[,5:8])
newdf[,5:8] <- as.numeric(unlist(newdf[,5:8]))

listOFBeds <- unique(newdf$Bed)

unique(newdf$Bed)

listOFreplacements <- c("02867-1358092141_Covered","IAD34847_Designed","OCP3_20140506_designed","OCP_20130724_designed",
                        "OCP3_20140506_designed","02867-1358092141_Covered","IAD46903_31_Designed","CHP2_designed_20120806",
                        "IAD41516_47_Designed","CCP","IAD79597_173_Designed","CCP","KH_IAD41516_47",
                        "SAT_Prostate_2015_1_Covered","CHP2_designed_20120806","None","OCP_20150630_designed", "CCP",
                        "OCP3_20140506_designed","OCP_20130724_designed","SAT_prostate_2014_1_Covered","TargetSeq-Exome-50Mb-hg19_revA",
                        "TargetSeq-Exome-50Mb-hg19_revA","OCP_20130724_designed","CHP2_designed_20120806","OCP3_20140506_designed")


namesRemapConsol <- data.frame(listOFBeds, listOFreplacements)
namesRemapConsol2425 <- namesRemapConsol[24:25,]
namesRemapConsol <- namesRemapConsol[-c(24,25),]

for(i in seq_along(namesRemapConsol$listOFBeds)){
  newdf$Bed <- gsub(pattern = namesRemapConsol$listOFBeds[i],
                    replacement = namesRemapConsol$listOFreplacements[i], x = newdf$Bed)
}

newdf$Bed[which(newdf$Bed == "OCP_20130724_designed.GC.geneName.Hayes_trimmed")] <- "OCP_20130724_designed"
newdf$Bed[which(newdf$Bed == "CHP2_designed_20120806.noTrack.GC.Hayes_trimmed")] <- "CHP2_designed_20120806"


newdf$Sample <- gsub("[[:punct:]]", "", newdf$Sample)
newdf$Sample <- sapply(newdf$Sample, toupper)
newdf$Sample <- gsub("DNA","", newdf$Sample)


simpasSamples <- unique(subSample$SAMPLE[grep("SPR", subSample$SAMPLE)])
simpasSamples.1 <- gsub("S","",simpasSamples)

simpasSamples2 <- cbind(simpasSamples.1,simpasSamples)
simpasSamples2 <- as.data.frame(simpasSamples2, stringsAsFactors = FALSE)
colnames(simpasSamples2) <- c("orig","new")

for(i in seq_along(simpasSamples2$orig)){
  newdf$Sample <- gsub(simpasSamples2$orig[i],simpasSamples2$new[i], newdf$Sample)
}


unique(newdf$Bed)
#newdf <- newdf[-which(newdf$Bed == "None"),]

numbers_only <- function(x) !grepl("\\D", x)

###changing all of the regular numbers into ion_express barcode
for(i in which(numbers_only(newdf$Barcode))){
  if(nchar(newdf$Barcode[i]) == 1){
    newdf$Barcode[i] <- paste0("IonXpress_00",newdf$Barcode[i])
  }
  if(nchar(newdf$Barcode[i]) == 2){
    newdf$Barcode[i] <- paste0("IonXpress_0",newdf$Barcode[i])
  }
}
unique(newdf$Barcode)

write.table(newdf, "/home/kevhu/data/20170808consolidatedTable.tsv", sep = '\t',row.names = FALSE)


listMatch <- unique(apply(X = newdf[,c("Sample","Barcode","Bed","Reportnum")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))
listMatch2 <- unique(apply(X = newdf[,c("Barcode","Bed","Reportnum")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))

listMatch3 <- unique(apply(X = newdf[,c("Sample","Bed")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))
listMatch4 <- unique(apply(X = newdf[,c("Sample","Barcode","Bed")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))

listMatch5 <- unique(apply(X = newdf[,c("Sample","Barcode","Bed","Report")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))
newdf$tmp <- apply(X = newdf[,c("Sample","Barcode","Bed","Report")],MARGIN = 1, FUN = function(x) paste(x, collapse = ","))


total3 <- apply(X = newdf[,c("Sample","Bed")],MARGIN = 1, FUN = function(x) paste(x, collapse = ","))
total4 <- apply(X = newdf[,c("Sample","Barcode","Bed")],MARGIN = 1, FUN = function(x) paste(x, collapse = ","))


length(which(listMatch2 %in% testNames2))
2861/5644
length(which(testNames2 %in% listMatch2))
2865/4574

length(which(listMatch %in% testNames1))
2814/4627
length(which(testNames1 %in% listMatch))
2814/5672


length(which(listMatch3 %in% testNames3))
3351/3868
length(which(testNames3 %in% listMatch3))
3351/3778

length(which(listMatch4 %in% testNames4))
3409/4278
length(which(testNames4 %in% listMatch4))
3409/3859

length(which(listMatch5 %in% testNames5))
3583/5670
length(which(testNames5 %in% listMatch5))
3583/(4589-100)


###another way to match inclusiveness: first I will match Sample, barcode, bed and report - if nothing is found then search Sample barcdoe and report

newdf.cc <- newdf
newdf.cc$tmp <- NULL
newdf.cc$tmp <- apply(X = newdf.cc[,c("Sample","Barcode","Bed","Report")],MARGIN = 1, FUN = function(x) paste(x, collapse = ","))
newdf.cc <- newdf.cc[!duplicated(newdf.cc$tmp),]
newdf.cc$tmp2 <- apply(X = newdf.cc[,c("Sample","Barcode","Report")],MARGIN = 1, FUN = function(x) paste(x, collapse = ","))


###what portion is/isn't found from mySQL to uniformity
in1 <- which(subSample.cc$tmp %in% newdf.cc$tmp)
in1.1 <- which(! subSample.cc$tmp %in% newdf.cc$tmp)
in2 <- which(subSample.cc$tmp2[in1.1] %in% newdf.cc$tmp2)
in3 <- c(in1, in2)
length(in3)/4589
out1 <- which(! subSample.cc$tmp2[in1.1] %in% newdf.cc$tmp2)

subSample.In <- subSample.cc[in3,]
subSample.Out <- subSample.cc[out1,]

subSample.In$tmp <- NULL
subSample.In$tmp2 <- NULL
subSample.Out$tmp <- NULL
subSample.Out$tmp2 <- NULL

write.table(subSample.In, "/home/kevhu/data/20170811mySQLmatched.tsv", sep = '\t',row.names = FALSE)
write.table(subSample.Out, "/home/kevhu/data/20170811mySQLnotMatched.tsv", sep = '\t',row.names = FALSE)


uniformIn <- which(newdf.cc$tmp %in% subSample.cc$tmp)
uniformIn.1 <- which(! newdf.cc$tmp %in% subSample.cc$tmp)
uniformIn2 <- which(newdf.cc$tmp2[uniformIn.1] %in% subSample.cc$tmp2)
uniformIn3 <- c(uniformIn, uniformIn2)
length(uniformIn3)/5670
uniformOut1 <- which(! newdf.cc$tmp2[uniformIn.1] %in% subSample.cc$tmp2)

newdf.In <- newdf[uniformIn3,]
newdf.Out <- newdf[uniformOut1,]

newdf.In$tmp <- NULL
newdf.Out$tmp <- NULL

write.table(newdf.In, "/home/kevhu/data/20170811consolTableMatched.tsv", sep = '\t',row.names = FALSE)
write.table(newdf.Out, "/home/kevhu/data/20170811consolTableNotMatched.tsv", sep = '\t',row.names = FALSE)


#testNotIn <- testNames5[which(!testNames5 %in% newdf$tmp)]
#testNotIn.1 <- NULL
#testIn <- testNames5[which(testNames5 %in% newdf$tmp)]
#testIn.1 <- NULL
#for(i in seq_along(testNotIn)){
#  testNotIn.1 <- rbind(testNotIn.1, unlist(strsplit(testNotIn[i],",")))
#}
#colnames(testNotIn.1) <- c("Sample","Barcode","Bed","Report")
#testNotIn.1 <- as.data.frame(testNotIn.1, stringsAsFactors = FALSE)
#
#for(i in seq_along(testIn)){
#  testIn.1 <- rbind(testIn.1, unlist(strsplit(testIn[i],",")))
#}
#colnames(testIn.1) <- c("Sample","Barcode","Bed","Report")
#testIn.1 <- as.data.frame(testIn.1, stringsAsFactors = FALSE)


#write.table(newdf.subsetIn, "/home/kevhu/data/20170810matchedIn.tsv", sep = '\t',row.names = FALSE)
#write.table(newdf.subseNotIn, "/home/kevhu/data/20170810matchedNotIn.tsv", sep = '\t',row.names = FALSE)



#list of things to exclude
upgradedDups <- system('sudo -kS ls /mnt/DATA3/zeus_temp_2017/',
       input="sat1840", intern = TRUE)


###finding the new matches

subSampleTotal.report <- unique(apply(X = subSample[,c("SAMPLE","BARCODE","BED","REPORT")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))
subSampleTotal.3ident <-  unique(apply(X = subSample[,c("SAMPLE","BARCODE","BED")],MARGIN = 1, FUN = function(x) paste(x, collapse = ",")))
#uniq var is the same as testNames5
#subSample$tmp <- apply(X = subSample[,c("SAMPLE","BARCODE","BED","REPORTNUM")],MARGIN = 1, FUN = function(x) paste(x, collapse = ","))
#subSampleTotal.reportnum <- subSample[!duplicated(subSample$tmp),]
#subSampleTotal.reportnum$tmp2 <- apply(X = subSampleTotal.reportnum[,c("SAMPLE","BARCODE","BED","REPORT")],MARGIN = 1, FUN = function(x) paste(x, collapse = ","))
tableOfMatches <- NULL

for(i in seq_along(subSampleTotal.3ident)){
  tableOfMatches[i] <- list(grep(subSampleTotal.3ident[i], subSampleTotal.report))
}

multipleMatches <- tableOfMatches[which(lengths(tableOfMatches) > 1)]

tableOfMatches2 <- NULL
for(i in seq_along(multipleMatches)){
  a <- NULL
  for(j in seq_along(multipleMatches[[i]])){
    a <- c(a, subSampleTotal.report[multipleMatches[[i]][j]])
  }
  tableOfMatches2[i] <- list(a)
}

save(tableOfMatches2, file = "/home/kevhu/data/listOfDuplicatedSampleRuns.Robj")


