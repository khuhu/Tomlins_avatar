#listOfDirectories <- c("/mnt/DATA3/Yoda_data_dump/Kelly_TPRKB_1_TranscriptomeReanalysis_430/",
#                       "/mnt/DATA3/Yoda_data_dump/Kelly_TPRKB_2_TranscriptomeReanalysis_431/")

listOfDirectories <- c("/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-173-Cho_mouse_20171130_362_439/plugin_out/coverageAnalysis_out.767/",
                       "/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-174-Cho_mouse_20171213_363_441/plugin_out/coverageAnalysis_out.769/")


#directory <-"/mnt/DATA3/Yoda_data_dump/LorenaUTReanalysisTxn_429/plugin_out/ampliSeqRNA_out.754"
#setwd(directory)
barcodes.list <- NULL
file.paths <- NULL
sample.Names <- NULL
for(i in seq_along(listOfDirectories)){
  setwd(listOfDirectories[i])
  dummy1 <- system('find . -name "*amplicon.cov.xls"', intern = TRUE)
  dummy1 <- gsub("\\./", listOfDirectories[i], dummy1)
  dummy2 <- system('find . -name "*bc_summary.xls"', intern = TRUE)
  dummy2 <- gsub("\\./", listOfDirectories[i], dummy2)
  file.paths <- c(file.paths, dummy1)
  barcodes <- system('find . -type d -name "IonXpress*"', intern = TRUE)
  #system('find . -type d -name "IonXpress*"')
  for(i in seq_along(barcodes)){
    barcodes[i] <- gsub(".*Ion","Ion",barcodes[i])
  }
  barcodes.list <- c(barcodes.list, barcodes)
  sample.Names <- c(sample.Names, dummy2)
}


#file.paths <- system('find . -name "*amplicon.cov.xls"', intern = TRUE)
file.paths2 <- file.paths[-c(grep("scraper", file.paths))]
sample.Names2 <- sample.Names[-c(grep("scraper", sample.Names))]
a <- read.table(file.paths2[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
a <- a[order(a$contig_id),]


### this loops is an artefact and I'm too lazy to change it - i know it's inefficient
listOfVarTabs <- NULL
for(i in seq_along(file.paths2)){
  dummyTable <- read.table(file.paths2[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  assign(paste0("VarTab",i), dummyTable[order(dummyTable$contig_id),])
  listOfVarTabs <- c(listOfVarTabs, paste0("VarTab",i))
}

agglo.table <- NULL
for(i in seq_along(listOfVarTabs)){
  dummyMatrix <- eval(as.name(listOfVarTabs[i]))
  dummyMatrix.2 <- apply(dummyMatrix[,c("fwd_e2e","rev_e2e")],1, sum)
  dummyMatrix.3 <- dummyMatrix[,"total_reads"]
  agglo.table <- cbind(agglo.table, dummyMatrix.2, dummyMatrix.3)
  colnames(agglo.table)[ncol(agglo.table)-1] <- paste0(barcodes.list[i],".total_cov")
  colnames(agglo.table)[ncol(agglo.table)] <- paste0(barcodes.list[i],".total_reads")
}


tableOfSamps <- NULL
for(i in seq_along(sample.Names2)){
  c <- read.table(paste(sample.Names2[i]), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  c.1 <- c[,1:2]
  tableOfSamps <- rbind(tableOfSamps, c.1)
}

###for Kelly's data in particular 
tableOfNames <- read.table("/mnt/DATA4/kevhu/KellyData/KellyNameList.csv", sep = ",", header = FALSE,
                           stringsAsFactors = FALSE)

###cool bit I learned from kelly's data about exact matching using anchors
for(i in 1:nrow(tableOfNames)){
  tableOfSamps$Sample.Name <- gsub(paste("^(",tableOfNames[i,1],")$", sep = ""),tableOfNames[i,2], tableOfSamps$Sample.Name)
}



a.1 <- a[,1:6]
finalTable <- agglo.table
for(i in c(6:1)){
  finalTable <- cbind(a[i], finalTable)
}


finalColNames <- colnames(finalTable)[7:ncol(finalTable)]
for(i in 1:nrow(tableOfSamps)){
  finalColNames <- gsub(tableOfSamps[i,1], tableOfSamps[i,2], finalColNames)
}



colnames(finalTable)[7:ncol(finalTable)] <- finalColNames

totReadCols <- grep("total_reads",colnames(finalTable))
finalTable2 <- finalTable[,-totReadCols]

brca1.subset <- finalTable2[which(grepl("Brca1",finalTable2$attributes)),1:ncol(finalTable2)]
brca1.subset <- brca1.subset[order(brca1.subset$contig_id, brca1.subset$contig_srt),]
par(cex.axis=0.5)
boxplot(brca1.subset, las = 2)
#finalTable <- finalTable[,order(colnames(finalTable))]
#write.table(finalTable, "/mnt/DATA4/kevhu/KellyData/20171204TRKB_Cov_Combined.Uniq.txt", quote = FALSE, sep = "\t", row.names = FALSE)






