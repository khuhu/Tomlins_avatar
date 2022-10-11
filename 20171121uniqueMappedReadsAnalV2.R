listOfDirectories <- c("/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-172-UT_transcriptome_359_424/plugin_out/coverageAnalysis_out.746/")


file.paths <- NULL
for(i in seq_along(listOfDirectories)){
  setwd(listOfDirectories[i])
  dummy1 <- system('find . -name "*amplicon.cov.xls"', intern = TRUE)
  dummy1 <- gsub("\\./", listOfDirectories[i], dummy1)
  file.paths <- c(file.paths, dummy1)
}
file.paths2 <- file.paths[-c(grep("scraper", file.paths))]
a <- read.table(file.paths2[1], header = TRUE, stringsAsFactors = FALSE, sep = "\t")
a <- a[order(a$contig_id),]
sumColsNames <- colnames(a)[7:12]
numRows <- nrow(read.table(file.paths2[1], header = TRUE))
listOfVarTabs <- NULL
for(i in seq_along(file.paths2)){
  dummyTable <- read.table(file.paths2[i], header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  assign(paste0("VarTab",i), dummyTable[order(dummyTable$contig_id),])
  listOfVarTabs <- c(listOfVarTabs, paste0("VarTab",i))
}
agglo.table <- NULL
for(i in 1:numRows){
  matrixForSumming <- NULL
  for(j in seq_along(listOfVarTabs)){
    matrixForSumming <- rbind(matrixForSumming,eval(as.name(listOfVarTabs[j]))[i, sumColsNames])
    }
  aggloRow <- colSums(matrixForSumming)
  agglo.table <- rbind(agglo.table, aggloRow)
}

finalTable <- a[,1:6]
for(i in 1:ncol(agglo.table)){
  finalTable <- cbind(finalTable, agglo.table[,i])
  colnames(finalTable)[i+6] <- colnames(agglo.table)[i]
}



write.table(finalTable, "/mnt/DATA4/kevhu/LorenaData/20171128CombinedCov.1cohort.nonUniq.txt", quote = FALSE, sep = "\t", row.names = FALSE)



###checking for just one gene, to see whether or not, it added correctly
varA.1 <- NULL
for(i in seq_along(listOfVarTabs)){
  varA <- eval(as.name(listOfVarTabs[i]))[grep("ADSL",eval(as.name(listOfVarTabs[i]))$attributes),"total_reads"]
  varA.1 <- c(varA.1, varA)
}

sum(varA.1)

