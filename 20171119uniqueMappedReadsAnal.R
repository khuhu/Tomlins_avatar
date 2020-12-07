listOfDirectories <- c("/mnt/DATA3/Zeus_data_dump/UT_Transcriptome_3Reanalysis_292/plugin_out/ampliSeqRNA_out.480/",
                       "mnt/DATA3/Yoda_data_dump/SATProton-172-UT_transcriptomeReAnalysis_432/plugin_out/ampliSeqRNA_out.757/",
                       "/mnt/DATA3/Zeus_data_dump/UT_Transcriptome_2Reanalysis_293/plugin_out/ampliSeqRNA_out.481")

listOfDirectories <- "/mnt/DATA3/Zeus_data_dump/plugin_out/ampliSeqRNA_out.480/"
barcodes.list <- NULL
file.paths <- NULL
> sample.Names <- NULL
> for(i in seq_along(listOfDirectories)){
  +   setwd(listOfDirectories[i])
  +   dummy1 <- system('find . -name "*amplicon.cov.xls"', intern = TRUE)
  +   dummy1 <- gsub("\\./", listOfDirectories[i], dummy1)
  +   dummy2 <- system('find . -name "*bc_summary.xls"', intern = TRUE)
  +   dummy2 <- gsub("\\./", listOfDirectories[i], dummy2)
  +   file.paths <- c(file.paths, dummy1)
  +   barcodes <- system('find . -type d -name "IonXpress*"', intern = TRUE)
  +   #system('find . -type d -name "IonXpress*"')
    +   for(i in seq_along(barcodes)){
      +     barcodes[i] <- gsub(".*Ion","Ion",barcodes[i])
      +   }
  +   barcodes.list <- c(barcodes.list, barcodes)
  +   sample.Names <- c(sample.Names, dummy2)
  + }
Error in setwd(listOfDirectories[i]) : cannot change working directory
> getwd()
[1] "/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-165-Kelly_TPRKB_1_Transcriptome_Human_Gene_Expression_351_408/plugin_out/coverageAnalysis_out.716"
> listOfDirectories <- "/mnt/DATA3/Zeus_data_dump/UT_Transcriptome_3Reanalysis_292/plugin_out/ampliSeqRNA_out.480/"
> barcodes.list <- NULL
> file.paths <- NULL
> sample.Names <- NULL
> for(i in seq_along(listOfDirectories)){
  +   setwd(listOfDirectories[i])
  +   dummy1 <- system('find . -name "*amplicon.cov.xls"', intern = TRUE)
  +   dummy1 <- gsub("\\./", listOfDirectories[i], dummy1)
  +   dummy2 <- system('find . -name "*bc_summary.xls"', intern = TRUE)
  +   dummy2 <- gsub("\\./", listOfDirectories[i], dummy2)
  +   file.paths <- c(file.paths, dummy1)
  +   barcodes <- system('find . -type d -name "IonXpress*"', intern = TRUE)
  +   #system('find . -type d -name "IonXpress*"')
    +   for(i in seq_along(barcodes)){
      +     barcodes[i] <- gsub(".*Ion","Ion",barcodes[i])
      +   }
  +   barcodes.list <- c(barcodes.list, barcodes)
  +   sample.Names <- c(sample.Names, dummy2)
  + }
> file.paths2 <- file.paths[-c(grep("scraper", file.paths))]
> sample.Names2 <- sample.Names[-c(grep("scraper", sample.Names))]
> a <- read.table(file.paths2[1], header = TRUE)
> sumColsNames <- colnames(a)[7:12]
> listOfVarTabs <- NULL
> for(i in seq_along(file.paths2)){
  +   a <- read.table(file.paths2[i], header = TRUE)
  +   assign(paste0("VarTab",i), a[,sumColsNames])
  +   listOfVarTabs <- c(listOfVarTabs, paste0("VarTab",i))
  + }
> agglo.table <- NULL
> for(i in seq_along(listOfVarTabs)){
  +     dummyMatrix <- eval(as.name(listOfVarTabs[i]))
  +     dummyMatrix.2 <- apply(dummyMatrix[,c("fwd_cov","rev_cov")],1, sum)
  +     dummyMatrix.3 <- dummyMatrix[,"total_reads"]
  +     agglo.table <- cbind(agglo.table, dummyMatrix.2, dummyMatrix.3)
  +     colnames(agglo.table)[ncol(agglo.table)-1] <- paste0(barcodes.list[i],".total_cov")
  +     colnames(agglo.table)[ncol(agglo.table)] <- paste0(barcodes.list[i],".total_reads")
  + }
Warning message:
  In .Method(..., deparse.level = deparse.level) :
  number of rows of result is not a multiple of vector length (arg 2)
> View(agglo.table)
> a.1 <- a[,1:6]
> finalTable <- agglo.table
> for(i in c(6:1)){
  +   finalTable <- cbind(a[i], finalTable)
  + }
> tableOfSamps <- NULL
> for(i in seq_along(sample.Names2)){
  +   c <- read.table(paste(sample.Names2[i]), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  +   c.1 <- c[,1:2]
  +   tableOfSamps <- rbind(tableOfSamps, c.1)
  + }
> finalColNames <- colnames(finalTable)
> for(i in 1:nrow(tableOfSamps)){
  +   finalColNames <- gsub(tableOfSamps[i,1], tableOfSamps[i,2], finalColNames)
  + }