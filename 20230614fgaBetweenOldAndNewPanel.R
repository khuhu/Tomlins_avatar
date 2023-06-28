### need to do fga comparison for samples sequenced between the two panels
allQcSamps <- read.table("/mnt/DATA5/tmp/kev/misc/20230614allQcSampsNewMousePanel.txt", sep = "\t",
                         stringsAsFactors = FALSE, header = TRUE)

yodaDir <- "/mnt/DATA6/Yoda_data_dump/"
mouseSampDirs <- c("Auto_user_SATProton-166-Cho_mouse_1_353_410", "Auto_user_SATProton-171-Cho_mouse_20170925_358_422",
                   "Auto_user_SATProton-173-Cho_mouse_20171130_362_439", "Auto_user_SATProton-174-Cho_mouse_20171213_363_441")

yodaDir2 <- paste0(yodaDir, mouseSampDirs)

oldMousedata <- NULL
for (i in yodaDir2) {
  setwd(i)
  tmp <- system("find . -name '*\\.bc_summary.xls'", intern = TRUE)
  if (length(tmp) > 1) {
    tmp <- tmp[-grep("link", tmp)]
  }
  tmpXls <- read.table(tmp, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  tmpXls$dir <- i
  oldMousedata  <- rbind(oldMousedata, tmpXls)
}

oldMousedata$sampleStripped <- str_remove(oldMousedata$Sample.Name, "\\_")
oldMousedata$sampleStripped <- str_remove(oldMousedata$sampleStripped, "\\-")
oldMousedata$sampleStripped <- str_remove_all(oldMousedata$sampleStripped, " ")
oldMousedata$sampleStripped <- tolower(oldMousedata$sampleStripped)
oldMousedata$sampleStripped <- str_remove_all(oldMousedata$sampleStripped, "o")

which(oldMousedata$sampleStripped %in% allQcSamps$strippedName2)
oldMousedata[which(oldMousedata$sampleStripped %in% allQcSamps$strippedName2),]
oldMousedata$sampleStripped[which(oldMousedata$sampleStripped %in% allQcSamps$strippedName2)]
allQcSamps$Sample.Name[which(allQcSamps$strippedName2 %in% oldMousedata$sampleStripped)]
### all 52 samples from cho-174 seem to be in between the two
### load in the samples and then do fraction of genome altered

allMNgsGeneMatrix <- read.table("/mnt/DATA5/tmp/kev/qcAllMousePaper/cnMatrix_gene.txt", sep = "\t",
                                header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

allMNgsGeneMatrix2 <- allMNgsGeneMatrix[, c(which(colnames(allMNgsGeneMatrix) %in% allQcSamps$Sample.Name[which(allQcSamps$strippedName2 %in% oldMousedata$sampleStripped)]))]

namesDf <- data.frame("newPanelNames" = colnames(allMNgsGeneMatrix2))
namesDf$oldPanelNames <- oldMousedata$sampleStripped[match(allQcSamps$strippedName2[which(allQcSamps$Sample.Name %in% namesDf$newPanelNames)],
                                                  oldMousedata$sampleStripped)]

oldMouseGeneMatrix <- read.table("/mnt/DATA5/tmp/kev/oldMouse/cnMatrix_gene.txt", sep = "\t",
                                 header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

colnames(oldMouseGeneMatrix) <- str_remove(colnames(oldMouseGeneMatrix), "\\_")
colnames(oldMouseGeneMatrix) <- str_remove(colnames(oldMouseGeneMatrix), "\\-")
colnames(oldMouseGeneMatrix) <- str_remove_all(colnames(oldMouseGeneMatrix), " ")
colnames(oldMouseGeneMatrix) <- tolower(colnames(oldMouseGeneMatrix))
colnames(oldMouseGeneMatrix) <- str_remove_all(colnames(oldMouseGeneMatrix), "o")


oldMouseGeneMatrix2 <- oldMouseGeneMatrix[, c(which(colnames(oldMouseGeneMatrix)  %in% allQcSamps$strippedName2[which(allQcSamps$strippedName2 %in% oldMousedata$sampleStripped)]))]
oldMouseGeneMatrix2 <- log2(oldMouseGeneMatrix2)
oldMouseGeneMatrix2[abs(oldMouseGeneMatrix2) < 0.2] <- 0

allMNgsGeneMatrix2 <- log2(allMNgsGeneMatrix2)
allMNgsGeneMatrix2[abs(allMNgsGeneMatrix2) < 0.2] <- 0

fga <- NULL
for (i in 1:nrow(namesDf)) {
  tmpOld <- length(which(oldMouseGeneMatrix2[ ,namesDf$oldPanelNames[i]] != 0))/35
  tmpNew <- length(which(allMNgsGeneMatrix2[ ,namesDf$newPanelNames[i]] != 0))/128
  fga <- rbind(fga, c(tmpOld, tmpNew))
}

fgaDf <- data.frame(fga)
colnames(fgaDf) <- c("oldPanel", "newPanel")
fgaDf <- fgaDf[order(fgaDf$oldPanel, decreasing = FALSE), ]

ggplot(fgaDf, aes(x = oldPanel, y = newPanel)) + geom_point() + xlab("35 gene panel fraction of genes altered") + 
  ylab("128 gene panel fraction of genes altered") + geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  ggtitle("Comparison of fraction of genes altered between small and large gene panel (R = 0.63; n = 51)") +
  theme(plot.title = element_text(hjust = 0.5))
  

