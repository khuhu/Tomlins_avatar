library(vcfR)
vcfDir <- "/home/kevhu/scripts/newMousePanelPipeline/vcfs/Reanalysis_AUS5-76-MG_test1_217/"

fileList <- system(paste("ls", vcfDir), intern = TRUE)
fileList <- fileList[grep("*\\.vcf.gz",fileList)]
fileList <- fileList[-grep("*tbi",fileList)]
vcfPaths <- paste0(vcfDir, fileList)

allVcfDf <- NULL
for (i in vcfPaths) {
  testVcf <- read.vcfR(i)
  testVcf_df <- data.frame(testVcf@fix, stringsAsFactors = FALSE)
  testVcf_df2 <- testVcf_df[grep("SNP", testVcf_df$ID),]
  snp_idx <- grep("SNP", testVcf_df$ID)
  
  testVcf_df_gt <- data.frame(testVcf@gt, stringsAsFactors = FALSE)
  testVcf_df_gt <- testVcf_df_gt[snp_idx,]
  testVcf_df_gt$FORMAT[1]
  
  sampNames <- colnames(testVcf_df_gt)[2]
  tmpGt <- NULL
  for (j in 1:nrow(testVcf_df_gt)) {
    print(i)
    tmpGtVec <- strsplit(x = testVcf_df_gt[,2][j], split = ":")
    tmp <- c(unlist(lapply(tmpGtVec, '[[', 1)), unlist(lapply(tmpGtVec, '[[', 2)), 
             unlist(lapply(tmpGtVec, '[[', 4)), unlist(lapply(tmpGtVec, '[[', 8)),
             unlist(lapply(tmpGtVec, '[[', 9)), unlist(lapply(tmpGtVec, '[[', 14)),
             unlist(lapply(tmpGtVec, '[[', 15)))
    tmpGt <- rbind(tmpGt, tmp)
  }
  
  colnames(tmpGt) <- c("GT", "GQ", "FDP", "FAO", "AF", "FSAR", "FSAF")
  testVcf_df3 <- cbind("Sample" = rep(sampNames, nrow(tmpGt)), testVcf_df2[,1:7], tmpGt)
  
  allVcfDf <- rbind(allVcfDf, testVcf_df3)
}

allVcfDf$GT <- as.character(allVcfDf$GT)

### filtering out SNP calls where indels are found - honestly, I think don't need to filter these SNPs, just make them no calls if odd
### so SNPS I'm comepletely giving all no calls to are ones where incidence rate of indel > 2/3
### eventually might 

allVcfDf2 <- allVcfDf
possibleIndels <- allVcfDf[which(nchar(allVcfDf$REF) > 1 | nchar(allVcfDf$ALT) > 1),]
badSNPs <- names(which(table(possibleIndels$ID) > (37 * 1/2)))
#badSNPs <- names(which(table(possibleIndels$ID) > 1))
allVcfDf2 <- allVcfDf2[-which(allVcfDf2$ID  %in% badSNPs),]
allVcfDf2$GT[which(nchar(allVcfDf2$REF) > 1 | nchar(allVcfDf2$ALT) > 1)] <- "./."

### get general idea of missing genotypes per sample
for (i in unique(allVcfDf2$Sample)) {
  tmpDf <- allVcfDf2[which(allVcfDf2$Sample == i),]
  print(i)
  print(length(which(tmpDf$GT == "./."))/length(unique(allVcfDf2$ID)))
}


### R admixture might be too complicated to get going at this point.
### going to create a bed file these unique regions and then use
### admixture on the cmd line

bedRegionsFile <- data.frame(str_remove(allVcfDf2$CHROM[1:449],"chr"), allVcfDf2$POS[1:449])
write.table(bedRegionsFile, "/mnt/DATA5/tmp/kev/misc/20210131snpsRegionFile.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
bedRegionsFile2 <- bedRegionsFile
bedRegionsFile2$str_remove.allVcfDf2.CHROM.1.449....chr.. <- paste0("chr", bedRegionsFile2$str_remove.allVcfDf2.CHROM.1.449....chr..)

write.table(bedRegionsFile2, "/mnt/DATA5/tmp/kev/misc/20210131snpsRegionFile_withchr.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

### structure plots
###

fullQTable <- read.table("/mnt/DATA5/tmp/kev/geno2/20210202myplinkSnps.5.Q",
                         sep = " ", stringsAsFactors = FALSE, header = FALSE)
sampNames <- read.table("/mnt/DATA5/tmp/kev/geno2/20210202myplinkSnps.fam",
                        stringsAsFactors = FALSE, sep = " ")

colnames(fullQTable) <- c("129", "6J", "BALB", "C3H", "FVB")

#fullQTable2 <- cbind(sampNames$V2, fullQTable,
#                     paste0(signif(fullQTable$`6J`, digits = 2), "/", signif(fullQTable$`129S1`, digits = 2)))
fullQTable2 <- data.frame(cbind(sampNames$V2, fullQTable), stringsAsFactors = FALSE)


fullQTable2 <- melt(fullQTable2)
fullQTable2$variable <- as.character(fullQTable2$variable)
#fullQTable2$`rownames(fullQTable)` <- fullQTable2$sampNames.V2
fullQTable2 <- fullQTable2[order(fullQTable2$value, decreasing = TRUE),]
fullQTable2$value <- signif(fullQTable2$value, digits = 2)
colnames(fullQTable2) <- c("SampleName", "label","Fraction")
fullQTable2$SampleName <- factor(fullQTable2$SampleName, levels = unique(fullQTable2$SampleName))
fullQTable2$label <- as.character(fullQTable2$label)
fullQTable2$label <- str_remove(fullQTable2$label, "X")
#tmp1 <- nrow(fullQTable2)/2 + 1
#tmp2 <- nrow(fullQTable2)
#fullQTable2$label[eval(tmp1:tmp2)] <- " "

colPalette <- c("#009E73","#D55E00", "#001a9e", "#9e007e", "#9e9600")

ggplot(fullQTable2, aes(x = SampleName, y = Fraction, fill=label)) +
  geom_bar(position="fill", stat="identity", color = "black") + 
  scale_fill_manual(values=colPalette) + theme_bw() + 
  theme(axis.text.x=element_text(angle=90,hjust=1))


