### use below to quickly read the 4.5m variants from cosmic

#library(data.table)
#cmc_table <- fread("/mnt/DATA5/tmp/kev/misc/cmc_export.v92.edited.tsv")
#cmc_tsg <- cmc_table[grep(cmc_table$`Mutation Description AA`,
#                          pattern = "Substitution - Missense"),]
#small_test <- cmc_tsg[sample(1:nrow(cmc_tsg), 300),]

source("/home/kevhu/scripts/20201020hotspotConversionFunctions.R")


### slight alteration for matching gene names
genomeToAA <- function(liftOverRes){
  ### make this work for any species to any species at a later date
  # tmpGeneName <- geneNameDf$mmusculus_homolog_associated_gene_name[which(geneNameDf$external_gene_name == liftOverRes$gene)]
  tmpGeneName <- geneNameDf$mmusculus_homolog_associated_gene_name[grep(paste0("^", liftOverRes$gene, "$"), geneNameDf$external_gene_name)]
  genomeCdsTmpTable <- mm10biomartTable[which(mm10biomartTable$external_gene_name == tmpGeneName),]
  genomeCdsTmpTable <- cdsConversion(genomeCdsTmpTable)
  
  cdsRangeBp <- IRanges(start = liftOverRes$start, end = liftOverRes$end)
  cdsExonRangeBp <- IRanges(start = genomeCdsTmpTable$exon_chrom_start, end = genomeCdsTmpTable$exon_chrom_end)
  cdsExon <- genomeCdsTmpTable[subjectHits(findOverlaps(cdsRangeBp, cdsExonRangeBp)),]
  
  ### the subtraction of end/start depends on strand
  tmpPos <- cdsExon$cds_start + abs(cdsExon$strand_cds_start_conv - mean(c(liftOverRes$start, liftOverRes$end))) - 1
  tmpAaPos <- ceiling(tmpPos/3)
  
  ### might need to return more
  return(tmpAaPos)
}

small_cosmic <- read.table("/home/kevhu/data/20190225mouseConvertedOncogeneHotspots.txt",
                           sep = "\t", stringsAsFactors = FALSE, header = TRUE)
hotspots <- unique(small_cosmic$Hotspot.Residue)

hotspotDf <- NULL
for (i in hotspots) {
  tmpVar <- strsplit(i, "\\|")
  tmpGene <- tmpVar[[1]][1]
  tmpVar2 <- tmpVar[[1]][2]
  tmpAa <- substr(tmpVar2, 1, 1)
  tmpPos <- substr(tmpVar2, 2, nchar(tmpVar2))
  tmpLine <- c(tmpGene, tmpAa, tmpPos)
  hotspotDf <- rbind(hotspotDf, tmpLine)
}

rownames(hotspotDf) <- NULL
colnames(hotspotDf) <- c("Gene", "Aa", "Position")
hotspotDf <- data.frame(hotspotDf, stringsAsFactors = FALSE)


hotspotDf <- hotspotDf[-which(hotspotDf == "GNAS"),] # different transcript ID
# hotspotDf <- hotspotDf[-which(hotspotDf == "GNA11"),]
hotspotDf <- hotspotDf[-which(hotspotDf == "GNAQ"),] # the name searching method is off -----------------------fixed
hotspotDf <- hotspotDf[-which(hotspotDf == "H3F3A"),] # Lower case is mouse name H3F3A, hg38 name is H3-3A
hotspotDf <- hotspotDf[-which(hotspotDf == "DNMT3A"),] # not canonical isoform? - liftOver empty
hotspotDf <- hotspotDf[-which(hotspotDf == "FOXL2"),] # msa error - no seq2 - missing mouse peptide
hotspotDf <- hotspotDf[-which(hotspotDf == "HIST1H3B"),] # not the gene name

### reworking naming system - lowercase search for then use the table - table makes
### weird/missing results i.e SPOP


### so to redesign portions. as mentioned before, first (1) the gene name change should be easy enough to implement
### (2) fix other issues like msa matching and what-not (3) the wrong transcript-isoform being useds
### the last portion is the hardest to deal with since, I need to get all the different isoforms ... and figure out
### an efficient way to query them - for the time being just skip ones like this




### 
###
###

resTbl <- NULL
for (i in 1:nrow(hotspotDf)) {
  print(i)
  testAaConversion <- tryCatch(aaToGenome(gene = hotspotDf$Gene[i],
                                 position = as.numeric(hotspotDf$Position[i])),
                               error = function(x) return(NULL))
  #print(testAaConversion)
  #if (testAaConversion == "empty") {
  #  next()
  #}
  if (is.null(testAaConversion)) {
    next()
  }
  
  liftOverOut <- cdsLiftOver(testAaConversion, chainFile = h38toMm10Chain)
  #print(liftOverOut)
  if (is.null(liftOverOut)) {
    next()
  }
  
  checkInput <- list(hotspotDf$Gene[i],
                     as.numeric(hotspotDf$Position[i]) ,
                     "liftOverOut" = liftOverOut)
  #print(checkInput)
  conversionCheckOut <- aaConversionCheck(checkInput)
  #print(conversionCheckOut)
  tmpSeq1 <- conversionCheckOut$seq1
  tmpSeq2 <- conversionCheckOut$seq2
  tmpConSeq <- conversionCheckOut$conSeq
  tmpCheck <- conversionCheckOut$check
  tmpScore <- mean(conversionCheckOut$conservationScore)
  resTbl <- rbind(resTbl, c(tmpSeq1, tmpSeq2, tmpConSeq, tmpCheck, tmpScore))
}




