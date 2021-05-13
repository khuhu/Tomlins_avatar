library(GenomicRanges)
library(rtracklayer)
library(stringr)

# loading data tables
# eventually get the exons (i.e include UTRs)

geneCdsHs <- read.table("/home/kevhu/data/20210203hg38KnownCanbiomartQuery.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

geneCdsMm <- read.table("/home/kevhu/data/20210203Mm10KnownCanbiomartQuery.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

cBioportalIds <- read.table("/mnt/DATA5/tmp/kev/misc/20210404cbioBiorMart_genes_tid.txt",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)


geneNamesDf <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/20200720genesDf.txt",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)


# used for converting gene name between speicies
firstUpper <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}


# cds conversion from cds to genomic position
# updated 20210512: simple and works as intended
cdsConversion <- function(df){
  if(df$strand[1] == 1){
    df$exon_chrom_start_strand <- df$exon_chrom_start
    df$exon_chrom_end_strand <- df$exon_chrom_end
  } else{
    df$exon_chrom_start_strand <- df$exon_chrom_end
    df$exon_chrom_end_strand <- df$exon_chrom_start
  }
  
  strand_cds_start_conv <- NULL
  strand_cds_end_conv <- NULL
  
  for (i in 1:nrow(df)) {
    strand_cds_start_conv <- c(strand_cds_start_conv,
                               df$exon_chrom_start_strand[i] + df$strand[i] * (df$cdna_coding_start[i] - 1))
    strand_cds_end_conv <- c(strand_cds_end_conv,
                             df$exon_chrom_start_strand[i] + df$strand[i] * (df$cdna_coding_end[i] - 1))
  }
  df$strand_cds_start_conv <- strand_cds_start_conv
  df$strand_cds_end_conv <- strand_cds_end_conv
  return(df)
}




# GenomicRanges grange object to Df
grangeToDf <- function(df){
  require(GenomicRanges)
  out <- data.frame(df@unlistData, mcols(df@unlistData), stringsAsFactors = FALSE)
  return(out)
}


# makes excel list of targets into list for downstream 
parseInputTable <- function(tbl){
  tblExons <- NULL
  tblHotspots <- NULL
  tblCopynumber <- NULL
  tblAllExons <- NULL
  for (i in 1:nrow(tbl)) {
    tmpRow <- tbl[i,]
    typeOfGene <- unlist(strsplit(tmpRow$Type, ","))
    if ("Exons" %in% typeOfGene) {
      tblExons <- c(paste(tmpRow$Gene, tmpRow$ExonInfo, sep = "|"), tblExons)
    }
    if ("Hotspot" %in% typeOfGene) {
      tblHotspots <- c(paste(tmpRow$Gene, tmpRow$HotspotInfo, sep = "|"), tblHotspots)
    }
    if ("Copynumber" %in% typeOfGene) {
      ### should be a two item list one with gene name, and other with amplicon number
      ### default 12, other column for custom number of amplicons if not empty
      numAmplicons <- ifelse (is.na(tmpRow$Other), 12, tmpRow$Other) 
      tblCopynumber <- c(paste(tmpRow$Gene, numAmplicons, sep = "|"), tblCopynumber)
    }
    if ("AllExons" %in% typeOfGene) {
      tblAllExons <- c(tmpRow$Gene, tblAllExons)
    }
  }
  listOfTbls <- list(tblExons, tblHotspots, tblCopynumber, tblAllExons)
}




exonConversionLOv2 <- function(geneName, exonList){
  require(GenomicRanges)
  require(rtracklayer)
  
  tmpExonDat <- geneExonResHs[which(geneExonResHs$external_gene_name == geneName),]
  finalTable <- NULL
  for (i in exonList) {
    tmpDf <- tmpExonDat[which(as.numeric(tmpExonDat$rank) == i),]
    tmpDf2 <- data.frame("chromosome" = paste0("chr", tmpDf$chromosome_name),
                         "start" = tmpDf$exon_chrom_start,
                         "end" = tmpDf$exon_chrom_end, stringsAsFactors = FALSE)
    
    tmpGrange <- makeGRangesFromDataFrame(tmpDf2)
    loOut <- liftOver(tmpGrange, ch)
    if (isEmpty(loOut)) {
      next()
    }
    loOutDf <- grangeToDf(loOut)
    loOutDf <- cbind(loOutDf, "original" = paste0(geneName,"exon", i), "gene" = geneName)
    finalTable <- rbind(finalTable, loOutDf)
  }
  
  finalTable
}



orthoPosAaLoNoMsaV2 <- function(geneName, aminoAcid , position, speciesStart){
  require(GenomicRanges)
  require(rtracklayer)
  
  
  tmpCdsDat <- geneCdsResHs[which(geneCdsResHs$external_gene_name == geneName),]
  if (nrow(tmpCdsDat) == 0) {
    print("gene not found")
    return(data.frame("seqnames" = NA, "start" = NA, "end" = NA, "width" = NA,
                      "strand" = NA, "original" = paste0(geneName, position, aminoAcid),
                      stringsAsFactors = FALSE, "gene" = geneName))
  }
  
  tmpMin <- min(tmpCdsDat$cdna_coding_start)
  tmpExonDat <- geneExonResHs[which(geneExonResHs$external_gene_name == geneName),]
  tmpExonDat$tmp <- tmpExonDat$exon_chrom_end - tmpExonDat$exon_chrom_start
  tmpExonDat <- tmpExonDat[order(tmpExonDat$rank),]
  tmpExonDat <- cdsConversion(tmpExonDat)
  
  #tmpExonDat2 <- geneExonMm[which(geneExonMm$external_gene_name == firstUpper(geneName)), ]
  
  ### the added variable is used to deal with strandedness 
  if (tmpExonDat$strand[1] > 0) {
    tmpExonDat$cdsRemapVar <- tmpExonDat$cdsStartConv
  }
  else if(tmpExonDat$strand[1] < 0){
    tmpExonDat$cdsRemapVar <- tmpExonDat$cdsEndConv
  }
  
  cdsPos <- position * 3 + tmpMin
  cdsSubject <- IRanges(cdsPos - 2, cdsPos)
  cdsRangeQuery <- IRanges(tmpExonDat$cdsStartConv, tmpExonDat$cdsEndConv)
  overlapRes <- GenomicRanges::findOverlaps(cdsRangeQuery,cdsSubject)
  
  tmpPos <- abs(cdsPos - tmpExonDat$cdsRemapVar[queryHits(overlapRes)]) + tmpExonDat$exon_chrom_start[queryHits(overlapRes)]  
  aaGenomeGrange <- tryCatch(makeGRangesFromDataFrame(data.frame("chromosome" = paste0("chr" ,tmpExonDat$chromosome_name[1]), "start" = tmpPos, "end" = tmpPos + 2)), 
                             error = function(x) NULL)
  
  if (is.null(aaGenomeGrange)) {
    return(data.frame("seqnames" = NA, "start" = NA, "end" = NA, "width" = NA,
                      "strand" = NA, "original" = paste0(geneName, position, aminoAcid),
                      stringsAsFactors = FALSE))
  }
  
  loOut <- liftOver(aaGenomeGrange, ch)
  loOutDf <- grangeToDf(loOut)
  
  if (nrow(loOutDf) == 0) {
    return(data.frame("seqnames" = NA, "start" = NA, "end" = NA, "width" = NA,
                      "strand" = NA, "original" = paste0(geneName, position, aminoAcid),
                      stringsAsFactors = FALSE, "gene" = geneName))
  }
  
  
  if (tmpExonDat$strand[1] > 0) {
    loOutDf$start <- loOutDf$start - 3
    loOutDf$end <- loOutDf$end - 3
    loOutDf <- cbind(loOutDf, "original" = paste0(firstUpper(geneName), position, aminoAcid), "gene" = geneName)
  }
  else if(tmpExonDat$strand[1] < 0){
    loOutDf <- cbind(loOutDf, "original" = paste0(firstUpper(geneName), position, aminoAcid), "gene" = geneName)
  }
  
  
  loOutDf
}



### used to LO in cn portion

cnLo <- function(df, ch) {
  res <- NULL
  tmpDf <- data.frame("chromosome" = paste0("chr", df$chromosome_name), "start" = df$exon_chrom_start, "end" = df$exon_chrom_end)
  for (i in 1:nrow(df)) {
    tmpGrange <- makeGRangesFromDataFrame(tmpDf[i,])
    tmpCnLoOut <- liftOver(tmpGrange, ch)
    tmpDf2 <- grangeToDf(tmpCnLoOut)
    if (length(is.infinite(tmpDf2$start)) == 0) {
      next()
    }
    res <- rbind(res, c("chromosome" = as.character(tmpDf2$seqnames[1]),
                        "start" = min(tmpDf2$start),"end" = max(tmpDf2$end), "rank" = df$rank[i]))
  } 
  res <- data.frame(res, stringsAsFactors = FALSE)
  res
}



# function to get amplicon targets for either cds/exon
# given x amount of targets wanted - no conversions

getAmpliconTargets <- function(geneList, type = "CDS", species = "human",
                               numAmps = 12, ampSize = 100){
  
  # this can be done better ... else if more species are added, I'll have 2^n species
  # booleans for checking 
  if (type == "CDS" & species == "human") {
    print("Using CDS for coordinates")
    tmpDb <- geneCdsHs
  } else{
    print("Using CDS + UTR for coordinates")
    # uncomment out once I get exon version
    #tmpDb <- geneExonhs
  }
  
  copyNumberOut <- NULL
  for (i in geneList) {
    tmpGeneName <- unlist(strsplit(i, "\\|"))[1]
    print(tmpGeneName)
    tmpEnstId <- cBioportalIds$Tid[which(cBioportalIds$Gene == tmpGeneName)]
    tmpGeneDat <- tmpDb[which(tmpDb$ensembl_transcript_id == tmpEnstId),]
    if (nrow(tmpGeneDat) == 0) {
      print(paste("no ensembl record", tmpGeneName))
      next()
    }
    
    # cds conversion gives us strand specific start and stop positions of the cds region
    tmpGeneDat2 <- cdsConversion(tmpGeneDat)
    tmpGeneDat2$tmp <- abs(tmpGeneDat2$strand_cds_start_conv - tmpGeneDat2$strand_cds_end_conv)
    tmpGeneDat2 <- tmpGeneDat2[order(tmpGeneDat2$rank),]
    
    numExons <- nrow(tmpGeneDat)
    # getting targets based on number of required targets vs number of exons
    # want to spread evenly if possible
    
    if(numAmps == numExons){
      print("# amp = # exons")
      tmpCoord <- NULL
      for (j in 1:nrow(tmpGeneDat2)) {
        tmpMidpoint <- mean(c(tmpGeneDat2$strand_cds_start_conv[j],
                              tmpGeneDat2$strand_cds_end_conv[j]))
        tmpStart <- tmpMidpoint - ampSize/2
        tmpEnd <- tmpMidpoint + ampSize/2
        tmpStart2 <- min(c(tmpStart, tmpEnd))
        tmpEnd2 <- max(c(tmpStart, tmpEnd))
        tmpDf <- data.frame("chromosome" = tmpGeneDat2$chromosome[j],
                            "start" = tmpStart2,
                            "end" =  tmpEnd2,
                            "position" = paste0(tmpGeneName, "exon", tmpGeneDat2$rank),
                               stringsAsFactors = FALSE)
        tmpCoord <- rbind(tmpCoord, tmpDf)
      }
    } else if(numExons > numAmps){
      # selects the n largest exons where n is number of desired amplicons
      print("# amp < # exons")
      tmpGeneDat3 <- tmpGeneDat2
      tmpGeneDat3 <- tmpGeneDat3[order(tmpGeneDat3$tmp, decreasing = TRUE),]
      tmpGeneDat3  <- tmpGeneDat3[1:numAmps,]
      tmpCoord <- NULL
      for (j in 1:nrow(tmpGeneDat3)) {
        tmpMidpoint <- mean(c(tmpGeneDat3$strand_cds_start_conv[j],
                              tmpGeneDat3$strand_cds_end_conv[j]))
        tmpStart <- tmpMidpoint - ampSize/2
        tmpEnd <- tmpMidpoint + ampSize/2
        tmpStart2 <- min(c(tmpStart, tmpEnd))
        tmpEnd2 <- max(c(tmpStart, tmpEnd))
        tmpDf <- data.frame("chromosome" = tmpGeneDat3$chromosome[j],
                            "start" = tmpStart2,
                            "end" =  tmpEnd2, "position" = paste0(tmpGeneName, "exon", tmpGeneDat3$rank[j]),
                            stringsAsFactors = FALSE)
        tmpCoord <- rbind(tmpCoord, tmpDf)
      }
    } else if(numAmps > numExons){
      ### only real nuance, evenly distribute extra amplicons onto the y(amps - exons) largest exons
      print("# amp > # exons")
      
      diffAmps <- numAmps - numExons
      tmpGeneDat3 <- tmpGeneDat2
      tmpGeneDat3$numberAmps <- 1
      tmpGeneDat3 <- tmpGeneDat3[order(tmpGeneDat3$tmp, decreasing = TRUE), ]
      overflowExons <- ifelse(numExons > diffAmps, diffAmps %% numExons, numExons %% diffAmps)
      overflowAmps <- ifelse(numExons > diffAmps, ceiling(diffAmps/numExons), floor(diffAmps/numExons))
      tmpGeneDat3$numberAmps[1:overflowExons] <- tmpGeneDat3$numberAmps[1:overflowExons] + overflowAmps
      
      ### special case where number of amps isn't divisible by number of exons and modulous of that is greater than 1
      ### i. example is 5 exons and 12 amps
      overflowExons2 <- numAmps %/% numExons
      if(overflowExons2 > 0 & overflowExons2 < numExons & numExons < diffAmps){
        tmpGeneDat3$numberAmps[1:overflowExons2] <- tmpGeneDat3$numberAmps[1:overflowExons2] + 1
      }
      
      
      # Below separates the above into single/multi amp targets
      # to process
      tmpDat_1amp <- tmpGeneDat3[which(tmpGeneDat3$numberAmps == 1), ]
      tmpDat_moreAmps <- tmpGeneDat3[which(tmpGeneDat3$numberAmps > 1), ]
      
      
      tmpCoords_1amp <- NULL
      for (j in 1:nrow(tmpDat_1amp)) {
        tmpMidpoint1Amp <- mean(c(tmpDat_1amp$strand_cds_start_conv[j], tmpDat_1amp$strand_cds_end_conv[j]))
        tmpStart <- tmpMidpoint1Amp - ampSize/2
        tmpEnd <- tmpMidpoint1Amp + ampSize/2
        tmpStart2 <- min(c(tmpStart, tmpEnd))
        tmpEnd2 <- max(c(tmpStart, tmpEnd))
        
        tmpDf <- tryCatch(data.frame("chromosome" = tmpDat_1amp$chromosome[j],
                   "start" = tmpStart2,
                   "end" = tmpEnd2,
                   "position" = paste0(tmpGeneName, "exon", tmpDat_1amp$rank[j]),
                   stringsAsFactors = FALSE),
                   error = function(x) return(NULL))
        tmpCoords_1amp <- rbind(tmpCoords_1amp, tmpDf)
      }
      
      
      
      
      # first I want to check if any of the proposed small cds genes can be detected with this method
      # 20210513: yeah I don't know what this is for
      #if (sum(tmpGeneDat3$tmp) < 1300) {
      #  smallCdsGenes <- c(smallCdsGenes, paste0(i, "|", sum(tmpGeneDat3$tmp)))
      #}
      
      #print(tmpExonDat_1amp)
      #print(tmpExonDat_moreAmps)
      
      
      tmpCoords_moreAmps <- NULL
      for (j in 1:nrow(tmpDat_moreAmps)) {
        strand <- NULL
        if (tmpDat_moreAmps$strand_cds_start_conv[j] > tmpDat_moreAmps$strand_cds_end_conv[j]) {
          strand <- 1
        } else{
          strand <- -1
        }
        
        tmpDf <- NULL
        tmpIncrement <- floor(tmpDat_moreAmps$tmp[j]/tmpDat_moreAmps$numberAmps[j])
        for (k in 1:tmpDat_moreAmps$numberAmps[j]) {
          tmpStart <- tmpDat_moreAmps$strand_cds_start_conv[j] + strand * (tmpIncrement * k - ampSize/2)
          tmpEnd <- tmpStart + strand  * ampSize
          tmpStart2 <- min(c(tmpStart, tmpEnd))
          tmpEnd2 <- max(c(tmpStart, tmpEnd))
          
          tmpDf <- rbind(tmpDf, data.frame("chromosome" = tmpDat_moreAmps$chromosome[j],
                                           "start" = tmpStart2, "end"  = tmpEnd2, 
                                           "position" = paste0(tmpGeneName, "exon", tmpDat_moreAmps$rank[j], "|", k) ,
                                           stringsAsFactors = FALSE))
        }
        tmpCoords_moreAmps <- rbind(tmpCoords_moreAmps, tmpDf)
      }
      # combine after processing these specific cases
      
      tmpCoord <- rbind(tmpCoords_1amp, tmpCoords_moreAmps)
    }
    tmpCoord$gene <- tmpGeneName
    copyNumberOut <- rbind(copyNumberOut, tmpCoord)
  }
  return(copyNumberOut)
}
