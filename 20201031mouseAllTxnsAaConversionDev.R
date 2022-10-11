### new dev script - for the multi-query 
### general idea is input species gene - I can first narrow down the subject txns by searching if the targeted amino acid is matched in its position
### this way I can get rid of non-canonical isoforms -> I can do the same for the homologous isoforms; search converted to see if any of them matchs
### the resulting pairs will then go under MSA, from this we figure out the best isoforms 



### notes: I need easier ways to exit out of these functions and produce next() on the outside loop 

# loadingTables -----------------------------------------------------------



library(liftOver)
library(GenomicRanges)
library(liftOver)
library(Biostrings)
library(msa)
library(stringr)

data("BLOSUM100")

hg38toMm10Chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/hg38ToMm10.over.chain")

hg38biomartTable <- read.table("/home/kevhu/data/20201030hg38KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

mm10biomartTable <- read.table("/home/kevhu/data/20201030Mm10KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

hg38Peptide <- read.table("/home/kevhu/data/20201030proteinHg.txt",  header = TRUE,
                          stringsAsFactors = FALSE, sep = "\t")

mm10Peptide <- read.table("/home/kevhu/data/20201030proteinMm.txt",  header = TRUE,
                          stringsAsFactors = FALSE, sep = "\t")


geneNameDf <- read.table("/home/kevhu/data/20201021geneNameBiomart.txt", header = TRUE,
                         stringsAsFactors = FALSE, sep = "\t")


small_cosmic <- read.table("/home/kevhu/data/20190225mouseConvertedOncogeneHotspots.txt",
                           sep = "\t", stringsAsFactors = FALSE, header = TRUE)


# firstUpper --------------------------------------------------------------


firstUpper <- function(gene){
  firstLetter <- toupper(substr(gene, start = 1, stop = 1))
  restOfGene <- tolower(substr(gene, start = 2, stop = nchar(gene)))
  res <- paste0(firstLetter, restOfGene)
  return(res)
}






# cdsConversion -----------------------------------------------------------

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
  lastExon <- max(df$rank)
  for (i in 1:nrow(df)) {
    if(df$rank[i] == 1){
      strand_cds_start_conv <- c(strand_cds_start_conv,
                                 df$exon_chrom_start_strand[i] - 1 + df$cdna_coding_start[i])
      strand_cds_end_conv <- c(strand_cds_end_conv, df$exon_chrom_end_strand[i])
    } else if(df$rank[i] == lastExon){
      strand_cds_start_conv <- c(strand_cds_start_conv, df$exon_chrom_start_strand[i])
      strand_cds_end_conv <- c(strand_cds_end_conv,
                               df$exon_chrom_start_strand[i] - 1 + df$cdna_coding_end[i])
    } else{
      strand_cds_start_conv <- c(strand_cds_start_conv, df$exon_chrom_start_strand[i])
      strand_cds_end_conv <- c(strand_cds_end_conv, df$exon_chrom_end_strand[i])
    }
  }
  df$strand_cds_start_conv<- strand_cds_start_conv
  df$strand_cds_end_conv <- strand_cds_end_conv
  return(df)
}



# multiSearchAaConversion -------------------------------------------------------------

aaToGenome <- function(aminoAcid, gene, position){
  aaCdsEnd <- position * 3
  aaCdsStart <- aaCdsEnd - 2
  aaRange <- IRanges(start = aaCdsStart, end = aaCdsEnd)
  
  tmpTable <- hg38biomartTable[which(hg38biomartTable$external_gene_name == gene),]
  if(nrow(tmpTable) == 0){
    print(paste0(gene," not in database"))
    return(NULL)
  }
  
  ### after pulling all isoforms, reduce them by checking for the aa position
  peptideDf <- hg38Peptide[which(hg38Peptide$ensembl_transcript_id %in% unique(tmpTable$ensembl_transcript_id)), ]
  tmpPepPositions <- substr(peptideDf$peptide, position, position)
  matchingTxnIds <- peptideDf$ensembl_transcript_id[which(tmpPepPositions == aminoAcid)]
  tmpTable2 <- tmpTable[which(tmpTable$ensembl_transcript_id %in% matchingTxnIds),]
  
  if (nrow(tmpTable2) == 0) {
    print("no correct isoforms")
    return(NULL)
  }
  
  isoformGenomeCoords <- NULL
  tmpEnsId <- unique(tmpTable2$ensembl_transcript_id)
  for (i in seq_along(tmpEnsId)) {
    tmpTable3 <- tmpTable2[which(tmpTable2$ensembl_transcript_id == tmpEnsId[i]),]
    tmpTable3 <- cdsConversion(tmpTable3)
    cdsRange <- IRanges(start = tmpTable3$cds_start, end = tmpTable3$cds_end)
    aaExon <- tmpTable3[subjectHits(findOverlaps(aaRange, cdsRange)),]
     
    if (nrow(aaExon) > 1) {
      #print(i)
      #print(aaExon)
      aaExon <- aaExon[1,]
    }
    if(aaExon$strand == 1){
      tmpPos <- aaExon$strand_cds_start_conv + (aaCdsStart - aaExon$cds_start) - 1
    } else{
      tmpPos <- aaExon$strand_cds_end_conv - (aaCdsEnd - aaExon$cds_end) - 1
    }
    
    genomeCoords <- data.frame("chrom" = aaExon$chromosome_name,"start" = tmpPos,
                         "end" = tmpPos + 2, "strand" = aaExon$strand, 
                         "ensembl" = aaExon$ensembl_transcript_id,
                         "gene" = aaExon$external_gene_name,
                         "original" = paste(gene, position, aminoAcid,sep = "|"),
                         stringsAsFactors = FALSE)
    
   isoformGenomeCoords <- rbind(isoformGenomeCoords, genomeCoords)
   }
  return(isoformGenomeCoords)
}


# liftOver ----------------------------------------------------------------

cdsLiftOver <- function(cdsToGenomeOut, chainFile){
  out <- liftOver(GRanges(seqnames = paste0("chr", unlist(cdsToGenomeOut$chrom)),
                          ranges = IRanges(start = unlist(cdsToGenomeOut$start), end = unlist(cdsToGenomeOut$end)),
                          strand = unlist(cdsToGenomeOut$strand)),
                  chain = chainFile)
  
  if (isEmpty(out)) {
    return(NULL)
  }
  out2 <- c("chrom" = as.character(out@unlistData@seqnames[1]),
            "start" = out@unlistData@ranges@start,
            "end" = out@unlistData@ranges@start + 2,
            "strand" = as.character(out@unlistData@strand@values),
            "gene" = unlist(cdsToGenomeOut$gene), 
            "ensembl" = cdsToGenomeOut$ensembl,
            "original" = cdsToGenomeOut$original)
  
  return(out2)
}




# genomeToAa --------------------------------------------------------------

genomeToAA <- function(liftOverRes){
  ### make this work for any species to any species at a later date
  tmpGeneName <- firstUpper(liftOverRes$gene)
  genomeCdsTmpTable <- mm10biomartTable[which(mm10biomartTable$external_gene_name %in% tmpGeneName),]
  if(nrow(genomeCdsTmpTable) == 0){
    tmpGeneName <- geneNameDf$mmusculus_homolog_associated_gene_name[which(geneNameDf$external_gene_name == liftOverRes$gene)]
    genomeCdsTmpTable <- mm10biomartTable[which(mm10biomartTable$external_gene_name %in% tmpGeneName),]
  }
  
  genomeCdsTmpTable <- cdsConversion(genomeCdsTmpTable)
  
  res <- NULL
  ### for this part I need to query each starting ens to every 
  for (i in seq_along(liftOverRes$ensembl)) {
    for (j in unique(genomeCdsTmpTable$ensembl)) {
      tmpCdsTable <- genomeCdsTmpTable[which(genomeCdsTmpTable$ensembl == j),]
      
      cdsRangeBp <- IRanges(start = liftOverRes$start[i], end = liftOverRes$end[i])
      cdsExonRangeBp <- IRanges(start = tmpCdsTable$exon_chrom_start, end = tmpCdsTable$exon_chrom_end)
      cdsExon <- tmpCdsTable[subjectHits(findOverlaps(cdsRangeBp, cdsExonRangeBp)),]
      tmpPos <- cdsExon$cds_start + abs(cdsExon$strand_cds_start_conv - mean(c(liftOverRes$start[i], liftOverRes$end[i])))
      #tmpPos <- cdsExon$cds_start + abs(cdsExon$strand_cds_start_conv - mean(c(liftOverRes$start[i], liftOverRes$end[i]))) - 1
      
      tmpAaPos <- ceiling(tmpPos/3)
      
      tmpRow <- c("gene" = liftOverRes$gene[i], "ensemblHg" = liftOverRes$ensembl[i],
                  "ensemblMm" = unique(tmpCdsTable$ensembl), "position" = tmpAaPos,
                  "original" = liftOverRes$original[i])
      
      if (length(tmpRow) != 5) {
        next()
      }
      res <- rbind(res, tmpRow)
    }
  }
  return(res)
}


# aaCheck -----------------------------------------------------------------

substrMiddlePep <- function(prt1, prt2, adjacentLength, aaPos1, aaPos2, aaArea){
  if (aaArea == "end") {
    sub_prt1 <- substr(prt1, aaPos1 - adjacentLength, aaPos1)
    sub_prt2 <- substr(prt2, aaPos2 - adjacentLength, aaPos2)
  } else if(aaArea == "start"){
    sub_prt1 <- substr(prt1, aaPos1, aaPos1 + adjacentLength)
    sub_prt2 <- substr(prt2, aaPos2, aaPos2 + adjacentLength)
  } else{
    sub_prt1 <- substr(prt1, aaPos1 - ceiling(adjacentLength/2), aaPos1 + ceiling(adjacentLength/2))
    sub_prt2 <- substr(prt2, aaPos2 - ceiling(adjacentLength/2), aaPos2 + ceiling(adjacentLength/2))
  }
  
  res <- list("sub_prt1" = sub_prt1, "sub_prt2" = sub_prt2)
  return(res)
}

aaConversionCheck <- function(singleRowDf){
  ### get proteins then subset based on the length of the shortest peptide - with a minimum of 10 probably
  convertedAaPos <- singleRowDf$position
  originalStrSplit <- unlist(str_split(singleRowDf$original, "\\|"))
  originalAaPos <- as.numeric(originalStrSplit[2])
  
  originalPeptide <- hg38Peptide$peptide[which(hg38Peptide$ensembl_transcript_id == singleRowDf$ensemblHg)]
  otherPeptide <- mm10Peptide$peptide[which(mm10Peptide$ensembl_transcript_id == singleRowDf$ensemblMm)]
  
  if (length(originalPeptide) < 1 | length(otherPeptide) < 1) {
    print(paste0("bad ", singleRowDf$original))
    return(NULL)
  }
  
  aaLength <- min(nchar(originalPeptide), nchar(otherPeptide))
  aaSurround <- aaLength * 0.025
  ### also need to check if they amino acids are on the ends/start of the peptide
  if (aaSurround < 10) {
    if((originalAaPos - 10) < 1){
      print("beginning < 10")
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10, 
                                     originalAaPos, convertedAaPos, "start")
    } else if ((originalAaPos + 10) < aaLength){
      print("end < 10")
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10, 
                                     originalAaPos, convertedAaPos, "end")
    } else{
      print("middle < 10")
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10,
                                     originalAaPos, convertedAaPos, "middle")
    }
  } else{
    if((originalAaPos - aaSurround) < 1){
      print("beginning > 10")
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                     originalAaPos, convertedAaPos, "start")
    } else if ((originalAaPos + aaSurround) < aaLength){
      print("end > 10")
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                     originalAaPos, convertedAaPos, "end")
    } else{
      print("middle > 10")
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                     originalAaPos, convertedAaPos, "middle")
    }
  }
  
  
  if (nchar(adjacentAas$sub_prt1[1]) == 0 | nchar(adjacentAas$sub_prt2[1]) == 0) {
    print(paste0("bad ", singleRowDf$original))
    return(NULL)
  }

  
  aaCheckOriginal <- substr(originalPeptide, originalAaPos, originalAaPos)
  aaCheckConverted <- substr(otherPeptide, convertedAaPos, convertedAaPos)
  matchingCheck <- ifelse(aaCheckConverted == aaCheckOriginal, "yes", "no")
  
  msaRes <- msa(AAStringSet(c(adjacentAas$sub_prt1[1], adjacentAas$sub_prt2[1])),
                type = "protein", substitutionMatrix = BLOSUM100);
  conScore <- mean(msaConservationScore(msaRes, substitutionMatrix = BLOSUM100))
  ConSeq <- msaConsensusSequence(msaRes)
  
  #print(msaRes)
  
  return(data.frame("seq1" = adjacentAas$sub_prt1, "seq2"=adjacentAas$sub_prt2, "conSeq"= ConSeq,
              "check" = matchingCheck, "conservationScore"=conScore, "original" = singleRowDf$original,
              "ensemblHg" = singleRowDf$ensemblHg, "ensemblMm"=singleRowDf$ensemblMm,
              stringsAsFactors = FALSE))
  
  #return(data.frame("seq1" = adjacentAas$sub_prt1, "seq2"=adjacentAas$sub_prt2, "conSeq"= ConSeq,
  #                  "check" = matchingCheck, stringsAsFactors = FALSE))
}



# testing -----------------------------------------------------------------



### (1) it should be simple enough to just loop everything from cds conversion to the the genomeAA
### (2) afeter the AA I should just do the isoform check I mentioned above
### (3) it'd be best to just do this for mouse at the same time


### everything below testing just the initial conversion

#test <- aaToGenome("R", "DNMT3A", 882)
#test <- data.frame(test, stringsAsFactors = FALSE)
#cdsLiftOver(test, hg38toMm10Chain)

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
hotspotDf <- data.frame(hotspotDf, stringsAsFactors = FALSE)

testDf <- NULL
for (i in 1:nrow(hotspotDf)) {
  tmpTest <- aaToGenome(hotspotDf$X2[i], hotspotDf$X1[i], as.numeric(hotspotDf$X3[i]))
  tmpTest <- data.frame(tmpTest, stringsAsFactors = FALSE)
  testDf <- rbind(testDf, tmpTest)
}
### should test everything up until the msa matching


### for msa matching and variant verfication - search all destination species isoforms 
### i.e genome to AA 


loOut <- NULL
for (i in unique(testDf$gene)) {
  tmpDf <- testDf[which(testDf$gene == i),]
  for (j in 1:nrow(tmpDf)) {
    tmpLoOut <- cdsLiftOver(tmpDf[j,], hg38toMm10Chain)
    if (is.null(tmpLoOut)) {
      print(paste0(i, " ",j))
      next()
    }
    loOut <- rbind(loOut, tmpLoOut)
  }
}

loOut <- data.frame(loOut, stringsAsFactors = FALSE)
loOut <- loOut[-which(!(loOut$strand == "+" | loOut$strand == "-")),]
loOut$start <- as.numeric(loOut$start)
loOut$end <- as.numeric(loOut$end)

genomeToAaDf <- NULL
for (i in unique(loOut$gene)) {
  tmpGeneDf <- loOut[which(loOut$gene == i),]
  tmpDf <- genomeToAA(tmpGeneDf)
  genomeToAaDf <- rbind(genomeToAaDf, tmpDf)
}

genomeToAaDf <- data.frame(genomeToAaDf, stringsAsFactors = FALSE)
genomeToAaDf$position <- as.numeric(genomeToAaDf$position)
rownames(genomeToAaDf) <- NULL


### below is for conversion check

checkDf <- NULL
for (i in 1:nrow(genomeToAaDf)) {
  tmpRow <- genomeToAaDf[i,]
  res <- aaConversionCheck(tmpRow)
  checkDf <- rbind(checkDf, res)
}


### percentage of original positions that were converted - this can be made into a graphic
yesAndNo <- unique(checkDf$original[which(checkDf$check == "yes")])
tmpCheckDf <- checkDf[-which(checkDf$original %in% yesAndNo),]
length(unique(tmpCheckDf$original))

length(unique(checkDf$original))/length(unique(small_cosmic$Hotspot.Residue))
(length(unique(checkDf$original)) - 16)/length(unique(small_cosmic$Hotspot.Residue))


### quick pie chart of the results or descriptive stats

length(unique(checkDf$original))
length(unique(small_cosmic$Hotspot.Residue))
slices <- c(13, 16, 120) # 13 no isoform with matching aa, 16 no conversion, 120 good
lbls <- c("no matching amino acid", "no conversion",
          "converted")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210311aaConversionStats.pdf", useDingbats = FALSE)
pie(slices, labels = lbls, main="Amino Acid Conversion Stats")
dev.off()



badCheck <- checkDf[which(checkDf$check == "no"), ]
summary(checkDf$conservationScore[which(checkDf$check == "yes")])
boxplot(checkDf$conservationScore[which(checkDf$check == "yes")])

genomeToAaDf_braf <- genomeToAaDf[which(genomeToAaDf$gene == "BRAF"),]
tmpCheck <- NULL
for (i in 1:nrow(genomeToAaDf_braf)) {
  tmpRow <- genomeToAaDf_braf[i,]
  res <- aaConversionCheck(tmpRow)
  tmpCheck <- rbind(tmpCheck, res)
}



###202100312 plot for conservation scores from converted

dfToGraph <- checkDf
dfToGraph$mLength <- nchar(checkDf$seq1)
### this simple plot shows there is no bias towards length of the amino acid matched
### and the conservation scores
plot(dfToGraph$mLength, dfToGraph$conservationScore)
###
dfToGraph$check <- str_replace(dfToGraph$check, "yes", "1")
dfToGraph$check <- str_replace(dfToGraph$check, "no", "0")
dfToGraph$check  <- as.numeric(dfToGraph$check)

dfToGraph <- dfToGraph[order(dfToGraph$check, dfToGraph$conservationScore,
                             decreasing = TRUE),]
dfToGraph2 <- dfToGraph[-which(duplicated(dfToGraph$original)),]

ggplot(data = dfToGraph2, aes(x = conservationScore, y = mLength,
                              color = relevel(factor(check), ref = "1"))) +
  geom_point(size = 1)

