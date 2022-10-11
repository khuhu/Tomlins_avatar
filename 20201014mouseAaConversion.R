
### parts of AA conversion algo. - (1) AA position to genomic (2) liftOVer 
### need to make it backwards compatible too i.e liftOver to human coordinates then 


# loading -----------------------------------------------------------------



library(liftOver)
library(GenomicRanges)
library(liftOver)
library(Biostrings)
library(msa)

data("BLOSUM100")

h38toMm10Chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/hg38ToMm10.over.chain")

hg38biomartTable <- read.table("/home/kevhu/data/20201020hg38KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

mm10biomartTable <- read.table("/home/kevhu/data/20201020Mm10KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

#hg38Cdna <- read.table("/home/kevhu/data/20201021HgUcscCdna.txt", header = TRUE,
#                       stringsAsFactors = FALSE, sep = "\t")

#mm10Cdna <- read.table("/home/kevhu/data/20201021MmUcscCdna.txt", header = TRUE,
#                       stringsAsFactors = FALSE, sep = "\t")

geneNameDf <- read.table("/home/kevhu/data/20201021geneNameBiomart.txt", header = TRUE,
                         stringsAsFactors = FALSE, sep = "\t")

hg38Peptide <- read.table("/home/kevhu/data/20201022proteinHg.txt",  header = TRUE,
                          stringsAsFactors = FALSE, sep = "\t")

mm10Peptide <- read.table("/home/kevhu/data/20201022proteinMm.txt",  header = TRUE,
                          stringsAsFactors = FALSE, sep = "\t")


### converts cds to genomic coordinates, takes strand into account
### normalizes the genomic coordinates to where the coding within the exon starts
### only really important for first and last exon


# firstUpper --------------------------------------------------------------


firstUpper <- function(gene){
  firstLetter <- toupper(substr(gene, start = 1, stop = 1))
  restOfGene <- tolower(substr(gene, start = 2, stop = nchar(gene)))
  res <- paste0(firstLetter, restOfGene)
  return(res)
}



# cdsConversion -----------------------------------------------------------


cdsConversion <- function(df){
  if(nrow(df) == 0){
    return("empty")
  }
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



# aaToGenome --------------------------------------------------------------
### converts position of amino acid for canonical transcript into genomic coordinates
aaToGenome <- function(gene, position){
  tmpTable <- hg38biomartTable[which(hg38biomartTable$external_gene_name == gene),]
  if(nrow(tmpTable) == 0){
    print(paste0(gene," not in database"))
    break()
  }
  tmpTable <- cdsConversion(tmpTable)
  aaCdsEnd <- position * 3
  aaCdsStart <- aaCdsEnd - 2
  aaRange <- IRanges(start = aaCdsStart, end = aaCdsEnd)
  cdsRange <- IRanges(start = tmpTable$cds_start, end = tmpTable$cds_end)
  aaExon <- tmpTable[subjectHits(findOverlaps(aaRange, cdsRange)),]
  if(aaExon$strand == 1){
    tmpPos <- aaExon$strand_cds_start_conv + (aaCdsStart - aaExon$cds_start) - 1
  } else{
    tmpPos <- aaExon$strand_cds_end_conv - (aaCdsEnd - aaExon$cds_end) - 1
  }
  
  genomeCoords <- list("chrom" = aaExon$chromosome_name,"start" = tmpPos,
                       "end" = tmpPos + 2, "strand" = aaExon$strand, 
                       "ensembl" = aaExon$ensembl_transcript_id,
                       "cdna_coding_start" = aaExon$cdna_coding_start,
                       "cdna_coding_end" = aaExon$cdna_coding_end,
                       "gene" = aaExon$external_gene_name)
  return(genomeCoords)
}



# genomeToAa --------------------------------------------------------------

### this should be simpler than cdsToGenome because I only care about position, and not strand
genomeToAA <- function(liftOverRes){
  ### make this work for any species to any species at a later date
  ### revised to search for lowercase first then table of genes
  tmpGeneName <- firstUpper(liftOverRes$gene)
  genomeCdsTmpTable <- mm10biomartTable[which(mm10biomartTable$external_gene_name == tmpGeneName),]
  if(nrow(genomeCdsTmpTable) == 0){
    tmpGeneName <- geneNameDf$mmusculus_homolog_associated_gene_name[which(geneNameDf$external_gene_name == liftOverRes$gene)]
    genomeCdsTmpTable <- mm10biomartTable[which(mm10biomartTable$external_gene_name == tmpGeneName),]
  }
  
  genomeCdsTmpTable <- cdsConversion(genomeCdsTmpTable)
  cdsRangeBp <- IRanges(start = liftOverRes$start, end = liftOverRes$end)
  cdsExonRangeBp <- IRanges(start = genomeCdsTmpTable$exon_chrom_start, end = genomeCdsTmpTable$exon_chrom_end)
  cdsExon <- genomeCdsTmpTable[subjectHits(findOverlaps(cdsRangeBp, cdsExonRangeBp)),]
  
  ### the subtraction of end/start depends on strand
  tmpPos <- cdsExon$cds_start + abs(cdsExon$strand_cds_start_conv - mean(c(liftOverRes$start, liftOverRes$end)))
  tmpAaPos <- ceiling(tmpPos/3)
  
  ### might need to return more
  return(tmpAaPos)
}




# cdsLiftOver -------------------------------------------------------------


### takes genomic coordinates of aa and liftovers to mouse coordinates
### output goes into function to check + final table
cdsLiftOver <- function(cdsToGenomeOut, chainFile){
  out <- liftOver(GRanges(seqnames = paste0("chr", unlist(cdsToGenomeOut$chrom)),
                          ranges = IRanges(start = unlist(cdsToGenomeOut$start), end = unlist(cdsToGenomeOut$end)),
                          strand = unlist(cdsToGenomeOut$strand)),
                  chain = chainFile)
  
  if (isEmpty(out)) {
    print("empty liftOver")
    return(NULL)
  }
  
  out2 <- list("chrom" = as.character(out@unlistData@seqnames[1]),
               "start" = out@unlistData@ranges@start,
               "end" = out@unlistData@ranges@start + 2,
               "strand" = as.character(out@unlistData@strand@values),
               "gene" = unlist(cdsToGenomeOut$gene))
  
  return(out2)
}



# aaCheck -----------------------------------------------------------------



substrMiddlePep <- function(prt1, prt2, adjacentLength, aaPos1, aaPos2, aaArea){
  if (aaArea == "end") {
    sub_prt1 <- substr(prt1, (aaPos1 - adjacentLength), aaPos1)
    sub_prt2 <- substr(prt2, (aaPos2 - adjacentLength), aaPos2)
  } else if(aaArea == "start"){
    sub_prt1 <- substr(prt1, aaPos1, (aaPos1 + adjacentLength))
    sub_prt2 <- substr(prt2, aaPos2, (aaPos2 + adjacentLength))
  } else if(aaArea == "middle"){
    sub_prt1 <- substr(prt1, (aaPos1 - adjacentLength/2), (aaPos1 + adjacentLength/2))
    sub_prt2 <- substr(prt2, (aaPos2 - adjacentLength/2), (aaPos2 + adjacentLength/2))
  }

  res <- list("sub_prt1" = sub_prt1, "sub_prt2" = sub_prt2)
  return(res)
}



### input is a list of two things:(1) original gene + aa position (2) the liftOver output
aaConversionCheck <- function(coordinatesList){
  ### get proteins then subset based on the length of the shortest peptide - with a minimum of 10 probably
  convertedAaPos <- genomeToAA(coordinatesList[[3]])
  originalAaPos <- coordinatesList[[2]]
  originalGene <- coordinatesList[[1]]
  
  newGeneName <- firstUpper(originalGene)
  otherPeptide <- mm10Peptide$peptide[which(mm10Peptide$external_gene_name == newGeneName)]
  if(length(otherPeptide) == 0){
    newGeneName <- geneNameDf$mmusculus_homolog_associated_gene_name[which(geneNameDf$external_gene_name == originalGene)]
    otherPeptide <- mm10Peptide$peptide[which(mm10Peptide$external_gene_name == newGeneName)]
  }
  
  originalPeptide <- hg38Peptide$peptide[which(hg38Peptide$external_gene_name == originalGene)]
  
  aaLength <- min(nchar(originalPeptide), nchar(otherPeptide))
  aaSurround <- floor(aaLength * 0.01)
  
  print(aaLength)
  print(originalAaPos)
  print(convertedAaPos)
  
  ### also need to check if they amino acids are on the ends/start of the peptide
  ### puts a minimum of 10 surrouding peptides
  if (aaSurround < 10) {
    if((originalAaPos - 10) < 1){
      print("start")
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10, 
                                     originalAaPos, convertedAaPos, "start")
    } else if ((originalAaPos + 10) > aaLength){
      print("end")
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10, 
                                     originalAaPos, convertedAaPos, "end")
    } else{
      print("middle")
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10,
                                     originalAaPos, convertedAaPos, "middle")
    }
  } else{
    if((originalAaPos - aaSurround) < 1){
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                     originalAaPos, convertedAaPos, "start")
    } else if ((originalAaPos + aaSurround) > aaLength){
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                     originalAaPos, convertedAaPos, "end")
    } else{
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                     originalAaPos, convertedAaPos, "middle")
    }
  }
  
  
  ### now do a quick msa + aa check. for msa use the strict table again
  ### do a quick dry run below, and then make new document to source functions
  ### this current document will serve as a Dev version
  
  aaCheckOriginal <- substr(originalPeptide, originalAaPos,originalAaPos)
  aaCheckConverted <- substr(otherPeptide, convertedAaPos, convertedAaPos)
  matchingCheck <- ifelse(aaCheckConverted == aaCheckOriginal, "yes", "no")
  
  print(aaCheckOriginal)
  print(aaCheckConverted)
  
  msaRes <- msa(AAStringSet(c(adjacentAas$sub_prt1,adjacentAas$sub_prt2)), type = "protein")
  conScore <- msaConservationScore(msaRes, substitutionMatrix = BLOSUM100)
  ConSeq <- msaConsensusSequence(msaRes)
  
  return(list("seq1" = adjacentAas$sub_prt1,"seq2"=adjacentAas$sub_prt1, "conSeq"=ConSeq,
              "check" = matchingCheck, "conservationScore"=conScore))
}


# randomTesting -----------------------------------------------------------



### check to see if cdsToGenome works + liftOver
#brafV600Conversion <- aaToGenome(gene = "BRAF", position = 600)
#tmpLiftOver <- liftOver(GRanges(seqnames = paste0("chr", unlist(brafV600Conversion$chrom)),
#                 ranges = IRanges(start = unlist(brafV600Conversion$start), end = unlist(brafV600Conversion$end)),
#                 strand = unlist(brafV600Conversion$strand)),
#         chain = h38toMm10Chain)

### test set for BRAF V600 E
### (1) convert aa position into cds position (2) interset the cds position with the gene's 
### (3) depending on the strand subtract from either start/end (4) get genomic coordinates
### (5) use liftover 

### 1 (1-3) 2 (4-6) 3 (7-9)

aaPos <- 882
tmpTable <- hg38biomartTable[which(hg38biomartTable$external_gene_name == "DNMT3A"),]
tmpTable <- cdsConversion(tmpTable)
aaCdsEnd <- aaPos * 3
aaCdsStart <- aaCdsEnd - 2
aaRange <- IRanges(start = aaCdsStart, end = aaCdsEnd)
cdsRange <- IRanges(start = tmpTable$cds_start, end = tmpTable$cds_end)
aaExon <- tmpTable[subjectHits(findOverlaps(aaRange, cdsRange)),]

### I guess doing it this way if strand is positive, add the difference of aaPos and cds_start to the genomic_start
### If the strand is negative subtract the difference of cds_end and aaPos from genomic_end
### idea is the tiling the gene forward or backward based on +/- strand

if(aaExon$strand == 1){
  tmpPos <- aaExon$strand_cds_start_conv + (aaCdsStart - aaExon$cds_start) - 1
} else{
  tmpPos <- aaExon$strand_cds_end_conv - (aaCdsEnd - aaExon$cds_end) - 1
}

tmpPos  


### testing out the aaMatching - this is odd becuase to do this properly I need have
### the cdna positions for both mouse and humans - area I would want to look at
### easier for the starting species, since I already have it - (look at cds_start/end)
### since I'm only looking at cdna and not the exons

### for mice ... I have the converted genomic coordinates, so I need to convert that into
### the cds_start/end - might be easiest to just use peptide themselves

tmpGeneName <- geneNameDf$mmusculus_homolog_associated_gene_name[which(geneNameDf$external_gene_name == liftOverOut$gene)]
genomeCdsTmpTable <- mm10biomartTable[which(mm10biomartTable$external_gene_name == tmpGeneName),]
genomeCdsTmpTable <- cdsConversion(genomeCdsTmpTable)
otherProtein <- mm10Peptide$peptide[which(mm10Peptide$external_gene_name == tmpGeneName)]

cdsRangeBp <- IRanges(start = liftOverOut$start, end = liftOverOut$end)
cdsExonRangeBp <- IRanges(start = genomeCdsTmpTable$exon_chrom_start, end = genomeCdsTmpTable$exon_chrom_end)
cdsExon <- genomeCdsTmpTable[subjectHits(findOverlaps(cdsRangeBp, cdsExonRangeBp)),]

### the subtraction of end/start depends on strand
tmpPos <- cdsExon$cds_start + abs(cdsExon$strand_cds_start_conv - mean(c(liftOverOut$start, liftOverOut$end))) - 1
tmpAaPos <- ceiling(tmpPos/3)
substr(otherProtein, tmpAaPos - 2, tmpAaPos + 2)

### aa-check
###
###

brafV600Conversion <- aaToGenome(gene = "DNMT3A", position = 882)
liftOverOut <- cdsLiftOver(brafV600Conversion, chainFile = h38toMm10Chain)
newGeneName <- geneNameDf$mmusculus_homolog_associated_gene_name[which(geneNameDf$external_gene_name == checkInput[[3]]$gene)]


### need to create a function to internally test the checkInput section
checkInput <- list("originalGene" = "GNAQ","originalPos" = 209, "liftOverOut" = liftOverOut)
convertedAaPos <- genomeToAA(checkInput[[3]])
originalAaPos <- checkInput[[2]]
originalGene <- checkInput[[1]]


originalPeptide <- hg38Peptide$peptide[which(hg38Peptide$external_gene_name == "PIK3CA")]
otherPeptide <- mm10Peptide$peptide[which(mm10Peptide$external_gene_name == "Pik3ca")]

aaLength <- min(nchar(originalPeptide), nchar(otherPeptide))
aaSurround <- aaLength * 0.05
### also need to check if they amino acids are on the ends/start of the peptide
if (aaSurround < 10) {
  if((originalAaPos - 10) < 1){
    adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10, 
                                   originalAaPos, convertedAaPos, "start")
  } else if ((originalAaPos + 10) < aaLength){
    adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10, 
                                   originalAaPos, convertedAaPos, "end")
  } else{
    adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, 10,
                                   originalAaPos, convertedAaPos, "middle")
  }
} else{
  if((originalAaPos - aaSurround) < 1){
    adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                   originalAaPos, convertedAaPos, "start")
  } else if ((originalAaPos + aaSurround) < aaLength){
    adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                   originalAaPos, convertedAaPos, "end")
  } else{
    adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide, aaSurround,
                                   originalAaPos, convertedAaPos, "middle")
  }
}

### now do a quick msa + aa check. for msa use the strict table again
### do a quick dry run below, and then make new document to source functions
### this current document will serve as a Dev version

aaCheckOriginal <- substr(originalPeptide, originalAaPos,originalAaPos)
aaCheckConverted <- substr(otherPeptide, convertedAaPos, convertedAaPos)
matchingCheck <- ifelse(aaCheckConverted == aaCheckOriginal, "yes", "no")

msaRes <- msa(AAStringSet(c(adjacentAas$sub_prt1,adjacentAas$sub_prt2)), type = "protein")
conScore <- msaConservationScore(msaRes, substitutionMatrix = BLOSUM100)
ConSeq <- msaConsensusSequence(msaRes)
