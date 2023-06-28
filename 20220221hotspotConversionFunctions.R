library(liftOver)
library(GenomicRanges)
library(liftOver)
library(Biostrings)
library(msa)

data("BLOSUM100")

vtml10Table <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/VTML/VTML_10.mat",skip = 10, sep = "" , header = TRUE,
                          na.strings ="", stringsAsFactors= FALSE, check.names = FALSE)

VTML10 <- as.matrix(vtml10Table)
### premade db files

h38toMm10Chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/hg38ToMm10.over.chain")

hg19toHg38Chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/hg19ToHg38.over.chain")


hg38biomartTable <- read.table("/home/kevhu/data/20230320hg38KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

mm10biomartTable <- read.table("/home/kevhu/data/20210203Mm10KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")


geneNameDf <- read.table("/home/kevhu/data/20230315geneNameBiomart.txt", header = TRUE,
                         stringsAsFactors = FALSE, sep = "\t")

hg38Peptide <- read.table("/mnt/DATA5/tmp/kev/misc/20230301allPeptideQuery.txt",  header = TRUE,
                          stringsAsFactors = FALSE, sep = "\t")

mm10Peptide <- read.table("/mnt/DATA5/tmp/kev/misc/20230301allPeptideQueryMm10.txt",  header = TRUE,
                          stringsAsFactors = FALSE, sep = "\t")



### how new gene name table was queried
# library(biomaRt)
# 
# hg38biomartTable <- read.table("/home/kevhu/data/20201030hg38KnownCanbiomartQuery.txt", header = TRUE,
#                                stringsAsFactors = FALSE, sep = "\t")
# 
# ensemblHg <- useMart("ensembl", dataset= "hsapiens_gene_ensembl")
# 
# geneName <- getBM(ensemblHg, attributes = attrHg_geneName, filters = "ensembl_transcript_id",
#                   values = hg38biomartTable$ensembl_transcript_id)
# 
# 
# write.table(geneName, "/home/kevhu/data/20220221geneNameBiomart.txt",
#             quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)

#df <- genomeCdsTmpTable[which(genomeCdsTmpTable$ensembl_transcript_id == "ENST00000275493"), ]

firstUpper <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

df <- tmpTable

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
    print(paste(i, length(strand_cds_start_conv), length(strand_cds_end_conv)))
  }
  df$strand_cds_start_conv<- strand_cds_start_conv
  df$strand_cds_end_conv <- strand_cds_end_conv
  return(df)
}

# cdsToGenomeOut <- testAaConversion[j, ]
# chainFile <-  h38toMm10Chain

cdsLiftOver <- function(cdsToGenomeOut, chainFile){
  out <- liftOver(GRanges(seqnames = paste0("chr", unlist(cdsToGenomeOut$chrom)),
                          ranges = IRanges(start = unlist(cdsToGenomeOut$start), end = unlist(cdsToGenomeOut$end)),
                          strand = unlist(cdsToGenomeOut$strand)),
                  chain = chainFile)
  
  if (isEmpty(out)) {
    return()
  }
  
  out2 <- data.frame("chrom" = as.character(out@unlistData@seqnames[1]),
               "start" = out@unlistData@ranges@start,
               "end" = out@unlistData@ranges@start + 2,
               "strand" = as.character(out@unlistData@strand@values),
               "gene" = unlist(cdsToGenomeOut$gene), stringsAsFactors = FALSE)
  
  return(out2)
}


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



### testing with BRAFV600 because multiple hits on ensids, can test NRAS after b/c only 1 transcript
### above functions unaffected
# i <- 4
# gene <- hotspotDf$Gene[i]
# position <- as.numeric(hotspotDf$Position[i])

aaToGenome <- function(gene, position){
  tmpTable <- hg38biomartTable[which(hg38biomartTable$external_gene_name == gene),]
  if (any(is.na(tmpTable$cds_start))) {
    tmpTable <- tmpTable[-which(is.na(tmpTable$cds_start)), ]
  }
  tmpTable <- cdsConversion(tmpTable)
  aaCdsEnd <- position * 3
  aaCdsStart <- aaCdsEnd - 2
  aaRange <- IRanges(start = aaCdsStart, end = aaCdsEnd)
  cdsRange <- IRanges(start = tmpTable$cds_start, end = tmpTable$cds_end)
  aaExon <- tmpTable[subjectHits(findOverlaps(aaRange, cdsRange)),]
  tmpPos <- NULL
  for (z in 1:nrow(aaExon)) {
    if(aaExon$strand[z] == 1){
      tmpPos <- c(tmpPos, aaExon$strand_cds_start_conv[z] + (aaCdsStart - aaExon$cds_start[z]) - 1)
    } else{
      tmpPos <- c(tmpPos, aaExon$strand_cds_end_conv[z] - (aaCdsEnd - aaExon$cds_end[z]) - 1)
    }
  }

  
  genomeCoords <- data.frame("chrom" = aaExon$chromosome_name,"start" = tmpPos,
                       "end" = tmpPos + 2, "strand" = aaExon$strand, 
                       "ensembl" = aaExon$ensembl_transcript_id,
                       "cdna_coding_start" = aaExon$cdna_coding_start,
                       "cdna_coding_end" = aaExon$cdna_coding_end,
                       "gene" = aaExon$external_gene_name, stringsAsFactors = FALSE)
  
  #genomeCoords <- GRanges(seqnames = paste0("chr", aaExon$chromosome_name),
  #                        ranges = IRanges(start = tmpPos, end = tmpPos + 2),
  #                        strand = aaExon$strand)
  
  return(genomeCoords)
}


### slight alteration for matching gene names

# liftOverRes <- coordinatesList[[3]]
genomeToAA <- function(liftOverRes){
  ### make this work for any species to any species at a later date
  tmpGeneName <- geneNameDf$mmusculus_homolog_associated_gene_name[grep(paste0("^", liftOverRes$gene, "$"), geneNameDf$external_gene_name)]
  ### stop gap for no match
  if (any(tmpGeneName == "")) {
    tmpGeneName <- unique(mm10biomartTable$external_gene_name)[grep(firstUpper(liftOverRes$gene), unique(mm10biomartTable$external_gene_name))]
    if (any(tmpGeneName == "") | isEmpty(tmpGeneName)){
      return(return(data.frame("ENSID" = NA, "Position" = NA, stringsAsFactors = FALSE)))
    } else{
      genomeCdsTmpTable <- mm10biomartTable[which(mm10biomartTable$external_gene_name %in% tmpGeneName),]
    }
  } else{
    genomeCdsTmpTable <- mm10biomartTable[which(mm10biomartTable$external_gene_name %in% tmpGeneName),]
  }
  genomeCdsTmpTable <- mm10biomartTable[which(mm10biomartTable$external_gene_name %in% tmpGeneName),]
  genomeCdsTmpTable <- cdsConversion(genomeCdsTmpTable)
  
  cdsRangeBp <- IRanges(start = liftOverRes$start, end = liftOverRes$end)
  cdsExonRangeBp <- IRanges(start = genomeCdsTmpTable$exon_chrom_start, end = genomeCdsTmpTable$exon_chrom_end)
  cdsExon <- genomeCdsTmpTable[subjectHits(findOverlaps(cdsRangeBp, cdsExonRangeBp)),]
  
  ### reducing to homolgous transcript of mane

  ### the subtraction of end/start depends on strand
  ### 20230313 I don't know why I subtracted 1 before, but when I do human back and forth, it off by 1
  ### test ran with bottom code using all EGFR to test 
  # tmpPos <- cdsExon$cds_start + abs(cdsExon$strand_cds_start_conv - mean(c(liftOverRes$start, liftOverRes$end))) - 1
  tmpPos <- cdsExon$cds_start + abs(cdsExon$strand_cds_start_conv - mean(c(liftOverRes$start, liftOverRes$end)))
  tmpAaPos <- ceiling(tmpPos/3)
  
  ### might need to return more
  return(data.frame("ENSID" = cdsExon$ensembl_transcript_id, "Position" = tmpAaPos, stringsAsFactors = FALSE))
}

### used to test backwartds compat with human itself
# liftOverRes <- testAaConversion[1,]
genomeToAAHg <- function(liftOverRes){
  ### make this work for any species to any species at a later date
  
  tmpGeneName <- liftOverRes$gene
  genomeCdsTmpTable <- hg38biomartTable[which(hg38biomartTable$external_gene_name %in% tmpGeneName),]
  genomeCdsTmpTable <- cdsConversion(genomeCdsTmpTable)
  
  cdsRangeBp <- IRanges(start = liftOverRes$start, end = liftOverRes$end)
  cdsExonRangeBp <- IRanges(start = genomeCdsTmpTable$exon_chrom_start, end = genomeCdsTmpTable$exon_chrom_end)
  cdsExon <- genomeCdsTmpTable[subjectHits(findOverlaps(cdsRangeBp, cdsExonRangeBp)),]
  cdsExon <- cdsExon[which(cdsExon$ensembl_transcript_id == liftOverRes$ensembl),]
  
  ### the subtraction of end/start depends on strand
  ### 20230313 I don't know why I subtracted 1 before, but when I do human back and forth, it off by 1
  # tmpPos <- cdsExon$cds_start + abs(cdsExon$strand_cds_start_conv - mean(c(liftOverRes$start, liftOverRes$end))) - 1
  tmpPos <- cdsExon$cds_start + abs(cdsExon$strand_cds_start_conv - mean(c(liftOverRes$start, liftOverRes$end)))
  tmpAaPos <- ceiling(tmpPos/3)
  
  ### might need to return more
  return(tmpAaPos)
}

# coordinatesList <- checkInput

aaConversionCheck <- function(coordinatesList){
  msaDf <- NULL
  ### get proteins then subset based on the length of the shortest peptide - with a minimum of 10 probably
  convertedAaPos <- genomeToAA(coordinatesList[[3]])
  if(any(is.na(convertedAaPos$ENSID))){
    convertedAaPos <- convertedAaPos[-which(is.na(convertedAaPos$ENSID)),]
  }
  
  ### for now won't implement this because it should return a list of two dfs - one for homologs and one for all search
  convertedAaPos <- tryCatch({convertedAaPos[which(convertedAaPos$ENSID == mm10Peptide$ensembl_transcript_id[which(mm10Peptide$ensembl_peptide_id == coordinatesList$liftOverOut$homo_pep)]),]},
                             warning = function(x){
                               return(NULL)
                             })
  
  
  if (is.null(convertedAaPos)) {
    msaDf <- rbind(msaDf,
                   data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                              "check" = NA, "conservationScore"= NA,
                              "originalEns" = originalEns, "ENSID" = NA,
                              "Hotspot" = NA,
                              "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                              "liftOverChrom" = coordinatesList$liftOverOut$chrom,
                              "liftOverStart" = coordinatesList$liftOverOut$start,
                              "liftOverEnd" = coordinatesList$liftOverOut$end,
                              "annoMane" = "", stringsAsFactors = FALSE))
    return(msaDf)
  } else if (nrow(convertedAaPos) == 0) {
    msaDf <- rbind(msaDf,
                   data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                              "check" = NA, "conservationScore"= NA,
                              "originalEns" = originalEns, "ENSID" = NA,
                              "Hotspot" = NA,
                              "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                              "liftOverChrom" = coordinatesList$liftOverOut$chrom,
                              "liftOverStart" = coordinatesList$liftOverOut$start,
                              "liftOverEnd" = coordinatesList$liftOverOut$end,
                              "annoMane" = "", stringsAsFactors = FALSE))
    return(msaDf)
  }
  originalAaPos <- coordinatesList[[2]]
  originalGene <- coordinatesList[[1]]
  
  newGeneName <- unique(geneNameDf$mmusculus_homolog_associated_gene_name[grep(paste0("^",checkInput[[3]]$gene, "$"),
                                                                               geneNameDf$external_gene_name)])

  if (newGeneName == "") {
    newGeneName <- unique(geneNameDf$external_gene_name[grep(paste0("^", toupper(checkInput[[3]]$gene), "$"),
                                                             geneNameDf$external_gene_name)])
  }

  
  originalEns <- coordinatesList$liftOverOut$ensembl
  originalPeptide <- hg38Peptide$peptide[which(hg38Peptide$ensembl_transcript_id == originalEns)]
  originalManeAnno <- hg38Peptide$transcript_mane_select[which(hg38Peptide$ensembl_transcript_id == originalEns)]
  if (length(originalPeptide) == 0) {
    return()
  }
  originalAa <- coordinatesList$liftOverOut$aa
  
  
  ### 20230302: b/c of multiple ens homoglous genes names ... used %in% instead of ==
  ### 20230315: switched it to match MANE protein's homologous peptide 
  
  # otherPeptide <- mm10Peptide$peptide[which(mm10Peptide$external_gene_name %in% newGeneName)]
  # ensPeptide <- mm10Peptide$ensembl_transcript_id[which(mm10Peptide$external_gene_name %in% newGeneName)]
  # 
  # otherPeptide <- otherPeptide[match(convertedAaPos$ENSID, ensPeptide)]
  # ensPeptide <- ensPeptide[match(convertedAaPos$ENSID, ensPeptide)]
  # 
  # if (any(otherPeptide == "Sequence unavailable")) {
  #   ensPeptide <- ensPeptide[-which(otherPeptide == "Sequence unavailable")]
  #   otherPeptide <- otherPeptide[-which(otherPeptide == "Sequence unavailable")]
  # }
  
  
  ### match it to the homolog peptide and if not available do above with all peptides i guess
  if (coordinatesList$liftOverOut$homo_pep == "") {
    
    otherPeptide <- mm10Peptide$peptide[which(mm10Peptide$external_gene_name %in% firstUpper(newGeneName))]
    ensPeptide <- mm10Peptide$ensembl_transcript_id[which(mm10Peptide$external_gene_name %in% firstUpper(newGeneName))]
    
    otherPeptide <- otherPeptide[match(convertedAaPos$ENSID, ensPeptide)]
    ensPeptide <- ensPeptide[match(convertedAaPos$ENSID, ensPeptide)]
    
    if (any(otherPeptide == "Sequence unavailable")) {
      ensPeptide <- ensPeptide[-which(otherPeptide == "Sequence unavailable")]
      otherPeptide <- otherPeptide[-which(otherPeptide == "Sequence unavailable")]
    }
    
  } else{
    
    ensPeptide <- mm10Peptide$ensembl_transcript_id[which(mm10Peptide$ensembl_peptide_id == coordinatesList$liftOverOut$homo_pep)]
    otherPeptide <- mm10Peptide$peptide[which(mm10Peptide$ensembl_peptide_id == coordinatesList$liftOverOut$homo_pep)]
    
  }
  

  
  if(length(ensPeptide) == 0){
    msaDf <- rbind(msaDf,
                   data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                              "check" = NA, "conservationScore"= NA,
                              "originalEns" = originalEns, "ENSID" = NA,
                              "Hotspot" = paste0(originalGene, originalAaPos, originalAa),
                              "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                              "liftOverChrom" = coordinatesList$liftOverOut$chrom,
                              "liftOverStart" = coordinatesList$liftOverOut$start,
                              "liftOverEnd" = coordinatesList$liftOverOut$end,
                              "annoMane" = "", stringsAsFactors = FALSE))
    return(msaDf)
  } else if (any(otherPeptide == "Sequence unavailable")) {
    ensPeptide <- ensPeptide[-which(otherPeptide == "Sequence unavailable")]
    otherPeptide <- otherPeptide[-which(otherPeptide == "Sequence unavailable")]
    if (length(ensPeptide) == 0) {
      msaDf <- rbind(msaDf,
                     data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                "check" = NA, "conservationScore"= NA,
                                "originalEns" = originalEns, "ENSID" = NA,
                                "Hotspot" = paste0(originalGene, originalAaPos, originalAa),
                                "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                "liftOverChrom" = coordinatesList$liftOverOut$chrom,
                                "liftOverStart" = coordinatesList$liftOverOut$start,
                                "liftOverEnd" = coordinatesList$liftOverOut$end,
                                "annoMane" = "", stringsAsFactors = FALSE))
      return(msaDf)
    }
  }


  
  for (z in 1:length(ensPeptide)) {
    print(z)
    
    aaLength <- min(nchar(originalPeptide), nchar(otherPeptide[z]))
    aaSurround <- aaLength * 0.05
    
    ### if aa isoform too short
    if (convertedAaPos$Position[z] > aaLength) {
      msaDf <- rbind(msaDf,
                     data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                "check" = NA, "conservationScore"= NA,
                                "originalEns" = originalEns, "ENSID" = ensPeptide[z],
                                "Hotspot" = paste0(originalGene, originalAaPos, originalAa),
                                "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                "liftOverChrom" = coordinatesList$liftOverOut$chrom,
                                "liftOverStart" = coordinatesList$liftOverOut$start,
                                "liftOverEnd" = coordinatesList$liftOverOut$end,
                                "annoMane" = "", stringsAsFactors = FALSE))
      next()
    }
    
    ### also need to check if they amino acids are on the ends/start of the peptide
    ### 20230313 trying only 15 amino acid
    # if (aaSurround < 15) {
    #   if((originalAaPos - 15) < 1){
    #     adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide[z], 15, 
    #                                    originalAaPos, convertedAaPos$Position[z], "start")
    #   } else if ((originalAaPos + 15) < aaLength){
    #     adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide[z], 15, 
    #                                    originalAaPos, convertedAaPos$Position[z], "end")
    #   } else{
    #     adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide[z], 15,
    #                                    originalAaPos, convertedAaPos$Position[z], "middle")
    #   }
    # } else{
    #   if((originalAaPos - aaSurround) < 1){
    #     adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide[z], aaSurround,
    #                                    originalAaPos, convertedAaPos$Position[z], "start")
    #   } else if ((originalAaPos + aaSurround) < aaLength){
    #     adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide[z], aaSurround,
    #                                    originalAaPos, convertedAaPos$Position[z], "end")
    #   } else{
    #     adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide[z], aaSurround,
    #                                    originalAaPos, convertedAaPos$Position[z], "middle")
    #   }
    # }
    
    aaSurround <- 15
    if((originalAaPos - aaSurround) < 1){
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide[z], aaSurround,
                                     originalAaPos, convertedAaPos$Position[z], "start")
    } else if ((originalAaPos + aaSurround) < aaLength){
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide[z], aaSurround,
                                     originalAaPos, convertedAaPos$Position[z], "end")
    } else{
      adjacentAas <- substrMiddlePep(originalPeptide, otherPeptide[z], aaSurround,
                                     originalAaPos, convertedAaPos$Position[z], "middle")
    }
    
    aaCheckOriginal <- substr(originalPeptide, originalAaPos, originalAaPos)
    ### the original peptide doesn't contain the human hotspot location
    if (nchar(aaCheckOriginal) == 0) {
      msaDf <- rbind(msaDf,
                     data.frame("seq1" = NA,"seq2"= NA, "conSeq"= NA,
                                "check" = NA, "conservationScore"= NA,
                                "originalEns" = originalEns, "ENSID" = ensPeptide[z],
                                "Hotspot" = paste0(originalGene, originalAaPos, originalAa),
                                "OriginalPosition" = NA, "ConvertedPosition" = NA, "MutEff" = NA,
                                "liftOverChrom" = coordinatesList$liftOverOut$chrom,
                                "liftOverStart" = coordinatesList$liftOverOut$start,
                                "liftOverEnd" = coordinatesList$liftOverOut$end,
                                "annoMane" = "", stringsAsFactors = FALSE))
      next()
    }
    aaCheckConverted <- substr(otherPeptide[z], convertedAaPos$Position[z], convertedAaPos$Position[z])
    matchingCheck <- ifelse(aaCheckConverted == aaCheckOriginal, "yes", "no")
    
    msaRes <- msaMuscle(AAStringSet(c(adjacentAas$sub_prt1,adjacentAas$sub_prt2)),
                        substitutionMatrix = VTML10)
    conScore <- msaConservationScore(msaRes, substitutionMatrix = VTML10)
    ConSeq <- msaConsensusSequence(msaRes)
    
    msaDf <- rbind(msaDf,
                    data.frame("seq1" = adjacentAas$sub_prt1,"seq2"=adjacentAas$sub_prt2, "conSeq"=ConSeq,
                               "check" = matchingCheck, "conservationScore"= mean(conScore),
                               "originalEns" = originalEns, "ENSID" = ensPeptide[z],
                               "Hotspot" = paste0(originalGene, originalAaPos, originalAa),
                               "OriginalPosition" = originalAaPos,"ConvertedPosition" = convertedAaPos$Position[z],
                               "MutEff" = coordinatesList[[4]],
                               "liftOverChrom" = coordinatesList$liftOverOut$chrom,
                               "liftOverStart" = coordinatesList$liftOverOut$start,
                               "liftOverEnd" = coordinatesList$liftOverOut$end,
                               "annoMane" = originalManeAnno, stringsAsFactors = FALSE))
    
  }
  return(msaDf)
}

### 20230313 this is example of quick check on cds conversion ran for human back and forth 
# 
# reverseConversionCheck <- hotspotDf
# for (i in 1:100) {
#   print(i)
#   testAaConversion <- tryCatch(aaToGenome(gene = hotspotDf$Gene[i],
#                                           position = as.numeric(hotspotDf$Position[i])),
#                                error = function(x) return(NULL))
#   if (is.null(testAaConversion)) {
#     reverseConversionCheck$PositionCheck[i] <- ""
#     next()
#   } else{
#     allChecks <- NULL
#     for (j in 1:nrow(testAaConversion)) {
#       print(j)
#       tmpCheck <- genomeToAAHg(testAaConversion[j, ])
#       allChecks <- c(allChecks, tmpCheck)
#     }
#     reverseConversionCheck$PositionCheck[i] <- paste(allChecks, collapse = "|")
#   }
#   
# }
