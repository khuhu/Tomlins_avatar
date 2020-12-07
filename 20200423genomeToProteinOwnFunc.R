

### protein position to genome position - assumes canonical protein for all of these

### for this function given the gene and the amino acid number

### amino acid function seems to be the easiest one to make - there is no mapping between cdna and positions
### there needs to be a check to see if exon is within protein 
### if not print a warning saying sequence not within cds - doing nucleic seq alignment instead

### after I finish initial version, probably try to make it into a Grange, b/c we can do coordinate intersect



firstUpper <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

cdsConversion <- function(df){
  print("cds conversion start")
  cdsStartConv <- NULL
  cdsEndConv <- NULL
  for (i in 1:nrow(df)) {
    if (i == 1) {
      cdsStartConv <- c(cdsStartConv, 1)
      cdsEndConv <- c(cdsEndConv, df$tmp[i] + 1)
      next()
    }
    else{
      cdsStartConv <- c(cdsStartConv, cdsEndConv[i - 1] + 1)
      cdsEndConv <- c(cdsEndConv, cdsEndConv[i - 1] + df$tmp[i] + 1)
    }
  }
  #print(paste(df, cdsStartConv, cdsEndConv))
  df <- cbind(df, cdsStartConv, cdsEndConv)
  #print("cds conversion end")
  return(df)
}


### multiple sequence alignment
### used for a specific amino acid (i)
###



msaMatching <- function(adjAA, hprt, mprt, df, i){
  msaOut <- msa(c(AAStringSet(hprt), AAStringSet(mprt)), order = "input");
  matchIdxPrt <- matchPattern(adjAA, msaOut@unmasked[[1]])
  
  ### if no initial match then look for match with 5, then 10 mismatches allowed
  if(isEmpty(matchIdxPrt)) {
    print("using 5 max.mismatch")
    matchIdxPrt <- matchPattern(adjAA, msaOut@unmasked[[1]], max.mismatch = 5)
    print('test okay')
    if(isEmpty(matchIdxPrt)){
      print("using 10 max.mismatch")
      matchIdxPrt <- matchPattern(adjAA, msaOut@unmasked[[1]], max.mismatch = 10)
    }
    tmpIdx <- which.min(abs(matchIdxPrt@ranges@start - df$position[i]))
    mouseRecip <- msaOut@unmasked[[2]][matchIdxPrt@ranges[tmpIdx]]
    out <- matchPattern(mouseRecip, mprt)
    if(isEmpty(out)){
      out <- NULL
      return(out)
    }
    else{
      return(list(out, matchIdxPrt)) 
    }
  }
  else{
    mouseRecip <- msaOut@unmasked[[2]][matchIdxPrt@ranges]
    out <- matchPattern(mouseRecip, mprt)
    if(isEmpty(out)){
      out <- NULL
      return(out)
    }
    return(list(out, matchIdxPrt))
  }
}



### matching for just 2 proteins
###
###



msaToMatchingV2 <- function(prt1, prt2){
  #print(c(AAStringSet(prt1), AAStringSet(prt2)))
  
  ### doing this to 
  #prt1 <- str_remove(prt1, "\\*")
  #prt2 <- str_remove(prt2, "\\*")
  
  msaOut <- msa(c(AAStringSet(prt1), AAStringSet(prt2)), order = "input");
  matchIdxPrt <- matchPattern(str_remove(prt1, "\\*"), msaOut@unmasked[[1]])
  
  print(paste("matchIdxPrt",matchIdxPrt))
  #print(msaOut, show = "complete")
  
  mismatchVar <- 0
  ### if no initial match then look for match with 5, then 10 mismatches allowed
  if(isEmpty(matchIdxPrt)) {
    print("using 1/5 max.mismatch")
    matchIdxPrt <- matchPattern(prt1, msaOut@unmasked[[1]], max.mismatch = nchar(prt1)/5)
    mismatchVar <- nchar(prt1)/5
    if(isEmpty(matchIdxPrt)){
      print("using 1/4 max.mismatch")
      matchIdxPrt <- matchPattern(prt1, msaOut@unmasked[[1]], max.mismatch = nchar(prt1)/4)
      mismatchVar <- nchar(prt1)/4
    }
    
    if(isEmpty(matchIdxPrt)){
      out <- NULL
      return(out)
    }
    
    recipPrt <- msaOut@unmasked[[2]][matchIdxPrt@ranges[1]]
    recipPrt <- str_remove_all(recipPrt,"-")
    out <- matchPattern(recipPrt, prt2, max.mismatch = mismatchVar)
    if(isEmpty(out)){
      out <- NULL
      return(out)
    }
    else{
      return(list(out, matchIdxPrt)) 
    }
  }
  else{
    recipPrt <- msaOut@unmasked[[2]][matchIdxPrt@ranges]
    print(paste("recipPrt", recipPrt))
    out <- matchPattern(recipPrt, prt2)
    if(isEmpty(out)){
      out <- NULL
      return(out)
    }
    return(list(out, matchIdxPrt))
  }
}


### converting amino acid matched to genomic coordinates

aaPositionToGenome <- function(msaMatchOut, species2Cds, species2Exon, aaCorStart, aaCorEnd){
  ### quick check then maps amino acid position to cds poistion
  if (is.null(msaMatchOut)) {
    print("no conversion: no msa result")
    exonConversions <- data.frame("chromosome" = "NA","genomic_start" = "NA", "genomic_end" = "NA", stringsAsFactors = FALSE)
    return(exonConversions)
  }
  
  
  ### convert matching positions into cds coordinates
  species2CdsStart <- min(as.numeric(unlist(strsplit(species2Cds$cdna_coding_start, ";"))))
  msaMatchCdsStart <- msaMatchOut[[1]]@ranges@start * 3 - 3 + species2CdsStart - aaCorStart
  msaMatchCdsEnd <- (msaMatchOut[[1]]@ranges@start + msaMatchOut[[1]]@ranges@width) * 3 - 4 + species2CdsStart - aaCorEnd
  
  ### mapping cds position back to genomic
  
  cdsSubject <- IRanges(msaMatchCdsStart, msaMatchCdsEnd)
  cdsRangeQuery <- IRanges(species2Exon$cdsStartConv, species2Exon$cdsEndConv)
  
  overlapRes <- GenomicRanges::findOverlaps(cdsRangeQuery, cdsSubject)
  
  #print(species2Exon)
  print(cdsSubject)
  print(cdsRangeQuery)
  print(overlapRes)
  print(msaMatchCdsStart)
  print(msaMatchCdsEnd)
  print(species2Exon$cdsRemapVar[queryHits(overlapRes)])
  
  ### three posibilities (1) only one match (2) two matches i.e spannign two exons (3) spanning more than 2
  if (length(queryHits(overlapRes)) == 1) {
    print("AA overlap 1")
    #tmpStart <- msaMatchCdsStart - species2Exon$cdsStartConv[queryHits(overlapRes)] + species2Exon$exon_chrom_start[queryHits(overlapRes)]  - 1
    #tmpEnd <- msaMatchCdsEnd - species2Exon$cdsStartConv[queryHits(overlapRes)] + species2Exon$exon_chrom_start[queryHits(overlapRes)]  - 1
    
    tmpStart <- abs(msaMatchCdsStart - species2Exon$cdsRemapVar[queryHits(overlapRes)]) + species2Exon$exon_chrom_start[queryHits(overlapRes)]  - 1
    tmpEnd <- abs(msaMatchCdsEnd - species2Exon$cdsRemapVar[queryHits(overlapRes)]) + species2Exon$exon_chrom_start[queryHits(overlapRes)]  - 1
  }
  else if(length(queryHits(overlapRes) == 2)){
    print("AA overlap 2")
    
    #tmpStartCoord1 <- msaMatchCdsStart - species2Exon$cdsStartConv[queryHits(overlapRes)[1]] + species2Exon$exon_chrom_start[queryHits(overlapRes)[1]]  -  1
    #tmpEndCoord1 <- species2Exon$exon_chrom_end[queryHits(overlapRes)[1]] 
    
    #tmpStartCoord2 <- species2Exon$exon_chrom_start[queryHits(overlapRes)[2]]
    #tmpEndCoord2 <- msaMatchCdsEnd - species2Exon$cdsStartConv[queryHits(overlapRes)[2]] + species2Exon$exon_chrom_start[queryHits(overlapRes)[2]]  - 1
    
    
    #tmpStartCoord1 <- abs(msaMatchCdsStart - species2Exon$cdsRemapVar[queryHits(overlapRes)[1]]) + species2Exon$exon_chrom_start[queryHits(overlapRes)[1]]  -  1
    #tmpEndCoord1 <- species2Exon$exon_chrom_end[queryHits(overlapRes)[1]]
    
    #tmpStartCoord2 <- species2Exon$exon_chrom_start[queryHits(overlapRes)[2]]
    #tmpEndCoord2 <- msaMatchCdsEnd - species2Exon$cdsRemapVar[queryHits(overlapRes)[2]] + species2Exon$exon_chrom_start[queryHits(overlapRes)[2]]  - 1
    
    
    if (species2Exon$strand[1] > 0) {
      tmpStartCoord1 <- abs(msaMatchCdsStart - species2Exon$cdsRemapVar[queryHits(overlapRes)[1]]) + species2Exon$exon_chrom_start[queryHits(overlapRes)[1]]  -  1
      tmpEndCoord1 <- species2Exon$exon_chrom_end[queryHits(overlapRes)[1]]
      
      tmpStartCoord2 <- species2Exon$exon_chrom_start[queryHits(overlapRes)[2]]
      tmpEndCoord2 <- abs(msaMatchCdsEnd - species2Exon$cdsRemapVar[queryHits(overlapRes)[2]]) + species2Exon$exon_chrom_start[queryHits(overlapRes)[2]]  - 1
      
    }
    
    else if(species2Exon$strand[1] < 0){
      tmpStartCoord1 <- species2Exon$exon_chrom_end[queryHits(overlapRes)[1]] - abs(msaMatchCdsStart - species2Exon$cdsRemapVar[queryHits(overlapRes)[1]]) - 1
      tmpEndCoord1 <- species2Exon$exon_chrom_end[queryHits(overlapRes)[1]]
      
      tmpStartCoord2 <- species2Exon$exon_chrom_end[queryHits(overlapRes)[2]]
      tmpEndCoord2 <- abs(msaMatchCdsEnd - species2Exon$cdsRemapVar[queryHits(overlapRes)[2]]) + species2Exon$exon_chrom_start[queryHits(overlapRes)[2]]  - 1
    }

    
    tmpStart <- c(tmpStartCoord1, tmpStartCoord2)
    tmpEnd <- c(tmpEndCoord1, tmpEndCoord2)
  }
  else if(length(queryHits(overlapRes)) > 2){
    print("AA overlap >2")
    tmpStartCoord1 <- msaMatchCdsStart - species2Exon$cdsStartConv[queryHits(overlapRes)[1]] + species2Exon$exon_chrom_start[queryHits(overlapRes)[1]]  -  1
    tmpEndCoord1 <- species2Exon$exon_chrom_end[queryHits(overlapRes)[1]] 
    
    ### this doesn't sufficiently obtain the middle overlaps, also the strand fix doesn't work for this yet, need to think how to properly do this from some examples
    ### down the line
    
    tmpStartCoord2 <- species2Exon$exon_chrom_start[2:(length(queryHits(overlapRes)) - 1)]
    tmpEndCoord2 <- species2Exon$exon_chrom_end[2:(length(queryHits(overlapRes)) - 1)]
    
    tmpStartCoord3 <- species2Exon$exon_chrom_start[queryHits(overlapRes)[length(queryHits(overlapRes))]]
    tmpEndCoord3 <- msaMatchCdsEnd - species2Exon$cdsStartConv[queryHits(overlapRes)[length(queryHits(overlapRes))]] + species2Exon$exon_chrom_start[queryHits(overlapRes)[length(queryHits(overlapRes))]]  - 1

    
    tmpStart <- c(tmpStartCoord1, tmpStartCoord2, tmpStartCoord3)
    tmpEnd <- c(tmpEndCoord1, tmpEndCoord2, tmpEndCoord3)
  }
  
  ### need chromosome - just re -query both exon tables with the chromosome numbers
  print(tmpStart)
  print(tmpEnd)
  exonConversions <- data.frame("chromosome" = as.numeric(species2Exon$chromosome_name[1]),"genomic_start" = as.numeric(tmpStart), "genomic_end" = as.numeric(tmpEnd))
  return(exonConversions)
}









### exon conversion 

exonConversion <- function(geneName, exonList, species1, species2){
  
  
  ### beginning to subset the large exons tables on a per gene basis - later could be query so loading not needed
  if(species1 == "human"){
    
    ### humans
    tmpTxnId <- queryResUCSCfiltHg38$ensembl_transcript_id_version[which(queryResUCSCfiltHg38$external_gene_name == toupper(geneName))]
    tmpExonDat <- geneExonResHs[which(geneExonResHs$ensembl_transcript_id_version == tmpTxnId),]
    tmpCdsDat <- geneCdsResHs[which(geneCdsResHs$ensembl_transcript_id_version == tmpTxnId),]
    tmpPeptideDat <- peptideHs[which(peptideHs$ensembl_transcript_id_version == tmpTxnId),]
    
    
    ### mouse reciprocal
    tmpTxnId2 <- queryResUcscMm10Filt$ensembl_transcript_id_version[which(queryResUcscMm10Filt$external_gene_name == firstUpper(geneName))]
    if (isEmpty(tmpTxnId2)) {
      tmpGeneName <- geneNamesDf$mmusculus_homolog_associated_gene_name[which(geneNamesDf == toupper(geneName))]
      tmpTxnId2 <- queryResUcscMm10Filt$ensembl_transcript_id_version[which(queryResUcscMm10Filt$external_gene_name == tmpGeneName)]
    }
    tmpCdsDat2 <- geneCdsResMm[which(geneCdsResMm$ensembl_transcript_id_version == tmpTxnId2),]
    tmpExonDat2 <- geneExonResMmV2[which(geneExonResMmV2$ensembl_transcript_id_version == tmpTxnId2),]
    prtSpecies2 <- peptideMm$peptide[which(peptideMm$ensembl_transcript_id_version == tmpTxnId2)]
  }
  
  ### just the opposite of the above
  else{
    if(species1 == "mouse"){
      tmpTxnId <- queryResUcscMm10Filt$ensembl_transcript_id_version[which(queryResUcscMm10Filt$external_gene_name == firstUpper(geneName))]
      tmpExonDat <- geneExonResMmV2[which(geneExonResMmV2$external_gene_name == tmpTxnId),]
      tmpCdsDat <- geneCdsResMm[which(geneCdsResMm$external_gene_name == tmpTxnId),]
      tmpPeptideDat <- peptideMm[which(peptideMm$ensembl_transcript_id_version == tmpTxnId),]
      
      ### currently no name change for usage of mouse to human but easily doable later on
      tmpTxnId2 <- queryResUCSCfiltHg38$ensembl_transcript_id_version[which(queryResUCSCfiltHg38$external_gene_name == toupper(geneName))]
      tmpExonDat2 <- geneExonResHs[which(geneExonResHs$ensembl_transcript_id_version == tmpTxnId2),]
      tmpCdsDat2 <- geneCdsResHs[which(geneCdsResHs$ensembl_transcript_id_version == tmpTxnId2),]
      prtSpecies2 <- peptideHs$peptide[which(peptideHs$ensembl_transcript_id_version == tmpTxnId2)]
    }
  }
  
  ### exon processing - down the line I can just create functions to process both exon and cds data properly - would make things a lot cleaner
  tmpExonDat$tmp <- tmpExonDat$exon_chrom_end - tmpExonDat$exon_chrom_start
  tmpExonDat <- tmpExonDat[order(tmpExonDat$rank),]
  tmpExonDat <- cdsConversion(tmpExonDat)
  
  
  tmpExonDat2$tmp <- tmpExonDat2$exon_chrom_end - tmpExonDat2$exon_chrom_start
  tmpExonDat2 <- tmpExonDat2[order(tmpExonDat2$rank),]
  tmpExonDat2 <- cdsConversion(tmpExonDat2)
  
  
  ### the added variable is used to deal with strandedness
  if (tmpExonDat2$strand[1] > 0) {
    tmpExonDat2$cdsRemapVar <- tmpExonDat2$cdsStartConv
  }
  else if(tmpExonDat2$strand[1] < 0){
    tmpExonDat2$cdsRemapVar <- tmpExonDat2$cdsEndConv
    #tmpStart <- tmpExonDat2$exon_chrom_start
  }
  
  
  
  ### cds processing
  cdsStart <- as.numeric(unlist(strsplit(tmpCdsDat$cdna_coding_start, ";")))
  cdsEnd <- as.numeric(unlist(strsplit(tmpCdsDat$cdna_coding_end, ";")))
  cdsRangeQuery <- IRanges(cdsStart, cdsEnd)
  print(paste("cdsRangeQuery", cdsRangeQuery))
  
  finalTable <- NULL
  for (i in seq_along(exonList)) {
    tmpDf <- tmpExonDat[which(as.numeric(tmpExonDat$rank) == i),]
    View(tmpDf)
    cdsSubject <- IRanges(tmpDf$cdsStartConv, tmpDf$cdsEndConv)
    # check to see if exon is within cds
    intersectRes <- GenomicRanges::intersect(cdsRangeQuery, cdsSubject)
    
    print(paste("cdsSubject",cdsSubject))
    print(GenomicRanges::intersect(cdsRangeQuery, cdsSubject))
    check <- "good"
    if (length(intersectRes) < 1) {
      ### for now I'll skip, later i can figure out sequence alignmnet with nucleic acid
      check <- "bad"
      print(paste("exon", i,  "not within cds"))
      next()
    }
    
    
    ### while loop is for checking translation/cds is correct. example NRAS exon 4
    #translateCheck <- "bad"
    tmpPrtStart <- intersectRes@start 
    tmpPrtEnd <- intersectRes@start + intersectRes@width - 1
    
    ### usage of counter assumes only offset by 1 bp
    
    print(paste(tmpPrtStart, tmpPrtEnd))
    counter <- 0
    aaCorrectionStart <- 0
    aaCorrectionEnd <- 0
    
    ### checking if cds are correct 
    ### NRAS is example of offset by 1, and PTEN is example of offset by 2
    
    while (1) {
      tmpTranslate <- Biostrings::translate(DNAString(substring(tmpCdsDat$cdna, tmpPrtStart, tmpPrtEnd)))
      translateCheck <- matchPattern(tmpTranslate, tmpPeptideDat$peptide)
      if (length(translateCheck) > 0) {
        break()
      }
      else{
        tmpPrtStart <- tmpPrtStart + 1
        aaCorrectionStart <- aaCorrectionStart + 1
      }
    }

    
    prtSpecies1  <- Biostrings::translate(DNAString(substring(tmpCdsDat$cdna, tmpPrtStart, tmpPrtEnd)))
    
    ### msa function - everything else looks fine
    print(prtSpecies1)
    msaOut <- msaToMatchingV2(prtSpecies1, prtSpecies2)
    
    ### take amino acid position and bring it back to genomic coordinates
    ### feed above into the aaToGenome
    
    genomeRes <- aaPositionToGenome(msaOut, tmpCdsDat2, tmpExonDat2, aaCorrectionStart, aaCorrectionEnd)
    genomeRes <- cbind(genomeRes, "original" = paste0(geneName,"_","exon","_", exonList[i]))
    finalTable <- rbind(finalTable, genomeRes)
  }
  return(finalTable)
}


### function for lookng at individual amino acid positions
###
###


orthoMutsV2 <- function(geneName, aminoAcid , position, species1, species2){
  
  ### subsetting necessary data - later on this can be a SQL query if needed
  if (species1 == "human") {
    tmpTxnId <- queryResUCSCfiltHg38$ensembl_transcript_id_version[which(queryResUCSCfiltHg38$external_gene_name == toupper(geneName))]
    prt1 <- peptideHs$peptide[which(peptideHs$ensembl_transcript_id_version == tmpTxnId)]
    
    
    tmpTxnId2 <- queryResUcscMm10Filt$ensembl_transcript_id_version[which(queryResUcscMm10Filt$external_gene_name == firstUpper(geneName))]
    if (isEmpty(tmpTxnId2)) {
      tmpGeneName <- geneNamesDf$mmusculus_homolog_associated_gene_name[which(geneNamesDf == toupper(geneName))]
      tmpTxnId2 <- queryResUcscMm10Filt$ensembl_transcript_id_version[which(queryResUcscMm10Filt$external_gene_name == tmpGeneName)]
      
    }
    prt2 <- peptideMm$peptide[which(peptideMm$ensembl_transcript_id_version == tmpTxnId2)]
    tmpCdsDat2 <- geneCdsResMm[which(geneCdsResMm$ensembl_transcript_id_version == tmpTxnId2),] 
    tmpExonDat2 <- geneExonResMmV2[which(geneExonResMmV2$ensembl_transcript_id_version == tmpTxnId2),]
    tmpExonDat2$tmp <- tmpExonDat2$exon_chrom_end - tmpExonDat2$exon_chrom_start
    tmpExonDat2 <- tmpExonDat2[order(tmpExonDat2$rank),]
    tmpExonDat2 <- cdsConversion(tmpExonDat2)
  }
  else if(species2 == "mouse"){
    tmpTxnId2 <- queryResUcscMm10Filt$ensembl_transcript_id_version[which(queryResUcscMm10Filt$external_gene_name == firstUpper(geneName))]
    prt1 <-peptideMm$peptide[which(peptideMm$external_gene_name == firstUpper(geneName))]
    prt2 <- peptideHs$peptide[which(peptideHs$external_gene_name == toupper(geneName))]
    
    tmpCdsDat2 <- geneCdsResHs[which(geneCdsResHs$external_gene_name == toupper(geneName)),]
    tmpExonDat2 <- geneExonResHs[which(geneExonResHs$external_gene_name == toupper(geneName)),]
    tmpExonDat2$tmp <- tmpExonDat2$exon_chrom_end - tmpExonDat2$exon_chrom_start
    tmpExonDat2 <- tmpExonDat2[order(tmpExonDat2$rank),]
    tmpExonDat2 <- cdsConversion(tmpExonDat2)
  }
  
  ### check
  aaCheck <- substr(prt1, position, position)
  if (aaCheck != aminoAcid) {
    print("warning: amino acid position does not match with ensembl canonical transcript/protein, so skipping")
    return(data.frame("chromosome" = "NA","genomic_start" = "NA", "genomic_end" = "NA", "original" = paste0(geneName, position, aminoAcid),
                      stringsAsFactors = FALSE))
  }
  
  
  ### the added variable is used to deal with strandedness
  if (tmpExonDat2$strand[1] > 0) {
    tmpExonDat2$cdsRemapVar <- tmpExonDat2$cdsStartConv
  }
  else if(tmpExonDat2$strand[1] < 0){
    tmpExonDat2$cdsRemapVar <- tmpExonDat2$cdsEndConv
    #tmpStart <- tmpExonDat2$exon_chrom_start
  }
  
  
  ### the next 3 if else statements are for positioning purposes since  I match surrounding areas
  
  ### starts
  if (floor(position - .005 * nchar(prt1) < 1)) {
    adjAA <- substr(prt1, start = position,
                    stop = ceiling(position + .01 * nchar(prt)))
    msaOut <- msaToMatchingV2(adjAA, prt2)
    
    ### since the adjAA is XNNNNNN make sure to use first coordinate of genome results
    genomeRes <- aaPositionToGenome(msaOut, tmpCdsDat2, tmpExonDat2, aaCorStart = 0 , aaCorEnd = 0)
    aaConversion <- data.frame("chromosome" = genomeRes$chromosome, "start" = genomeRes$genomic_start,
                               "end" = genomeRes$genomic_start + 2, "original" = paste0(geneName, position, aminoAcid),stringsAsFactors = FALSE)
    print("start")
  }
  
  ### end
  else if(ceiling(position + .005 * nchar(prt1)) > nchar(prt1)){
    adjAA <- substr(prt1, start = floor(position - .01 * nchar(prt1)),
                    stop = position)
    msaOut <- msaToMatchingV2(adjAA, prt2)
    
    ### since the adjAA is NNNNNNX make sure to use the last coordinate of genome results
    genomeRes <- aaPositionToGenome(msaOut, tmpCdsDat2, tmpExonDat2, aaCorStart = 0 , aaCorEnd = 0)
    aaConversion <- data.frame("chromosome" = genomeRes$chromosome, "start" = genomeRes$genomic_start - 2,
                               "end" = genomeRes$genomic_start, "original" = paste0(geneName, position, aminoAcid),stringsAsFactors = FALSE)
    print("end")
  }
  
  ### everything in the middle
  else {
    adjAA <- substr(prt1, start = floor(position - .005 * nchar(prt1)),
                    stop = ceiling(position + .005 * nchar(prt1)))
    msaOut <- msaToMatchingV2(adjAA, prt2)
    
    ### since the adjAA is NNNXNNN make sure to use the middle (mean) coordinate of genome results
    genomeRes <- aaPositionToGenome(msaOut, tmpCdsDat2, tmpExonDat2, aaCorStart = 0 , aaCorEnd = 0)
    
    print(genomeRes)
    
    ### this middle conversion should just be the mean - 1 and mean + , but only if the length of the genomic coordinate has length 1, i.e one contiguous region
    ### assumes at most it spans 2 exons
    if (nrow(genomeRes) == 2) {
      print(2)
      cdsLength1 <- genomeRes$genomic_end[1] - genomeRes$genomic_start[1]
      cdsLength2 <- genomeRes$genomic_end[2] - genomeRes$genomic_start[2]
      meanCds <- floor(mean(c(cdsLength1, cdsLength2)))
      if (cdsLength1 > cdsLength2) {
        aaConversion <- data.frame("chromosome" = genomeRes$chromosome[1], "start" = genomeRes$genomic_start[1] + meanCds,
                                   "end" = genomeRes$genomic_start[1] + meanCds + 3, "original" = paste0(geneName, position, aminoAcid),stringsAsFactors = FALSE)
      }
      else if(cdsLength2 > cdsLength1){
        aaConversion <- data.frame("chromosome" = genomeRes$chromosome[2], "start" = genomeRes$genomic_end[2] - meanCds,
                                   "end" = genomeRes$genomic_end[2] - meanCds + 3, "original" = paste0(geneName, position, aminoAcid),stringsAsFactors = FALSE)
      }
    }
    else if(nrow(genomeRes) == 1){
      print(1)
      aaConversion <- data.frame("chromosome" = genomeRes$chromosome, "start" = mean(c(genomeRes$genomic_start, genomeRes$genomic_end)) - 1,
                                 "end" = mean(c(genomeRes$genomic_start, genomeRes$genomic_end)) + 2, "original" = paste0(geneName, position, aminoAcid),stringsAsFactors = FALSE)
    }
    
    print("else")
  }
  ###apparently faster to return since it's a non-local jump
  #return(aaConversion)
  aaConversion
}


















### make below into function .. can be used to map back coordinates too


geneExonRes <- geneExonResHs[which(geneExonResHs$external_gene_name == "PTEN"),]
cdsStartConv <- NULL
cdsEndConv <- NULL
tmpExonDat <- geneExonRes
tmpExonDat$tmp <- tmpExonDat$exon_chrom_end - tmpExonDat$exon_chrom_start
tmpExonDat <- tmpExonDat[order(geneExonRes$rank),]
for (i in 1:nrow(tmpExonDat)) {
  if (i == 1) {
    cdsStartConv <- c(cdsStartConv, 1)
    cdsEndConv <- c(cdsEndConv, tmpExonDat$tmp[i] + 1)
    
    next()
  }
  else{
    cdsStartConv <- c(cdsStartConv, cdsEndConv[i - 1] + 1)
    cdsEndConv <- c(cdsEndConv, cdsEndConv[i - 1] + tmpExonDat$tmp[i] + 1)
  }
}
tmpExonDat <- cbind(tmpExonDat, cdsStartConv, cdsEndConv)
tmpExonDat$cdsRevStart <- 1
tmpExonDat$cdsReVEnd <- tmpExonDat$tmp + 1

hsNrasPep <- peptideHs$peptide[which(peptideHs$external_gene_name == "NRAS")]
test <- msaToMatchingV2("TIQLIQNHFVD", hsNrasPep)
test[[1]]@ranges@start


tmpQuery <- IRanges(100, 130)
tmpSub <- IRanges(tmpExonDat$cdsStartConv, tmpExonDat$cdsEndConv)
findOverlaps(tmpQuery, tmpSub)
115 + 114716050 - 1 - 114 


testTryCatch <- tryCatch(Biostrings::translate(DNAString(substring(geneCdsResHs$cdna[9394], 846, 924))), warning = function(x) return("bad"))
testTryCatch

Biostrings::translate(DNAString(substring(geneCdsResHs$cdna[9394], 925, 1009)))

testNras <- Biostrings::translate(DNAString(substring(geneCdsResHs$cdna[1669], 132, 242)))
matchPattern(testNras, peptideHs$peptide[1674])
testNras <- Biostrings::translate(DNAString(substring(geneCdsResHs$cdna[1669], 243, 421)))
matchPattern(testNras, peptideHs$peptide[1674])
testNras <- Biostrings::translate(DNAString(substring(geneCdsResHs$cdna[1669], 422 + 1, 581)))
matchPattern(testNras, peptideHs$peptide[1674])


testPten <- Biostrings::translate(DNAString(substring(geneCdsResHs$cdna[9394], 846, 924)))
matchPattern(testPten, peptideHs$peptide[9495])
testPten <- Biostrings::translate(DNAString(substring(geneCdsResHs$cdna[9394], 925 + 1, 1009)))
testPtenOut <- matchPattern(testPten, peptideHs$peptide[9495])
testPten <- Biostrings::translate(DNAString(substring(geneCdsResHs$cdna[9394], 925 + 2, 1009)))
matchPattern(testPten, peptideHs$peptide[9495])
testPten <- Biostrings::translate(DNAString(substring(geneCdsResHs$cdna[9394], 1010 , 1054)))
matchPattern(testPten, peptideHs$peptide[9495])
testPten <- Biostrings::translate(DNAString(substring(geneCdsResHs$cdna[9394], 1010 +1, 1054)))
matchPattern(testPten, peptideHs$peptide[9495])



testAA1 <- "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIE"
testAA2 <- peptideMm$peptide[2841]
testMsaOut <- msa(c(AAStringSet(testAA1), AAStringSet(testAA2)), order = "input")
testMsaOut
testMsaOut@unmasked[[1]]
matchPattern(testAA1, testMsaOut@unmasked[[1]], max.mismatch = nchar(testAA1)/3)


testAA1 <- "LGPISLNWFEELSSEAPPYNSEPAEESEHKNNNYEPNLFKTPQRKPSYNQLASTPIIFKEQGLTLPLYQSPVKELDKFKLDL"
testAA2 <- peptideMm$peptide[5453]
testMsaOut <- msa(c(AAStringSet(testAA1), AAStringSet(testAA2)), order = "input")
testMsaOut


testMsaOut@unmasked[[1]]
testMatchIdx <- matchPattern(str_remove(testAA1, "\\*"), testMsaOut@unmasked[[1]])
#testMatchIdx <- matchPattern(testAA1, testMsaOut@unmasked[[1]], max.mismatch = nchar(testAA1)/5)
tmpRecipPrt <- testMsaOut@unmasked[[2]][testMatchIdx@ranges[1]]
tmpRecipPrt <- str_remove_all(tmpRecipPrt, "-")
#matchPattern(tmpRecipPrt, testAA2, max.mismatch = nchar(testAA1)/5)
matchPattern(tmpRecipPrt, testAA2)


testMsaOut@unmasked@ranges[1]@start
testMsaOut[[1]]@ranges@start + testMsaOut[[1]]@ranges@width



### testing exon conversion. the ends are positionally (genomically) off by a couple base pairs 

tmpExonList <- c(1:7)
exonOut <- exonConversion("NRAS", exonList = tmpExonList, species1 = "human", species2 = "mouse")
tmpExonList <- c(1:9)
exonOut <- exonConversion("PTEN", exonList = tmpExonList, species1 = "human", species2 = "mouse")
tmpExonList <- c(1:18)
exonOut <- exonConversion("BRAF", exonList = tmpExonList, species1 = "human", species2 = "mouse")
tmpExonList <- c(1:10)
exonOut <- exonConversion("BRCA2", exonList = tmpExonList, species1 = "human", species2 = "mouse")




### they all look good. can run script to run majority of genes exons later to give a decent statistic - only weird thing is the second matched one for BRAF exon 18




### testing the amino acid translator/converter
### BRAF V600E, G469A, L597
### BRCA1 A1708E, A1789T, M1775R  - presents problem of when the adjAA spands more than one exon
### TP53 V173M, C242S, L330R - this test mainly to see if genes with slightly different names still work



### if strand is negative, then i want to subtract using cdsEndConv, whereas if strand is positive, 
### I want to subtract the cdsStartConv, this should all be within abs for absolute diffrence, then add genomic start


### figured out this "works" only to a degree where the position is off by like -4 bps (should be upstream)
testOrth <- orthoMutsV2("BRAF", "V", 600, species1 = "human", species2 = "mouse")
testOrth <- orthoMutsV2("BRAF", "G", 469, species1 = "human", species2 = "mouse")
testOrth <- orthoMutsV2("BRAF", "L", 597, species1 = "human", species2 = "mouse")


testOrth <- orthoMutsV2("BRCA1", "A", 1708, species1 = "human", species2 = "mouse")
testOrth <- orthoMutsV2("BRCA1", "A", 1789, species1 = "human", species2 = "mouse")
testOrth <- orthoMutsV2("BRCA1", "M", 1775, species1 = "human", species2 = "mouse")

testOrth <- orthoMutsV2("TP53", "V", 1703, species1 = "human", species2 = "mouse")
testOrth <- orthoMutsV2("TP53", "C", 242, species1 = "human", species2 = "mouse")
testOrth <- orthoMutsV2("TP53", "L", 330, species1 = "human", species2 = "mouse")

testAA1 <- "GLATVKSRW"
testAA2 <- peptideHs$peptide[6956]
testMsaOut <- msa(c(AAStringSet(testAA1), AAStringSet(testAA2)), order = "input")
testMatchIdxPrt <- matchPattern(str_remove(testAA1, "\\*"), testMsaOut@unmasked[[1]])
testRecipPrt <- testMsaOut@unmasked[[2]][testMatchIdxPrt]
testAAinput <- matchPattern(testRecipPrt, peptideMm$peptide[6070])
testAAinput

testSpecies2CdsStart <- min(as.numeric(unlist(strsplit(geneCdsResMm$cdna_coding_start[6080], ";"))))
testMsaMatchCdsStart <- testAAinput@ranges@start * 3 - 3 + testSpecies2CdsStart
testMsaMatchCdsEnd <- (testAAinput@ranges@start + testAAinput@ranges@width) * 3 - 4 + testSpecies2CdsStart

tmpIdtest <- "ENSMUST00000002487.14"
tmpExonDattest  <- geneExonResMmV2[which(geneExonResMmV2$ensembl_transcript_id_version ==  tmpIdtest),]
tmpExonDattest$tmp <- tmpExonDattest$exon_chrom_end - tmpExonDattest$exon_chrom_start
tmpExonDattest <- tmpExonDattest[order(tmpExonDattest$rank),]
tmpExonDattest <- cdsConversion(tmpExonDattest)


tmpId <- "ENST00000646891.1"
tmpExonDat  <- geneExonResHs[which(geneExonResHs$ensembl_transcript_id_version ==  tmpId),]
tmpExonDat$tmp <- tmpExonDat$exon_chrom_end - tmpExonDat$exon_chrom_start
tmpExonDat <- tmpExonDat[order(tmpExonDat$rank),]
tmpExonDat <- cdsConversion(tmpExonDat)




msa::print(testMsaOut, show="complete")

tmpMatchOut <- matchPattern("GLATVKSRW", peptideMm$peptide[6070])





  