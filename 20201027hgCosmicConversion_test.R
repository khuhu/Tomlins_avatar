### creating mouse cosmic convserion

library(stringr)
library(data.table)
library(GenomicRanges)
library(foreach)
library(doParallel)

cmc_table <- fread("/mnt/DATA5/tmp/kev/misc/cmc_export.v92.edited.tsv")
cmc_tsg <- cmc_table[grep(cmc_table$`Mutation Description AA`,
                          pattern = "Substitution - Missense"),]

h38toMm10Chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/hg38ToMm10.over.chain")

small_test <- cmc_tsg[sample(1:nrow(cmc_tsg), 300),]

### making function for the foreach test
convertingSnps <- function(i){
  tmpCoords <- small_test$`Mutation genome position GRCh38`[i]
  badRowOut <- c("hgGene" = NA, "ENST_ID_VER" = NA,
                 "hgCdsMutation" = NA, "hgAaMutation" = NA,
                 "CosmicID" = NA, "EXAC_AF" = NA,
                 "GNOMAD_AF" = NA,
                 "Gr38_position" = NA,
                 "mm10chr" = NA, 
                 "mm10start" = NA,
                 "mm10end" = NA)
  
  if (tmpCoords == ":-") {
    return(badRowOut)
  }
  tmpChr <- str_split(small_test$`Mutation genome position GRCh38`[i], "\\:")[[1]][1]
  tmpRange <- str_split(small_test$`Mutation genome position GRCh38`[i], "\\:")[[1]][2]
  tmpStart <- strsplit(tmpRange, "\\-")[[1]][1]
  tmpEnd <- strsplit(tmpRange, "\\-")[[1]][2]
  
  if (tmpChr == "23") {
    return(badRowOut)
  }
  
  convertedSNV <- liftOver(GRanges(seqnames = paste0("chr", tmpChr),
                                   ranges = IRanges(start = as.numeric(tmpStart),
                                                    end = as.numeric(tmpEnd))),
                           chain = h38toMm10Chain)
  
  if(isEmpty(convertedSNV)){
    return(badRowOut)
  }
  tmpRow <- c("hgGene" = small_test$GENE_NAME[i], "ENST_ID_VER" = small_test$ACCESSION_NUMBER[i],
              "hgCdsMutation" = small_test$`Mutation CDS`[i], "hgAaMutation" = small_test$`Mutation AA`[i],
              "CosmicID" = small_test$GENOMIC_MUTATION_ID[i], "EXAC_AF" = small_test$EXAC_AF[i],
              "GNOMAD_AF" = small_test$GNOMAD_EXOMES_AF[i],
              "Gr38_position" = small_test$`Mutation genome position GRCh38`[i],
              "mm10chr" = as.character(convertedSNV@unlistData@seqnames[1]), 
              "mm10start" = convertedSNV@unlistData@ranges@start,
              "mm10end" = convertedSNV@unlistData@ranges@start + convertedSNV@unlistData@ranges@width - 1)
  return(tmpRow)
}


finalTable <- NULL
registerDoParallel(16)
finalTable <- foreach(i=1:nrow(small_test),.combine=rbind) %dopar% {
  convertingSnps(i)
}
stopImplicitCluster()
finalTable <- data.frame(finalTable, stringsAsFactors = FALSE)
finalTable <- finalTable[-which(is.na(finalTable)),]


### quick test for 100 positions before 4.5m
small_test_converted_table <- NULL
timestamp()
#### above tries to parallelize it
for (i in 1:nrow(small_test)) {
  print(i)
  tmpCoords <- small_test$`Mutation genome position GRCh38`[i]
  if (tmpCoords == ":-") {
    next()
  }
  tmpChr <- str_split(small_test$`Mutation genome position GRCh38`[i], "\\:")[[1]][1]
  tmpRange <- str_split(small_test$`Mutation genome position GRCh38`[i], "\\:")[[1]][2]
  tmpStart <- strsplit(tmpRange, "\\-")[[1]][1]
  tmpEnd <- strsplit(tmpRange, "\\-")[[1]][2]
  
  if (tmpChr == "23") {
    next()
  }
  
  convertedSNV <- liftOver(GRanges(seqnames = paste0("chr", tmpChr),
                   ranges = IRanges(start = as.numeric(tmpStart),
                                    end = as.numeric(tmpEnd))),
           chain = h38toMm10Chain)
  
  if(isEmpty(convertedSNV)){
    print("skip")
    next()
  }
  tmpRow <- c("hgGene" = small_test$GENE_NAME[i], "ENST_ID_VER" = small_test$ACCESSION_NUMBER[i],
              "hgCdsMutation" = small_test$`Mutation CDS`[i], "hgAaMutation" = small_test$`Mutation AA`[i],
              "CosmicID" = small_test$GENOMIC_MUTATION_ID[i], "EXAC_AF" = small_test$EXAC_AF[i],
              "GNOMAD_AF" = small_test$GNOMAD_EXOMES_AF[i],
              "Gr38_position" = small_test$`Mutation genome position GRCh38`[i],
              "mm10chr" = as.character(convertedSNV@unlistData@seqnames[1]), 
              "mm10start" = convertedSNV@unlistData@ranges@start,
              "mm10end" = convertedSNV@unlistData@ranges@start + convertedSNV@unlistData@ranges@width - 1)
  small_test_converted_table <- rbind(small_test_converted_table, tmpRow)
}
timestamp()

small_test_converted_table <- NULL
for (i in 1:nrow(small_test)) {
  print(i)
  tmpRow <- convertingSnps(i)
  if(is.null(tmpRow)){
    next()
  }
  small_test_converted_table <- rbind(small_test_converted_table, tmpRow)
}


