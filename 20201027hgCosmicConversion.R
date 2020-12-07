#!/usr/bin/env Rscript
.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))

library(stringr)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(foreach)
library(doParallel)

#cmc_table <- fread("/mnt/DATA5/tmp/kev/misc/cmc_export.v92.edited.tsv")
h38toMm10Chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/hg38ToMm10.over.chain")


convertingSnps <- function(i){
  tmpCoords <- cmc_table$`Mutation genome position GRCh38`[i]
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
  tmpChr <- str_split(cmc_table$`Mutation genome position GRCh38`[i], "\\:")[[1]][1]
  tmpRange <- str_split(cmc_table$`Mutation genome position GRCh38`[i], "\\:")[[1]][2]
  tmpStart <- strsplit(tmpRange, "\\-")[[1]][1]
  tmpEnd <- strsplit(tmpRange, "\\-")[[1]][2]
  
  
  ### add a search parameter for chromosome 24
  
  if (tmpChr == "23") {
    return(badRowOut)
  }
  
  if (tmpChr == "24") {
    return(badRowOut)
  }
  
  convertedSNV <- liftOver(GRanges(seqnames = paste0("chr", tmpChr),
                                   ranges = IRanges(start = as.numeric(tmpStart),
                                                    end = as.numeric(tmpEnd))),
                           chain = h38toMm10Chain)
  
  if(isEmpty(convertedSNV)){
    return(badRowOut)
  }
  tmpRow <- c("hgGene" = cmc_table$GENE_NAME[i], "ENST_ID_VER" = cmc_table$ACCESSION_NUMBER[i],
              "hgCdsMutation" = cmc_table$`Mutation CDS`[i], "hgAaMutation" = cmc_table$`Mutation AA`[i],
              "CosmicID" = cmc_table$GENOMIC_MUTATION_ID[i], "EXAC_AF" = cmc_table$EXAC_AF[i],
              "GNOMAD_AF" = cmc_table$GNOMAD_EXOMES_AF[i],
              "Gr38_position" = cmc_table$`Mutation genome position GRCh38`[i],
              "mm10chr" = as.character(convertedSNV@unlistData@seqnames[1]), 
              "mm10start" = convertedSNV@unlistData@ranges@start,
              "mm10end" = convertedSNV@unlistData@ranges@start + convertedSNV@unlistData@ranges@width - 1)
  return(tmpRow)
}


setwd("/mnt/DATA5/tmp/kev/misc/")
fileNames <- system("ls *cosmic*", intern = TRUE)

for (j in seq_along(fileNames)) {
  print(j)
  cmc_table <- fread(fileNames[j])
  finalTable <- NULL
  registerDoParallel(24)
  finalTable <- foreach(i=1:nrow(cmc_table),.combine=rbind) %dopar% {
    convertingSnps(i)
  }
  stopImplicitCluster()
  
  finalTable <- data.frame(finalTable, stringsAsFactors = FALSE)
  finalTable <- finalTable[-which(is.na(finalTable)),]
  if (j == 1) {
    write.table(finalTable ,"20201103convertedCosmic_par.txt", sep = "\t", quote = FALSE,
                col.names = TRUE, row.names = FALSE)
  } else{
    write.table(finalTable ,"20201103convertedCosmic_par.txt", sep = "\t", quote = FALSE,
                col.names = TRUE, row.names = FALSE, append = TRUE)
  }
  print(j)
}

#finalTable <- NULL
#registerDoParallel(24)
#finalTable <- foreach(i=1:nrow(cmc_table),.combine=rbind) %dopar% {
#  convertingSnps(i)
#}
#stopImplicitCluster()

#finalTable <- data.frame(finalTable, stringsAsFactors = FALSE)
#finalTable <- finalTable[-which(is.na(finalTable)),]

#write.table(finalTable ,"/home/kevhu/data/20201031convertedCosmic_par.txt", sep = "\t", quote = FALSE,
#            col.names = TRUE, row.names = FALSE)