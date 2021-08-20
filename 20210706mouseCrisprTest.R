listOfGenes <- c("Brca1", "Trp53", "Rb1", "Nf1")

tmpVcf <- read.table("/mnt/DATA5/tmp/kev/CRISPRtest/kc_ovc7.mm10_multianno.txt", sep = "\t", 
                     stringsAsFactors = FALSE, fill = TRUE, skip = 1)
tmpVcf2 <- tmpVcf[which(tmpVcf$V18 == "PASS" & tmpVcf$V6 == "exonic"),]
af <- as.numeric(unlist(lapply(strsplit(tmpVcf2$V21, ":"), '[[', 3)))
ad <- unlist(lapply(strsplit(tmpVcf2$V21, ":"), '[[', 2))
refCount <- as.numeric(unlist(lapply(strsplit(ad, ","), '[[', 1)))
altCount <- as.numeric(unlist(lapply(strsplit(ad, ","), '[[', 2)))
totalFiltCount <- refCount + altCount


tmpVcf3 <- tmpVcf2[which(af > .10 & totalFiltCount > 20),]
tmpVcf3 <- tmpVcf3[grep(paste0(listOfGenes, collapse = "|"), tmpVcf3$V7),]

fileDir <- "/mnt/DATA5/tmp/kev/CRISPRtest/"

setwd(fileDir)
listOfFiles <- system("ls *mm10_multianno.txt",  intern = TRUE)

finalTable <- NULL
for (i in listOfFiles) {
  fname <- paste0(fileDir, i)
  tmpVcf <- read.table(fname, sep = "\t", 
                       stringsAsFactors = FALSE, fill = TRUE, skip = 1)
  sampleName <- str_remove(i, "\\.mm10.*")
  tmpVcf2 <- tmpVcf[which(tmpVcf$V18 == "PASS" & tmpVcf$V6 == "exonic"),]
  #tmpVcf2 <- tmpVcf[which(tmpVcf$V6 == "exonic"),]
  af <- as.numeric(unlist(lapply(strsplit(tmpVcf2$V21, ":"), '[[', 3)))
  ad <- unlist(lapply(strsplit(tmpVcf2$V21, ":"), '[[', 2))
  refCount <- as.numeric(unlist(lapply(strsplit(ad, ","), '[[', 1)))
  altCount <- as.numeric(unlist(lapply(strsplit(ad, ","), '[[', 2)))
  totalFiltCount <- refCount + altCount
  tmpVcf2$af <- af
  tmpVcf2$ref <- refCount
  tmpVcf2$alt <- altCount
  
  #tmpVcf3 <- tmpVcf2[which(af > .05 & totalFiltCount > 20),]
  tmpVcf3 <- tmpVcf2[which(altCount >= 10 & totalFiltCount > 20),]
  tmpVcf3 <- tmpVcf3[grep(paste0(listOfGenes, collapse = "|"), tmpVcf3$V7),]
  res <- cbind("Sample" = rep(sampleName,  nrow(tmpVcf3)),
               tmpVcf3)
  res <- res[which(res$V11 == ""),]
  
  finalTable <- rbind(finalTable, res[,c(1:11,23:25)])
}

### finds p53 in 1784-LT where the ICE analysis doens't detect. otherwise all consistent
