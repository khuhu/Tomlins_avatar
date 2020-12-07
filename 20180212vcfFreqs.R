#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")



### things to note is the hotspot file used for calling should have the exact same number of lines as a merged of nonHotspot vcf, but it has 13 less;
### rerun the the tvcutils to create a hotspot file with new merged unison and see if it matches, othertwise tvc util is getting rid of some

### List of things to do  ... first test only two sample merged and see if everything adds up
### For the chromosomes do a grep -v # just for the two merged samp and use code below for chromLabel
### Afterwards, do it for entire data set


library(VariantAnnotation)
vcf <- readVcf("/mnt/DATA4/kevhu/choLab/vcfs/hotspot.merged.vcf", "mm10")

str(vcf)

###fixed is what I'll iterate through
vcf@fixed
rowRanges(vcf)
genomePos <- rownames(vcf)

vcf@rowRanges@elementMetadata
###can just grab the ratios like this
vcf@info["FAO"][[1]]/vcf@info["FDP"][[1]]
vcf@info["OPOS"]
#vcf@info[,39]
vcf@info[,which(colnames(vcf@info) == "OREF")][2]
vcf@info[,which(colnames(vcf@info) == "OALT")][2]

vcf@info["FAO"][[1]][2]
vcf@info["FDP"][[1]][2]

alleFreqTab <- NULL
alleFreqTab <- c("OPOS","FAO","FDP","freq")
for(i in 1:nrow(vcf@info)){
  OPOS <- unlist(vcf@info["OPOS"][[1]][i])
  for(j in seq_along(OPOS)){
    FAO <- unlist(vcf@info["FAO"][[1]][i])
    FAO.1 <- FAO[j]
    #print(FAO.1)
    FDP <- unlist(vcf@info["FDP"][[1]][i])
    #print(FDP.1)
    freq <- FAO.1/FDP.1
    #print(freq)
    rowOfDat <- list(OPOS[j], FAO.1, FDP, freq)
    print(rowOfDat)
    alleFreqTab <- rbind(alleFreqTab, rowOfDat)
  }
}

rownames(alleFreqTab) <- NULL
colnames(alleFreqTab) <- c("OPOS","FAO","FDP","freq")
alleFreqTab <- alleFreqTab[-1,]
alleFreqTab <- data.frame(alleFreqTab, stringsAsFactors = FALSE, row.names = FALSE)


###below file can only be made from combined of all vars
tabOfCoords <- read.table("/mnt/DATA4/kevhu/choLab/vcfs/hotspotIdentifiers.txt",sep = "\t", stringsAsFactors = FALSE)
chromLab <- NULL
for(i in 1:nrow(tabOfCoords)){
  b <- length(unlist(strsplit(tabOfCoords[,5][i],",")))
  chr <- tabOfCoords[,1][i]
  chrom <- rep(chr, b)
  chromLab <- c(chromLab, chrom)
}

