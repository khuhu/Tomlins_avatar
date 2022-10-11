hgsc_gistic <- read.table("/mnt/DATA5/tmp/kev/programs/GISTIC2/mouse_hgsc/all_lesions.conf_99.txt", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)




adenoma_gistic <- read.table("/mnt/DATA5/tmp/kev/programs/GISTIC2/mouse_adenoma/all_lesions.conf_99.txt", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)

numberMarkers <- NULL
meanSeg <- NULL
varSeg <- NULL
zeroMeanCount <- NULL
posCount <- NULL
negCount <- NULL

for (i in 1:nrow(hgsc_gistic)) {
  tmpString <- str_remove(str_remove(hgsc_gistic$Peak.Limits[i], ".*probes "), "\\).*")
  tmpString2 <- as.numeric(unlist(str_split(tmpString, ":")))
  tmpMarks <- round((tmpString2[2] - tmpString2[1])/2) + 1
  numberMarkers <- c(numberMarkers, tmpMarks)
  
  print(tmpString2)
  print(tmpMarks)
  
  tmpSegs  <- as.numeric(hgsc_gistic[i,10:112])
  tmpSegs[which(abs(tmpSegs) < 0.2)] <- 0
  
  tmpMean <- mean(tmpSegs)
  meanSeg <- c(meanSeg, tmpMean)
  
  tmpVar <- var(tmpSegs)
  varSeg <- c(varSeg, tmpVar)
  
  tmpZC <- length(which(tmpSegs == 0))
  zeroMeanCount <- c(zeroMeanCount, tmpZC)
  
  tmpNeg <- length(which(tmpSegs < 0))
  negCount <- c(negCount, tmpNeg)
  
  tmpPos <- length(which(tmpSegs > 0))
  posCount <- c(posCount, tmpPos)
}

hgsc_gistic2 <- hgsc_gistic[,c(1:4,6)]

hgsc_gistic2$numMarks <- numberMarkers
hgsc_gistic2$segMean <- meanSeg
hgsc_gistic2$segVar <- varSeg
hgsc_gistic2$zeroCount <- zeroMeanCount
hgsc_gistic2$posCount <- posCount
hgsc_gistic2$negCount <- negCount

hgsc_gistic2 <- hgsc_gistic2[grep("CN", hgsc_gistic2$Unique.Name),]
hgsc_gistic2 <- hgsc_gistic2[order(as.numeric(hgsc_gistic2$q.values), decreasing = FALSE),]

rownames(hgsc_gistic2) <- NULL
### self anno - for text labels on plot

peakNames <- c("Trp53Del", "Brca1Del", "NND", "Brca2Amp", "MycAmp", "Rictor/RAD1/DroshaAmp",
  "Ccnd1/HrasAmp", "Cdkn1bAmp", "Ccne2Amp", "EgfrAmp", "GnasAmp", "TertDel",
  "NND","Sox9Amp","Kmt5cAmp", "NoSig", "NoSig")
hgsc_gistic2$peakNames <- peakNames





crc_gistic <- read.table("/mnt/DATA5/tmp/kev/programs/GISTIC2/mouse_crc/all_lesions.conf_99.txt", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)




tm <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220113cholabAffymetrixListAnno.xlsx")
tmp2 <- tm[-which(tm$Preference %in%  c("Low", "Bad")),]


write.table(tmp2, "/mnt/DATA5/tmp/kev/misc/20220113cholabAffymetrixListAnnoFilt.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


load("/mnt/DATA5/tmp/kev/misc/clusters.gigamuga.Rdata")



library(devtools)

## allow R to look for pacakges in both CRAN and Bioconductor
# setRepositories(ind = 1:2)

## install from Github source
install_github("andrewparkermorgan/argyle", build = TRUE, build_manual = TRUE,
               build_vignettes = TRUE, force = TRUE)
vignette("package" = "argyle")

library(argyle)
library(plyr)

data(ex)
summary(ex)
intensity(ex)$x[1:10]
intensity(ex)$y[1:10]
