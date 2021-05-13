library(ggplot2)

cyto_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_cytoDf.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = FALSE)
mouse_cyto <- cyto_hg38mm10[grep("m_chr",cyto_hg38mm10$V1),]
mouse_cyto$V6 <- paste0(str_remove_all(mouse_cyto$V1, "m_chr"), mouse_cyto$V4)
mouse_cyto$V7 <- mouse_cyto$V3 - mouse_cyto$V2

brca_brca1del <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/benDavidAnalysis/BenDavidGemmSupp10.xls",
                                  sheet = 9) 

brca_trp53del <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/benDavidAnalysis/BenDavidGemmSupp10.xls",
                                   sheet = 3) 

brca_sv40tag <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/benDavidAnalysis/BenDavidGemmSupp10.xls",
                                   sheet = 5) 

brca_mycdel <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/benDavidAnalysis/BenDavidGemmSupp10.xls",
                                   sheet = 2) 

# so the definition of LST with visualizaiton:
# chromosomal break between adjacent regions of at least 10 Mb
# some groups need another min a distance between them not larger than 3Mb.
# https://github.com/sztup/scarHRD - not using that
# checks segment is at least > 10Mb

lstCheck <- function(cnDf){
  cnLength <- sum(mouse_cyto$V7[which(mouse_cyto$V6 %in% cnDf$Band)])
  # length check
  if (cnLength > 10000000) {
    return(1)
  } else{
    return(0)
  }
}



lstCount <- function(x){
  lstCountTable <- NULL
  for (i in 2:ncol(x)) {
    sampleDf <- x[,c(1,i)]
    # first reduce the matrix to anything that as a CNA
    cna_idx <- which(sampleDf[,2] != 0)
    # finds consecutive regions ... function is like magic
    adjObj <- split(cna_idx, cumsum(c(1, diff(cna_idx) != 1)))
    lstCount <- 0
    for (j in 1:length(adjObj)) {
      tmpDf <- sampleDf[adjObj[[j]],]
      #check that consecutive cnas have same copy-number
      valChecker <- c(1,1+which(diff(tmpDf[,2])!=0))
      valChecker <- c(valChecker, nrow(tmpDf))
      if (length(valChecker) > 2) {
        # instance of more than one in region
        for (j in 1:(length(valChecker)-1)) {
          tmpDf2 <- tmpDf[k:k+1, 2]
          lstCount <- lstCount + lstCheck(tmpDf2)
        }
      } else{
        # instance of one in region
        lstCount <- lstCount + lstCheck(tmpDf)
      }
    }
    lstCountTable <- rbind(lstCountTable,c(colnames(sampleDf[,2]), lstCount))
  }
  lstCountTable <- data.frame(lstCountTable, stringsAsFactors = FALSE)
  colnames(lstCountTable) <- c("Samples", "Count")
  lstCountTable$Count <- as.numeric(lstCountTable$Count)
  return(lstCountTable)
}

# testing function idea 

lst_brca1 <- lstCount(brca_brca1del)
lst_brca1$type <- "brca1"
lst_trp53 <- lstCount(brca_trp53del)
lst_trp53$type <- "trp53"
lst_sv40tag <- lstCount(brca_sv40tag)
lst_sv40tag$type <- "sv40tag"
lst_myc <- lstCount(brca_mycdel)
lst_myc$type <- "myc"

lst_obj <- rbind(lst_brca1[,2:3], lst_trp53[,2:3],
                 lst_sv40tag[,2:3], lst_myc[,2:3])

ggplot(lst_obj) + geom_boxplot(aes(x = type,y = Count, fill = type)) + xlab("Genotype") +
  ylab("Number of LSTs")
# be low could possibly explain
# since the majority of cytobands larger than 10Mb
# LSTs could be undercounted

hist(mouse_cyto$V7/100000)

