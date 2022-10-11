### quick script to get get percent identity for proteins of my genes of interest, only need to do this TSGs - expected to have higher levels of 
### conservation, and then find number of exons ... 

library(ggplot2)
data("BLOSUM100")
listOfAllGenes <- readxl::read_xlsx("/home/kevhu/data/20200803mouseInput.xlsx", sheet = 1)
TSGs <- listOfAllGenes$Gene[grep("AllExons", listOfAllGenes$Type)]
oncogenes <- listOfAllGenes$Gene[grep("Copynumber", listOfAllGenes$Type)]

firstUpper <- function(gene){
  firstLetter <- toupper(substr(gene, start = 1, stop = 1))
  restOfGene <- tolower(substr(gene, start = 2, stop = nchar(gene)))
  res <- paste0(firstLetter, restOfGene)
  return(res)
}


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

testDf <- NULL
for (i in TSGs) {
  tmpList <- hg38Peptide$peptide[which(hg38Peptide$external_gene_name == i)]
  
  tmpHgExon <- max(hg38biomartTable$rank[which(hg38biomartTable$external_gene_name == i)])
  
  mmGene <- firstUpper(i)
  tmpList2 <- mm10Peptide$peptide[which(mm10Peptide$external_gene_name == mmGene)]
  if(length(tmpList2) == 0){
    mmGene <- geneNameDf$mmusculus_homolog_associated_gene_name[which(geneNameDf$external_gene_name == i)]
    tmpList2 <- mm10Peptide$peptide[which(mm10Peptide$external_gene_name == mmGene)]
  }
  
  tmpMmExon <- max(mm10biomartTable$rank[which(mm10biomartTable$external_gene_name == mmGene)])
  
  tmpPid_list <- NULL
  for (j in 1:length(tmpList)) {
    if (tmpList[j] %in% c("Sequence unavailable", "NA")) {
      next()
    }
    for (k in 1:length(tmpList2)) {
      if (tmpList2[k] %in% c("Sequence unavailable", "NA")) {
        next()
      }
      tmpPid <- pid(pairwiseAlignment(tmpList[j], tmpList2[k], substitutionMatrix = BLOSUM100), type = "PID2")
      tmpPid_list <- c(tmpPid_list, tmpPid)
    }
  }
  testDf <- rbind(testDf, c("gene" = i, "median" = median(tmpPid_list), "mean" = mean(tmpPid_list),
                  "max" = max(tmpPid_list), "hg_Exons" = tmpHgExon, "mm_Exons" = tmpMmExon))
}

testDf <- data.frame(testDf, stringsAsFactors = FALSE)
testDf[,2:6] <- lapply(testDf[,2:6], as.numeric)
mean(testDf$max)
min(testDf$max)
max(testDf$max)


testDf_onco <- NULL
for (i in oncogenes) {
  tmpList <- hg38Peptide$peptide[which(hg38Peptide$external_gene_name == i)]
  
  tmpHgExon <- max(hg38biomartTable$rank[which(hg38biomartTable$external_gene_name == i)])
  
  mmGene <- firstUpper(i)
  tmpList2 <- mm10Peptide$peptide[which(mm10Peptide$external_gene_name == mmGene)]
  if(length(tmpList2) == 0){
    mmGene <- geneNameDf$mmusculus_homolog_associated_gene_name[which(geneNameDf$external_gene_name == i)]
    tmpList2 <- mm10Peptide$peptide[which(mm10Peptide$external_gene_name == mmGene)]
  }
  
  tmpMmExon <- max(mm10biomartTable$rank[which(mm10biomartTable$external_gene_name == mmGene)])
  
  tmpPid_list <- NULL
  for (j in 1:length(tmpList)) {
    if (tmpList[j] %in% c("Sequence unavailable", "NA")) {
      next()
    }
    for (k in 1:length(tmpList2)) {
      if (tmpList2[k] %in% c("Sequence unavailable", "NA")) {
        next()
      }
      tmpPid <- pid(pairwiseAlignment(tmpList[j], tmpList2[k], substitutionMatrix = BLOSUM100), type = "PID2")
      tmpPid_list <- c(tmpPid_list, tmpPid)
    }
  }
  testDf_onco <- rbind(testDf_onco, c("gene" = i, "median" = median(tmpPid_list), "mean" = mean(tmpPid_list),
                            "max" = max(tmpPid_list), "hg_Exons" = tmpHgExon, "mm_Exons" = tmpMmExon))
}

testDf_onco <- data.frame(testDf_onco, stringsAsFactors = FALSE)
testDf_onco[,2:6] <- lapply(testDf_onco[,2:6], as.numeric)
mean(testDf_onco$max)
min(testDf_onco$max)
max(testDf_onco$max)






### do one for TSGs, and oncogenes
dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210311_TSGs_PID.pdf", useDingbats = FALSE)
ggplot(data = testDf, aes(x = hg_Exons, y = mm_Exons, label = gene)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(aes(colour = max), size = 1) + scale_color_gradient(low = "#00008b",
                                                       high = "#ff0000",
                                                       name = "PID") +
  ggrepel::geom_text_repel(size = 2) + ggtitle("TSGs (fully-tiling genes)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

### onco
dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210311_Onco_PID.pdf", useDingbats = FALSE)
ggplot(data = testDf_onco, aes(x = hg_Exons, y = mm_Exons, label = gene)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(aes(colour = max), size = 1.0) + scale_color_gradient(low = "#00008b",
                                                       high = "#ff0000",
                                                       name = "PID") + 
  ggtitle("Oncogenes (12 amplicons)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()


### from this we have the overrall conservation stats, use range of max and median of max. 
### after pick one of a a gene with different amounts of exons, then can see if any of them are split
### like the one I want to show off


### max in hg is only conserved in the first 3 exons


