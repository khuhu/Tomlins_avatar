### the enrichment analysis
### eventually I need to combine this 
### with the synteny analysis

library(lsa)
library(ggplot2)
library(philentropy)

### main steps (1) get gene set list and locations on both mouse and human i.e table


### from there I can use the bed files from the synteny plots in order to do some type 
### of synteny enrichment calculation
hg38biomartTable <- read.table("/home/kevhu/data/20201030hg38KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

mm10biomartTable <- read.table("/home/kevhu/data/20201030Mm10KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

geneNameDf <- read.table("/home/kevhu/data/20201021geneNameBiomart.txt", header = TRUE,
                         stringsAsFactors = FALSE, sep = "\t")


human_hallmarks <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/MSigDb/h.all.v7.4.symbols.gmt",
                              sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

human_oncogenic <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/MSigDb/c6.all.v7.4.symbols.gmt",
                              sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

human_immunologic <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/MSigDb/c7.all.v7.4.symbols.gmt",
                              sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

allSynTable <- read.table("/mnt/DATA5/tmp/kev/misc/20210621fearonAllSyn.bed", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

h_exon_boundaries <- read.table("/mnt/DATA5/tmp/kev/misc/20210617hg38ExonBoundaries.txt", sep = "\t",
                                stringsAsFactors = FALSE, header = TRUE)

m_exon_boundaries <- read.table("/mnt/DATA5/tmp/kev/misc/20210617mm10ExonBoundaries.txt", sep = "\t",
                                stringsAsFactors = FALSE, header = TRUE)

segZFilt <- read.table("/mnt/DATA5/tmp/kev/misc/20210618testFearonSegFilt.txt", sep = "\t",
                       header = TRUE, stringsAsFactors = FALSE)



### getting coordinates for all genes human and mouse
# h_exon_boundaries <- NULL
# for (i in unique(hg38biomartTable$external_gene_name)) {
#   tmpTable  <- hg38biomartTable[which(hg38biomartTable$external_gene_name == i),]
#   tmpChr <- unique(tmpTable$chromosome_name)
#   tmpStart <- min(tmpTable$exon_chrom_start)
#   tmpEnd <- max(tmpTable$exon_chrom_end)
#   tmpVector <- c("gene" = i, "chrom" = tmpChr, "start" = tmpStart, "end" = tmpEnd)
#   h_exon_boundaries <- rbind(h_exon_boundaries, tmpVector)
# }
# 
# h_exon_boundaries <- data.frame(h_exon_boundaries, stringsAsFactors = FALSE)
# write.table(h_exon_boundaries,"/mnt/DATA5/tmp/kev/misc/20210617hg38ExonBoundaries.txt", sep = "\t",
#             quote = FALSE, row.names = FALSE, col.names = TRUE)
# 
# m_exon_boundaries <- NULL
# for (i in unique(mm10biomartTable$external_gene_name)) {
#   tmpTable  <- mm10biomartTable[which(mm10biomartTable$external_gene_name == i),]
#   tmpChr <- unique(tmpTable$chromosome_name)
#   tmpStart <- min(tmpTable$exon_chrom_start)
#   tmpEnd <- max(tmpTable$exon_chrom_end)
#   tmpVector <- c("gene" = i, "chrom" = tmpChr, "start" = tmpStart, "end" = tmpEnd)
#   m_exon_boundaries <- rbind(m_exon_boundaries, tmpVector)
# }
# 
# m_exon_boundaries <- data.frame(m_exon_boundaries, stringsAsFactors = FALSE)
# write.table(m_exon_boundaries,"/mnt/DATA5/tmp/kev/misc/20210617mm10ExonBoundaries.txt", sep = "\t",
#             quote = FALSE, row.names = FALSE, col.names = TRUE)





# the general idea is to pull one gene list per syntenic block or change from mouse
# potential problems later on: current process works great for single sample
# where it's easy to disingtuish which boundary one would want to find synteny for

# potential solution is to summarize frequency data in order to fine regions
# to perform the look ups on i.e from mouse to human


hGrange <- GRanges(seqnames = h_exon_boundaries$chrom,
                   IRanges(start = h_exon_boundaries$start,
                           end = h_exon_boundaries$end))
mGrange <- GRanges(seqnames = m_exon_boundaries$chrom,
                   IRanges(start = m_exon_boundaries$start,
                           end = m_exon_boundaries$end))
colnames(allSynTable) <- c("h_chr", "h_start", "h_end","m_chr",
                           "m_start", "m_end", "m_freq","h_freq")

#instead of doing it per-sample we'll do per mouse chromosome
#i.e find contigous regions within a the mouse chromosome 

geneList <- NULL
for (i in unique(allSynTable$m_chr)) {
  
  tmpSynDf <- allSynTable[which(allSynTable$m_chr == i),]
  
  for (j in unique(tmpSynDf$m_freq)) {
    tmpSegDf <- tmpSynDf[which(tmpSynDf$m_freq == j),]
    synGrangeM <- GRanges(seqnames =  str_remove(tmpSegDf$m_chr[1],"m_chr"),
                          IRanges(start = min(tmpSegDf$m_start),
                                  end = max(tmpSegDf$m_end)))
    synGrangeH <- GRanges(seqnames = str_remove(tmpSegDf$h_chr, "h_chr"),
                          IRanges(start = tmpSegDf$h_start,
                                  end = tmpSegDf$h_end))
    
    hGenes <- h_exon_boundaries$gene[subjectHits(findOverlaps(synGrangeH, hGrange))]
    mGenes <- m_exon_boundaries$gene[subjectHits(findOverlaps(synGrangeM, mGrange))]
    hFreq <- tmpSegDf$h_freq[queryHits(findOverlaps(synGrangeH, hGrange))]
    mFreq <- tmpSegDf$m_freq[queryHits(findOverlaps(synGrangeM, mGrange))]
    
    tmpGeneLists <- c("m_chr" = tmpSegDf$m_chr[1], "m_start" = min(tmpSegDf$m_start),
                      "m_end" = max(tmpSegDf$m_end),
                      "h_gene" = paste(hGenes, collapse = ","),
                      "m_gene" = paste(mGenes, collapse = ","),
                      "h_freq" = paste(hFreq, collapse = ","),
                      "m_freq" = paste(mFreq, collapse = ","))
    geneList <- rbind(geneList, tmpGeneLists)
  }
}



rownames(geneList) <- NULL
geneList <- data.frame(geneList, stringsAsFactors = FALSE)

# convert mouse gene names to human for enrichment analysis

geneList$convert <- "empty"
for (i in 1:nrow(geneList)) {
  tmpMGene <- unlist(str_split(geneList$m_gene[i], ","))
  tmpConvert <- NULL
  for (j in tmpMGene) {
    tmpVar <- geneNameDf$external_gene_name[which(geneNameDf$mmusculus_homolog_associated_gene_name == j)]
    if (length(tmpVar) == 0) {
      tmpVar <- toupper(j)
    } else if(length(tmpVar) > 1){
      tmpVar <- tmpVar[1]
    }
    tmpConvert <- c(tmpConvert, tmpVar)
  }
  geneList$convert[i] <- paste(tmpConvert, collapse = ",")
}


human_hallmarks_list <- NULL
for (i in 1:nrow(human_hallmarks)) {
  tmpList <- list(human_hallmarks[i, 3:ncol(human_hallmarks)])
  names(tmpList) <- human_hallmarks[i, 1]
  human_hallmarks_list[[i]] <- tmpList
}

# easiest way is to get list of all genes and their respective value
# and do enrichment from that. from that graphing it will only need 
# to show each pathway once

all_h_gene_symbol <- unlist(str_split(paste(geneList$h_gene, collapse = ","), ","))
all_m_gene_symbol <- unlist(str_split(paste(geneList$convert, collapse = ","), ","))

all_h_value <- as.numeric(unlist(str_split(paste(geneList$h_freq, collapse = ","), ",")))
all_m_value <- as.numeric(unlist(str_split(paste(geneList$m_freq, collapse = ","),",")))

pathwayRes <- NULL
for(i in 1:length(human_hallmarks_list)){
  tmpHallList <- unlist(human_hallmarks_list[[i]])
  # idea is how to deal with uneven sets and what if genes are not present?
  # simple answer is to treat it as zero so i get equal vectors
  tmpMIdx <- which(all_m_gene_symbol %in% tmpHallList)
  tmpHIdx <- which(all_h_gene_symbol %in% tmpHallList)
  
  tmpMGene <- all_m_value[tmpMIdx]
  names(tmpMGene) <- all_m_gene_symbol[tmpMIdx]
  tmpMGene <- tmpMGene[-which(duplicated(names(tmpMGene)))]
  
  tmpHGene <- all_h_value[tmpHIdx]
  names(tmpHGene) <- all_h_gene_symbol[tmpHIdx]
  tmpHGene <- tmpHGene[-which(duplicated(names(tmpHGene)))]
  
  # using union to help create ordered empty vectors to compare
  tmpPathway <- data.frame("gene" = union(all_m_gene_symbol[tmpMIdx], all_h_gene_symbol[tmpHIdx]),
                           "h" = 0,
                           "m" = 0)
  tmpPathway$h[match(names(tmpHGene), tmpPathway$gene)] <- tmpHGene
  tmpPathway$m[match(names(tmpMGene), tmpPathway$gene)] <- tmpMGene
  tmpPathway$h <- ifelse(abs(tmpPathway$h) > 0.1,  tmpPathway$h, 0)
  tmpPathway$m <- ifelse(abs(tmpPathway$m) > 0.1,  tmpPathway$m, 0)
  
  tmpPathway2 <- tmpPathway
  tmpPathway2$h <- ifelse(tmpPathway2$h > 0,  1, tmpPathway2$h)
  tmpPathway2$h <- ifelse(tmpPathway2$h < 0,  -1, tmpPathway2$h)
  tmpPathway2$m <- ifelse(tmpPathway2$m  > 0, 1, tmpPathway2$m)
  tmpPathway2$m <- ifelse(tmpPathway2$m  < 0, -1, tmpPathway2$m)
  
  exactMat <- matrix(ncol = 2, nrow = 2)
  exactMat[1,1] <- table(tmpPathway2$h)[3]
  exactMat[2,1] <- table(tmpPathway2$h)[1]
  exactMat[1,2] <- table(tmpPathway2$m)[3]
  exactMat[2,2] <- table(tmpPathway2$m)[1]
  
  exactMat[which(is.na(exactMat))] <- 0
  
  fisherRes <- fisher.test(exactMat)
  
  corTest <- cor.test(tmpPathway$h, tmpPathway$m)
  #jac <- jaccard(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
  #tani <- tanimoto(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
  cosi <- cosine_dist(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
  ruzi <- ruzicka(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
  euc <- euclidean(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
  pathw <- names(human_hallmarks_list[[i]])
  pathwayRes <- rbind(pathwayRes, c("pathway" = pathw, corTest$estimate,
                                    "pval" = corTest$p.value, "cosine" = cosi,
                                    "ruzicka" = ruzi, "euclidean" = euc,
                                    "fisher" = fisherRes$p.value))
}

pathwayRes <- data.frame(pathwayRes, stringsAsFactors = FALSE)
pathwayRes$cor <- as.numeric(pathwayRes$cor)
pathwayRes$pval <- as.numeric(pathwayRes$pval)
pathwayRes_filt <- pathwayRes[which(abs(pathwayRes$cor) > .1),]
pathwayRes_filt <- pathwayRes_filt[order(pathwayRes_filt$cor),]
pathwayRes_filt$fill <- ifelse(pathwayRes_filt$cor > 0, "firebrick1", "lightblue")
pathwayRes_filt$pathway <- factor(pathwayRes_filt$pathway, levels = unique(pathwayRes_filt$pathway))

pdf("/mnt/DATA5/tmp/kev/misc/20210622syntenypathway.pdf", useDingbats = FALSE)
ggplot(pathwayRes_filt, aes(x = pathway, y = cor)) + geom_bar(stat = "identity", color = pathwayRes_filt$fill) + 
  ylim(c(-0.4, 0.4)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Pathway correlations based on syntenic aneuploidy") + xlab("Hallmakrs pathways") + ylab("Pearson corr")
dev.off()

# below shows duplicated genes seems to be the same sign, but slightly different magnitudes
# for now take random one - shouldn't change much
# tmp <- all_m_value[tmpMIdx]
# names(tmp) <- all_m_gene_symbol[tmpMIdx]
# tmp2 <- tmp[which(names(tmp) %in% names(tmp)[which(duplicated(names(tmp)))])]
# tmp2[order(names(tmp2))]


### function for pathwayRes - change as tests are added or remove
getGeneList <- function(allSynTable){
  # extracts gene list from synteny table
  colnames(allSynTable) <- c("h_chr", "h_start", "h_end","m_chr",
                             "m_start", "m_end", "m_freq","h_freq")
  geneList <- NULL
  for (i in unique(allSynTable$m_chr)) {
    
    tmpSynDf <- allSynTable[which(allSynTable$m_chr == i),]
    
    for (j in unique(tmpSynDf$m_freq)) {
      tmpSegDf <- tmpSynDf[which(tmpSynDf$m_freq == j),]
      synGrangeM <- GRanges(seqnames =  str_remove(tmpSegDf$m_chr[1],"m_chr"),
                            IRanges(start = min(tmpSegDf$m_start),
                                    end = max(tmpSegDf$m_end)))
      synGrangeH <- GRanges(seqnames = str_remove(tmpSegDf$h_chr, "h_chr"),
                            IRanges(start = tmpSegDf$h_start,
                                    end = tmpSegDf$h_end))
      
      hGenes <- h_exon_boundaries$gene[subjectHits(findOverlaps(synGrangeH, hGrange))]
      mGenes <- m_exon_boundaries$gene[subjectHits(findOverlaps(synGrangeM, mGrange))]
      hFreq <- tmpSegDf$h_freq[queryHits(findOverlaps(synGrangeH, hGrange))]
      mFreq <- tmpSegDf$m_freq[queryHits(findOverlaps(synGrangeM, mGrange))]
      
      tmpGeneLists <- c("m_chr" = tmpSegDf$m_chr[1], "m_start" = min(tmpSegDf$m_start),
                        "m_end" = max(tmpSegDf$m_end),
                        "h_gene" = paste(hGenes, collapse = ","),
                        "m_gene" = paste(mGenes, collapse = ","),
                        "h_freq" = paste(hFreq, collapse = ","),
                        "m_freq" = paste(mFreq, collapse = ","))
      geneList <- rbind(geneList, tmpGeneLists)
    }
  }
  
  
  rownames(geneList) <- NULL
  geneList <- data.frame(geneList, stringsAsFactors = FALSE)
  
  # convert mouse gene names to human for enrichment analysis
  
  geneList$convert <- "empty"
  for (i in 1:nrow(geneList)) {
    tmpMGene <- unlist(str_split(geneList$m_gene[i], ","))
    tmpConvert <- NULL
    for (j in tmpMGene) {
      tmpVar <- geneNameDf$external_gene_name[which(geneNameDf$mmusculus_homolog_associated_gene_name == j)]
      if (length(tmpVar) == 0) {
        tmpVar <- toupper(j)
      } else if(length(tmpVar) > 1){
        tmpVar <- tmpVar[1]
      }
      tmpConvert <- c(tmpConvert, tmpVar)
    }
    geneList$convert[i] <- paste(tmpConvert, collapse = ",")
  }
  return(geneList)
}


enrichmentStats <- function(geneList, pathwayList){
  # pathway list should be processed 
  pathwayRes <- NULL
  
  all_h_gene_symbol <- unlist(str_split(paste(geneList$h_gene, collapse = ","), ","))
  all_m_gene_symbol <- unlist(str_split(paste(geneList$convert, collapse = ","), ","))
  all_h_value <- as.numeric(unlist(str_split(paste(geneList$h_freq, collapse = ","), ",")))
  all_m_value <- as.numeric(unlist(str_split(paste(geneList$m_freq, collapse = ","),",")))
  
  
  
  for(i in 1:length(pathwayList)){
    
    tmpHallList <- unlist(pathwayList[[i]])
    tmpMIdx <- which(all_m_gene_symbol %in% tmpHallList)
    tmpHIdx <- which(all_h_gene_symbol %in% tmpHallList)
    
    tmpMGene <- all_m_value[tmpMIdx]
    names(tmpMGene) <- all_m_gene_symbol[tmpMIdx]
    tmpMGene <- tmpMGene[-which(duplicated(names(tmpMGene)))]
    
    tmpHGene <- all_h_value[tmpHIdx]
    names(tmpHGene) <- all_h_gene_symbol[tmpHIdx]
    tmpHGene <- tmpHGene[-which(duplicated(names(tmpHGene)))]
    
    # using union to help create ordered empty vectors to compare
    tmpPathway <- data.frame("gene" = union(all_m_gene_symbol[tmpMIdx], all_h_gene_symbol[tmpHIdx]),
                             "h" = 0,
                             "m" = 0)
    tmpPathway$h[match(names(tmpHGene), tmpPathway$gene)] <- tmpHGene
    tmpPathway$m[match(names(tmpMGene), tmpPathway$gene)] <- tmpMGene
    tmpPathway$h <- ifelse(abs(tmpPathway$h) > 10,  tmpPathway$h, 0)
    tmpPathway$m <- ifelse(abs(tmpPathway$m) > 10,  tmpPathway$m, 0)
    
    tmpPathway2 <- tmpPathway
    tmpPathway2$h <- ifelse(tmpPathway2$h > 0,  1, tmpPathway2$h)
    tmpPathway2$h <- ifelse(tmpPathway2$h < 0,  -1, tmpPathway2$h)
    tmpPathway2$m <- ifelse(tmpPathway2$m  > 0, 1, tmpPathway2$m)
    tmpPathway2$m <- ifelse(tmpPathway2$m  < 0, -1, tmpPathway2$m)
    
    exactMat <- matrix(ncol = 2, nrow = 2)
    exactMat[1,1] <- table(tmpPathway2$h)[3]
    exactMat[2,1] <- table(tmpPathway2$h)[1]
    exactMat[1,2] <- table(tmpPathway2$m)[3]
    exactMat[2,2] <- table(tmpPathway2$m)[1]
    
    exactMat[which(is.na(exactMat))] <- 0
    
    fisherRes <- fisher.test(exactMat)
    
    corTest <- cor.test(tmpPathway$h, tmpPathway$m)
    #jac <- jaccard(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
    #tani <- tanimoto(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
    cosi <- cosine_dist(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
    ruzi <- ruzicka(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
    euc <- euclidean(P = tmpPathway$h, Q = tmpPathway$m, testNA = TRUE)
    pathw <- names(human_hallmarks_list[[i]])
    pathwayRes <- rbind(pathwayRes, c("pathway" = pathw, corTest$estimate,
                                      "pval" = corTest$p.value, "cosine" = cosi,
                                      "ruzicka" = ruzi, "euclidean" = euc,
                                      "fisher" = fisherRes$p.value))
    print(paste0(tmpPathway$h, collapse = ","))
    print(paste0(tmpPathway$m, collapse = ","))
  }
  
  pathwayRes <- data.frame(pathwayRes, stringsAsFactors = FALSE)
  pathwayRes$cor <- as.numeric(pathwayRes$cor)
  pathwayRes$pval <- as.numeric(pathwayRes$pval)
  # pathwayRes_filt <- pathwayRes[which(abs(pathwayRes$cor) > .1),]
  # pathwayRes_filt <- pathwayRes_filt[order(pathwayRes_filt$cor),]
  # pathwayRes_filt$fill <- ifelse(pathwayRes_filt$cor > 0, "firebrick1", "lightblue")
  # pathwayRes_filt$pathway <- factor(pathwayRes_filt$pathway, levels = unique(pathwayRes_filt$pathway))
  # 
  
  
  return(pathwayRes)
}


### make sure to zero out things less than 10% to test it out

geneListCoadArm <- getGeneList(coad_arm_allSynTable)
pathwayResCoadArm <- enrichmentStats(geneListCoadArm, human_hallmarks_list)

geneListCoadCna <- getGeneList(coad_cna_allSynTable)
pathwayResCoadCna <- enrichmentStats(geneListCoadCna, human_hallmarks_list)

