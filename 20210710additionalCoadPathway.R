finalTable <- read.table("/mnt/DATA5/tmp/kev/misc/20210628panCancerAllSegReduced_hg19.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t")

finalTable$length <- finalTable$End - finalTable$Start
finalTable2 <- finalTable[which(finalTable$Chromosome %in% c(1:22)),]
armRes <- separateSegments_tcga(finalTable2, tcga_ploidy)
armTable <- armRes[[1]]
cnaTable <- armRes[[2]]


tcga_arm <- armTable
tcga_cna <- cnaTable

coad_arm <- tcga_arm[which(tcga_arm$Sample %in% coad$Sample.ID),]
coad_cna <- tcga_cna[which(tcga_cna$Sample %in% coad$Sample.ID),]

kap_arm <- tcga_arm[grep(paste(coad_KAP$samples, collapse = "|"), tcga_arm$Sample),]
kap_cna <- tcga_cna[grep(paste(coad_KAP$samples, collapse = "|"), tcga_cna$Sample),]


coad_arm <- coad_arm[-grep(paste(coad_KAP$samples, collapse = "|"), coad_arm$Sample),]
coad_cna <- coad_cna[-grep(paste(coad_KAP$samples, collapse = "|"), coad_cna$Sample),]


coad_arm2 <- data.frame("sampleID" = coad_arm$Sample, "chrom" = coad_arm$Chromosome,
                        "start.pos" = coad_arm$Start, "end.pos" = coad_arm$End,
                        "n.probes" = NA, "mean" = coad_arm$newTotalCN, stringsAsFactors = FALSE)

coad_cna2 <- data.frame("sampleID" = coad_cna$Sample, "chrom" = coad_cna$Chromosome,
                        "start.pos" = coad_cna$Start, "end.pos" = coad_cna$End,
                        "n.probes" = NA, "mean" = coad_cna$newTotalCN, stringsAsFactors = FALSE)

kap_arm2 <- data.frame("sampleID" = kap_arm$Sample, "chrom" = kap_arm$Chromosome,
                       "start.pos" = kap_arm$Start, "end.pos" = kap_arm$End,
                       "n.probes" = NA, "mean" = kap_arm$newTotalCN, stringsAsFactors = FALSE)
kap_cna2 <- data.frame("sampleID" = kap_cna$Sample, "chrom" = kap_cna$Chromosome,
                       "start.pos" = kap_cna$Start, "end.pos" = kap_cna$End,
                       "n.probes" = NA, "mean" = kap_cna$newTotalCN, stringsAsFactors = FALSE)



coad_arm_freq_out <- getFreqData(coad_arm2)
ampDels_coad_arm <- ampsDels(coad_arm_freq_out)
coad_arm_amp_bed <- reducingFreqBed(ampDels_coad_arm[[1]], ampDels_coad_arm[[2]])
coad_arm_del_bed <- reducingFreqBed(ampDels_coad_arm[[3]], ampDels_coad_arm[[4]])

coad_cna_freq_out <- getFreqData(coad_cna2)
ampDels_coad_cna <- ampsDels(coad_cna_freq_out)
coad_cna_amp_bed <- reducingFreqBed(ampDels_coad_cna[[1]], ampDels_coad_cna[[2]])
coad_cna_del_bed <- reducingFreqBed(ampDels_coad_cna[[3]], ampDels_coad_cna[[4]])


kap_arm_freq_out <- getFreqData(kap_arm2)
ampDels_kap_arm <- ampsDels(kap_arm_freq_out)
kap_arm_amp_bed <- reducingFreqBed(ampDels_kap_arm[[1]], ampDels_kap_arm[[2]])
kap_arm_del_bed <- reducingFreqBed(ampDels_kap_arm[[3]], ampDels_kap_arm[[4]])

kap_cna_freq_out <- getFreqData(kap_cna2)
ampDels_kap_cna <- ampsDels(kap_cna_freq_out)
kap_cna_amp_bed <- reducingFreqBed(ampDels_kap_cna[[1]], ampDels_kap_cna[[2]])
kap_cna_del_bed <- reducingFreqBed(ampDels_kap_cna[[3]], ampDels_kap_cna[[4]])




arm_amp_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210712_1copy_fearon_arm_amp.bed", sep = "\t",
                          stringsAsFactors = FALSE, header = TRUE)
arm_del_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210712_1copy_fearon_arm_del.bed", sep = "\t",
                          stringsAsFactors = FALSE, header = TRUE)
tcga_coad_arm_amp_bed <- coad_arm_amp_bed
tcga_coad_arm_del_bed <- coad_arm_del_bed

coad_arm_allSynTable <- circosFreq(arm_amp_bed, arm_del_bed, tcga_coad_arm_amp_bed, tcga_coad_arm_del_bed,
                                   filename = "20210712fearonCoad_arm")


cna_amp_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210712_1copy_fearon_arm_amp.bed", sep = "\t",
                          stringsAsFactors = FALSE, header = TRUE)
cna_del_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210712_1copy_fearon_cna_del.bed", sep = "\t",
                          stringsAsFactors = FALSE, header = TRUE)

tcga_coad_cna_amp_bed <- coad_cna_amp_bed
tcga_coad_cna_del_bed <- coad_cna_del_bed

coad_cna_allSynTable <- circosFreq(cna_amp_bed, cna_del_bed, tcga_coad_cna_amp_bed, tcga_coad_cna_del_bed,
                                   filename = "20210712fearonCoad_cna")


kap_coad_cna_amp_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210629_kapCnaAmp.bed", sep = "\t",
                                   stringsAsFactors = FALSE, header = TRUE)
kap_coad_cna_del_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210629_kapCnaDel.bed", sep = "\t",
                                   stringsAsFactors = FALSE, header = TRUE)
kap_coad_arm_amp_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210629_kapArmAmp.bed", sep = "\t",
                                   stringsAsFactors = FALSE, header = TRUE)
kap_coad_arm_del_bed <- read.table("/mnt/DATA5/tmp/kev/misc/20210629_kapCnaDel.bed", sep = "\t",
                                   stringsAsFactors = FALSE, header = TRUE)

kap_arm_allSynTable <- circosFreq(arm_amp_bed, arm_del_bed, kap_coad_arm_amp_bed, kap_coad_arm_del_bed,
                                  filename = "20210712fearonKap_arm")

kap_cna_allSynTable <- circosFreq(cna_amp_bed, cna_del_bed, kap_coad_cna_amp_bed, kap_coad_cna_del_bed,
                                  filename = "20210712fearonKap_cna")

h_exon_boundaries <- read.table("/mnt/DATA5/tmp/kev/misc/20210617hg38ExonBoundaries.txt", sep = "\t",
                                stringsAsFactors = FALSE, header = TRUE)

m_exon_boundaries <- read.table("/mnt/DATA5/tmp/kev/misc/20210617mm10ExonBoundaries.txt", sep = "\t",
                                stringsAsFactors = FALSE, header = TRUE)



hGrange <- GRanges(seqnames = h_exon_boundaries$chrom,
                   IRanges(start = as.numeric(h_exon_boundaries$start),
                           end = as.numeric(h_exon_boundaries$end)))
mGrange <- GRanges(seqnames = m_exon_boundaries$chrom,
                   IRanges(start = as.numeric(m_exon_boundaries$start),
                           end = as.numeric(m_exon_boundaries$end)))


human_hallmarks_list <- NULL
for (i in 1:nrow(human_hallmarks)) {
  tmpList <- list(human_hallmarks[i, 3:ncol(human_hallmarks)])
  names(tmpList) <- human_hallmarks[i, 1]
  human_hallmarks_list[[i]] <- tmpList
}


geneListCoadArm <- getGeneList(coad_arm_allSynTable)
pathwayResCoadArm <- enrichmentStats(geneListCoadArm, human_hallmarks_list)

geneListCoadCna <- getGeneList(coad_cna_allSynTable)
pathwayResCoadCna <- enrichmentStats(geneListCoadCna, human_hallmarks_list)

geneListKapArm <- getGeneList(kap_arm_allSynTable)
pathwayResKapArm <- enrichmentStats(geneListKapArm, human_hallmarks_list)

geneListKapCna <- getGeneList(kap_cna_allSynTable)
pathwayResKapCna <- enrichmentStats(geneListKapCna, human_hallmarks_list)



pathwayResCoadArm_filt <- pathwayResCoadArm[which(abs(pathwayResCoadArm$cor) > .1),]
pathwayResCoadArm_filt$pathway <- factor(pathwayResCoadArm_filt $pathway, levels = unique(pathwayResCoadArm_filt$pathway))
#pathwayResCoadArm_filt$fill <- ifelse(pathwayResCoadArm_filt$cor > 0, "firebrick1", "lightblue")
#pathwayResCoadArm_filt$fill <- ifelse(pathwayResCoadArm_filt$fisher < 0.05/50, pathwayResCoadArm_filt$fill, "white")
# a <- ggplot(pathwayResCoadArm_filt, aes(x = pathway, y = cor)) + geom_bar(stat = "identity", color = pathwayResCoadArm_filt$fill) +
#   ylim(c(-1, 1)) + theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4), plot.title = element_text(hjust = 0.5)) +
#   ggtitle("Coad Arm") + xlab("") + ylab("Pearson corr")
a <- ggplot(pathwayResCoadArm_filt, aes(x = pathway, y = cor)) + geom_bar(stat = "identity") +
  ylim(c(-1, 1)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Coad Arm") + xlab("") + ylab("Pearson corr")




pathwayResCoadCna_filt <- pathwayResCoadCna[which(abs(pathwayResCoadCna$cor) > .1),]
pathwayResCoadCna_filt$pathway <- factor(pathwayResCoadCna_filt $pathway, levels = unique(pathwayResCoadCna_filt$pathway))
#pathwayResCoadCna_filt$fill <- ifelse(pathwayResCoadCna_filt$cor > 0, "firebrick1", "lightblue")
#pathwayResCoadCna_filt$fill <- ifelse(pathwayResCoadCna_filt$fisher < 0.05/50, pathwayResCoadCna_filt$fill, "white")
# b <- ggplot(pathwayResCoadCna_filt, aes(x = pathway, y = cor)) + geom_bar(stat = "identity", color = pathwayResCoadCna_filt$fill) + 
#   ylim(c(-1, 1)) + theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4), plot.title = element_text(hjust = 0.5)) +
#   ggtitle("Coad CNA") + xlab("") + ylab("")

b <- ggplot(pathwayResCoadCna_filt, aes(x = pathway, y = cor)) + geom_bar(stat = "identity") + 
  ylim(c(-1, 1)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Coad CNA") + xlab("") + ylab("")

pathwayResKapArm_filt <- pathwayResKapArm[which(abs(pathwayResKapArm$cor) > .1),]
pathwayResKapArm_filt$pathway <- factor(pathwayResKapArm_filt $pathway, levels = unique(pathwayResKapArm_filt$pathway))
# pathwayResKapArm_filt$fill <- ifelse(pathwayResKapArm_filt$cor > 0, "firebrick1", "lightblue")
# pathwayResKapArm_filt$fill <- ifelse(pathwayResKapArm_filt$fisher < 0.05/50, pathwayResKapArm_filt$fill, "white")
# c <- ggplot(pathwayResKapArm_filt, aes(x = pathway, y = cor)) + geom_bar(stat = "identity", color = pathwayResKapArm_filt$fill) + 
#   ylim(c(-1, 1)) + theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4), plot.title = element_text(hjust = 0.5)) +
#   ggtitle("Kap Arm") + xlab("Hallmark pathways") + ylab("Pearson corr")

c <- ggplot(pathwayResKapArm_filt, aes(x = pathway, y = cor)) + geom_bar(stat = "identity") + 
  ylim(c(-1, 1)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Kap Arm") + xlab("Hallmark pathways") + ylab("Pearson corr")

pathwayResKapCna_filt <- pathwayResKapCna[which(abs(pathwayResKapCna$cor) > .1),]
pathwayResKapCna_filt$pathway <- factor(pathwayResKapCna_filt $pathway, levels = unique(pathwayResKapCna_filt$pathway))
# pathwayResKapCna_filt$fill <- ifelse(pathwayResKapCna_filt$cor > 0, "firebrick1", "lightblue")
# pathwayResKapCna_filt$fill <- ifelse(pathwayResKapCna_filt$fisher < 0.05/50, pathwayResKapCna_filt$fill, "white")
# d <- ggplot(pathwayResKapCna_filt, aes(x = pathway, y = cor)) + geom_bar(stat = "identity", color = pathwayResKapCna_filt$fill) + 
#   ylim(c(-1, 1)) + theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4), plot.title = element_text(hjust = 0.5)) +
#   ggtitle("Kap CNA") + xlab("Hallmark pathways") + ylab("")

d <- ggplot(pathwayResKapCna_filt, aes(x = pathway, y = cor)) + geom_bar(stat = "identity") + 
  ylim(c(-1, 1)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4), plot.title = element_text(hjust = 0.5)) +
  ggtitle("Kap CNA") + xlab("Hallmark pathways") + ylab("")


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210712syntenypathway.pdf", useDingbats = FALSE)
grid.arrange(a,b,c,d, ncol = 2)
dev.off()





