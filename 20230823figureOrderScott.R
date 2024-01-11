### order of the figure progression according to scott
### (1) freq plots do for both cna and aneuploidy


a <- freqPlotv2(armGisticTp53_amp_bed, armGisticTp53_del_bed,
                main = "hgsoc n = 368", speciesType = "human")
b <- freqPlotv2(allMouseAneu_bprn_amp_bed, allMouseAneu_bprn_del_bed,
                main = "bprn n = 78", speciesType = "mouse")
c <- freqPlotv2(allMouseAneu_bpn_amp_bed, allMouseAneu_bpn_del_bed,
                main = "bpn n = 7", speciesType = "mouse")
d <- freqPlotv2(allMouseAneu_bpp_amp_bed, allMouseAneu_bpp_del_bed,
                main = "bpp n = 11", speciesType = "mouse")

grid.arrange(a, b, c, d, nrow = 4)

e <- freqPlotv2(armGisticCoadApc_amp_bed, armGisticCoadApc_del_bed,
                main = "coadread n = 349", speciesType = "human")
g <- freqPlotv2(allMouseAneu_coadAdenoCar_amp_bed, allMouseAneu_coadAdenoCar_del_bed,
                main = "mouse coadread n = 27", speciesType = "mouse")
f <- freqPlotv2(allMouseAneu_coadAdeno_amp_bed, allMouseAneu_coadAdeno_del_bed,
                main = "mouse adenoma n = 11", speciesType = "mouse")

grid.arrange(e, f, g, nrow = 3)

pdf(file = "/mnt/DATA5/tmp/kev/misc/20230908humanCancersFreq.pdf", useDingbats = TRUE,
    width = 12, height = 5)
grid.arrange(a, e, nrow = 2)
dev.off()
pdf(file = "/mnt/DATA5/tmp/kev/misc/20230908mouseCancersFreq.pdf", useDingbats = TRUE,
    width = 12, height = 10)
grid.arrange(b,c,d, f, g, nrow = 5)
dev.off()

pdf(file = "/mnt/DATA5/tmp/kev/misc/20231003allCoadreadFreq.pdf", useDingbats = TRUE,
    width = 10, height = 8)
grid.arrange(e, g, f, nrow = 3)
dev.off()




a <- freqPlotv2(gisticOvCnaamp_bed, gisticOvCnadel_bed,
                main = "hgsoc n = 368", speciesType = "human")
b <- freqPlotv2(allMouseCna_bprn_amp_bed, allMouseCna_bprn_del_bed,
                main = "bprn n = 78", speciesType = "mouse")
c <- freqPlotv2(allMouseCna_bpn_amp_bed, allMouseCna_bpn_del_bed,
                main = "bpn n = 7", speciesType = "mouse")
d <- freqPlotv2(allMouseCna_bpp_amp_bed, allMouseCna_bpp_del_bed,
                main = "bpp n = 11", speciesType = "mouse")

grid.arrange(a, b, c, d, nrow = 4)


### 20231024: creating dummy frequency plot for methodology figure

dummyFreqAmp <- armGisticCoadApc_amp_bed
dummyFreqAmp$Freq[-which(dummyFreqAmp$Chr == "13")] <- 0
dummyFreqDel <- armGisticCoadApc_del_bed
dummyFreqDel$Freq[-which(dummyFreqDel$Chr == "18")] <- 0

pdf(file = "/mnt/DATA5/tmp/kev/misc/20231024exampleFreq.pdf", useDingbats = TRUE,
    width = 10, height = 4)
freqPlotv2(dummyFreqAmp, dummyFreqDel,
           main = "example", speciesType = "human")
dev.off()

### 20231024: single chromosome for method
circosFreqSingleChr(allMouseAneu_coadAdenoCar_amp_bed, allMouseAneu_coadAdenoCar_del_bed, armGisticCoadApc_amp_bed,
                    armGisticCoadApc_del_bed, filename = "20231024coadreadSingleChr13", ref = "human",
                    chromosome = "h_chr13")


### 20231024: human and mouse freq histogram for method

human_freq_10000 <- rnorm(n = 10000, mean = 0.5, sd = 0.1)
mouse_freq_10000 <- rnorm(n = 10000, mean = 0.3, sd = 0.1)

human_freq_10000 <- c(human_freq_10000, 0, 1)
mouse_freq_10000 <- c(mouse_freq_10000, 0, 1)


pdf(file = "/mnt/DATA5/tmp/kev/misc/20231024exampleHumanPerm.pdf", useDingbats = TRUE,
    width = 10, height = 4)
ggplot() + geom_histogram(aes(human_freq_10000), bins = 100, color="darkgreen", fill="white") +
  geom_vline(xintercept = c(0.6), linetype = "dashed") + 
  geom_vline(xintercept = c(0.75), linetype = "dashed", color = "darkred") + scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(breaks = seq(0, 500, 100)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")) + 
  xlab("Human Frequency 10,000 Permutations") + ylab("Frequency Count")
dev.off()

pdf(file = "/mnt/DATA5/tmp/kev/misc/20231024exampleMousePerm.pdf", useDingbats = TRUE,
    width = 10, height = 4)
ggplot() + geom_histogram(aes(mouse_freq_10000), bins = 100, color="darkblue", fill="white") +
  geom_vline(xintercept = c(0.7), linetype = "dashed") +
  geom_vline(xintercept = c(0.6), linetype = "dashed", color = "darkred") +scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(breaks = seq(0, 500, 100)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")) + 
  xlab("Mouse Freqeuncy 10,000 Permutations") + ylab("Frequency Count")
dev.off()

### fix the cna later
### do fga plots (b)

### need to do parallel loop to calculate all fga for graph

dir <- "/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/"
setwd(dir)
allDirs <- list.dirs()
allDirs <- allDirs[grep("TP", allDirs)]

listOfFiles <- allDir

allCans <- c("BLCA", "BRCA", "CESC", "COADREAD", "ESCA", "GBM",
             "GBMLGG", "HNSC", "KIPAN", "KIRC", "LGG", "LIHC",
             "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD",
             "SARC", "STES", "TGCT", "THYM")


allCancerFga <- NULL
for (j in seq_along(listOfFiles)) {
  setwd(dir)
  setwd(listOfFiles[j])
  
  tmpBroad <- read.table("broad_values_by_arm.txt", sep = "\t", stringsAsFactors = FALSE,
                         header = TRUE, check.names = FALSE)
  colnames(tmpBroad)[2:ncol(tmpBroad)] <- gsub("(^.*?-.{3}?)-.*", "\\1",
                                               colnames(tmpBroad)[2:ncol(tmpBroad)])
  if (allCans[j] == "COADREAD") {
    tmpBroad <- tmpBroad[, c(1, grep(paste0(annotationTable$Sample.ID, collapse = "|"), colnames(tmpBroad)))] 
  }
  tmpBroad <- tmpBroad[-grep("X", tmpBroad$`Chromosome Arm`), ]
  tmpBroadMat <- tmpBroad[, 2:ncol(tmpBroad)]
  tmpBroad2 <- cbind(hg19ArmLocations[, c("chrStripped", "start", "end")], tmpBroadMat)
  tmpBroadMelt <- melt(tmpBroad2, id.vars = c("chrStripped", "start", "end"))
  tmpBroadMelt2 <- data.frame("sampleID" = tmpBroadMelt$variable, "chrom" = tmpBroadMelt$chrStripped,
                              "start.pos" = tmpBroadMelt$start, "end.pos" = tmpBroadMelt$end, 
                              "n.probes" = NA, "mean" = tmpBroadMelt$value)
  tmpBroadMelt2$str <- paste0(tmpBroadMelt2$sampleID, tmpBroadMelt2$chrom, 
                              tmpBroadMelt2$start.pos, tmpBroadMelt2$end.pos)
  tmpBroadMelt2$length <- tmpBroadMelt2$end.pos - tmpBroadMelt2$start.pos
  tmpRes <- fgaCalculator_tcga_aneuploidy(tmpBroadMelt2)
  tmpRes$Cancer <- allCans[j]
  allCancerFga <- rbind(allCancerFga, tmpRes)
}

adenoFga <- fgaCalculator_amp_aneuploidy(allMouseAneu_coadAdeno)
adenoCarFga <- fgaCalculator_amp_aneuploidy(allMouseAneu_coadAdenoCar)
adenoFga$Cancer <- "mm_adeno"
adenoCarFga$Cancer <- "mm_coadread"

allFga <- rbind(allCancerFga, adenoFga, adenoCarFga)
allFga$Cancer <- reorder(allFga$Cancer, allFga$fga, median)
allFga$Cancer <- factor(allFga$Cancer, rev(levels(allFga$Cancer)))
fgaColorVector <- rep("#000000", length(unique(allFga$Cancer)))
fgaColorVector[which(levels(allFga$Cancer) %in% c("COADREAD", "mm_coadread", "mm_adeno" ))] <- "#8B0000"




pdf(file = "/mnt/DATA5/tmp/kev/misc/20231005fgaComp.pdf", useDingbats = TRUE,
    width = 7, height = 4)
ggplot(allFga, aes(x = Cancer, y = fga)) + geom_boxplot(color = fgaColorVector) + 
  ylab("fraction of genome altered (aneuploidy)") + xlab("cancer type/model") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



### genomic plots highlighting the cancer census genes (c)
### should feed in a list of genes to highlight i.e label the important ones
bprnGenes <- c("NRAS", "PIK3CA", "SOX2", "FGFR1", "MYC", "BRCA1", "ERBB2", "NF2")
pdf(file = "/mnt/DATA5/tmp/kev/misc/20230825bprnSimpleGenes.pdf", useDingbats = TRUE,
    width = 10, height = 2)
genePlotSynteny(hgsc_bprn_arm_allSynTable_simple,
                geneList = bprnGenes)
dev.off()

bprnGenes <- c("NRAS", "PIK3CA", "SOX2", "FGFR1", "MYC", "BRCA1", "ERBB2", "NF2")
pdf(file = "/mnt/DATA5/tmp/kev/misc/20230825bprnSimpleGenes.pdf", useDingbats = TRUE,
    width = 10, height = 2)
genePlotSynteny(hgsc_bprn_arm_allSynTable_simple,
                geneList = bprnGenes)
dev.off()


coadGenes <- c("PIK3R1", "APC", "SMAD2", "SMAD4")
pdf(file = "/mnt/DATA5/tmp/kev/misc/20230825coadSimpleGenes.pdf", useDingbats = TRUE,
    width = 10, height = 2)
genePlotSynteny(coad_adenoCar_arm_allSynTable_simple,
                geneList = coadGenes)
dev.off()

### synteny overlap graphs

pdf(file = "/mnt/DATA5/tmp/kev/misc/20230908bprnSyntenyColors.pdf", useDingbats = TRUE,
    width = 10, height = 2)
syntenyOverlap(allMouseAneu_bprn_amp_bed, allMouseAneu_bprn_del_bed,
               armGisticTp53_amp_bed, armGisticTp53_del_bed)
dev.off()

pdf(file = "/mnt/DATA5/tmp/kev/misc/20230908coadSyntenyColors.pdf", useDingbats = TRUE,
    width = 10, height = 2)
syntenyOverlap(allMouseAneu_coadAdenoCar_amp_bed, allMouseAneu_coadAdenoCar_del_bed,
               armGisticCoadApc_amp_bed, armGisticCoadApc_del_bed)
dev.off()

### empty syntenyPlot i.e no freq

hgsc_bprn_arm_allSynTable <- circosFreq2(allMouseAneu_bprn_amp_bed, allMouseAneu_bprn_del_bed, armGisticTp53_amp_bed,
                                         armGisticTp53_del_bed, filename = "20230825emptyCircos",
                                         ref = "human", empty = TRUE)


### scott's metric of synteny - he wants a percent of synteny vs fga

bprnSynFraction <- syntenyFraction(allMouseAneu_bprn_amp_bed, allMouseAneu_bprn_del_bed,
                                   armGisticTp53_amp_bed, armGisticTp53_del_bed,
                                   allMouseAneu_bprn)
bppSynFraction <- syntenyFraction(allMouseAneu_bpp_amp_bed, allMouseAneu_bpp_del_bed,
                                   armGisticTp53_amp_bed, armGisticTp53_del_bed,
                                   allMouseAneu_bpp)
bpnSynFraction <- syntenyFraction(allMouseAneu_bpn_amp_bed, allMouseAneu_bpn_del_bed,
                                   armGisticTp53_amp_bed, armGisticTp53_del_bed,
                                   allMouseAneu_bpn)

coadAdenoSynFraction <- syntenyFraction(allMouseAneu_coadAdeno_amp_bed, allMouseAneu_coadAdeno_del_bed,
                                  armGisticCoadApc_amp_bed, armGisticCoadApc_del_bed,
                                  allMouseAneu_coadAdeno)
coadAdenoCarSynFraction <- syntenyFraction(allMouseAneu_coadAdenoCar_amp_bed, allMouseAneu_coadAdenoCar_del_bed,
                                           armGisticCoadApc_amp_bed, armGisticCoadApc_del_bed,
                                           allMouseAneu_coadAdenoCar)



# save.image(file = "/mnt/DATA4/test_nextflow/20230829allSyntenyEnivorn.RData")
# load(file = "/mnt/DATA4/test_nextflow/20230829allSyntenyEnivorn.RData")

### below shows intersting results b/c in general there are higher fractions of conserved synteny
### between the mouse hgsc and mm coad in terms of recurrent changes >= 10$ of samples
### however synteny fraction is a function of fga, still makes sense b.c our mouse model
### shows fewer CNA events in coad. conclusion is some aneuploidies have function 
### while others don't .... makes some sense, but we already knew that
### good b/c this leads to the idea of which are important in the next set of figures

allSynFractionDf <- rbind(data.frame("type" = "bprn_hgsoc", "fraction" = bprnSynFraction),
                          data.frame("type" = "bpp_hgsoc", "fraction" = bppSynFraction),
                          data.frame("type" = "bpn_hgsoc", "fraction" = bpnSynFraction),
                          data.frame("type" = "mm_adeno", "fraction" = coadAdenoSynFraction),
                          data.frame("type" = "mm_coad", "fraction" = coadAdenoCarSynFraction))
allSynFractionDf2 <- allSynFractionDf
allSynFractionDf$type <- factor(allSynFractionDf$type, levels = c("bpn_hgsoc", "bprn_hgsoc", "bpp_hgsoc", "mm_adeno", "mm_coad"))



pdf(file = "/mnt/DATA5/tmp/kev/misc/20230830fractionSyntenyCon.pdf", useDingbats = TRUE,
    width = 7, height = 4)
ggplot(data = allSynFractionDf) + geom_boxplot(aes(x = type, y = fraction)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("fraction of conserved synteny")
dev.off()


pdf(file = "/mnt/DATA5/tmp/kev/misc/20230830fractionSyntenyCon.pdf", useDingbats = TRUE,
    width = 7, height = 4)
allFraction <- rbind(bprnFga, bppFga, bpnFga, adenoFga, adenoCarFga)
allFraction$synFrac <- allSynFractionDf2$fraction
ggplot(allFraction, aes(x = fga, y = synFrac)) + geom_point() + facet_wrap(~type) + 
  geom_smooth(method = lm, se = FALSE) + 
  scale_x_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1)) + 
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1)) + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"))
dev.off()


cor(bprnFga$fga, bprnSynFraction)
cor(adenoCarFga$fga, coadAdenoCarSynFraction)

### second set of figures would be 





### 20230901 notes from 20230901 meeting

# look at conservation of mouse to human when calculating because less events
# more informative with more splits b/c it helps narrows down which areas or genes to look at
### across cancer types look at aneuploidies without focal changes - as  a primary driver versus a second hit - highly recurrent
### for simple analysis this is the pattern of aneuploidues we expect if these human aneuploidies were important
### stacked barplot of expected chromsoome changes i.e for each human change which mouse chromsomes do you expect to see changed



### didn't look at colorectal b/c there are few gains i.e low level ~13% in pool of 23 samples that were gains
### everything else were losses



armGisticTp53_amp_bed2 <- armGisticTp53_amp_bed
armGisticTp53_del_bed2 <- armGisticTp53_del_bed
armGisticTp53_amp_bed2$type <- "gain"
armGisticTp53_del_bed2$type <- "loss"
hgscFreqsDf <- rbind(armGisticTp53_amp_bed2, armGisticTp53_del_bed2)
hgscFreqsDf$Zscore <- sapply(hgscFreqsDf$Freq, function(x) (x - mean(hgscFreqsDf$Freq))/sd(hgscFreqsDf$Freq))


### by chromosome arms I look at the most frequent changes
armGisticCoadApc_amp_bed2 <- armGisticCoadApc_amp_bed
armGisticCoadApc_del_bed2 <- armGisticCoadApc_del_bed
armGisticCoadApc_amp_bed2$type <- "gain"
armGisticCoadApc_del_bed2$type <- "loss"
coadFreqsDf <- rbind(armGisticCoadApc_amp_bed2, armGisticCoadApc_del_bed2)
coadFreqsDf$Zscore <- sapply(coadFreqsDf$Freq, function(x) (x - mean(coadFreqsDf$Freq))/sd(coadFreqsDf$Freq))

### stacked barplot 
barplotFraction <- NULL
for (i in unique(coad_adenoCar_arm_allSynTable$h_chr)) {
  # i <- unique(hgsc_bprn_arm_allSynTable$h_chr)[1]
  tmpDf <- coad_adenoCar_arm_allSynTable[which(coad_adenoCar_arm_allSynTable$h_chr == i), ]
  tmpDf$length <- tmpDf$h_end - tmpDf$h_start
  for (j in unique(tmpDf$m_chr)) {
    tmpRes <- sum(tmpDf$length[which(tmpDf$m_chr == j)])/sum(tmpDf$length)
    barplotFraction <- rbind(barplotFraction, data.frame("h_chr" = i, "m_chr" = j, "fraction" = tmpRes))
  }
}

barplotFractionFilt <- barplotFraction[which(barplotFraction$fraction > 0.05), ]
barplotFractionFilt$h_chr <- factor(barplotFractionFilt$h_chr, levels = paste0("h_chr", 1:22))
barplotFractionFilt$m_chr <- factor(barplotFractionFilt$m_chr, levels = paste0("m_chr", 1:19))

barplotFractionFiltCoad <- barplotFractionFilt
barplotFractionFiltCoad$color <- "#FFFFFF"
barplotFractionFiltCoad$color[which(barplotFractionFiltCoad$h_chr == "h_chr20")] <- "#8B0000"
barplotFractionFiltCoad$color[which(barplotFractionFiltCoad$h_chr == "h_chr7")] <- "#8B0000"
barplotFractionFiltCoad$color[which(barplotFractionFiltCoad$h_chr == "h_chr13")] <- "#8B0000"
barplotFractionFiltCoad$color[which(barplotFractionFiltCoad$h_chr == "h_chr8")] <- "#8B0000"
barplotFractionFiltCoad$color[which(barplotFractionFiltCoad$h_chr == "h_chr18")] <- "#00008B"
barplotFractionFiltCoad$color[which(barplotFractionFiltCoad$h_chr == "h_chr17")] <- "#00008B"


pdf(file = "/mnt/DATA5/tmp/kev/misc/20230908fractionSyntenyCon.pdf", useDingbats = TRUE)
ggplot(barplotFractionFiltCoad, aes(fill= m_chr, y= fraction, x = h_chr)) + 
  geom_bar(position="stack", stat="identity", color = barplotFractionFiltCoad$color) + 
  scale_fill_manual(values = colorVector) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + 
  ylab("mouse chromosome synteny fraction") + xlab("human chromsome")
dev.off()


coad_adenoCar_arm_allSynTable$tmpString <- paste(coad_adenoCar_arm_allSynTable$h_chr, coad_adenoCar_arm_allSynTable$m_chr)
coadSynChroms <- c("h_chr7 m_chr5", "h_chr7 m_chr12", 
                   "h_chr8 m_chr15", "h_chr8 m_chr8", "h_chr8 m_chr3", 
                   "h_chr13 m_chr5", "h_chr13 m_chr8","h_chr13 m_chr3",
                   "h_chr20 m_chr2")
coad_adenoCar_arm_allSynTable2 <- coad_adenoCar_arm_allSynTable[which(coad_adenoCar_arm_allSynTable$tmpString %in% coadSynChroms), ]
coadZscoreCompGr <- GRanges(seqnames = str_remove(coad_adenoCar_arm_allSynTable2$h_chr, "h\\_chr"),
                            IRanges(start = coad_adenoCar_arm_allSynTable2$h_start,
                                    end = coad_adenoCar_arm_allSynTable2$h_end))

cancerGeneCensus[subjectHits(findOverlaps(coadZscoreCompGr,cancerGeneCensusGr)),]

### looking at same thing but with mouse as base
### need synteny plot of other way around
coad_adenoCar_arm_allSynTableMouseRef <- circosFreq2(allMouseAneu_coadAdenoCar_amp_bed, allMouseAneu_coadAdenoCar_del_bed, armGisticCoadApc_amp_bed,
                                                     armGisticCoadApc_del_bed, filename = "tmp", ref = "mouse", plot = FALSE)

coad_adenoCar_arm_allSynTableMouseRef$str <- paste(coad_adenoCar_arm_allSynTableMouseRef$h_chr, coad_adenoCar_arm_allSynTableMouseRef$h_start, coad_adenoCar_arm_allSynTableMouseRef$h_end,
                                                   coad_adenoCar_arm_allSynTableMouseRef$m_chr, coad_adenoCar_arm_allSynTableMouseRef$m_start, coad_adenoCar_arm_allSynTableMouseRef$m_end)
coad_adenoCar_arm_allSynTableMouseRef <- coad_adenoCar_arm_allSynTableMouseRef[-which(duplicated(coad_adenoCar_arm_allSynTableMouseRef$str)), ]

barplotFractionMouse <- NULL
i <- unique(coad_adenoCar_arm_allSynTableMouseRef$m_chr)[1]
for (i in unique(coad_adenoCar_arm_allSynTableMouseRef$m_chr)) {
  tmpDf <- coad_adenoCar_arm_allSynTableMouseRef[which(coad_adenoCar_arm_allSynTableMouseRef$m_chr == i), ]
  tmpDf$length <- tmpDf$m_end - tmpDf$m_start
  for (j in unique(tmpDf$h_chr)) {
    tmpRes <- sum(tmpDf$length[which(tmpDf$h_chr == j)])/sum(tmpDf$length)
    barplotFractionMouse <- rbind(barplotFractionMouse, data.frame("h_chr" = j, "m_chr" = i,"fraction" = tmpRes,
                                                                   "numberRegions" = length(which(tmpDf$h_chr == j))))
  }
}


# barplotFractionFiltMouse <- barplotFractionMouse[which(barplotFractionMouse$fraction > 0.01), ]
barplotFractionFiltMouse <- barplotFractionMouse
barplotFractionFiltMouse$h_chr <- factor(barplotFractionFiltMouse$h_chr, levels = paste0("h_chr", 1:22))
barplotFractionFiltMouse$m_chr <- factor(barplotFractionFiltMouse$m_chr, levels = paste0("m_chr", 1:19))

barplotFractionFiltMouseCoad <- barplotFractionFiltMouse
barplotFractionFiltMouseCoad$color <- "#FFFFFF"
barplotFractionFiltMouseCoad$color[which(barplotFractionFiltMouseCoad$m_chr == "m_chr5")] <- "#8B0000"
barplotFractionFiltMouseCoad$color[which(barplotFractionFiltMouseCoad$m_chr == "m_chr12")] <- "#8B0000"
barplotFractionFiltMouseCoad$color[which(barplotFractionFiltMouseCoad$m_chr == "m_chr18")] <- "#00008B"
barplotFractionFiltMouseCoad$color[which(barplotFractionFiltMouseCoad$m_chr == "m_chr14")] <- "#00008B"
barplotFractionFiltMouseCoad$color[which(barplotFractionFiltMouseCoad$m_chr == "m_chr9")] <- "#00008B"
barplotFractionFiltMouseCoad$color[which(barplotFractionFiltMouseCoad$m_chr == "m_chr7")] <- "#00008B"
barplotFractionFiltMouseCoad$color[which(barplotFractionFiltMouseCoad$m_chr == "m_chr4")] <- "#00008B"
barplotFractionFiltMouseCoad$alpha <- ifelse(barplotFractionFiltMouseCoad$color == "#FFFFFF", 0.1, 1)

pdf(file = "/mnt/DATA5/tmp/kev/misc/20231011fractionSyntenyConMouseRef.pdf", useDingbats = TRUE,
    width = 10, height = 5)
ggplot(barplotFractionFiltMouseCoad, aes(fill= h_chr, y= fraction, x = m_chr)) + 
  geom_bar(position="stack", stat="identity", color = barplotFractionFiltMouseCoad$color,
           alpha =  barplotFractionFiltMouseCoad$alpha) + 
  scale_fill_manual(values = colorVector) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + 
  ylab("human chromosome synteny fraction") + xlab("mouse chromsome")
dev.off()


barplotFractionFiltMouseCoad2 <- barplotFractionFiltMouseCoad
barplotFractionFiltMouseCoad2$str <- paste(barplotFractionFiltMouseCoad2$h_chr, barplotFractionFiltMouseCoad2$m_chr)
barplotFractionFiltMouseCoad2$alpha <- 0.1
barplotFractionFiltMouseCoad2$alpha[which(barplotFractionFiltMouseCoad2$str == "h_chr13 m_chr5" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr7 m_chr5" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr7 m_chr12" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr15 m_chr9" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr15 m_chr7" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr14 m_chr14" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr8 m_chr14" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr1 m_chr4" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr18 m_chr18")] <- 1
barplotFractionFiltMouseCoad2$color <- "#FFFFFF"
barplotFractionFiltMouseCoad2$color[which(barplotFractionFiltMouseCoad2$str == "h_chr13 m_chr5" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr7 m_chr5" | 
                                            barplotFractionFiltMouseCoad2$str == "h_chr7 m_chr12")] <- "#8B0000"

barplotFractionFiltMouseCoad2$color[which(barplotFractionFiltMouseCoad2$str == "h_chr18 m_chr18" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr15 m_chr9" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr15 m_chr7" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr14 m_chr14" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr8 m_chr14" |
                                            barplotFractionFiltMouseCoad2$str == "h_chr1 m_chr4")] <- "#00008B"



pdf(file = "/mnt/DATA5/tmp/kev/misc/20231011fractionSyntenyConMouseRefPerm.pdf", useDingbats = TRUE,
    width = 10, height = 5)
ggplot(barplotFractionFiltMouseCoad2, aes(fill= h_chr, y= fraction, x = m_chr)) + 
  geom_bar(position="stack", stat="identity", color = barplotFractionFiltMouseCoad2$color,
           alpha =  barplotFractionFiltMouseCoad2$alpha) + 
  scale_fill_manual(values = colorVector) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + 
  ylab("human chromosome synteny fraction") + xlab("mouse chromsome")
dev.off()
