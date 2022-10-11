
syntenyFreq <- function(segmentRange, syntenyRange, freqTable){
  tmpVar <- queryHits(findOverlaps(segmentRange, syntenyRange))
  if (isEmpty(tmpVar)) {
    tmpVar <- 0
  } else if(length(tmpVar) > 1){
    # need to make sure weighted length isn't inflated if start
    # start position is before or end position is after the given block
    tmpStart <-  freqTable$Start[tmpVar]
    tmpEnd <-  freqTable$End[tmpVar]
    tmpStart[tmpStart < synteny_hg38mm10$ref_start_pos[i] ] <- synteny_hg38mm10$ref_start_pos[i]
    tmpEnd[tmpEnd > synteny_hg38mm10$ref_end_pos[i]] <- synteny_hg38mm10$ref_end_pos[i]
    tmpVar <- weighted.mean(freqTable$Freq[tmpVar], tmpEnd - tmpStart)
  } else{
    tmpVar <- freqTable$Freq[queryHits(findOverlaps(segmentRange, syntenyRange))]
  }
  return(tmpVar)
}


synteny_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_syntenyDf.txt",sep = "\t", stringsAsFactors = FALSE, header = TRUE)

mouseBedFile <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.gc.bed", stringsAsFactors = FALSE, sep = "\t",
                           header = FALSE)

probes_1 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/cnAmplicon_matrix.txt",
                                       header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
probes_2 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/cnAmplicon_matrix.txt",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
probes_3 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/cnAmplicon_matrix.txt",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
probes_4 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-76-MG_test1_255_185/cnAmplicon_matrix.txt",
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)

# (1) easiest method is to treat each block like a gene - biggest assumption is these syntenic blocks
# uniformly change i.e output may be inaccurate if part of block is gained, while other is lost
# i.e mean of 0 can just be high gains and losses

# (2) alternatively, I can trust segmentation and look at overall gains and losses of syntenic blocks
# done by taking segmentation calls and overlapping them with syntenic blocks 

# have to go with alternate method b/c I would need the probe probe files for the respective cancers
# plus with the probe files ploidy isn't accounted for


# mouseAmplicons <- probes_3
# extraCols <- c("AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene", "NumGC",
#                "Length", "GC", "TotalPool", "Weights")
# mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
# mouseBed2 <- data.frame("AmpliconId" = mouseBedIdx, "ChromNum" = mouseBedFile$V1[mouseBedIdx],
#                         "Start" = mouseBedFile$V2[mouseBedIdx], "End" = mouseBedFile$V3[mouseBedIdx])
# mouseBed2$ChromNum <- as.numeric(str_remove(str_replace(mouseBed2$ChromNum, "chrX", "chr23"), "chr"))
# mouseAmplicons2 <- mouseAmplicons[,-which(colnames(mouseAmplicons) %in% extraCols)]
# mouseAmplicons2[,2:(ncol(mouseAmplicons2))] <- log2(mouseAmplicons2[,2:(ncol(mouseAmplicons2))])
# mouseBed2$ChromNum <- str_replace(mouseBed2$ChromNum, "23", "20")


# i need to calculate these by sample - create a sample by syntenic block matrix
# resMatrix <- matrix(0, ncol = nrow(synteny_hg38mm10), nrow = ncol(mouseAmplicons2) - 1)
# colnames(resMatrix) <- synteny_hg38mm10$symbol
synteny_hg38mm10$comp_chr <- str_replace(synteny_hg38mm10$comp_chr, "X", "20")
syntenicGrange <- GRanges(seqnames = synteny_hg38mm10$comp_chr, 
                          IRanges(start = synteny_hg38mm10$comp_start_pos, end = synteny_hg38mm10$comp_end_pos))

hgsc_arm_amp_bed_filt10 <- hgsc_arm_amp_bed[-which(hgsc_arm_amp_bed$Freq < 10),]
hgsc_arm_del_bed_filt10 <- hgsc_arm_del_bed[-which(hgsc_arm_del_bed$Freq < 10),]
bprn_arm_amp_bed_filt10 <- bprn_arm_amp_bed[-which(bprn_arm_amp_bed$Freq < 10),]
bprn_arm_del_bed_filt10 <- bprn_arm_del_bed[-which(bprn_arm_del_bed$Freq < 10),]

hgsc_cna_amp_bed_filt10 <- hgsc_cna_amp_bed[-which(hgsc_cna_amp_bed$Freq < 10),]
hgsc_cna_del_bed_filt10 <- hgsc_cna_del_bed[-which(hgsc_cna_del_bed$Freq < 10),]
bprn_cna_amp_bed_filt10 <- bprn_cna_amp_bed[-which(bprn_cna_amp_bed$Freq < 10),]
bprn_cna_del_bed_filt10 <- bprn_cna_del_bed[-which(bprn_cna_del_bed$Freq < 10),]

segGrangeArmGain_h <- GRanges(seqnames = hgsc_arm_amp_bed_filt10$Chr,
                             IRanges(start = hgsc_arm_amp_bed_filt10$Start,
                                     end = hgsc_arm_amp_bed_filt10$End))
segGrangeArmLoss_h <- GRanges(seqnames = hgsc_arm_del_bed_filt10$Chr,
                             IRanges(start = hgsc_arm_del_bed_filt10$Start,
                                     end = hgsc_arm_del_bed_filt10$End))
segGrangeArmGain_m <- GRanges(seqnames = bprn_arm_amp_bed_filt10$Chr,
                             IRanges(start = bprn_arm_amp_bed_filt10$Start,
                                     end = bprn_arm_amp_bed_filt10$End))
segGrangeArmDel_m <- GRanges(seqnames = bprn_arm_del_bed_filt10$Chr,
                              IRanges(start = bprn_arm_del_bed_filt10$Start,
                                      end = bprn_arm_del_bed_filt10$End))

segGrangeCnaGain_h <- GRanges(seqnames = hgsc_cna_amp_bed_filt10$Chr,
                              IRanges(start = hgsc_cna_amp_bed_filt10$Start,
                                      end = hgsc_cna_amp_bed_filt10$End))
segGrangeCnaLoss_h <- GRanges(seqnames = hgsc_cna_del_bed_filt10$Chr,
                              IRanges(start = hgsc_cna_del_bed_filt10$Start,
                                      end = hgsc_cna_del_bed_filt10$End))
segGrangeCnaGain_m <- GRanges(seqnames = bprn_cna_amp_bed_filt10$Chr,
                              IRanges(start = bprn_cna_amp_bed_filt10$Start,
                                      end = bprn_cna_amp_bed_filt10$End))
segGrangeCnaDel_m <- GRanges(seqnames = bprn_cna_del_bed_filt10$Chr,
                             IRanges(start = bprn_cna_del_bed_filt10$Start,
                                     end = bprn_cna_del_bed_filt10$End))


ptm <- proc.time()
resMatrix <- matrix(0, ncol = nrow(synteny_hg38mm10), 4)
colnames(resMatrix) <- synteny_hg38mm10$symbol
for (i in 1:nrow(synteny_hg38mm10)) {
  tmpGrange_h <- GRanges(seqnames = synteny_hg38mm10$ref_chr[i],
                       IRanges(start = synteny_hg38mm10$ref_start_pos[i],
                               end = synteny_hg38mm10$ref_end_pos[i]))
  
  tmpGrange_m <- GRanges(seqnames = synteny_hg38mm10$comp_chr[i],
                         IRanges(start = synteny_hg38mm10$comp_start_pos[i],
                                 end = synteny_hg38mm10$comp_end_pos[i]))
  
  hgGainVar <- syntenyFreq(segGrangeArmGain_h, tmpGrange_h, hgsc_arm_amp_bed_filt10)
  hgLossVar <- syntenyFreq(segGrangeArmLoss_h, tmpGrange_h, hgsc_arm_del_bed_filt10)
  mmGainVar <- syntenyFreq(segGrangeArmGain_m, tmpGrange_m, bprn_arm_amp_bed_filt10)
  mmLossVar <- syntenyFreq(segGrangeArmDel_m, tmpGrange_m, bprn_arm_del_bed_filt10)
  
  resMatrix[1,i] <- hgGainVar
  resMatrix[2,i] <- hgLossVar
  resMatrix[3,i] <- mmGainVar
  resMatrix[4,i] <- mmLossVar
  
}

proc.time() - ptm

rownames(resMatrix) <- c("HgArmGain", "HgArmLoss", "MmArmGain", "MmArmLoss")

### for the graphs - I can separate them into changes 
resMatrix_filt <- resMatrix[,-which(apply(resMatrix,2, sum) == 0)]
resMatrix_melt <- reshape2::melt(resMatrix_filt)
resMatrix_melt$value[grep("Loss",resMatrix_melt$Var1)] <- resMatrix_melt$value[grep("Loss",resMatrix_melt$Var1)] * -1

ggplot(resMatrix_melt, aes(fill=Var1, y=value, x=Var2)) + 
  geom_bar(position="stack", stat="identity")


### create pie chart stats i.e how man have no changes
armSynCountTable <- resMatrix_melt[-which(resMatrix_melt$value == 0),]
singularChangesNames <- names(table(armSynCountTable$Var2))[which(table(armSynCountTable$Var2) == 1)]
armSynSingular <- armSynCountTable[which(armSynCountTable$Var2 %in% singularChangesNames),]
armSynMultiple <- armSynCountTable[-which(armSynCountTable$Var2 %in% singularChangesNames),]

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210823hgscBprnSynBlock.pdf", width = 14,useDingbats = FALSE)
ggplot(armSynMultiple, aes(fill=Var1, y=value, x=Var2)) + 
  geom_bar(position="stack", stat="identity") + theme_bw() + 
  theme(axis.text.x=element_text(angle=90,hjust=1)) + 
  ylab("Frequency of Alteration * gain(+)/loss(-)") +
  xlab("Syntenic Block") + guides(fill=guide_legend(title="Type of Aneuploidy & Species"))
dev.off()

# stats total 913 syntenic blocks mean of ~2.8Mb and median ~0.25Mb; 
# 544 have some change; 84 have changes in both species - filtered by "frequent changes" >= 10%
# should scale the dots by number of genes -  create a graphic 
dim(resMatrix)
length(unique(as.character(resMatrix_melt$Var2)))
summary((synteny_hg38mm10$ref_end_pos- synteny_hg38mm10$ref_start_pos)/1e6)
summary((synteny_hg38mm10$comp_end_pos- synteny_hg38mm10$comp_start_pos)/1e6)


### need some statistic for these concordant changes happnning more than chance

# (1) permute 
# (2) count concordant events
# (3) is my original number of concordant events more than by chance? not specific to syntenic block
# (4) for a specific event could I just 

resMatrix_bin <- resMatrix
resMatrix_bin[resMatrix_bin > 0] <- 1
gains <- apply(resMatrix_bin[c(1,3),], 2, sum)
losses <- apply(resMatrix_bin[c(2,4),], 2, sum)
mix1 <- apply(resMatrix_bin[c(1,4),], 2, sum)
mix2 <- apply(resMatrix_bin[c(2,3),], 2, sum)
table(gains)[3]
table(losses)[3]
table(mix1)[3] + table(mix2)[3]

### (1) find genes within concordant changes - wonder if they overlap with known pathways?
### (2) are the size distributions the same for concordant changes vs discordant?
synteny_hg38mm10$mLength <- synteny_hg38mm10$comp_end_pos - synteny_hg38mm10$comp_start_pos
synteny_hg38mm10$hLength <- synteny_hg38mm10$ref_end_pos - synteny_hg38mm10$ref_start_pos

# looked at correlation to make sure general length of these blocks make sense
# genic content though? should be the same as nature of these blocks are built

cor(synteny_hg38mm10$mLength, synteny_hg38mm10$hLength)

gainNames <- names(which(gains == 2))
lossNames <- names(which(losses == 2))
disNames <- c(names(which(mix1 == 2)), names(which(mix2 == 2)))

gainLength <- (synteny_hg38mm10$hLength/1e6)[which(synteny_hg38mm10$symbol %in% gainNames)]
lossLength <- (synteny_hg38mm10$hLength/1e6)[which(synteny_hg38mm10$symbol %in% lossNames)]
disLength <- (synteny_hg38mm10$hLength/1e6)[which(synteny_hg38mm10$symbol %in% disNames)]
armLengthDf <- data.frame("Type" = c(rep("armGain (39)", length(gainLength)), rep("armLoss (5)", length(lossNames)),
                                     rep("armDisc (42)", length(disLength)), rep("all", nrow(synteny_hg38mm10))),
                          "Length" = c(gainLength, lossLength, disLength, synteny_hg38mm10$hLength/1e6), stringsAsFactors = FALSE)

ggplot(armLengthDf, aes(x = Type, y = Length)) + geom_boxplot()

### calculate the gains, losses and sum of occurrences of mix1 + mix2
### then take difference for what we observed above
### to make it reproducible, need to randomize which seed is set in each loop and append it to list of seeds
### uhhh maybe K-S test with permutation?

gainsList <- NULL
lossesList <- NULL
discordantList <- NULL
i <- 0
while (i < 10000) {
  resampSynRes <- rbind(sample(resMatrix_bin[1,], ncol(resMatrix_bin)),
                        sample(resMatrix_bin[2,], ncol(resMatrix_bin)),
                        sample(resMatrix_bin[3,], ncol(resMatrix_bin)),
                        sample(resMatrix_bin[4,], ncol(resMatrix_bin)))
  gains <- apply(resampSynRes[c(1,3),], 2, sum)
  losses <- apply(resampSynRes[c(2,4),], 2, sum)
  mix1 <- apply(resampSynRes[c(1,4),], 2, sum)
  mix2 <- apply(resampSynRes[c(2,3),], 2, sum)
  
  gainsList <- c(gainsList, table(gains)[3])
  lossesList <- c(lossesList, table(losses)[3])
  discordantList <- c(discordantList, table(mix1)[3] + table(mix2)[3])
  i <- i + 1
}

length(which(gainsList > 39))/1e4 * 2
length(which(lossesList < 5))/1e4 * 2
length(which(discordantList < 42))/1e4 * 2

a <- ggplot(data = as.data.frame(gainsList), aes(x = gainsList)) + 
  geom_histogram(binwidth = 1, color="black", fill="grey") + theme_bw() +
  geom_vline(xintercept = 39, color = "blue") + xlab("Number of gain") + 
  ggtitle("Number of concordant gains - perm (n=10,000); 39")

b <- ggplot(data = as.data.frame(lossesList), aes(x = lossesList)) + 
  geom_histogram(binwidth = 1, color="black", fill="grey") + theme_bw() +
  geom_vline(xintercept = 5, color = "blue") +
  ggtitle("Number of concordant losses - perm (n=10,000); 5")

c <- ggplot(data = as.data.frame(discordantList), aes(x = discordantList)) + 
  geom_histogram(binwidth = 1, color="black", fill="grey") + theme_bw() +
  geom_vline(xintercept = 42, color = "blue") + 
  ggtitle("Number of discordant changes - perm (n=10,000); 42")


pdf("/mnt/DATA5/tmp/kev/misc/20210828permTestArmHgscSynBlocks.pdf", 
    useDingbats = FALSE, height = 7, width = 15)
grid.arrange(a,b,c, ncol = 3)
dev.off()




### CNA


resMatrixCna <- matrix(0, ncol = nrow(synteny_hg38mm10), 4)
colnames(resMatrixCna) <- synteny_hg38mm10$symbol
for (i in 1:nrow(synteny_hg38mm10)) {
  tmpGrange_h <- GRanges(seqnames = synteny_hg38mm10$ref_chr[i],
                         IRanges(start = synteny_hg38mm10$ref_start_pos[i],
                                 end = synteny_hg38mm10$ref_end_pos[i]))
  
  tmpGrange_m <- GRanges(seqnames = synteny_hg38mm10$comp_chr[i],
                         IRanges(start = synteny_hg38mm10$comp_start_pos[i],
                                 end = synteny_hg38mm10$comp_end_pos[i]))
  
  hgGainVar <- syntenyFreq(segGrangeCnaGain_h, tmpGrange_h, hgsc_cna_amp_bed_filt10)
  hgLossVar <- syntenyFreq(segGrangeCnaLoss_h, tmpGrange_h, hgsc_cna_del_bed_filt10)
  mmGainVar <- syntenyFreq(segGrangeCnaGain_m, tmpGrange_m, bprn_cna_amp_bed_filt10)
  mmLossVar <- syntenyFreq(segGrangeCnaDel_m, tmpGrange_m, bprn_cna_del_bed_filt10)
  
  resMatrixCna[1,i] <- hgGainVar
  resMatrixCna[2,i] <- hgLossVar
  resMatrixCna[3,i] <- mmGainVar
  resMatrixCna[4,i] <- mmLossVar
  
}

rownames(resMatrixCna) <- c("HgCnaGain", "HgCnaLoss", "MmCnaGain", "MmCnaLoss")
resMatrixCna_filt <- resMatrixCna[,-which(apply(resMatrixCna,2, sum) == 0)]
resMatrixCna_melt <- reshape2::melt(resMatrixCna_filt)
resMatrixCna_melt$value[grep("Loss",resMatrixCna_melt$Var1)] <- resMatrixCna_melt$value[grep("Loss",resMatrixCna_melt$Var1)] * -1

cnaSynCountTable <- resMatrixCna_melt[-which(resMatrixCna_melt$value == 0),]
singularChangesNamesCna <- names(table(cnaSynCountTable$Var2))[which(table(cnaSynCountTable$Var2) == 1)]
cnaSynSingular <- cnaSynCountTable[which(cnaSynCountTable$Var2 %in% singularChangesNamesCna),]
cnaSynMultiple <- cnaSynCountTable[-which(cnaSynCountTable$Var2 %in% singularChangesNamesCna),]


ggplot(cnaSynMultiple, aes(fill=Var1, y=value, x=Var2)) + 
  geom_bar(position="stack", stat="identity") + theme_bw() + 
  theme(axis.text.x=element_text(angle=90,hjust=1)) + 
  ylab("Frequency of Alteration * gain(+)/loss(-)") +
  xlab("Syntenic Block") + guides(fill=guide_legend(title="Type of Aneuploidy & Species"))


resMatrixCna_bin <- resMatrixCna
resMatrixCna_bin[resMatrixCna_bin > 0] <- 1
gainsCna <- apply(resMatrixCna_bin[c(1,3),], 2, sum)
lossesCna <- apply(resMatrixCna_bin[c(2,4),], 2, sum)
mix1Cna <- apply(resMatrixCna_bin[c(1,4),], 2, sum)
mix2Cna <- apply(resMatrixCna_bin[c(2,3),], 2, sum)
table(gainsCna)[3]
table(lossesCna)[3]
table(mix1Cna)[3] + table(mix2Cna)[3]

gainNamesCna <- names(which(gainsCna == 2))
lossNamesCna <- names(which(lossesCna == 2))
disNamesCna <- c(names(which(mix1Cna == 2)), names(which(mix2Cna == 2)))

gainLengthCna <- (synteny_hg38mm10$hLength/1e6)[which(synteny_hg38mm10$symbol %in% gainNamesCna)]
lossLengthCna <- (synteny_hg38mm10$hLength/1e6)[which(synteny_hg38mm10$symbol %in% lossNamesCna)]
disLengthCna <- (synteny_hg38mm10$hLength/1e6)[which(synteny_hg38mm10$symbol %in% disNamesCna)]
cnaLengthDf <- data.frame("Type" = c(rep("cnaGain (78)", length(gainLengthCna)), rep("cnaLoss (36)", length(lossNamesCna)),
                                     rep("cnaDisc (98)", length(disLengthCna)), rep("all", nrow(synteny_hg38mm10))),
                          "Length" = c(gainLengthCna, lossLengthCna, disLengthCna, synteny_hg38mm10$hLength/1e6), stringsAsFactors = FALSE)

ggplot(cnaLengthDf, aes(x = Type, y = Length)) + geom_boxplot()

armCnaLength <- rbind(armLengthDf[-which(armLengthDf$Type == "all"),],
                      cnaLengthDf)

pdf("/mnt/DATA5/tmp/kev/misc/20210828lengthDisMatchingSyn.pdf", useDingbats = FALSE)
ggplot(armCnaLength, aes(x = Type, y = Length)) + geom_boxplot() + xlab("Type of matches between species") +
  theme_bw() + ggtitle("Length distribution in matching changes (syn. blocks)")
dev.off()

lengthTTest <- pairwise.t.test(armCnaLength$Length, armCnaLength$Type, p.adjust.method = "BH")
lengthTTest <- data.frame(lengthTTest$p.value, stringsAsFactors = FALSE)
lengthTTest <- data.frame(lapply(lengthTTest, function(x) signif(x, digits = 2)))
colnames(lengthTTest) <- c("all", "armDisc(42)", "armGain(39)", "armLoss(5)", 
                           "cnaDisc (98)", "cnaGain(78)")
rownames(lengthTTest) <- c("armDisc(42)", "armGain(39)", "armLoss(5)", 
                           "cnaDisc (98)", "cnaGain(78)", "cnaLoss(36)")

pdf("/mnt/DATA5/tmp/kev/misc/20210828pairwiseTTestLength.pdf", width = 12, useDingbats = FALSE)
grid.table(lengthTTest)
dev.off()
### perm test CNA


gainsListCna <- NULL
lossesListCna <- NULL
discordantListCna <- NULL
i <- 0
while (i < 10000) {
  resampSynRes <- rbind(sample(resMatrixCna_bin[1,], ncol(resMatrixCna_bin)),
                        sample(resMatrixCna_bin[2,], ncol(resMatrixCna_bin)),
                        sample(resMatrixCna_bin[3,], ncol(resMatrixCna_bin)),
                        sample(resMatrixCna_bin[4,], ncol(resMatrixCna_bin)))
  gains <- apply(resampSynRes[c(1,3),], 2, sum)
  losses <- apply(resampSynRes[c(2,4),], 2, sum)
  mix1 <- apply(resampSynRes[c(1,4),], 2, sum)
  mix2 <- apply(resampSynRes[c(2,3),], 2, sum)
  
  gainsListCna <- c(gainsListCna, table(gains)[3])
  lossesListCna <- c(lossesListCna, table(losses)[3])
  discordantListCna <- c(discordantListCna, table(mix1)[3] + table(mix2)[3])
  i <- i + 1
}

length(which(gainsListCna > 78))/1e4 * 2
length(which(lossesListCna > 36))/1e4 * 2
length(which(discordantListCna < 99))/1e4 * 2

d <- ggplot(data = as.data.frame(gainsListCna), aes(x = gainsListCna)) + 
  geom_histogram(binwidth = 1, color="black", fill="grey") + theme_bw() +
  geom_vline(xintercept = 78, color = "blue") + xlab("Number of gain") + 
  ggtitle("Number of concordant gains - perm (n=10,000); 53")

e <- ggplot(data = as.data.frame(lossesListCna), aes(x = lossesListCna)) + 
  geom_histogram(binwidth = 1, color="black", fill="grey") + theme_bw() +
  geom_vline(xintercept = 36, color = "blue") +
  ggtitle("Number of concordant losses - perm (n=10,000); 26")

f <- ggplot(data = as.data.frame(discordantListCna), aes(x = discordantListCna)) + 
  geom_histogram(binwidth = 1, color="black", fill="grey") + theme_bw() +
  geom_vline(xintercept = 99, color = "blue") + 
  ggtitle("Number of discordant changes - perm (n=10,000); 79")

grid.arrange(d, e, f, ncol = 3)

pdf("/mnt/DATA5/tmp/kev/misc/20210828permTestCna5MbHgscSynBlocks.pdf", 
    useDingbats = FALSE, height = 7, width = 15)
grid.arrange(d, e, f, ncol = 3)
dev.off()


### see if any similar blocks between cna and arms,
### then genes from these -> pathways for importance (?)

intersect(gainNames, gainNamesCna)

# ptch1 shared loss and fancc
intersect(lossNames, lossNamesCna)
intersect(gainNames, gainNamesCna)

# the types I need to look at are arm gains, and all losses
gainNamesBlocks <- synteny_hg38mm10[which(synteny_hg38mm10$symbol %in% gainNames),]
lossNamesBlocks <- synteny_hg38mm10[which(synteny_hg38mm10$symbol %in% lossNames),]
lossNamesCnaBlocks <- synteny_hg38mm10[which(synteny_hg38mm10$symbol %in% lossNamesCna),]

armGainGlDf <- getGeneList_synblock(gainNamesBlocks)
armGainGl <- unlist(str_split(paste(armGainGlDf$h_gene, collapse = ","), ","))

armLossGlDf <- getGeneList_synblock(lossNamesBlocks)
armLossGl <- unlist(str_split(paste(armLossGlDf$h_gene, collapse = ","), ","))

cnaLossGlDf <- getGeneList_synblock(lossNamesCnaBlocks)
cnaLossGl <- unlist(str_split(paste(cnaLossGlDf$h_gene, collapse = ","), ","))


tmpHallList <- NULL
for (i in 1:length(human_hallmarks_list)) {
  tmpHallList <- c(tmpHallList, unlist(human_hallmarks_list[[i]]))
}

tmpHallList <- tmpHallList[-which(tmpHallList == "")]
tmpHallList <- unique(tmpHallList)

enrichmentStats_synblock(armGainGl, human_hallmarks_list, tmpHallList)
enrichmentStats_synblock(armLossGl, human_hallmarks_list, tmpHallList)
cnaLoss_pathwayDf <- data.frame(enrichmentStats_synblock(cnaLossGl, human_hallmarks_list, tmpHallList), stringsAsFactors = FALSE)
cnaLoss_pathwayDf$qval <- p.adjust(cnaLoss_pathwayDf$hyper.pval, method = "BH")




