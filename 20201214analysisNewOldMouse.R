library(ggplot2)
library(ggrepel)

newPanelCalls <- read.table("/mnt/DATA5/tmp/kev/newMouse/combinedCalls.txt", sep = "\t", header = TRUE,
                            stringsAsFactors = FALSE)

newPanelCalls$Gene <- tolower(newPanelCalls$Gene)

oldPanelCalls <- read.table("/mnt/DATA5/tmp/kev/testMouse/combinedCalls.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
oldPanelCalls$Gene <- tolower(oldPanelCalls$Gene)


interGene <- intersect(newPanelCalls$Gene, oldPanelCalls$Gene)
newPanelCalls <- newPanelCalls[which(newPanelCalls$Gene %in% interGene),]
oldPanelCalls <- oldPanelCalls[which(oldPanelCalls$Gene %in% interGene),]

crossTableIdx <- readxl::read_xlsx("/home/kevhu/data/20201207annotations.xlsx")
crossTableIdx <- crossTableIdx[,1:3]
crossTableIdx$old_name2 <- tolower(crossTableIdx$old_name)
crossTableIdx$old_name2 <- str_remove(crossTableIdx$old_name2 , "[[:punct:]]")



oldPanelCalls_names <- str_remove(tolower(oldPanelCalls$Sample), "[[:punct:]]")
newPanelCalls_names <- as.numeric(str_remove(str_remove(newPanelCalls$Sample, "MG_"), "X.*"))
newPanelCalls_names2 <- newPanelCalls_names

for (i in unique(newPanelCalls_names)) {
  newPanelCalls_names2[which(newPanelCalls_names %in% i)] <- crossTableIdx$old_name2[which(crossTableIdx$mg_id %in% i)]
}


oldPanelCalls$Sample <- oldPanelCalls_names
newPanelCalls$Sample <- newPanelCalls_names2

oldPanelCalls_subset <- oldPanelCalls[which(oldPanelCalls$Sample %in% crossTableIdx$old_name2), ]
newPanelCalls_subset <- newPanelCalls[which(newPanelCalls$Sample %in% oldPanelCalls$Sample), ]


oldPanelCalls_subset <- oldPanelCalls_subset[order(oldPanelCalls_subset$Sample, oldPanelCalls_subset$Gene), ]
newPanelCalls_subset <- newPanelCalls_subset[order(newPanelCalls_subset$Sample,newPanelCalls_subset$Gene), ]
genesToRemove <- c("trp53", "nf1", "rb1", "brca1", "pten")

oldPanelCalls_subset_2 <- oldPanelCalls_subset[-which(oldPanelCalls_subset$Gene %in% genesToRemove),] 
newPanelCalls_subset_2 <- newPanelCalls_subset[-which(newPanelCalls_subset$Gene %in% genesToRemove),]


oldPanel_subset2_calls <- intersect(which(oldPanelCalls_subset_2$Log10QValue < -1.30103), which(abs(log2(oldPanelCalls_subset_2$CopyNumberRatio)) > 0.2))

cor(log2(oldPanelCalls_subset_2$CopyNumberRatio[oldPanel_subset2_calls]), log2(newPanelCalls_subset_2$CopyNumberRatio[oldPanel_subset2_calls]))
plot(log2(oldPanelCalls_subset_2$CopyNumberRatio[oldPanel_subset2_calls]), log2(newPanelCalls_subset_2$CopyNumberRatio[oldPanel_subset2_calls]))


newColors <- rep("black", length(oldPanel_subset2_calls))
tmpDf <- newPanelCalls_subset_2[oldPanel_subset2_calls,]
for (i in 1:length(newColors)) {
  if(tmpDf$Log10QValue[i] < -1.3){
    newColors[i] <- "good"
  } else{
    newColors[i] <- "bad"
  }
}


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20201214somaticCN_point.pdf",useDingbats = FALSE)
ggplot(data = data.frame("oldPanelCnr" = log2(oldPanelCalls_subset_2$CopyNumberRatio[oldPanel_subset2_calls]),
                         "newPanelCnr" = log2(newPanelCalls_subset_2$CopyNumberRatio[oldPanel_subset2_calls]),
                         "label" = oldPanelCalls_subset_2$Gene[oldPanel_subset2_calls],
                         "match" = newColors),
       aes(oldPanelCnr, newPanelCnr, label = label)) +
  geom_point(aes(x = oldPanelCnr, y = newPanelCnr, color = match)) +
  scale_colour_manual(values = c("good" = "darkgreen", "bad" = "darkred")) + 
  geom_text_repel() + ylab("New Panel log2(Cnr)") + xlab("Old Panel log2(Cnr)") + 
  ggtitle("Somatic gains/losses: 54 alterations, 12 genes, 12 samples; R = 0.91") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


### same as above but for Del genes
newPanelDel <- read.table("/mnt/DATA5/tmp/kev/newMouse2/combinedCalls.txt", sep = "\t", header = TRUE,
                            stringsAsFactors = FALSE)
oldPanelDel <- read.table("/mnt/DATA5/tmp/kev/testMouse2/combinedCalls.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

newPanelDel_subset <- newPanelDel[which(grepl("Del", newPanelDel$Gene)),]
oldPanelDel_subset <- oldPanelDel[which(grepl("Del", oldPanelDel$Gene)),]

oldPanelDel_names <- str_remove(tolower(oldPanelDel_subset$Sample), "[[:punct:]]")
newPanelDel_names <- as.numeric(str_remove(str_remove(newPanelDel_subset$Sample, "MG_"), "X.*"))
newPanelDel_names2 <- newPanelDel_names

for (i in unique(newPanelDel_names)) {
  newPanelDel_names2[which(newPanelDel_names %in% i)] <- crossTableIdx$old_name2[which(crossTableIdx$mg_id %in% i)]
}

oldPanelDel_subset$Sample <- oldPanelDel_names
newPanelDel_subset$Sample <- newPanelDel_names2

oldPanelDel_subset2 <- oldPanelDel_subset[which(oldPanelDel_subset$Sample %in% crossTableIdx$old_name2), ]
newPanelDel_subset2 <- newPanelDel_subset[which(newPanelDel_subset$Sample %in% oldPanelDel_subset$Sample), ]

oldPanelDel_subset2 <- oldPanelDel_subset2[order(oldPanelDel_subset2$Sample, oldPanelDel_subset2$Gene), ]
newPanelDel_subset2 <- newPanelDel_subset2[order(newPanelDel_subset2$Sample,newPanelDel_subset2$Gene), ]

oldPanel_subset2_Del <- intersect(which(oldPanelDel_subset2$Log10QValue < -1.30103), which(abs(log2(oldPanelDel_subset2$CopyNumberRatio)) > 0.2))

cor(log2(oldPanelDel_subset2$CopyNumberRatio[oldPanel_subset2_Del]), log2(newPanelDel_subset2$CopyNumberRatio[oldPanel_subset2_Del]))




newColors2 <- rep("black", length(oldPanel_subset2_Del))
tmpDf2 <- newPanelDel_subset2[oldPanel_subset2_Del,]
for (i in 1:length(newColors2)) {
  if(tmpDf2$Log10QValue[i] < -1.3){
    newColors2[i] <- "good"
  } else{
    newColors2[i] <- "bad"
  }
}


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20201214creloxCN_point.pdf",useDingbats = FALSE)
ggplot(data = data.frame("oldPanelCnr" = log2(oldPanelDel_subset2$CopyNumberRatio[oldPanel_subset2_Del]),
                         "newPanelCnr" = log2(newPanelDel_subset2$CopyNumberRatio[oldPanel_subset2_Del]),
                         "label" = oldPanelDel_subset2$Gene[oldPanel_subset2_Del],
                         "match" = newColors2),
       aes(oldPanelCnr, newPanelCnr, label = label)) +
  geom_point(aes(x = oldPanelCnr, y = newPanelCnr, color = match)) +
  scale_colour_manual(values = c("good" = "darkgreen", "bad" = "darkred")) + 
  geom_text_repel() + ylab("New Panel log2(Cnr)") + xlab("Old Panel log2(Cnr)") + 
  ggtitle("Cre-lox sub-gene losses: 39 alterations, 5 genes, 16 samples; R = 0.97") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


### looking at new info called from same samples sequenced

newPanelCalls <- read.table("/mnt/DATA5/tmp/kev/newMouse/combinedCalls.txt", sep = "\t", header = TRUE,
                            stringsAsFactors = FALSE)

newPanelCalls$Gene <- tolower(newPanelCalls$Gene)

oldPanelCalls <- read.table("/mnt/DATA5/tmp/kev/testMouse/combinedCalls.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
oldPanelCalls$Gene <- tolower(oldPanelCalls$Gene)

crossTableIdx <- readxl::read_xlsx("/home/kevhu/data/20201207annotations.xlsx")
crossTableIdx <- crossTableIdx[,1:3]

crossTableIdx$old_name2 <- tolower(crossTableIdx$old_name)
crossTableIdx$old_name2 <- str_remove(crossTableIdx$old_name2 , "[[:punct:]]")


oldPanelCalls_names <- str_remove(tolower(oldPanelCalls$Sample), "[[:punct:]]")
newPanelCalls_names <- as.numeric(str_remove(str_remove(newPanelCalls$Sample, "MG_"), "X.*"))
newPanelCalls_names2 <- newPanelCalls_names

for (i in unique(newPanelCalls_names)) {
  newPanelCalls_names2[which(newPanelCalls_names %in% i)] <- crossTableIdx$old_name2[which(crossTableIdx$mg_id %in% i)]
}


oldPanelCalls$Sample <- oldPanelCalls_names
newPanelCalls$Sample <- newPanelCalls_names2

oldPanelCalls_subset <- oldPanelCalls[which(oldPanelCalls$Sample %in% crossTableIdx$old_name2), ]
newPanelCalls_subset <- newPanelCalls[which(newPanelCalls$Sample %in% oldPanelCalls$Sample), ]

oldPanelCalls_subset <- oldPanelCalls_subset[order(oldPanelCalls_subset$Sample, oldPanelCalls_subset$Gene), ]
newPanelCalls_subset <- newPanelCalls_subset[order(newPanelCalls_subset$Sample,newPanelCalls_subset$Gene), ]

normalSamps <- c("13604n", "14104t", "14154n", "14433n", "2405n", "2519n", "2796n", "3867n","8234n", "2611n")

oldPanelCalls_subset <- oldPanelCalls_subset[-which(oldPanelCalls_subset$Sample %in% normalSamps),]
newPanelCalls_subset <- newPanelCalls_subset[-which(newPanelCalls_subset$Sample %in% normalSamps),]

oldPanelCalls_subset$Color <- "old"
newPanelCalls_subset$Color <- "new"

dfForBoxplot <- NULL
for (i in unique(oldPanelCalls_subset$Sample)) {
  print(i)
  tmpDfOld <- oldPanelCalls_subset[which(oldPanelCalls_subset$Sample == i),]
  tmpDfNew <- newPanelCalls_subset[which(newPanelCalls_subset$Sample == i),]
  tmpDfOld$Sample <- paste0(tmpDfOld$Sample,"_v1")
  tmpDfNew$Sample <- paste0(tmpDfNew$Sample,"_v2")
  
  tmpDfOld <- tmpDfOld[intersect(which(tmpDfOld$Log10QValue < -1.301303),
                                 which(abs(log2(tmpDfOld$CopyNumberRatio)) > 0.2)),]
  tmpDfNew <- tmpDfNew[intersect(which(tmpDfNew$Log10QValue < -1.301303),
                                 which(abs(log2(tmpDfNew$CopyNumberRatio)) > 0.2)),]
  
  if (nrow(tmpDfOld) == 0) {
    tmpDfOld2 <- c(paste0(i, "_v1"), 1, "old")
  }
  else{
    tmpDfOld2 <- tmpDfOld[,c("Sample", "CopyNumberRatio", "Color")]
  }
  
  if (nrow(tmpDfNew) == 0) {
    tmpDfNew2 <- c(paste0(i, "_v2"), 1, "new")
  }
  else{
    tmpDfNew2 <- tmpDfNew[,c("Sample", "CopyNumberRatio", "Color")]
  }
  dfForBoxplot <- rbind(dfForBoxplot, tmpDfOld2, tmpDfNew2)
}

dfForBoxplot$CopyNumberRatio <- log2(as.numeric(dfForBoxplot$CopyNumberRatio))


ggplot(data = dfForBoxplot, aes(Sample, CopyNumberRatio, Color)) + geom_boxplot(aes(fill = Color)) + geom_jitter(size = 1.0) + 
  scale_fill_brewer(palette="Spectral") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### the graphs above don't give me too much to go off ... in general 
library(DNAcopy)

bedFile <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.gc.bed",
                      sep = "\t", header = FALSE)
new_amplicon <- read.table("/mnt/DATA5/tmp/kev/newMouse/cnAmplicon_matrix.txt", sep = "\t", stringsAsFactors = FALSE,
                           header = TRUE)
bedIdx <- as.numeric(str_remove(new_amplicon$AmpliconId, "AMP_"))
new_amplicon2 <- cbind(bedFile$V1[bedIdx], bedFile$V2[bedIdx], new_amplicon)
new_amplicon2[,6:42] <- log2(new_amplicon2[,6:42])

CNA.object <- CNA(cbind(new_amplicon2$MG_5X37),
                  new_amplicon2$ChromNum, new_amplicon2$`bedFile$V2[bedIdx]`,
                  data.type="logratio",sampleid="2163lot_new")

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
View(segments.p(segment.smoothed.CNA.object))


plot(segment.smoothed.CNA.object, plot.type="s")



