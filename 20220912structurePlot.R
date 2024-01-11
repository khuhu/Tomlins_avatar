require(ggplot2)

### structure plots seem to go in the format where they list the reference populations in the first few bars
### showing them with there color, then you have the admixed samples with example in link below:
### https://www.genetics.org/content/210/2/719.figures-only


### testing below

### for our Q table the columns needed for the correct format will be SampleName, Ref Pop, value
### so that for k number of ref poputations, there will be k entries of SampleName for ref pop and their estimated ancestry
### for that ref pop


### best results I've had so far is the 100 reps using 500 SNP markers. gets more than half of the 8th and some 9th

sampleMap <- read.table("/br_z1/kevin_storage/advancedGenomicsCore/tmp/Sample_Map.txt", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)

# qtable <- read.table("/br_z1/kevin_storage/ADMIXTURE/merge48ForAdmixture.5.Q", sep = " ",
#                      header = FALSE, stringsAsFactors = FALSE)


qtable <- read.table("/br_z1/kevin_storage/ADMIXTURE/20220912merge48Autosomes.5.Q", sep = " ",
                     header = FALSE, stringsAsFactors = FALSE)

colnames(qtable) <- c("129", "BALB", "C3H", "C57BL/J", "FVB001")

# qtable_small <- qtable[c(71:76, 23),]
# qtable_small <- apply(qtable_small, 2, function(x) signif(x, 2))

qtable_small <- qtable
qtable_small <- apply(qtable_small, 2, function(x) signif(x, 2))

qtable_small_melt  <- reshape2::melt(qtable_small )
qtable_small_melt$Var2 <- as.character(qtable_small_melt$Var2)
qtable_small_melt <- qtable_small_melt[,2:3]
qtable_small_melt$SampleName <- paste0("5516-KH-", 1:5)
colnames(qtable_small_melt)[1:2] <- c("Strain", "Fraction")

fullQTable2 <- qtable_small_melt


fullQTable2$SampleName <- sampleMap$ID[match(fullQTable2$SampleName, sampleMap$Name)]

# fullQTable2 <- cbind(rownames(fullQTable), fullQTable,
#                      paste0(signif(fullQTable$`6J`, digits = 2), "/", signif(fullQTable$`129S1`, digits = 2)))
# fullQTable2 <- melt(fullQTable2)
# fullQTable2$variable <- as.character(fullQTable2$variable)
# fullQTable2$`rownames(fullQTable)` <- as.character(fullQTable2$`rownames(fullQTable)`)
# fullQTable2 <- fullQTable2[order(fullQTable2$value, decreasing = TRUE),]
# fullQTable2$value <- signif(fullQTable2$value, digits = 2)
# colnames(fullQTable2) <- c("SampleName", "label","Strain", "Fraction")
# fullQTable2$SampleName <- factor(fullQTable2$SampleName, levels = unique(fullQTable2$SampleName))
# fullQTable2$label <- as.character(fullQTable2$label)
# tmp1 <- nrow(fullQTable2)/2 + 1
# tmp2 <- nrow(fullQTable2)
# fullQTable2$label[eval(tmp1:tmp2)] <- " "

colPalette <- c("#009E73","#D55E00", "#8b0000", "#800080", "#ffff00")

# ggplot(fullQTable2, aes(x = SampleName, y = Fraction, fill=Strain)) +
#   geom_bar(position="fill", stat="identity", color = "black") + 
#   geom_text(aes(label = label, y = 1.01), size = 2) + 
#   scale_fill_manual(values=colPalette) + theme_bw() + 
#   theme(axis.text.x=element_text(angle=90,hjust=1))

### for more sampels

#ggplot(fullQTable2, aes(x = SampleName, y = Fraction, fill=Strain)) +
#  geom_bar(position="fill", stat="identity") + 
#  scale_fill_manual(values=colPalette) + theme_bw() + 
#  theme(axis.text.x=element_text(angle=90,hjust=1))



ggplot(fullQTable2, aes(x = SampleName, y = Fraction, fill=Strain)) +
  geom_bar(position="fill", stat="identity", color = "black") + 
  scale_fill_manual(values=colPalette) + theme_bw() + 
  theme(axis.text.x=element_text(angle=90,hjust=1))

