library(ggplot2)

### creating structure plots for previous admixture subsmapling examples

fullQTable <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/vcfs/20190624/1633/plink.5.Q",
                         sep = " ", stringsAsFactors = FALSE, header = FALSE)
sampNames <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/vcfs/20190624/1633/plink.fam",
                        stringsAsFactors = FALSE, sep = " ")

colnames(fullQTable) <- c("129", "BALB", "C3H", "FVB", "6J")

#fullQTable2 <- cbind(sampNames$V2, fullQTable,
#                     paste0(signif(fullQTable$`6J`, digits = 2), "/", signif(fullQTable$`129S1`, digits = 2)))
fullQTable2 <- data.frame(cbind(sampNames$V2, fullQTable), stringsAsFactors = FALSE)


fullQTable2 <- melt(fullQTable2)
fullQTable2$variable <- as.character(fullQTable2$variable)
#fullQTable2$`rownames(fullQTable)` <- fullQTable2$sampNames.V2
fullQTable2 <- fullQTable2[order(fullQTable2$value, decreasing = TRUE),]
fullQTable2$value <- signif(fullQTable2$value, digits = 2)
colnames(fullQTable2) <- c("SampleName", "label","Fraction")
fullQTable2$SampleName <- factor(fullQTable2$SampleName, levels = unique(fullQTable2$SampleName))
fullQTable2$label <- as.character(fullQTable2$label)
fullQTable2$label <- str_remove(fullQTable2$label, "X")
#tmp1 <- nrow(fullQTable2)/2 + 1
#tmp2 <- nrow(fullQTable2)
#fullQTable2$label[eval(tmp1:tmp2)] <- " "

colPalette <- c("#009E73","#D55E00", "#001a9e", "#9e007e", "#9e9600")

ggplot(fullQTable2, aes(x = SampleName, y = Fraction, fill=label)) +
  geom_bar(position="fill", stat="identity", color = "black") + 
  scale_fill_manual(values=colPalette) + theme_bw() + 
  theme(axis.text.x=element_text(angle=90,hjust=1))




### structure plot for only looking at 300

fullQTable <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/vcfs/20190624/300/myplinkBed.5.Q",
                         sep = " ", stringsAsFactors = FALSE, header = FALSE)
sampNames <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/vcfs/20190624/300/myplinkBed.fam",
                        stringsAsFactors = FALSE, sep = " ")

colnames(fullQTable) <- c("129", "BALB", "C3H", "FVB", "6J")

#fullQTable2 <- cbind(sampNames$V2, fullQTable,
#                     paste0(signif(fullQTable$`6J`, digits = 2), "/", signif(fullQTable$`129S1`, digits = 2)))
fullQTable2 <- data.frame(cbind(sampNames$V2, fullQTable), stringsAsFactors = FALSE)


fullQTable2 <- melt(fullQTable2)
fullQTable2$variable <- as.character(fullQTable2$variable)
#fullQTable2$`rownames(fullQTable)` <- fullQTable2$sampNames.V2
fullQTable2 <- fullQTable2[order(fullQTable2$value, decreasing = TRUE),]
fullQTable2$value <- signif(fullQTable2$value, digits = 2)
colnames(fullQTable2) <- c("SampleName", "label","Fraction")
fullQTable2$SampleName <- factor(fullQTable2$SampleName, levels = unique(fullQTable2$SampleName))
fullQTable2$label <- as.character(fullQTable2$label)
fullQTable2$label <- str_remove(fullQTable2$label, "X")
#tmp1 <- nrow(fullQTable2)/2 + 1
#tmp2 <- nrow(fullQTable2)
#fullQTable2$label[eval(tmp1:tmp2)] <- " "

colPalette <- c("#009E73","#D55E00", "#001a9e", "#9e007e", "#9e9600")

ggplot(fullQTable2, aes(x = SampleName, y = Fraction, fill=label)) +
  geom_bar(position="fill", stat="identity", color = "black") + 
  scale_fill_manual(values=colPalette) + theme_bw() + 
  theme(axis.text.x=element_text(angle=90,hjust=1))



### structure plot for only looking at 150

fullQTable <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/vcfs/20190624/150/myplinkBed.5.Q",
                         sep = " ", stringsAsFactors = FALSE, header = FALSE)
sampNames <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/vcfs/20190624/150/myplinkBed.fam",
                        stringsAsFactors = FALSE, sep = " ")

colnames(fullQTable) <- c("129", "BALB", "C3H", "FVB", "6J")

#fullQTable2 <- cbind(sampNames$V2, fullQTable,
#                     paste0(signif(fullQTable$`6J`, digits = 2), "/", signif(fullQTable$`129S1`, digits = 2)))
fullQTable2 <- data.frame(cbind(sampNames$V2, fullQTable), stringsAsFactors = FALSE)


fullQTable2 <- melt(fullQTable2)
fullQTable2$variable <- as.character(fullQTable2$variable)
#fullQTable2$`rownames(fullQTable)` <- fullQTable2$sampNames.V2
fullQTable2 <- fullQTable2[order(fullQTable2$value, decreasing = TRUE),]
fullQTable2$value <- signif(fullQTable2$value, digits = 2)
colnames(fullQTable2) <- c("SampleName", "label","Fraction")
fullQTable2$SampleName <- factor(fullQTable2$SampleName, levels = unique(fullQTable2$SampleName))
fullQTable2$label <- as.character(fullQTable2$label)
fullQTable2$label <- str_remove(fullQTable2$label, "X")
#tmp1 <- nrow(fullQTable2)/2 + 1
#tmp2 <- nrow(fullQTable2)
#fullQTable2$label[eval(tmp1:tmp2)] <- " "

colPalette <- c("#009E73","#D55E00", "#001a9e", "#9e007e", "#9e9600")

ggplot(fullQTable2, aes(x = SampleName, y = Fraction, fill=label)) +
  geom_bar(position="fill", stat="identity", color = "black") + 
  scale_fill_manual(values=colPalette) + theme_bw() + 
  theme(axis.text.x=element_text(angle=90,hjust=1))





