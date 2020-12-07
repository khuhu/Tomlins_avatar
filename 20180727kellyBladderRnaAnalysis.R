library(edgeR)
library(ggplot2)
library(ggrepel)

tableOfDat <- read.table(file = "/mnt/DATA4/kevhu/KellyData/201807bladderRna1.txt", sep = "\t", header = TRUE)
tableOfDat <- tableOfDat[-1,]
rownames(tableOfDat) <- tableOfDat$Gene
tableOfDat <- tableOfDat[,4:ncol(tableOfDat)]
tableOfDat <- tableOfDat[,order(colnames(tableOfDat))]

dgeObj <- DGEList(counts = tableOfDat)
bioExpressFilt <- c(unlist(which(rowSums(cpm(dgeObj) > 5) >= 2)))

dgeObj <- dgeObj[bioExpressFilt, , keep.lib.sizes = FALSE]
dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples

MorP <- factor(c(rep(c("M","P"),7), "M"))
design <- model.matrix(~0 + MorP, data=dgeObj$samples)
design

dgeObj <- estimateDisp(dgeObj, design)
fit <- glmFit(dgeObj, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
topTags(lrt)
lrt.table <- lrt$table
lrtPvals <- lrt$table$PValue
lrtPVals.ad <- p.adjust(lrtPvals, method = c("BH"))
lrt$table$QVal <- lrtPVals.ad
lrt.tab <- lrt$table
lrt.tab$Gene <- rownames(lrt.tab)
lrt.tab$color <- NULL


for(i in seq_along(lrt.tab$logFC)){
  if(abs(lrt.tab$logFC[i]) > 1 & lrt.tab$QVal[i] < 0.05){
    lrt.tab$color[i]  <- "red"
  }
  else{
    lrt.tab$color[i] <- "blue"
  }
}

plot2 <- ggplot(data = lrt.tab,aes(logFC, -log10(QVal)))
plot2 <- plot2 +  geom_point(pch= 20, color = lrt.tab$color, size = 1.25) + geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.25)
plot2 <-  plot2 + ggtitle("P vs M")
plot2 <- plot2 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot2 <- plot2 + geom_text_repel(data = lrt.tab[which(lrt.tab$color == "red"),], aes(label = Gene, vjust =-1.0), size = 3, alpha = 1.0)

setwd("/mnt/DATA4/kevhu/KellyData/")
sig.genes <- lrt.tab[which(lrt.tab$color == "red"),]
if(nrow(sig.genes) == 0){
  sig.genes <- topTags(lrt)
}
gNames <- rownames(sig.genes)
sig.genes <- cbind(gNames,sig.genes)
sig.genes <- sig.genes[,c(1,2,3,5,6)]
colnames(sig.genes)[1] <- "Genes"
write.table(sig.genes,"20180730.PvM.txt", quote = FALSE, row.names = FALSE, sep = "\t")


pdf(file = "20180730.PvM.pdf",onefile = TRUE, width = 5)
plot2
dev.off()





### trying above with paired analysis
Paired <- factor(c("1","1","2","2","3","3","4","4","none","5","5","none","6","6","none"))
design2 <- model.matrix(~0 + MorP + Paired, data=dgeObj$samples)
design2



dgeObj2 <- estimateDisp(dgeObj, design2)
fit2 <- glmFit(dgeObj2, design2)
lrt2<- glmLRT(fit2, contrast=c(-1,1,rep(0,6)))
topTags(lrt2)
lrt.table2 <- lrt2$table
lrtPvals2 <- lrt2$table$PValue
lrtPVals.ad2 <- p.adjust(lrtPvals2, method = c("BH"))
lrt2$table$QVal <- lrtPVals.ad2
lrt.tab2 <- lrt2$table
lrt.tab2$Gene <- rownames(lrt.tab2)
lrt.tab2$color <- NULL


for(i in seq_along(lrt.tab2$logFC)){
  if(abs(lrt.tab2$logFC[i]) > 1 & lrt.tab2$QVal[i] < 0.05){
    lrt.tab2$color[i]  <- "red"
  }
  else{
    lrt.tab2$color[i] <- "blue"
  }
}

plot3 <- ggplot(data = lrt.tab2,aes(logFC, -log10(QVal)))
plot3 <- plot3 +  geom_point(pch= 20, color = lrt.tab2$color, size = 1.25) + geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.25)
plot3 <-  plot3 + ggtitle("P vs M: corrected for paired effect" )
plot3 <- plot3 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot3 <- plot3 +  geom_text_repel(data = lrt.tab2[which(lrt.tab2$color == "red"),], aes(label = Gene, vjust =-1.0), size = 3, alpha = 1.0)

sig.genes2 <- lrt.tab2[which(lrt.tab2$color == "red"),]
if(nrow(sig.genes) == 0){
  sig.genes2 <- topTags(lrt2)
}
gNames <- rownames(sig.genes2)
sig.genes2 <- cbind(gNames,sig.genes2)
sig.genes2 <- sig.genes2[,c(1,2,3,5,6)]
colnames(sig.genes2)[1] <- "Genes"

write.table(sig.genes2,"20180730.PvMAdjusted.txt", quote = FALSE, row.names = FALSE, sep = "\t")


pdf(file = "20180730.PvMAdjusted.pdf",onefile = TRUE, width = 5)
plot3
dev.off()



###this portion removes part of the data i.e ones that are unpaired and bad couple
###
###
###
###
###

tableOfDat2 <- tableOfDat
tableOfDat2 <- tableOfDat2[,c(1:8,10,11)]

dgeObj <- DGEList(counts = tableOfDat2)
bioExpressFilt <- c(unlist(which(rowSums(cpm(dgeObj) > 5) >= 2)))

dgeObj <- dgeObj[bioExpressFilt, , keep.lib.sizes = FALSE]
dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples

MorP <- factor(c(rep(c("M","P"),5)))
design <- model.matrix(~0 + MorP, data=dgeObj$samples)
design

dgeObj <- estimateDisp(dgeObj, design)
fit <- glmFit(dgeObj, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
topTags(lrt)
lrt.table <- lrt$table
lrtPvals <- lrt$table$PValue
lrtPVals.ad <- p.adjust(lrtPvals, method = c("BH"))
lrt$table$QVal <- lrtPVals.ad
lrt.tab <- lrt$table
lrt.tab$Gene <- rownames(lrt.tab)
lrt.tab$color <- NULL


for(i in seq_along(lrt.tab$logFC)){
  if(abs(lrt.tab$logFC[i]) > 1 & lrt.tab$QVal[i] < 0.05){
    lrt.tab$color[i]  <- "red"
  }
  else{
    lrt.tab$color[i] <- "blue"
  }
}


for(i in seq_along(lrt.tab$logFC)){
  if(abs(lrt.tab$logFC[i]) > 5 | lrt.tab$QVal[i] < 1e-10){
    lrt.tab$secondThres[i]  <- "Yes"
  }
  else{
    lrt.tab$secondThres[i] <- "No"
  }
}


plot2 <- ggplot(data = lrt.tab,aes(logFC, -log10(QVal)))
plot2 <- plot2 +  geom_point(pch= 20, color = lrt.tab$color, size = 1.25) + geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.25)
plot2 <-  plot2 + ggtitle("P vs M")
plot2 <- plot2 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot2 <- plot2 + geom_text_repel(data = lrt.tab[which(lrt.tab$secondThres == "Yes"),], aes(label = Gene, vjust =-1.0), size = 3, alpha = 1.0)
plot2

setwd("/mnt/DATA4/kevhu/KellyData/")
sig.genes <- lrt.tab[which(lrt.tab$color == "red"),]
if(nrow(sig.genes) == 0){
  sig.genes <- topTags(lrt)
}
gNames <- rownames(sig.genes)
sig.genes <- cbind(gNames,sig.genes)
sig.genes <- sig.genes[,c(1,2,3,5,6)]
colnames(sig.genes)[1] <- "Genes"
write.table(sig.genes,"20180731.PvM.reduced.txt", quote = FALSE, row.names = FALSE, sep = "\t")


pdf(file = "20180731.PvM.reduced.pdf",onefile = TRUE, width = 5)
plot2
dev.off()





### trying above with paired analysis
Paired <- factor(c("1","1","2","2","3","3","4","4","5","5"))
design2 <- model.matrix(~0 + MorP + Paired, data=dgeObj$samples)
design2



dgeObj2 <- estimateDisp(dgeObj, design2)
fit2 <- glmFit(dgeObj2, design2)
lrt2<- glmLRT(fit2, contrast=c(-1,1,rep(0,4)))
topTags(lrt2)
lrt.table2 <- lrt2$table
lrtPvals2 <- lrt2$table$PValue
lrtPVals.ad2 <- p.adjust(lrtPvals2, method = c("BH"))
lrt2$table$QVal <- lrtPVals.ad2
lrt.tab2 <- lrt2$table
lrt.tab2$Gene <- rownames(lrt.tab2)
lrt.tab2$color <- NULL
lrt.tab2$secondThres <- NULL

for(i in seq_along(lrt.tab2$logFC)){
  if(abs(lrt.tab2$logFC[i]) > 1 & lrt.tab2$QVal[i] < 0.05){
    lrt.tab2$color[i]  <- "red"
  }
  else{
    lrt.tab2$color[i] <- "blue"
  }
}

for(i in seq_along(lrt.tab2$logFC)){
  if(abs(lrt.tab2$logFC[i]) > 5 | lrt.tab2$QVal[i] < 1e-10){
    lrt.tab2$secondThres[i]  <- "Yes"
  }
  else{
    lrt.tab2$secondThres[i] <- "No"
  }
}


plot3 <- ggplot(data = lrt.tab2,aes(logFC, -log10(QVal)))
plot3 <- plot3 +  geom_point(pch= 20, color = lrt.tab2$color, size = 1.25) + geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.25)
plot3 <-  plot3 + ggtitle("P vs M: corrected for paired effect" )
plot3 <- plot3 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot3 <- plot3 +  geom_text_repel(data = lrt.tab2[which(lrt.tab2$secondThres == "Yes"),], aes(label = Gene, vjust =-1.0), size = 3, alpha = 1.0)
plot3

sig.genes2 <- lrt.tab2[which(lrt.tab2$color == "red"),]
if(nrow(sig.genes) == 0){
  sig.genes2 <- topTags(lrt2)
}
gNames <- rownames(sig.genes2)
sig.genes2 <- cbind(gNames,sig.genes2)
sig.genes2 <- sig.genes2[,c(1,2,3,5,6)]
colnames(sig.genes2)[1] <- "Genes"

write.table(sig.genes2,"20180731.PvMAdjusted.reduced.txt", quote = FALSE, row.names = FALSE, sep = "\t")


pdf(file = "20180731.PvMAdjusted.reduced.pdf",onefile = TRUE, width = 5)
plot3
dev.off()

### there were some weird results and I should've checked this sooner, but PCA plot for Kelly's data - code taken from lorena PCA plot
### done on local testingDESeq.R 

dgeObj$counts
pca.out <- prcomp(t(log2(dgeObj$counts + 1)))
pca.graph.dat <- data.frame("PC1" = pca.out$x[,1], "PC2"=pca.out$x[,2], "class" = factor(c(rep(c("M","P"),4), "P",rep(c("P","M"),3))),
                            "pairing" = factor(c("1","1","2","2","3","3","4","4","none","5","5","none","6","6","none")),
                            "batch" = factor(c(rep(1,6),2,2,2,1,1,rep(2,4))))

screeplot(pca.out)


ggplot(data = pca.graph.dat,aes(PC1,PC2,colour=class, shape=batch, label=rownames(pca.graph.dat))) +
  geom_point(size = 2.5) + geom_text(nudge_y = 4) +
  ggtitle(label = "PCA biplot by samples PC1 vs PC2") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

#### below is just a guideline and code I used for lorena's stuff
ggplot(data = pca.graph.data) +
  geom_point(aes(PC1,PC2,colour=class, shape=pairing),size = 2.5) +
  geom_point(data = subset(pca.graph.data, CTBNN1 == "Mut"),
             aes(PC1,PC2,shape=AdenoSquam),colour="white",size = 1.0,inherit.aes = FALSE) +
  scale_colour_brewer(palette = "Paired", labels= a, name = c("Samples")) +
  #geom_text(aes(label=pca.graph.data$sampleName), size = 1.5, nudge_y = 5) +
  ggtitle(label = "PCA biplot by samples PC1 vs PC2") + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

### next section will be the newer compariosns by group

newParied <- c("BL464","BL465","BL466","BL467","BL468","BL469","BL472","BL473")
newUnpaired <- c("BL463","BL464","BL465","BL466","BL467","BL468","BL469","BL471","BL472","BL473","BL476","BL478")
newUnpaired2 <- c("BL463","BL464","BL465","BL466","BL467","BL468","BL469","BL472","BL473","BL476","BL478")


tableOfDat.paired <- tableOfDat[, which(colnames(tableOfDat) %in% newParied)]
tableOfDat.unpaired <- tableOfDat[, which(colnames(tableOfDat) %in% newUnpaired)]
tableOfDat.unpaired2 <- tableOfDat[, which(colnames(tableOfDat) %in% newUnpaired2)]

dgeObj.paired <- DGEList(counts = tableOfDat.paired)
bioExpressFilt.paired <- c(unlist(which(rowSums(cpm(dgeObj.paired) > 5) >= 2)))

dgeObj.paired <- dgeObj.paired[bioExpressFilt.paired, , keep.lib.sizes = FALSE]
dgeObj.paired <- calcNormFactors(dgeObj.paired)
dgeObj.paired$samples

MorP.paired <- factor(c(rep(c("M","P"),4)))
newPaired.classes <- factor(c(1,1,2,2,3,3,4,4))

design.paired <- model.matrix(~0 + MorP.paired + newPaired.classes, data=dgeObj.paired$samples)
design.paired

dgeObj.paired <- estimateDisp(dgeObj.paired, design.paired)
fit.paired <- glmFit(dgeObj.paired, design.paired)
lrt.paired <- glmLRT(fit.paired, contrast=c(-1,1,0,0,0))
topTags(lrt.paired)
lrt.table <- lrt.paired$table
lrtPvals <- lrt.paired$table$PValue
lrtPVals.ad <- p.adjust(lrtPvals, method = c("BH"))
lrt.paired$table$QVal <- lrtPVals.ad
lrt.tab <- lrt.paired$table
lrt.tab$Gene <- rownames(lrt.tab)
lrt.tab$color <- NULL


for(i in seq_along(lrt.tab$logFC)){
  if(abs(lrt.tab$logFC[i]) > 1 & lrt.tab$QVal[i] < 0.05){
    lrt.tab$color[i]  <- "red"
  }
  else{
    lrt.tab$color[i] <- "blue"
  }
}


for(i in seq_along(lrt.tab$logFC)){
  if(abs(lrt.tab$logFC[i]) > 5 | lrt.tab$QVal[i] < 1e-10){
    lrt.tab$secondThres[i]  <- "Yes"
  }
  else{
    lrt.tab$secondThres[i] <- "No"
  }
}


plot2 <- ggplot(data = lrt.tab,aes(logFC, -log10(QVal)))
plot2 <- plot2 +  geom_point(pch= 20, color = lrt.tab$color, size = 1.25) + geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.25)
plot2 <-  plot2 + ggtitle("P vs M paired analysis")
plot2 <- plot2 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot2 <- plot2 + geom_text_repel(data = lrt.tab[which(lrt.tab$secondThres == "Yes"),], aes(label = Gene, vjust =-1.0), size = 3, alpha = 1.0)
plot2

setwd("/mnt/DATA4/kevhu/KellyData/")
sig.genes <- lrt.tab[which(lrt.tab$color == "red"),]
if(nrow(sig.genes) == 0){
  sig.genes <- topTags(lrt)
}
gNames <- rownames(sig.genes)
sig.genes <- cbind(gNames,sig.genes)
sig.genes <- sig.genes[,c(1,2,3,5,6)]
colnames(sig.genes)[1] <- "Genes"
write.table(sig.genes,"20180802.PvM.paired.txt", quote = FALSE, row.names = FALSE, sep = "\t")


pdf(file = "20180802.PvM.paired.pdf",onefile = TRUE, width = 5)
plot2
dev.off()


###

dgeObj.unpaired <- DGEList(counts = tableOfDat.unpaired)
#dgeObj.unpaired <- DGEList(counts = tableOfDat.unpaired2)
bioExpressFilt.ununpaired <- c(unlist(which(rowSums(cpm(dgeObj.unpaired) > 5) >= 2)))

dgeObj.unpaired <- dgeObj.unpaired[bioExpressFilt.unpaired, , keep.lib.sizes = FALSE]
dgeObj.unpaired <- calcNormFactors(dgeObj.unpaired)
dgeObj.unpaired$samples

MorP.unpaired <- factor(c("P","M","P","M","P","M","P","P","M","P","M","M"))
#MorP.unpaired <- factor(c("P","M","P","M","P","M","P","M","P","M","M"))

design.unpaired <- model.matrix(~0 + MorP.unpaired, data=dgeObj.unpaired$samples)
design.unpaired

dgeObj.unpaired <- estimateDisp(dgeObj.unpaired, design.unpaired)
fit.unpaired <- glmFit(dgeObj.unpaired, design.unpaired)
lrt.unpaired <- glmLRT(fit.unpaired, contrast=c(-1,1))
topTags(lrt.unpaired)
lrt.table <- lrt.unpaired$table
lrtPvals <- lrt.unpaired$table$PValue
lrtPVals.ad <- p.adjust(lrtPvals, method = c("BH"))
lrt.unpaired$table$QVal <- lrtPVals.ad
lrt.tab <- lrt.unpaired$table
lrt.tab$Gene <- rownames(lrt.tab)
lrt.tab$color <- NULL


for(i in seq_along(lrt.tab$logFC)){
  if(abs(lrt.tab$logFC[i]) > 1 & lrt.tab$QVal[i] < 0.05){
    lrt.tab$color[i]  <- "red"
  }
  else{
    lrt.tab$color[i] <- "blue"
  }
}


for(i in seq_along(lrt.tab$logFC)){
  if(abs(lrt.tab$logFC[i]) > 5 | lrt.tab$QVal[i] < 1e-10){
    lrt.tab$secondThres[i]  <- "Yes"
  }
  else{
    lrt.tab$secondThres[i] <- "No"
  }
}


plot2 <- ggplot(data = lrt.tab,aes(logFC, -log10(QVal)))
plot2 <- plot2 +  geom_point(pch= 20, color = lrt.tab$color, size = 1.25) + geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.25)
plot2 <-  plot2 + ggtitle("P vs M unpaired analysis")
plot2 <- plot2 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot2 <- plot2 + geom_text_repel(data = lrt.tab[which(lrt.tab$secondThres == "Yes"),], aes(label = Gene, vjust =-1.0), size = 3, alpha = 1.0)
plot2

setwd("/mnt/DATA4/kevhu/KellyData/")
sig.genes <- lrt.tab[which(lrt.tab$color == "red"),]
if(nrow(sig.genes) == 0){
  sig.genes <- topTags(lrt)
}
gNames <- rownames(sig.genes)
sig.genes <- cbind(gNames,sig.genes)
sig.genes <- sig.genes[,c(1,2,3,5,6)]
colnames(sig.genes)[1] <- "Genes"
write.table(sig.genes,"20180802.PvM.unpairedNo71.txt", quote = FALSE, row.names = FALSE, sep = "\t")


pdf(file = "20180802.PvM.unpaired.pdf",onefile = TRUE, width = 5)
plot2
dev.off()









