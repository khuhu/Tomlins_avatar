##### Start of code that works
####
####
samp <- allNames[1]

library(ggplot2)
df$copynumber <- log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]])
#df.subset <- df[which(df$ChromNum < 5),]
library(ggplot2)
geneSpacing <- 5;
chromSpacing <- 40;
pointSpacing <- 40;
maxGenes <- 7;
dummyList <- c(1:7);
dummyGeneList <- unique(df$Gene)

#ONLY USE THE IS.EVEN FUNCTION BELOW IF YOU'RE LOOKING AT SEPARATE AMPLICONS
for(i in seq_along(df$GeneNum)){
  if(df$GeneNum[i] %% 2 == 0){
    df$GeneNum[i] <- df$GeneNum[i] - 1
  }
}

for(i in seq_along(df$GeneNum)){
  if(df$GeneNum[i] == 1){
    next()
  }
  else if(df$GeneNum[i] %% 2 != 0){
    df$GeneNum[i] <- df$GeneNum[i] - 1
  }
}


df$xCoord <- pointSpacing*(geneSpacing*df$GeneNum + df$AmpliconIndex)
#df$xCoord <- (geneSpacing*df$GeneNum + df$AmpliconIndex)

#mediansPerGene <- NULL
#xmin <- NULL
#xmax <- NULL
#Group <- NULL
geneBreaks <- NULL
for(i in 1:ceiling(length(unique(df$Gene))/maxGenes)){
  dummyGenes <- dummyGeneList[((i-1)*7) + dummyList]
  for(j in seq_along(dummyGenes))
  {
    df[which(df$Gene == dummyGenes[j]),"Group"] <- i 
    df[which(df$Gene == dummyGenes[j]),"geneEst"] <- log2(geneEst[which(geneEst$Gene == dummyGenes[j]), samp])
    #df[which(df$Gene == dummyGenes[j]),"medianCN"] <- median(df[which(df$Gene == dummyGenes[j]),"copynumber"])
    #xmin <- c(xmin,min(df[which(df$Gene == dummyGenes[j]),"xCoord"]))
    #xmax <- c(xmax, max(df[which(df$Gene == dummyGenes[j]),"xCoord"]))
    #Group <- c(Group, i)
    geneBreaks <- c(geneBreaks,mean(df[which(df$Gene == dummyGenes[j]),"xCoord"]))
  }
}


#df$Group <- factor(df$Group)
#mediansPerGene <- cbind(mediansPerGene, mediansPerGene, xmin, xmax, Group)
#mediansPerGene[is.infinite(mediansPerGene)] <- NA
#colnames(mediansPerGene) <- c("y","yend","x","xend","Group")
#mediansPerGene <- data.frame(mediansPerGene, stringsAsFactors = FALSE)
#dfMed <- df[,c("xCoord","medianCN","Gene")]

library(ggplot2)

testggplot <- ggplot(data = df, aes(x = df$xCoord, y = df$copynumber, colour = factor(df$Color), group = factor(df$Gene)))
testggplot <- testggplot + geom_point(size = 0.75, pch = 1) + theme_bw() + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + facet_wrap(~ Group, scales = "free_x", ncol = 2)
#testggplot <- testggplot + geom_point(size = 0.75, pch = 1) + geom_line(size = 0.1) + theme_bw() + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + facet_wrap(~ Group, scales = "free_x", ncol = 2)
#testggplot <- testggplot + geom_point(aes(x = df$xCoord, y = df$medianCN, group = factor(df$Gene)), colour = "black", size = 0.05) + geom_line(size = 0.05) + geom_point(size = 0.75, pch = 1) + geom_line(size = 0.1) + theme_bw() + theme(legend.position="none") + facet_wrap(~ Group, scales = "free")
#testggplot <- testggplot + geom_line(aes(x = df$xCoord, y = df$medianCN), colour = "black", size = 0.05)
testggplot <- testggplot + geom_line(aes(x = df$xCoord, y = df$geneEst), colour = "black", size = 0.05)
testggplot <- testggplot + ylab(label = "log2(copyNumber)") + xlab(label = "") 


geneLabels = unique(df$Gene) 
geneBreaks <- geneBreaks[1:length(geneLabels)]
chromLabels <- unique(df$ChromNum)
chromBreaks <- NULL
posAblines <- NULL
xpos <- NULL
chromDat <- NULL

for(i in seq_along(unique(df$ChromNum))){
  xpos <- mean(df[which(df$ChromNum == chromLabels[i]),"xCoord"])
  xpos <- c(xpos, unique(df[which(df$ChromNum == chromLabels[i]),"Group"]))
  chromBreaks <- rbind(chromBreaks, xpos)
  start <- min(df[which(df$ChromNum == chromLabels[i]),"xCoord"]) - 2*chromSpacing
  start <- c(start, unique(df[which(df$ChromNum == chromLabels[i]),"Group"]))
  #stop <- min(df[which(df$ChromNum == chromLabels[i]),"xCoord"]) + chromSpacing
  #stop <- c(stop, unique(df[which(df$ChromNum == chromLabels[i]),"Group"]))
  #posAblines <- rbind(posAblines, start, stop)
  posAblines <- rbind(posAblines, start)
}

chromDat <- data.frame("chromLabels" = chromLabels, "xpos" = chromBreaks[,1], "y" = rep(round(min(df$copynumber)),length(chromBreaks)/2), "Group" = factor(chromBreaks[,2]))
dummy1 <- data.frame(x = posAblines[,1], Group = factor(posAblines[,2]))
testggplot <- testggplot + scale_x_continuous(breaks = geneBreaks, labels = geneLabels) + theme(axis.text.x=element_text(angle = 90, hjust = 0))
testggplot <- testggplot + geom_vline(data = dummy1, aes(xintercept = x), colour = "gray")
testggplot <- testggplot + geom_text(data = chromDat, aes(x = chromDat$xpos, y = chromDat$y, label = chromDat$chromLabels), inherit.aes = FALSE)

#testggplot <- testggplot + geom_text(data = chromDat, aes(x = chromDat$chromBreaks, y = chromDat$y,label = chromDat$chromLabels))
#testggplot <- testggplot + stat_summary(fun.y=min,aes(label=paste0('N=',geneLabels)),geom='text',col='blue',cex=5)
#testggplot <- testggplot + scale_x_discrete(labels=c("0.5" = "Dose 0.5", "1" = "Dose 1","2" = "Dose 2"))
#testggplot + geom_hline(yintercept = mediansPerGene$y, inherit.aes + FALSE)

#for(i in 1:nrow(mediansPerGene)){
#  testggplot <- testggplot + geom_segment(data = mediansPerGene,aes_string(x = mediansPerGene$x[i],y = mediansPerGene$y[i],
#                                              xend = mediansPerGene$xend[i],yend = mediansPerGene$yend[i], group = mediansPerGene$Group[i]), colour = "black")
#}

testggplot

pdf(file = "/home/kevhu/data/15622TP1Quant.pdf",onefile = TRUE, height = 5*ceiling(length(unique(df$Gene))/maxGenes) , width = 30.00)
testggplot 
dev.off()




###actual implementation - going to put it in the for loop at the end


longpsamp <- function(samp, .....){
  geneSpacing <- 5;
  chromSpacing <- 40;
  pointSpacing <- 40;
  maxGenes <- 7;
  dummyList <- c(1:7);
  dummyGeneList <- unique(df$Gene)
  geneBreaks <- NULL
  df$copynumber <- log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]])
  df$xCoord <- pointSpacing*(geneSpacing*df$GeneNum + df$AmpliconIndex)
  for(i in 1:ceiling(length(unique(df$Gene))/maxGenes)){
    dummyGenes <- dummyGeneList[((i-1)*7) + dummyList]
    for(j in seq_along(dummyGenes))
    {
      df[which(df$Gene == dummyGenes[j]),"Group"] <- i 
      df[which(df$Gene == dummyGenes[j]),"geneEst"] <- log2(geneEst[which(geneEst$Gene == dummyGenes[j]), samp])
      geneBreaks <- c(geneBreaks,mean(df[which(df$Gene == dummyGenes[j]),"xCoord"]))
    }
  }
  testggplot <- ggplot(data = df, aes(x = df$xCoord, y = df$copynumber, colour = factor(df$Color), group = factor(df$Gene)))
  testggplot <- testggplot + geom_point(size = 0.75, pch = 1) + theme_bw() + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + facet_wrap(~ Group, scales = "free_x", ncol = 2)
  testggplot <- testggplot + geom_line(aes(x = df$xCoord, y = df$geneEst), colour = "black", size = 0.05)
  testggplot <- testggplot + ylab(label = "log2(copyNumber)") + xlab(label = "") 
  testggplot <- testggplot + geom_text(data = df[df$AmpliconId %in% outlierList,], aes(xCoord, copynumber,label = AmpliconId))
  geneLabels = unique(df$Gene) 
  geneBreaks <- geneBreaks[1:length(geneLabels)]
  chromLabels <- unique(df$ChromNum)
  chromBreaks <- NULL
  posAblines <- NULL
  xpos <- NULL
  chromDat <- NULL
  for(i in seq_along(unique(df$ChromNum))){
    xpos <- mean(df[which(df$ChromNum == chromLabels[i]),"xCoord"])
    xpos <- c(xpos, unique(df[which(df$ChromNum == chromLabels[i]),"Group"]))
    chromBreaks <- rbind(chromBreaks, xpos)
    start <- min(df[which(df$ChromNum == chromLabels[i]),"xCoord"]) - 2*chromSpacing
    start <- c(start, unique(df[which(df$ChromNum == chromLabels[i]),"Group"]))
    posAblines <- rbind(posAblines, start)
  }
  chromDat <- data.frame("chromLabels" = chromLabels, "xpos" = chromBreaks[,1], "y" = rep(round(min(df$copynumber)),length(chromBreaks)/2), "Group" = factor(chromBreaks[,2]))
  dummy1 <- data.frame(x = posAblines[,1], Group = factor(posAblines[,2]))
  testggplot <- testggplot + scale_x_continuous(breaks = geneBreaks, labels = geneLabels) + theme(axis.text.x=element_text(angle = 90, hjust = 0))
  testggplot <- testggplot + geom_vline(data = dummy1, aes(xintercept = x), colour = "gray")
  testggplot <- testggplot + geom_text(data = chromDat, aes(x = chromDat$xpos, y = chromDat$y, label = chromDat$chromLabels), inherit.aes = FALSE)
  print(testggplot)
  #pdf(file = paste("long.",samp, ".pdf",sep = ""),onefile = TRUE, height = 5*ceiling(length(unique(df$Gene))/maxGenes) , width = 30.00)
  #testggplot 
  #dev.off()
}



###removed dollar signs from aes



longpsamp <- function(samp, .....){
  geneSpacing <- 5;
  chromSpacing <- 40;
  pointSpacing <- 40;
  maxGenes <- 7;
  dummyList <- c(1:7);
  dummyGeneList <- unique(df$Gene)
  geneBreaks <- NULL
  df$copynumber <- log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]])
  df$xCoord <- pointSpacing*(geneSpacing*df$GeneNum + df$AmpliconIndex)
  for(i in 1:ceiling(length(unique(df$Gene))/maxGenes)){
    dummyGenes <- dummyGeneList[((i-1)*7) + dummyList]
    for(j in seq_along(dummyGenes))
    {
      df[which(df$Gene == dummyGenes[j]),"Group"] <- i 
      df[which(df$Gene == dummyGenes[j]),"geneEst"] <- log2(geneEst[which(geneEst$Gene == dummyGenes[j]), samp])
      geneBreaks <- c(geneBreaks,mean(df[which(df$Gene == dummyGenes[j]),"xCoord"]))
    }
  }
  testggplot <- ggplot(data = df, aes(x = xCoord, y = copynumber, colour = factor(Color), group = factor(Gene)))
  testggplot <- testggplot + geom_point(size = 0.75, pch = 1) + theme_bw() + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + facet_wrap(~ Group, scales = "free_x", ncol = 2)
  testggplot <- testggplot + geom_line(aes(x = xCoord, y = geneEst), colour = "black", size = 0.05)
  testggplot <- testggplot + ylab(label = "log2(copyNumber)") + xlab(label = "") 
  testggplot <- testggplot + geom_text(data = df[df$AmpliconId %in% outlierList,], aes(xCoord, copynumber,label = AmpliconId, color = "#000000"), size = 1)
  geneLabels = unique(df$Gene) 
  geneBreaks <- geneBreaks[1:length(geneLabels)]
  chromLabels <- unique(df$ChromNum)
  chromBreaks <- NULL
  posAblines <- NULL
  xpos <- NULL
  chromDat <- NULL
  for(i in seq_along(unique(df$ChromNum))){
    xpos <- mean(df[which(df$ChromNum == chromLabels[i]),"xCoord"])
    xpos <- c(xpos, unique(df[which(df$ChromNum == chromLabels[i]),"Group"]))
    chromBreaks <- rbind(chromBreaks, xpos)
    start <- min(df[which(df$ChromNum == chromLabels[i]),"xCoord"]) - 2*chromSpacing
    start <- c(start, unique(df[which(df$ChromNum == chromLabels[i]),"Group"]))
    posAblines <- rbind(posAblines, start)
  }
  chromDat <- data.frame("chromLabels" = chromLabels, "xpos" = chromBreaks[,1], "y" = rep(round(min(df$copynumber)),length(chromBreaks)/2), "Group" = factor(chromBreaks[,2]))
  dummy1 <- data.frame(x = posAblines[,1], Group = factor(posAblines[,2]))
  testggplot <- testggplot + scale_x_continuous(breaks = geneBreaks, labels = geneLabels) + theme(axis.text.x=element_text(angle = 90, hjust = 0))
  testggplot <- testggplot + geom_vline(data = dummy1, aes(xintercept = x), colour = "gray")
  testggplot <- testggplot + geom_text(data = chromDat, aes(x = xpos, y = y, label = chromLabels), inherit.aes = FALSE)
  print(testggplot)
  #pdf(file = paste("long.",samp, ".pdf",sep = ""),onefile = TRUE, height = 5*ceiling(length(unique(df$Gene))/maxGenes) , width = 30.00)
  #testggplot 
  #dev.off()
}




