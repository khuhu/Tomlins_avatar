
###editing psamp so that it makes graph that plots amplicons that tile across gene in respective genomci order

psamp.ggplot <- function(samp, ylim=NULL, cex.axis=0.75, ax=TRUE, main=NA, analysisResults=NULL) {
  geneSpacing <- 40;
  chromSpacing <- 0;
  if (is.na(main)) main=samp;
  
  xlim <- c(0, nrow(ampliconInfo) + 1 + geneSpacing*(1 + nrow(geneInfo)) + 23*chromSpacing);
  
  ## plot(y=log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]]), x=(geneSpacing*df$GeneNum + chromSpacing*(df$ChromNum-1) + df$AmpliconIndex), col=df$Color, cex=20.0*sqrt(df$Weights), xaxt="n", main=main, ylab="Log2(CN Ratio)", xlim=NULL, xlab="", ylim=ylim)
  
  ## c(0,nrow(ampliconInfo)+geneSpacing*nrow(geneInfo)+chromSpacing*23), xlab="", ylim=ylim);
  
  plot(y=log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]]), x=(geneSpacing*df$GeneNum + chromSpacing*(df$ChromNum-1) + df$AmpliconIndex), col=df$Color, cex=0.3, xaxt="n", main=main, ylab="Log2(CN Ratio)", xlim=xlim, xlab="", ylim=ylim, xaxs="i", cex.axis=0.9, cex.lab=0.9);
  
  for (gene in geneInfo$Gene) {
    lines(x=geneInfo[gene, c("MinIndex", "MaxIndex")] + geneSpacing*c(-0.5,0.5) + geneSpacing*geneInfo[gene, "GeneNum"] + chromSpacing*geneInfo[gene, "ChromNum"], y=rep(log2(geneEst[gene, samp]),2));
  }
  
  if (ax && cex.axis > 0) {
    if (!is.null(analysisResults)) {
      vals <- analysisResults[match(geneInfo$Gene, analysisResults$Gene),c("Log10QValue", "ZScore", "Call")];
      geneAxisLabels <- geneInfo$Gene;
      cols <- rep("", length(vals$Call));
      for (i in 1:length(vals$Call)) {
        v <- (1.0*vals[i,"Log10QValue"]);
        call <- vals[i,"Call"];
        ## vv[i] <- if (v <= -3.0) "(<0.001)" else { if (v <= -2.0) "(0.001-0.01)" else { if (v <= -1.0) "(0.01-0.1)" else "" } };
        cols[i] <- if (!colorSigGenes || (10.0**v) > colorSigGenesMaxQ) "black" else { if (call=="GAIN") "red" else "blue" };
      }
      if (showQ) {
        tmp <- paste("(Q=", format(10.0**(pmax(vals$Log10QValue, -16.00)), scientific=TRUE, digits=2), ")", sep="");
        if (labelSigGenesOnly) { tmp[cols=="black"] <- ""; }
        geneAxisLabels <- paste(geneAxisLabels, tmp);
      }
      if (showZ) {
        tmp <- paste("(Z=", format(vals$ZScore, scientific=FALSE, digits=2), ")", sep="");
        if (labelSigGenesOnly) { tmp[cols=="black"] <- ""; }
        geneAxisLabels <- paste(geneAxisLabels, tmp);
      }
      if (showZ && showQ) { cex.axis <- 0.8*cex.axis; }
    }
    labelPositions = 0.5*(geneInfo$MinIndex+geneInfo$MaxIndex)+geneSpacing*geneInfo$GeneNum+chromSpacing*(geneInfo$ChromNum-1);
    Map(function(x,y,z) axis(1, at=x, col.axis=y, labels=z, lwd=0, cex.axis=cex.axis, mgp=c(1,0.5,0), las=2), labelPositions, cols, geneAxisLabels);
    # axis(1, at=labelPositions, labels=geneAxisLabels, las=2, cex.axis=cex.axis, mgp=c(1,0.5,0), tick=FALSE);
  }
  
  chromMin <- rep(0, 24);
  chromMax <- rep(0, 24);
  chromProbes <- rep(0, 24);
  for (i in 1:24) {
    dfChrom <- df[df$ChromNum==i,];
    chromProbes[i] <- length(dfChrom$GeneNum);
    if (chromProbes[i] > 0) {
      chromMax[i] <- max(geneSpacing*dfChrom$GeneNum + chromSpacing*(dfChrom$ChromNum-1) + dfChrom$AmpliconIndex);
      chromMin[i] <- min(geneSpacing*dfChrom$GeneNum + chromSpacing*(dfChrom$ChromNum-1) + dfChrom$AmpliconIndex);
    }
  }
  
  for (i in 1:23) {
    if (chromProbes[i] > 0) {
      if (i>1) {
        lines(x=rep(chromMin[i] - 0.5*geneSpacing - 0.5*chromSpacing - 0.5, 2), y=c(-10,10), col="gray80");
      }
      if (chromProbes[i+1] == 0 && i<20) {
        lines(x=rep(chromMax[i] + 0.5*geneSpacing + 0.5*chromSpacing + 0.5, 2), y=c(-10,10), col="gray80");
      }
      text(x=(chromMin[i]+chromMax[i])*0.5, y=ylim[1], col="black", labels=paste(i, sep=""), cex=0.5, srt=90);
    }
  }
}



###testing bits and pieces
samp <- allNames[1]
geneSpacing <- 5;
chromSpacing <- 40;
pointSpacing <- 1500;
ylim <- NULL
main <- ""
xlim <- NULL;
plot(y=log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]]),
     x=(geneSpacing*df$GeneNum + chromSpacing*(df$ChromNum-1) + df$AmpliconIndex), col=df$Color,
     cex=0.3, xaxt="n", main=main, ylab="Log2(CN Ratio)", xlim=xlim, xlab="", ylim=ylim, xaxs="i", cex.axis=0.9, cex.lab=0.9);



pointSpacing <- 40;
par(mfrow=c(5,1))
for(i in 1:ceiling(length(unique(df$Gene))/maxGenes)){
  dummyGenes <- dummyGeneList[((i-1)*7) + dummyList]
  dummyDf <- df[which(df$Gene == dummyGenes),]
  plot(y = dummyDf$copynumber,
     x=(pointSpacing*(geneSpacing*dummyDf$GeneNum + dummyDf$AmpliconIndex)), col=df$Color,
     cex=0.5, xaxt="n", main=main, ylab="Log2(CN Ratio)", xlim=xlim, xlab="", ylim=ylim, xaxs="i", cex.axis=0.9, cex.lab=0.9, pch = 15);
}

### this basic test - copynumber is for samp 
df$copynumber <- log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]])
df.subset <- df[which(df$ChromNum < 5),]
library(ggplot2)

###from the code below: ~ 7 genes per each part of the graph makes sense
pdf(file = "/home/kevhu/data/testplot.pdf", height = 5.00, width = 20.00)

testggplot <- ggplot(data = df.subset, aes(pointSpacing*(geneSpacing*df.subset$GeneNum + df.subset$AmpliconIndex), df.subset$copynumber))
testggplot + geom_point(aes(colour = df.subset$Color), size = 0.75) + geom_line(aes(group = df.subset$Gene, color = df.subset$Color), size = 0.1) + theme_bw()

dev.off()

###now I will try and iteratively create graph containing entire dataset

geneSpacing <- 5;
chromSpacing <- 40;
pointSpacing <- 1500;
maxGenes <- 7;
dummyList <- c(1:7);
dummyGeneList <- unique(df$Gene)
listOfPlots <- NULL
rownames(df) <- NULL

for(i in 1:ceiling(length(unique(df$Gene))/maxGenes)){
  dummyGenes <- dummyGeneList[((i-1)*7) + dummyList]
  dummyDf <- df[which(df$Gene == dummyGenes),]
  #assign(paste("testggplot",i, sep = ""),
  #ggplot(aes(colour = df$Color, group = df$Gene)) + geom_point(aes(x = pointSpacing*(geneSpacing*df[which(df$Gene == dummyGenes),]$GeneNum + df[which(df$Gene == dummyGenes),]$AmpliconIndex), y = df[which(df$Gene == dummyGenes),]$copynumber), size = 0.75) + geom_line(size = 0.1) + theme_bw() + theme(legend.position="none"))
  testggplot <- ggplot(inherit.aes = TRUE,data = dummyDf, aes(x = pointSpacing*(geneSpacing*dummyDf$GeneNum + dummyDf$AmpliconIndex), y = dummyDf$copynumber, colour = factor(dummyDf$Color), group = factor(dummyDf$Gene)))
  testggplot <- testggplot + geom_point(size = 0.75) + geom_line(size = 0.1) + theme_bw() + theme(legend.position="none")
  listOfPlots[[i]] <- testggplot
  #print(testggplot)
}



library(gridExtra)

pdf(file = "/home/kevhu/data/testplot.pdf",onefile = TRUE, height = 2.5*ceiling(length(unique(df$Gene))/maxGenes) , width = 20.00)
grid.arrange(testggplot1,testggplot2)
dev.off()





testggplot <- ggplot(data = df, aes(x = pointSpacing*(geneSpacing*df$GeneNum + df$AmpliconIndex), y = df$copynumber, colour = df$Color, group = df$Gene))
testggplot + geom_point(size = 2.5) + geom_line(size = 1) + theme_bw()

pdf(file = "/home/kevhu/data/testplot.pdf", height = 5.00, width = 10.00)

###idea is to make a new column/s for the x-axis for graph that Scott wants ... my idea will be to try to create new start and end positions based on min start pos.

###for now I think easiest way is to just to put them in order i.e just use the row number b/c they all should be in order



###code below is used to ammend a group, which I will later utilize to do facet plots, and get y-coords for median lines

##### Start of code that works
####
####
samp <- allNames[1]
df$copynumber <- log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]])
df.subset <- df[which(df$ChromNum < 5),]
library(ggplot2)
geneSpacing <- 5;
chromSpacing <- 40;
pointSpacing <- 40;
maxGenes <- 7;
dummyList <- c(1:7);
dummyGeneList <- unique(df$Gene)
#df$xCoord <- pointSpacing*(geneSpacing*df$GeneNum + df$AmpliconIndex)
df$xCoord <- (geneSpacing*df$GeneNum + df$AmpliconIndex)

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

pdf(file = "/home/kevhu/data/testplot.pdf",onefile = TRUE, height = 5*ceiling(length(unique(df$Gene))/maxGenes) , width = 30.00)
testggplot 
dev.off()
















