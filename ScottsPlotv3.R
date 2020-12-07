###editing long samp

##creating second color scheme grouping for per chromosome

Color2 <- NULL
dfColorList <- c("firebrick1","darkolivegreen3","goldenrod2","dodgerblue4","darkorchid4")
k <- 0;
for(i in seq_along(unique(df$Gene))){
  k <- k + 1;
  Color2[which(df$Gene == unique(df$Gene)[i])] <- dfColorList[k]
  if(k == 5){
    k <- 0;
  }
}

df$Color2 <- Color2


newpsamp <- function(samp, .....){
  #samp <- allNames[1]
  geneSpacing <- 50;
  geneBreaks <- NULL
  df$copynumber <- log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]])
  df$xCoord <- (df$AmpliconIndex + df$GeneNum *geneSpacing)
  #below was used to previously create gene names
  for(i in seq_along(unique(df$Gene))){
    a <- df[which(df$Gene == unique(df$Gene)[i]), c("xCoord")]
    gPos <- NULL
    gPos <- mean(min(a), max(a))
    geneBreaks <- c(geneBreaks, gPos)
  }
  #geneBreaks <- as.data.frame(geneBreaks)
  
  for(i in seq_along(df$Gene)){
    df$geneEst[i] <- log2(geneEst[[samp]][which(geneEst$Gene == df$Gene[i])])
  }
  
  
  testggplot <- ggplot(data = df) + geom_point(data = df, aes(x = xCoord, y = copynumber),color = df$Color2,size = 1.5, alpha=0.5) +
    theme_bw() + theme(legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank()) +
    geom_line(data = df,aes(x=xCoord, y=geneEst, group=factor(Gene)), color="black", size=1.5, alpha=1.0 , inherit.aes = FALSE) +
    ylab(label = "log2(copyNumber)") + xlab(label = "") + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2), axis.ticks = element_line(colour = "black", size = 2))
  geneLabels = unique(df$Gene)
  chromLabels <- unique(df$ChromNum)
  chromLabPos <- NULL
  chromDat <- NULL
  chromLine <- NULL
  for(i in seq_along(chromLabels)){
    xmin <- min(df[which(df$ChromNum == chromLabels[i]),"xCoord"])
    xmax <- max(df[which(df$ChromNum == chromLabels[i]),"xCoord"])
    chromLabPos <- rbind(chromLabPos,  mean(c(xmin, xmax)))
    linepos <- xmax + 30
    if(i == length(chromLabels)){
      linepos <- NULL
    }
    chromLine <- c(chromLine, linepos)
  }
  chromLine <- as.data.frame(chromLine)
  chromLabels <- gsub("23","X", chromLabels)
  
  ### switched chrom data because andi's data doesn't work perferctly for code below
  #chromDat <- data.frame("chromLabels" = chromLabels, "y" = rep(round(max(df$copynumber)),length(chromLabels)), "x" = chromLabPos, "col"= rep(c("white","black"),round(length(chromLabels)/2)))
  chromDat <- data.frame("chromLabels" = chromLabels, "y" = rep(3.25,length(chromLabels)), "x" = chromLabPos, "col"= c(rep(c("white","black"),10), "white"))
  
  
  ###
  panelSqs <- NULL
  x1 <- NULL
  x2 <- NULL
  for(i in seq_along((chromLabels))){
    if(i < length(chromLabels)){
      xminn <- min(df[which(df$ChromNum == chromLabels[i]),"xCoord"]) - 20
      xmaxx <- chromLine[i,"chromLine"]
    }
    else{
      xminn <- max(df[which(df$ChromNum == chromLabels[i -1 ]),"xCoord"]) + 20
      xmaxx <- Inf
    }
    if(i == 1){
      xminn <- min(df[which(df$ChromNum == chromLabels[i]),"xCoord"]) - 50
      xmaxx <- chromLine[i,"chromLine"]
    }
    x1 <- c(x1, xminn)
    x2 <- c(x2, xmaxx)
  }
  
  ### same with chromDat but had to switch it up for panelSqs since Andi's number of chromosomes are uneven ... need to find a more robust way of doing these things
  
  #panelSqs <- data.frame("xminn"=x1,"xmax"=x2,"yminn"=c(0.90*round(max(df$copynumber))),
  #                       "ymaxx"=c(Inf), "col" = rep(c("black","white"),round(length(x1)/2)))
  
  panelSqs <- data.frame("xminn"=x1,"xmax"=x2,"col" = c(rep(c("black","white"),round(length(x1)/2)),"black"))
  
  #yminn =c(0.90*round(max(df$copynumber)))
  yminn = 3
  ymaxx =c(Inf)
  ###
  testggplot <- testggplot + geom_vline(data = chromLine, aes(xintercept = chromLine), colour = "gray", linetype = 2) + geom_hline(yintercept = 0, size = 1.15)
  ###gets rid of the space between borders and first and last points
  #testggplot <- testggplot + scale_x_continuous(limits = c(min(df[,"xCoord"]) - 2*geneSpacing,max(df[,"xCoord"])) + geneSpacing, expand = c(0,0))
  testggplot <- testggplot + geom_rect(data = panelSqs, aes(xmin=xminn, xmax=xmaxx, ymin=yminn, ymax=ymaxx), fill = panelSqs$col, inherit.aes = FALSE) +
    geom_text(data = chromDat, aes(x = x, y = y, label = chromLabels), inherit.aes = FALSE, size = 5, color = chromDat$col)
  testggplot <- testggplot + geom_hline(yintercept = yminn, size = 1.5) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 10))
  
  testggplot <- testggplot + scale_x_continuous(breaks = geneBreaks, labels = geneLabels,limits = c(min(df[,"xCoord"]) - 2*geneSpacing,max(df[,"xCoord"])) + geneSpacing, expand = c(0,0)) + theme(axis.text.x=element_text(angle = 90, hjust = 0.2)) +
    ylim(c(-3,3.5))
  
  geneBreaks <- unlist(sapply(geneBreaks, round))
  print(testggplot)
  #pdf(file = paste("long.",samp, ".pdf",sep = ""),onefile = TRUE, height = 0.5*ceiling(length(unique(df$Gene))/maxGenes) , width = 30.00)
  #testggplot 
  #dev.off()
}


###For facetsamp, need to rescale such that every new chrom starts with xcoord 0. This might screw up the gene est lines
###need to redo it by gene spacing too


AmpliconIndexZeroNormed <- NULL
GeneNumZeroNormed <- NULL
for(i in seq_along(unique(df$ChromNum))){
  AmpliconIndexZeroNormed[which(df$ChromNum == unique(df$ChromNum)[i])] <- df[which(df$ChromNum == unique(df$ChromNum)[i]), "AmpliconIndex"] - min(df[which(df$ChromNum == unique(df$ChromNum)[i]), "AmpliconIndex"]) + 1
  GeneNumZeroNormed[which(df$ChromNum == unique(df$ChromNum)[i])] <- df[which(df$ChromNum == unique(df$ChromNum)[i]), "GeneNum"] - min(df[which(df$ChromNum == unique(df$ChromNum)[i]), "GeneNum"]) + 1
}


df$AmpliconIndex2 <- AmpliconIndexZeroNormed
df$GeneNum2 <- GeneNumZeroNormed

facetsamp <- function(samp, .....){
  geneSpacing <- 20;
  chromSpacing <- 200;
  #dummyList <- c(1:7);
  #dummyGeneList <- unique(df$Gene)
  geneBreaks <- NULL
  df$copynumber <- log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]])
  df$xCoord2 <- (df$AmpliconIndex2 + (df$GeneNum2 -1) *geneSpacing + chromSpacing)
  for(i in seq_along(unique(df$Gene))){
    gPos <- NULL
    gPos <- mean(df[which(df$Gene == unique(df$Gene)[i]),]$xCoord2)
    df[which(df$Gene == unique(df$Gene)[i]),"geneEst"] <- log2(geneEst[which(geneEst$Gene ==unique(df$Gene)[i]), samp])
    geneBreaks <- c(geneBreaks, gPos)
  }
  geneBreaks <- as.data.frame(geneBreaks)
  testggplot <- ggplot(data = df, aes(x = xCoord2, y = copynumber))
  #testggplot <- testggplot + geom_point(size = 0.75, pch = 1) + theme_bw() + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + facet_wrap(~ ChromNum, scales = "free_x", nrow = 1)
  testggplot <- testggplot + geom_point(data = df, aes(color=factor(Color2)),size = 1.5)
  testggplot <- testggplot + theme_bw() + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + facet_grid(~ ChromNum, scales = "free_x", space = "free_x")
  testggplot <- testggplot + theme(panel.spacing = unit(0,"points"), panel.border = element_rect(linetype = "dashed", fill = NA)) 
  testggplot <- testggplot + geom_line(data = df, aes(x = xCoord2, y = geneEst, group = factor(Gene)), color = "black", size = 0.05, inherit.aes = FALSE)
  testggplot <- testggplot + ylab(label = "log2(copyNumber)") + xlab(label = "") + scale_x_continuous(breaks = NULL)
  geneLabels = unique(df$Gene)
  chromLabels <- unique(df$ChromNum)
  chromLabPos <- NULL
  chromDat <- NULL
  for(i in seq_along(unique(df$ChromNum))){
    xcurr <- max(df[which(df$ChromNum == chromLabels[i]),"xCoord2"])
    chromLabPos <- rbind(chromLabPos, xcurr)
  }
  chromLabPos <- as.data.frame(chromLabPos)
  ###need to calculate the chromBreaks similar to how I did gene Breaks
  ChromDat <- data.frame("chromLabels" = chromLabels, "x" = chromLabPos)
  testggplot <- testggplot + geom_hline(yintercept = 0)
  print(testggplot)
}

