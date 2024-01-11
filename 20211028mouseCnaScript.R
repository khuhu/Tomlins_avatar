library(stringr)
library(pracma)

### 20211028 pracma and numpy_digitize added & pool bed
numpy_digitize <- function(vec, breaks){
  res <- rep(0, length(vec))
  for (i in seq_along(breaks)) {
    if (i < length(breaks)) {
      res[which(vec >= breaks[i] & vec < breaks[i + 1])] <- i
    } else {
      res[which(vec == breaks[i])] <- i
    }
  }
  return(res)
}

poolBeds <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_pools.bed", header = TRUE, 
                       stringsAsFactors = FALSE, sep = "\t")


args <- commandArgs(trailingOnly=TRUE);


### comment out after testing using test args:

# args <- c("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349_median/amplicon.GCinput.txt",
#           "/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349_median/sampleInfo.input.txt",
#           "/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349_median/amplicon.combinedCoverage.input.txt",
#           "--min-amplicons-per-gene=3")


showQ <- TRUE;
showZ <- FALSE;
colorSigGenes <- TRUE;
colorSigGenesMaxQ <- 0.01;
labelSigGenesOnly <- FALSE;
minAmpliconsPerGene <- 4;
#noisyAmpliconQuantile <- 0.05;

args <- as.list(args);

if (length(args)==0) {
  args <- c("dummy");
}

parseFlag <- function(prefix, flag, castFn) {
  flag <- as.character(flag);
  m <- regmatches(flag, regexec(paste(prefix, "=(.+)", sep=""), flag))[[1]];
  ret <- if (length(m)<2) { NA } else { castFn(m[2]) };
  if (is.na(ret)) {
    cat("# Couldn't parse flag: ", flag, "\n");
    quit();
  }
  return(ret);
}

j <- 0;
for (i in 1:length(args)) {
  if (pmatch("--color-sig-genes=", args[i-j], nomatch=FALSE)) {
    val <- parseFlag("--color-sig-genes", args[i-j], as.logical);
    if (!is.na(val)) colorSigGenes <- val;
    args[i-j] <- NULL;
    j <- j+1;
  } else if (pmatch("--show-z=", args[i-j], nomatch=FALSE)) {
    val <- parseFlag("--show-z", args[i-j], as.logical);
    if (!is.na(val)) showZ <- val;
    args[i-j] <- NULL;
    j <- j+1;
  } else if (pmatch("--show-q=", args[i-j], nomatch=FALSE)) {
    val <- parseFlag("--show-q", args[i-j], as.logical);
    if (!is.na(val)) showQ <- val;
    args[i-j] <- NULL;
    j <- j+1;
  } else if (pmatch("--label-sig-genes-only=", args[i-j], nomatch=FALSE)) {
    val <- parseFlag("--label-sig-genes-only", args[i-j], as.logical);
    if (!is.na(val)) labelSigGenesOnly <- val;
    args[i-j] <- NULL;
    j <- j+1;
  } else if (pmatch("--color-sig-genes-max-q=", args[i-j], nomatch=FALSE)) {
    val <- parseFlag("--color-sig-genes-max-q", args[i-j], as.double);
    if (!is.na(val)) colorSigGenesMaxQ <- val;
    args[i-j] <- NULL;
    j <- j+1;
  } else if (pmatch("--min-amplicons-per-gene=", args[i-j], nomatch=FALSE)) {
    val <- parseFlag("--min-amplicons-per-gene", args[i-j], as.integer);
    if (!is.na(val)) minAmpliconsPerGene <- val;
    args[i-j] <- NULL;
    j <- j+1;
  } else if (pmatch("--noisy-amplicon-quantile=", args[i-j], nomatch=FALSE)) {
    val <- parseFlag("--noisy-amplicon-quantile", args[i-j], as.double);
    if (!is.na(val)) noisyAmpliconQuantile <- val;
    args[i-j] <- NULL;
    j <- j+1;
  } else if (args[i-j]!="--" && pmatch("--", args[i-j], nomatch=FALSE)) {
    cat("# Skipping unknown flag:", args[[i-j]], "\n");
    args[i-j] <- NULL;
    j <- j+1;
  }
}

args <- unlist(args);

if (length(args) < 3) {
  cat("\n");
  cat("Usage: [Rscript] AMPLICON_INFO SAMPLE_INFO READ_COUNTS [SNP_FILE] [OPTIONS]\n");
  cat("Supported options are:\n");
  cat("  --show-q=T/F           Add Q-values to gene labels? (default is true)\n");
  cat("  --show-z=T/F           Add Z-scores to gene labels? (default is false)\n");
  cat("  --color-sig-genes=T/F  Colored gene labels for sig. gains and losses (default is true)\n");
  cat("  --color-sig-genes-max-q=VAL    Q-value cutoff for coloring (default is 0.01)\n");
  cat("  --label-sig-genes-only=T/F     Show Q & Z only for sig. gains/losses (default is false)\n");
  cat("  --noisy-amplicon-quantile=VAL  Fraction of lowest-coverage amplicons to drop (default is 0.05)\n");
  cat("  --min-amplicons-per-gene=VAL   Minimum number of amplicons in analyzed genes (default 4)\n");
  cat("\n");
  quit();
}

ampliconInfoFilename <- args[1];
sampleInfoFilename <- args[2];
readCountsFilename <- args[3];
variantCallsFilename <- if (length(args) >= 4) args[4] else "(none)";

cat("## Color significant genes?:", colorSigGenes, "\n");
if (colorSigGenes) {
  cat("## Color significant genes with Q <= ", colorSigGenesMaxQ, "\n");
}
cat("## Label w/ Q-value?:", showQ, "\n");
cat("## Label w/ Z-score?:", showZ, "\n");
cat("## Label sig. genes only?:", labelSigGenesOnly, "\n");
cat("## Min. amplicons per gene:", minAmpliconsPerGene, "\n");
cat("## Amplicon info filename:", ampliconInfoFilename, "\n");
cat("## Sample info filename:", sampleInfoFilename, "\n");
cat("## Read counts filename:", readCountsFilename, "\n");
cat("## Variant calls filename:", variantCallsFilename, "\n");
cat("\n");

getGeneInfo <- function(df, cols=c("red", "green", "orange", "blue", "purple")) {
  genes <- unique(df$Gene);
  ## Fold gene information into data matrix
  geneInfo <- data.frame(Gene=genes, 
                         MinIndex = sapply(genes, function(g) min(df$AmpliconIndex[df$Gene==g])),
                         MaxIndex = sapply(genes, function(g) max(df$AmpliconIndex[df$Gene==g])),
                         ChromNum = as.numeric(sapply(genes, function(g) min(df$ChromNum[df$Gene==g]))),
                         NumProbes = sapply(genes, function(g) sum(1*(df$Gene==g))),
                         Label = rep("", length(genes)),
                         GeneNum = 1:length(genes),
                         Color = rep("red", length(genes)));
  rownames(geneInfo) <- geneInfo$Gene;
  geneInfo$Label <- paste(geneInfo$Gene, " (", geneInfo$NumProbes, ")", sep="");
  geneInfo$Color <- rep(cols, length(genes))[1:length(genes)];
  return(geneInfo);
}


## Read mapping from SNP sample names to read-count sample names
cat("Reading sample info...");
bcNameMapping <- read.table(sampleInfoFilename, header=TRUE, sep="\t", stringsAsFactors = FALSE);
cat("...read information on", nrow(bcNameMapping), "samples.\n");

## Read in the amplicon information
cat("Reading amplicon info...");
ampliconInfo <- read.table(ampliconInfoFilename, header=TRUE, sep="\t",stringsAsFactors=FALSE);
cat("...read information on", nrow(ampliconInfo), "amplicons.\n");
rownames(ampliconInfo) <- ampliconInfo$AmpliconIndex;

## Read in the matrix of (replicate-pooled) read counts
cat("Reading read counts...");
dfRaw <- read.table(readCountsFilename, header=TRUE, sep=",", stringsAsFactors = FALSE, check.names = FALSE);
if (length(which(duplicated(colnames(dfRaw)))) > 1) {
  dfRaw <- dfRaw[,-which(duplicated(colnames(dfRaw)))]
}


cat("...read table with", nrow(dfRaw), "rows.\n");

tumorNames <- as.vector(bcNameMapping$Sample[bcNameMapping$SampleClass == "Tumor" & bcNameMapping$Sample %in% colnames(dfRaw)]);
normalNames <- as.vector(bcNameMapping$Sample[bcNameMapping$SampleClass == "Normal" & bcNameMapping$Sample %in% colnames(dfRaw)]);
cat("Tumor names (", length(tumorNames), "):", tumorNames, "\n");
cat("Normal names (", length(normalNames), "):", normalNames, "\n");

if (length(which(duplicated(normalNames))) > 1) {
  normalNames <- normalNames[-which(duplicated(normalNames))]
}

allNames <- append(tumorNames, normalNames);


### gets rid of bad samples, and almost zero count samples
badSamps1 <- colnames(dfRaw)[which(is.na(dfRaw[1,]))]
if (length(badSamps1 > 0)) {
  dfRaw <- dfRaw[,-which(is.na(dfRaw[1,]))]
}

libCountCheck <- apply(dfRaw[,2:ncol(dfRaw)], 2, sum)
tmpStat <- which(libCountCheck/median(libCountCheck) < 0.001)

badSamps2 <- colnames(dfRaw)[tmpStat + 1]
if (length(badSamps2 > 0)) {
  dfRaw <- dfRaw[,-(tmpStat + 1)]
}

badSamps3 <- names(which(apply(dfRaw[,2:ncol(dfRaw)], 2, median) == 0))

if(length(c(badSamps1, badSamps2, badSamps3)) > 0){
  allNames <- allNames[-which(allNames %in% c(badSamps1, badSamps2, badSamps3))]
}

## Merge the amplicon information into the read count table
df <- dfRaw;
df <- merge(df, ampliconInfo, by="AmpliconId")
df <- df[order(df$AmpliconIndex),];


poolNames <- normalNames;
if (length(poolNames)==0) {
  cat("!!!!  No normal samples found.  Using median of tumor samples as baseline.\n");
  poolNames <- tumorNames;
  normalPool <- FALSE;  
} else {
  normalPool <- TRUE;
}




## Create geneInfo table
cat("Creating gene information table...");
genes <- unique(df$Gene);
geneInfo <- getGeneInfo(df);
cat("...analyzing", nrow(geneInfo), "genes.\n");

## Merge geneInfo into the read count table
df <- merge(df, geneInfo, by=c("Gene", "ChromNum"));
df <- df[order(df$AmpliconIndex),];
df$ChromNum <- as.numeric(df$ChromNum)
rownames(df) <- df$AmpliconId;


linBreaksGc <- linspace(min(ampliconInfo$GC), max(ampliconInfo$GC), 5);
linBreaksLen <- linspace(min(ampliconInfo$Length), max(ampliconInfo$Length), 5);
ampliconInfo$gcbin <- numpy_digitize(ampliconInfo$GC, linBreaksGc);
ampliconInfo$lenbin <- numpy_digitize(ampliconInfo$Length, linBreaksLen);
ampliconInfo$fiveByFive <- paste0(ampliconInfo$gcbin, ampliconInfo$lenbin);

##  new normalization methods
for (i in allNames) {
  # print(i)
  tmpVector <- df[[i]];
  tmpVector[which(poolBeds$pool == 1)] <- tmpVector[which(poolBeds$pool == 1)]/median(tmpVector[which(poolBeds$pool == 1)]);
  tmpVector[which(poolBeds$pool == 2)] <- tmpVector[which(poolBeds$pool == 2)]/median(tmpVector[which(poolBeds$pool == 2)]);
  for (j in unique(ampliconInfo$fiveByFive)) {
    # print(j)
    ### for some fringe amplicons on good samples created NAs or NaNs, bad samples like 
    ### 6943ROT_X63 in 142 will probabl turn out bad and be obvious - filt downstream
    if (median(tmpVector[which(ampliconInfo$fiveByFive == j)]) == 0 | is.na(median(tmpVector[which(ampliconInfo$fiveByFive == j)]))) {
      next()
    } else{
      tmpVector[which(ampliconInfo$fiveByFive == j)] <- tmpVector[which(ampliconInfo$fiveByFive == j)]/median(tmpVector[which(ampliconInfo$fiveByFive == j)]);
    }
  }
  df[[i]] <- tmpVector;
}


if (length(poolNames)>1) {
  normalAmpMedian <- apply(df[, poolNames], 1, median);
} else {
  normalAmpMedian <- df[[poolNames[1]]];
}


df <- df[-which(normalAmpMedian == 0),]
normalAmpMedian <- normalAmpMedian[-which(normalAmpMedian == 0)]
df2 <- df;
normalAmpMedian2 <- normalAmpMedian;

keepGenes <- geneInfo[geneInfo$NumProbes >= minAmpliconsPerGene,]$Gene;
normalAmpMedian <- normalAmpMedian[df$Gene %in% keepGenes]
df <- df[df$Gene %in% keepGenes,];

dfCorrected <- df;
dfCorrected2 <- df2;


geneEst <- data.frame(Gene=unique(df$Gene));
rownames(geneEst) <- unique(df$Gene);

### the correction itself in df sometimes gives NA values ..

for (x in allNames) {
  
  geneEst[[x]] <- rep(0, length(unique(df$Gene)));
  dfCorrected[[x]] <- dfCorrected[[x]]/normalAmpMedian;
  dfCorrected2[[x]] <- dfCorrected2[[x]]/normalAmpMedian2;
  
  for (gene in geneEst$Gene) {
    ix <- (df$Gene==gene);
    geneEst[gene, x] <- median(dfCorrected[ix, x]);
  }
}


psamp <- function(samp, ylim=NULL, cex.axis=0.25, ax=TRUE, main=NA, analysisResults=NULL) {
  geneSpacing <- 80;
  chromSpacing <- 0;
  if (is.na(main)) main=samp;
  
  xlim <- c(0, nrow(ampliconInfo) + 1 + geneSpacing*(1 + nrow(geneInfo)) + 23*chromSpacing);
  
  ## plot(y=log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]]), x=(geneSpacing*df$GeneNum + chromSpacing*(df$ChromNum-1) + df$AmpliconIndex), col=df$Color, cex=20.0*sqrt(df$Weights), xaxt="n", main=main, ylab="Log2(CN Ratio)", xlim=NULL, xlab="", ylim=ylim)
  
  ## c(0,nrow(ampliconInfo)+geneSpacing*nrow(geneInfo)+chromSpacing*23), xlab="", ylim=ylim);
  
  plot(y=log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]]), x=(geneSpacing*df$GeneNum + chromSpacing*(as.numeric(df$ChromNum)-1) + df$AmpliconIndex), col=df$Color, cex=0.3, xaxt="n", main=main, ylab="Log2(CN Ratio)", xlim=xlim, xlab="", ylim=ylim, xaxs="i", cex.axis=0.9, cex.lab=0.9);
  
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
      if (showZ && showQ) { cex.axis <- 1.0*cex.axis; }
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




dirToWrite <- str_remove(ampliconInfoFilename, "amplicon.GCinput.txt")
#print(dirToWrite)
combinedTable <- NULL
genes <- unique(df$Gene)
for (samp in allNames) {
      analysisResults <- data.frame(Gene=genes, Sample=rep(samp, length(genes)), NumProbes=geneInfo$NumProbes[which(geneInfo$Gene %in% keepGenes)], CopyNumberRatio=geneEst[,samp],
                                    ZScore=rep(0.0, length(genes)), Call=rep(0, length(genes)), Sig=rep(0, length(genes)), Log10PValue=rep(0.0, length(genes)));
  rownames(analysisResults) <- genes;
  sampleList <- c(samp, normalNames)
  for (gene in genes) {
    row <- analysisResults[gene,];
    est <- row$CopyNumberRatio;
    z <- 0.6745 * (est - apply(geneEst[gene, sampleList], 1, median))/apply(geneEst[gene, sampleList], 1, stats::mad);
   
    analysisResults[gene, "ZScore"] <- z;
    analysisResults[gene, "Log10PValue"] <- if (is.na(z)) { 0.0 } else { log10(2.0*pnorm(-abs(z))) };
    analysisResults[gene, "Call"] <- if (est > 1.0) "GAIN" else { if (est < 1.0) "LOSS" else "NONE" };
    analysisResults[gene, "Sig"] <- if (is.na(z)) "(??)" else { if (z >= 2 && z < 3) "(*)" else { if (z>=3 && z<5) "(**)" else { if (z>=5) "(***)" else "()" } } };
  }
  analysisResults$PValueRank <- analysisResults$Log10PValue;
  tmp <- order(analysisResults$PValue);
  for (i in 1:length(genes)) {
    analysisResults[tmp[i], "PValueRank"] <- i;
  }
  analysisResults <- analysisResults[tmp,];
  analysisResults$Log10QValue <- analysisResults$Log10PValue;
  for (i in 1:length(genes)) {
    # Benjamini-Hochberg FDR equation
    analysisResults[i, "Log10QValue"] <- analysisResults[i, "Log10PValue"] + log10(length(genes) / analysisResults[i, "PValueRank"]);
  }
  for (i in 1:length(genes)) {
    # monotonicity correction
    analysisResults[i, "Log10QValue"] = min(analysisResults[i:length(genes), "Log10QValue"]);
  }
  setwd(dirToWrite)
  pdf(paste("out.", samp, ".pdf", sep=""), width=7 * 2, height=5 * 2, useDingbats=FALSE);
  dev.off();
  write.table(format(analysisResults, scientific=FALSE, digits=4), file=paste("calls.", samp, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE);
  combinedTable <- rbind(combinedTable, format(analysisResults, scientific=FALSE, digits=4))
}



write.table(combinedTable, file="combinedCalls.txt", quote=FALSE, sep="\t", row.names=FALSE);
write.table(geneEst, file="cnMatrix_gene.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names = TRUE);
write.table(dfCorrected2, file = "cnAmplicon_matrix.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names = TRUE);
write.table(df2, file = "gcCorrectedCounts_matrix.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names = TRUE);


warnings();
