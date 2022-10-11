library(stringr)
args <- commandArgs(trailingOnly=TRUE);

showQ <- TRUE;
showZ <- FALSE;
colorSigGenes <- TRUE;
colorSigGenesMaxQ <- 0.01;
labelSigGenesOnly <- FALSE;
minAmpliconsPerGene <- 4;
noisyAmpliconQuantile <- 0.05;

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
cat("## Noisy amplicon quantile:", noisyAmpliconQuantile, "\n");
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

## Define a lowess-correction vector
lowessCorrect <- function(dat, gc, alpha=0.25) {
  lo <- lowess(y=dat, x=gc, alpha);
  loVec <- lo$y[match(gc, lo$x)];
  return(loVec);
}

##
## Read all raw data
##

## Read mapping from SNP sample names to read-count sample names
cat("Reading sample info...");
# bcNameMapping <- read.table(file.choose(), header=TRUE, sep="\t");
# 2014/02/26: DHH added 'stringsAsFactors' argument. Was treating gene column as factor, which was kicking 
# 						out an error when in psamp function (plot statement)
#bcNameMapping <- read.table(sampleInfoFilename, header=TRUE, sep="\t",stringsAsFactors=FALSE);
bcNameMapping <- read.table(sampleInfoFilename, header=TRUE, sep="\t", stringsAsFactors = FALSE);
cat("...read information on", nrow(bcNameMapping), "samples.\n");

## Read in the amplicon information
cat("Reading amplicon info...");
# ampliconInfo <- read.table(file.choose(), header=TRUE, sep="\t");
# 2014/02/26: DHH added 'stringsAsFactors' argument. Was treating gene column as factor, which was kicking 
# 						out an error when creating initial gene table
ampliconInfo <- read.table(ampliconInfoFilename, header=TRUE, sep="\t",stringsAsFactors=FALSE);
cat("...read information on", nrow(ampliconInfo), "amplicons.\n");
rownames(ampliconInfo) <- ampliconInfo$AmpliconIndex;

## Read in the matrix of (replicate-pooled) read counts
cat("Reading read counts...");
# dfRaw <- read.table(file.choose(), header=TRUE, sep=",");
# 2014/02/26: DHH added 'stringsAsFactors' argument. Was treating gene column as factor, which was kicking 
# 						out an error when in psamp function (plot statement)
#dfRaw <- read.table(readCountsFilename, header=TRUE, sep=",",stringsAsFactors=FALSE);
#2017/08/31: Kevinaded check names portion, b/c colnames got weird 
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

## Read in variant data
if (variantCallsFilename != "(none)") {
  cat("Reading variant calls...");
  # bcRaw <- read.table(file.choose(), header=TRUE, sep="\t");
  bcRaw <- read.table(variantCallsFilename, header=TRUE, sep="\t");
  cat("...read table with", nrow(bcRaw), "rows.\n");
  bc <- merge(bcRaw, bcNameMapping, by="DNA");
  variantData <- TRUE;
  cat("...reduced to", nrow(bc), "rows after restricting samples.\n");
} else {
  variantData <- FALSE;
  cat("No variant calls read.\n");
}

##
## Process the variant information
##

if (variantData) {
  
  ## Create table of distinct variants
  variants <- unique(bc[,c("Chrom", "Position", "Type", "Ref", "Variant", "AmpliconId")]);
  variants$VariantId <- 1:nrow(variants);
  variants$SymbolType <- rep(1, nrow(variants));
  variants$SymbolType[variants$Type=="INS"] = 2;
  variants$SymbolType[variants$Type=="DEL"] = 6;
  rownames(variants) <- variants$VariantId;
  bc <- merge(variants, bc);
  cat("...found", nrow(variants), "unique variants.\n");
  
  ##
  ## For each (variant, sample) pair, sum the total reference and variant coverage and compute
  ## the minor allele frequency.
  ##
  variantCounts <- unique(bc[,c("VariantId", "Sample")]);
  variantCounts$Ref.Cov <- rep(0, nrow(variantCounts));
  variantCounts$Var.Cov <- rep(0, nrow(variantCounts));
  for (i in 1:nrow(variantCounts)) {
    varId <- variantCounts[i, "VariantId"];
    samp <- variantCounts[i, "Sample"];
    tmp <- apply(bc[bc$VariantId==varId & bc$Sample==samp, c("Ref.Cov", "Var.Cov")], 2, sum);
    variantCounts[i, "Ref.Cov"] <- tmp["Ref.Cov"];
    variantCounts[i, "Var.Cov"] <- tmp["Var.Cov"];
  }
  variantCounts$Total.Cov <- variantCounts$Ref.Cov + variantCounts$Var.Cov;
  variantCounts$Minor.Allele.Frac <- pmin(variantCounts$Var.Cov, variantCounts$Ref.Cov) / variantCounts$Total.Cov;
  cat("...generated variant counts for", nrow(variantCounts), "sample/variant pairs.\n");
  
}


##
##
##

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

if(length(c(badSamps1, badSamps2)) > 0){
  allNames <- allNames[-which(allNames %in% c(badSamps1, badSamps2))]
}

## Merge the amplicon information into the read count table
df <- dfRaw;
###my added one-liner
#df[is.na(df)] <- 0

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

## Add a column for the total pool counts
if (length(poolNames)>1) {
  df$TotalPool <- apply(df[,poolNames], 1, sum);
} else {
  df$TotalPool <- df[[poolNames[1]]];
}

# genes <- unique(df$Gene);
# geneInfo <- getGeneInfo(df);
#
# boxplot(df$GC ~ df$Gene, las=2, cex.axis=0.6, at=order(genes), col=(geneInfo$Color)[order(genes)], names=(geneInfo$Label)[order(genes)], ylab="GC Fraction")

## Create initial geneInfo table before dropping any probes
cat("Creating initial gene information table...\n");
genes <- unique(df$Gene);
geneInfo <- getGeneInfo(df);


# df2 - used to keep SNP markers
df2 <- df

## Drop the noisiest amplicons
cat("Dropping", (noisyAmpliconQuantile*100), "percent of amplicons with least pool counts...\n");
minPoolCount <- quantile(df$TotalPool, noisyAmpliconQuantile);
df <- df[df$TotalPool >= minPoolCount,];
cat("...reduced to", nrow(df), "amplicons; dropped all with less than", minPoolCount, "total pool reads.\n");

## Drop amplicons in genes with a low total probe count
keepGenes <- geneInfo[geneInfo$NumProbes >= minAmpliconsPerGene,]$Gene;
df <- df[df$Gene %in% keepGenes,];
cat("...reduced to", nrow(df), "amplicons; dropped all genes with less than", minAmpliconsPerGene, "total probes.\n");

if (normalPool) {
  df$Weights <- pmax(1, df$TotalPool) / sum(df$TotalPool);
} else {
  tmpFrac <- df[,poolNames];
  for (poolName in poolNames) {
    #seems it usually expects more than one samples ... so a quick fix would be to make the data into a df
    #tmpFrac <- data.frame(tmpFrac, stringsAsFactors = FALSE)
    #names(tmpFrac) <- poolNames
    tmpFrac[[poolName]] <- tmpFrac[[poolName]] / sum(tmpFrac[[poolName]]);
  }
  df$Weights <- apply(tmpFrac[,poolNames], 1, median);
}

####what are these weights for .... so in my mind they are the (amplicon counts for segment/total number of amplicon segments)



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
#write.table(df,"df_amps.txt",quote=FALSE,row.names=FALSE)




## For each sample, calculate expected counts
##
dfExpected <- df;
for (x in allNames) {
  dfExpected[[x]] <- df$Weights * sum(df[[x]]);
}


##basically df$Weights is the chance of it happening * (the total number of occurences?) i.e expectation
## another way to look at it is to find a weighted average i.e probability of it occuring times the value it represents
## sum(df[[x]]) finds each individual value


## Calculate ratio between actual counts and expected counts
##
dfResidRatio <- df;
for (x in allNames) {
  dfResidRatio[[x]] <- pmax(1, df[[x]]) / dfExpected[[x]];
}

###dfResidRatio is equal to number of actual counts/expected counts -> can be used for likelihood ratio which
### is found to be ratio of actual counts/ expected counts


#write.table(dfResidRatio,"dfResidRatio.txt",quote=FALSE,row.names=FALSE)

## Apply lowess correction to (log of) residual ratio
##
dfCorrectionRatio <- df;
for (x in allNames) {
  dfCorrectionRatio[[x]] <- exp(lowessCorrect(log(dfResidRatio[[x]]), df$GC));
}


### Do a lowess corecction for the log likelihood - why use log likelihood - via wikiepdia stated below
###Finding the maximum of a function often involves taking the derivative of a function and solving for 
###the parameter being maximized, and this is often easier when the function being maximized is a
###log-likelihood rather than the aaron likelihood function.
###side-note: the logarithm of a function achieves its maximum value at the same points as the function itself,
###and hence the log-likelihood can be used in place of the likelihood in maximum likelihood estimation and related techniques


###In this case though you log and then exponentiate it ... but why? is it just to make a monotic transformation so that
### lowess regresssion has an easier run time?


#write.table(dfCorrectionRatio,"dfCorrectionRatio.txt",quote=FALSE,row.names=FALSE)

#for corrections I need to do 



geneEst <- data.frame(Gene=genes);
rownames(geneEst) <- geneEst$Gene;
rawGeneEst <- geneEst;
geneEstRelErr <- geneEst;
geneEstSampleRelErr <- geneEst;
sampleMedianGeneEst <- dfResidRatio;
#write.table(df,"df_master.txt",quote=FALSE)


# 20210601: KH
# for the gc corrected counts for all amplicons minus the bottom fifth percentile in terms of pooled reads
# Only thing I did was create a parallel process for df2 similar to df 1 without dropping features with amplicons less than 3
# 20211103: KH kept all amplicons b/c the segmentation process should smooth outliers anyways

if (normalPool) {
  df2$Weights <- pmax(1, df2$TotalPool) / sum(df2$TotalPool);
} else {
  tmpFrac <- df2[,poolNames];
  for (poolName in poolNames) {
    tmpFrac[[poolName]] <- tmpFrac[[poolName]] / sum(tmpFrac[[poolName]]);
  }
  df2$Weights <- apply(tmpFrac[,poolNames], 1, median);
}

dfExpected2 <- df2;
for (x in allNames) {
  dfExpected2[[x]] <- df2$Weights * sum(df2[[x]]);
}

dfResidRatio2 <- df2;
for (x in allNames) {
  dfResidRatio2[[x]] <- pmax(1, df2[[x]]) / dfExpected2[[x]];
}

dfCorrectionRatio2 <- df2;
for (x in allNames) {
  dfCorrectionRatio2[[x]] <- exp(lowessCorrect(log(dfResidRatio2[[x]]), df2$GC));
}

dfCorrectedCounts <- dfCorrectionRatio2;

for (x in allNames) {
  wts <- df$Weights;
  rawCounts <- df[[x]];
  uncorrectedCounts <- dfResidRatio[[x]];
  correctedCounts <- dfResidRatio[[x]] / dfCorrectionRatio[[x]];
  
  
  ### new for segmentation - the var2's are based off of amp-unfiltered df2
  correctedCounts2 <- dfResidRatio2[[x]] / dfCorrectionRatio2[[x]];
  dfCorrectedCounts[[x]] <- correctedCounts2;
  
  ### divide the likelihood ratio by corrected likelihood ratio from lowess (why do you divide?)
  ### isn't lowess just a regression techinique, so shouldn't you be subtracting the residuals
  ### resid ratio is eqaul to (actual counts/expected count) / (corrected actual/ corrected expected) (???)
  rawGeneEst[[x]] <- rep(0, length(genes));
  geneEst[[x]] <- rep(0, length(genes));
  geneEstRelErr[[x]] <- rep(0, length(genes));
  geneEstSampleRelErr[[x]] <- rep(0, length(genes));
  for (gene in geneEst$Gene) {
    ix <- (df$Gene==gene);
    #write.table(ix,paste("correctCounts.",x,".txt",sep=""),quote=FALSE)
    geneEst[gene, x] <- weighted.mean(correctedCounts[ix], wts[ix]);
    rawGeneEst[gene, x] <- weighted.mean(uncorrectedCounts[ix], wts[ix]);
    geneEstRelErr[gene, x] <- 1.0 / sqrt(sum(rawCounts[ix]));
    
    ### The relative error is often used to compare approximations of numbers of widely differing size
    
    # geneEstSampleRelErr[gene, x] <- sd(correctedCounts[ix]) / sqrt(length(correctedCounts[ix])) / geneEst[gene, x];
    geneEstSampleRelErr[gene, x] <- sqrt(weighted.mean((correctedCounts[ix] - geneEst[gene, x])**2, wts[ix]) / (length(correctedCounts[ix]) - 1)) / geneEst[gene, x];
    
    ###I guess this does the same, but it calculates per indivudal sample - definitely more complicated, but formula should be somehwere online im guessing
    
  }
  sampleMedianGeneEst[[x]] <- rep(median(geneEst[[x]]), length(sampleMedianGeneEst[[x]]));
  geneEst[[x]] <- geneEst[[x]] / median(geneEst[[x]]);
  rawGeneEst[[x]] <- rawGeneEst[[x]] / median(rawGeneEst[[x]]);
  ###This is where I should run the SVA + regression for correction 
  
  
  ### I'm guessing you do this b/c you want to find how many times more than the baseline (which in our base in the median)
  ### gene copy numbers there are
}

### these counts actually work for segmentation - idea is to make a copy-number raitos and denom is mean of normals (after correctedcounts and gc-correction)


### pretty sure df$weights is used to create normal ratios already
# if (length(poolNames)>1) {
#   normalAmpMean <- apply(dfCorrectedCounts[, poolNames], 1, mean);
# } else {
#   normalAmpMean <- dfCorrectedCounts[[poolNames[1]]];
# }
# 
# 
# 
# dfCorrectedCounts2 <- dfCorrectedCounts
# for (x in allNames) {
#   dfCorrectedCounts2[[x]] <- dfCorrectedCounts[[x]]/normalAmpMean
# }
# 
# 
# dfZscoresCounts <- df2
# for (x in allNames) {
#   dfZscoresCounts[[x]] <- dfResidRatio2[[x]]/dfCorrectionRatio2[[x]]
# }





## Output plotted amplicon-level CN ratio - DHH, 2015/12/10
dfPlotValue <- df;
for (x in allNames) {
  dfPlotValue[[x]] <- dfResidRatio[[x]] / dfCorrectionRatio[[x]] / sampleMedianGeneEst[[x]];
  ### I think this is jsut corrected counts / sampleMedianEst  to give ratio of copy-number
  
}
dfPlotValue$GenomeIndex <- seq(1,nrow(dfPlotValue),1)
#write.table(dfPlotValue,"ampLevel.allSamps.plottedValues.txt",quote=FALSE,row.names=FALSE,sep="\t")
#write.table(dfPlotValue,"/home/kevhu/tmp/ampLevel.allSamps.plottedValues.txt",quote=FALSE,row.names=FALSE,sep="\t")



if (length(poolNames)>1) {
  ## genePoolSD <- data.frame(Gene=genes, SD=apply(geneEst[,poolNames], 1, function(x) 1.4826 * mad(as.numeric(x))));
  genePoolSD <- data.frame(Gene=genes, SD=apply(geneEst[,poolNames], 1, function(x) sd(as.numeric(x))));
} else {
  genePoolSD <- data.frame(Gene=genes, SD=rep(NA, length(genes)));
}

psampMultiple <- function(samps, ylim=NULL, cex.axis=0.75, ax=TRUE, main=NA) {
  geneSpacing <- 0;
  chromSpacing <- 0;
  if (is.na(main)) main="";
  
  xlim <- c(0, nrow(ampliconInfo) + 1 + geneSpacing*(1 + nrow(geneInfo)) + 23*chromSpacing);
  colors = rep(c("black", "purple4", "purple", "blue", "darkturquoise", "green4", "green", "goldenrod1", "firebrick1"), 10);
  
  samp <- samps[length(samps)];
  plot(y=log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]]), x=(geneSpacing*df$GeneNum + chromSpacing*(df$ChromNum-1) + df$AmpliconIndex), col=colors[length(samps)], cex=0.08, xaxt="n", main=main, ylab="Log2(CN Ratio)", xlim=xlim, xlab="", ylim=ylim, xaxs="i", cex.axis=0.75, cex.lab=0.75);
  
  for (i in (length(samps)-1):1) {
    samp <- samps[i];
    points(y=log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]]), x=(geneSpacing*df$GeneNum + chromSpacing*(df$ChromNum-1) + df$AmpliconIndex), col=colors[i], cex=0.08);
  }
  
  for (i in 1:length(samps)) {
    samp <- samps[i];
    for (gene in c("ERBB2", "DDR2", "RB1")) {
      lines(x=geneInfo[gene, c("MinIndex", "MaxIndex")] + geneSpacing*c(-0.5,0.5) + geneSpacing*geneInfo[gene, "GeneNum"] + chromSpacing*geneInfo[gene, "ChromNum"], y=rep(log2(geneEst[gene, samp]),2), col=colors[i]);
    }
  }
  
  legend(x=30, y=5.0, fill=colors[9:1], legend=paste(c("8", "16", "24", "32", "40", "56", "64", "72", "80"), "%")[9:1], cex=0.5, border="white");
  
  if (ax && cex.axis > 0) {
    axis(1, at=0.5*(geneInfo$MinIndex+geneInfo$MaxIndex)+geneSpacing*geneInfo$GeneNum+chromSpacing*(geneInfo$ChromNum-1), labels=geneInfo$Gene, las=2, cex.axis=cex.axis, mgp=c(1,0.5,0), tick=FALSE);
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
  
  for (i in 1:20) {
    if (chromProbes[i] > 0) {
      text(x=(chromMin[i]+chromMax[i])*0.5, y=ylim[1], col="black", labels=paste(i, sep=""), cex=0.5, srt=90);
    }
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

psampUncorrected <- function(samp, ylim=NULL, cex.axis=0.75) {
  plot(y=log2(dfResidRatio[[samp]] / median(dfResidRatio[[samp]])), x=(10*df$GeneNum + df$AmpliconIndex), col=df$Color, cex=20.0*sqrt(df$Weights), xaxt="n", main=samp, ylab="Log2(Raw CN Ratio)", xlab="", ylim=ylim);
  
  for (gene in geneInfo$Gene) {
    lines(x=geneInfo[gene, c("MinIndex", "MaxIndex")] + c(-5,5) + 10*geneInfo[gene, "GeneNum"], y=rep(log2(rawGeneEst[gene, samp]),2));
  }
  
  if (cex.axis > 0) {
    axis(1, at=0.5*(geneInfo$MinIndex+geneInfo$MaxIndex)+10*geneInfo$GeneNum, labels=geneInfo$Gene, las=2, cex.axis=cex.axis, mgp=c(1,0.5,0), tick=FALSE);
  }
}

pmultsnps <- function(samps, ax=TRUE, shiftByRatio=1.0) {
  
  tmp <- merge(merge(merge(variantCounts[variantCounts$Sample %in% samps,], variants), ampliconInfo, by="AmpliconId"), geneInfo, by="Gene")
  
  plot(y=tmp$Minor.Allele.Frac, x=tmp$AmpliconIndex + 10*tmp$GeneNum, col=tmp$Color, cex=0.05*sqrt(tmp$Total.Cov), xaxt="n", xlim=c(0,nrow(ampliconInfo)+10*nrow(geneInfo)), ylim=c(0,0.5), pch=tmp$SymbolType, ylab="Minor Allele Fraction", xlab="");
  
  ###point of having minor allele freqs. (?)
  
  if (ax) {
    axis(1, at=0.5*(geneInfo$MinIndex+geneInfo$MaxIndex)+10*geneInfo$GeneNum, labels=geneInfo$Gene, las=2, cex.axis=0.75, mgp=c(1,0.5,0), tick=FALSE);
  }
  
  for (gene in geneInfo$Gene) {
    lines(x=rep(geneInfo[gene, "MinIndex"] + 10*geneInfo[gene,"GeneNum"] - 5, 2), y=c(0,0.5), col="grey");
  }
}


psnps <- function(samp, ax=TRUE, shiftByRatio=1.0) {
  
  tmp <- merge(merge(merge(variantCounts[variantCounts$Sample == samp,], variants), ampliconInfo, by="AmpliconId"), geneInfo, by="Gene")
  
  geneSpacing <- 40;
  chromSpacing <- 0;
  
  xlim <- c(0, nrow(ampliconInfo) + 1 + geneSpacing*(1 + nrow(geneInfo)) + 23*chromSpacing);
  
  plot(y=tmp$Minor.Allele.Frac, x=tmp$AmpliconIndex + geneSpacing*tmp$GeneNum, col=tmp$Color, cex=0.05*sqrt(tmp$Total.Cov), xaxt="n", xlim=xlim, ylim=c(0,0.5), pch=tmp$SymbolType, ylab="Minor Allele Fraction", xlab="", xaxs="i");
  
  if (ax) {
    labelPositions = 0.5*(geneInfo$MinIndex+geneInfo$MaxIndex)+geneSpacing*geneInfo$GeneNum+chromSpacing*(geneInfo$ChromNum-1);
    axis(1, at=labelPositions, labels=geneInfo$Gene, las=2, cex.axis=0.75, mgp=c(1,0.5,0), tick=FALSE);
  }
  
  for (gene in geneInfo$Gene) {
    lines(x=rep(geneInfo[gene, "MinIndex"] + geneSpacing*geneInfo[gene,"GeneNum"] - 0.5*geneSpacing, 2), y=c(0,0.5), col="grey");
  }
  
  alleleFreqFromCN <- function(cn) {
    if (cn>=1.0) { return(1.0 / (2.0*cn)); }
    return((2.0*cn - 1) / (2.0*cn));
  }
  
  for (gene in geneInfo$Gene) {
    lines(x=geneInfo[gene, c("MinIndex", "MaxIndex")] + geneSpacing*c(-0.5,0.5) + geneSpacing*geneInfo[gene, "GeneNum"], y=rep(alleleFreqFromCN(geneEst[gene, samp]*shiftByRatio), 2), col="black", lty=(if (geneEst[gene, samp]*shiftByRatio>1.0) "dotted" else "solid"));
  }
}

pboth <- function(samp, ylim=NULL, shiftByRatio=1.0, cex.axis=0.75, analysisResults=NULL) {
  def.par <- par(no.readonly = TRUE);
  layout(matrix(c(1,1,2),3,1,byrow=TRUE));
  par(mar=c(4,5,3,2));
  psamp(samp, ylim=ylim, analysisResults=analysisResults, cex.axis=cex.axis);
  par(mar=c(2,5,0,2));
  psnps(samp, ax=FALSE, shiftByRatio=shiftByRatio);
  par(def.par);
}

pCorrectedAndUncorrected <- function(samp, ylim=NULL) {
  def.par <- par(no.readonly = TRUE);
  layout(matrix(c(1,1,1,1,2,2,2),7,1,byrow=TRUE));
  par(mar=c(4,5,3,2));
  psampUncorrected(samp, ylim=ylim);
  par(mar=c(2,5,0,2));
  psamp(samp, ax=FALSE, main="", ylim=ylim);
  par(def.par);
}

binomPvals <- function(samp, genes) {
  
  tmp <- merge(merge(merge(variantCounts[variantCounts$Sample == samp,], variants), ampliconInfo, by="AmpliconId"), geneInfo, by="Gene");
  tmp <- tmp[tmp$Gene %in% genes, c("Minor.Allele.Frac", "Total.Cov")];
  tmp$Minor.Cov <- round(tmp$Total.Cov * tmp$Minor.Allele.Frac);
  # tmp <- tmp[tmp$Minor.Allele.Frac > 0.03,];
  return(tmp);
}


dirToWrite <- str_remove(ampliconInfoFilename, "amplicon.GCinput.txt")
#print(dirToWrite)

### 20220922: with pangu graphs large
tmpWidth <- round(nrow(geneEst)/20);
tmpHeight <- round(nrow(geneEst)/30);

combinedTable <- NULL
for (samp in allNames) {
  
  analysisResults <- data.frame(Gene=genes, Sample=rep(samp, length(genes)), NumProbes=geneInfo$NumProbes, CopyNumberRatio=geneEst[,samp], RawCopyNumberRatio=rawGeneEst[,samp], ProbeError=geneEstSampleRelErr[,samp], PoolError=genePoolSD$SD, ZScore=rep(0.0, length(genes)), Call=rep(0, length(genes)), Sig=rep(0, length(genes)), Log10PValue=rep(0.0, length(genes)));
  rownames(analysisResults) <- genes;
  for (gene in genes) {
    row <- analysisResults[gene,];
    est <- row$CopyNumberRatio;
    z1 <- (est - 1.) / (row$PoolError);
    z2 <- if (row$NumProbes > 1) { (est - 1.) / (row$ProbeError * est) } else NA;
    z <- min(abs(z1), abs(z2));
    if (is.na(z1)) {
      z <- abs(z2);
    } else if (is.na(z2)) {
      z <- abs(z1);
    }
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
  # pdf(paste("out.", samp, ".pdf", sep=""), width=7 * 2, height=5 * 2, useDingbats=FALSE);
  pdf(paste("out.", samp, ".pdf", sep=""), width=tmpWidth * 2, height=5 * 2, useDingbats=FALSE);
  if (variantData) pboth(samp, cex.axis=0.5, ylim=c(-2,4), analysisResults=analysisResults) else psamp(samp, cex.axis=0.5, ylim=c(-4,4), analysisResults=analysisResults);
  dev.off();
  write.table(format(analysisResults, scientific=FALSE, digits=4), file=paste("calls.", samp, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE);
  analysisResults$Sample <- str_remove_all(analysisResults$Sample, " ")
  combinedTable <- rbind(combinedTable, format(analysisResults, scientific=FALSE, digits=4))
}


colnames(geneEst)  <- str_remove_all(colnames(geneEst) , " ")
write.table(combinedTable, file="combinedCalls.txt", quote=FALSE, sep="\t", row.names=FALSE);
write.table(geneEst, file="cnMatrix_gene.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names = TRUE);
write.table(dfCorrectedCounts, file = "cnAmplicon_matrix.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names = TRUE);
write.table(dfCorrectedCounts, file = "gcCorrectedCounts_matrix.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names = TRUE);


warnings();
