
args <- commandArgs(trailingOnly=TRUE);

#args <- c("/mnt/DATA4/kevhu/tmp2/amplicon.GCinput.txt", "/mnt/DATA4/kevhu/tmp2/sampleInfo.input.txt", "/mnt/DATA4/kevhu/tmp2/amplicon.combinedCoverage.input.txt",
#          "--min-amplicons-per-gene=3")

args <- c("/mnt/DATA4/kevhu/tmp3/amplicon.GCinput.txt", "/mnt/DATA4/kevhu/tmp3/sampleInfo.input.txt", "/mnt/DATA4/kevhu/tmp3/amplicon.combinedCoverage.input.txt",
          "--min-amplicons-per-gene=3")

args <- c("/mnt/DATA4/kevhu/choLab/20191104mouseSamps/cnCalls/amplicon.GCinput.txt", "/mnt/DATA4/kevhu/choLab/20191104mouseSamps/cnCalls/sampleInfo.input.txt", "/mnt/DATA4/kevhu/choLab/20191104mouseSamps/cnCalls/amplicon.combinedCoverage.input.txt",
          "--min-amplicons-per-gene=3")


args <- c("/mnt/DATA4/kevhu/choLab/20200301analysis/amplicon.GCinput.txt", "/mnt/DATA4/kevhu/choLab/20200301analysis/sampleInfo.input.txt",
          "/mnt/DATA4/kevhu/choLab/20200301analysis/amplicon.combinedCoverage.input.txt",
          "--min-amplicons-per-gene=3")


args <- c("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-144-NEN_OCAP_2_361_357/amplicon.GCinput.txt",
          "/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-144-NEN_OCAP_2_361_357/sampleInfo.input.txt",
          "/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-144-NEN_OCAP_2_361_357/amplicon.combinedCoverage.input.txt",
          "--min-amplicons-per-gene=3")

showQ <- TRUE;
showZ <- FALSE;
colorSigGenes <- TRUE;
colorSigGenesMaxQ <- 0.01;
labelSigGenesOnly <- FALSE;
minAmpliconsPerGene <- 4;
noisyAmpliconQuantile <- 0.05;
#noisyAmpliconQuantile <- 0.00;

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



####commented out the part dropping the noisiest amplicon so I can keep the Rb1 deletion as its own gene

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

dfCorrectedCounts <- dfCorrectionRatio;
tmp <- dfCorrectionRatio;

dfNormals <- df[, poolNames]

for (x in allNames) {
  wts <- df$Weights;
  rawCounts <- df[[x]];
  uncorrectedCounts <- dfResidRatio[[x]];
  correctedCounts <- dfResidRatio[[x]] / dfCorrectionRatio[[x]];
  
  ### new for segmentation
  dfCorrectedCounts[[x]] <- correctedCounts;
  
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
  #if (x == "MG_21X53") {
  #  print(geneEst[[x]])
  #  print(median(geneEst[[x]]))
  #}
  geneEst[[x]] <- geneEst[[x]] / median(geneEst[[x]]);
  rawGeneEst[[x]] <- rawGeneEst[[x]] / median(rawGeneEst[[x]]);
  ###This is where I should run the SVA + regression for correction 
  
  
  ### I'm guessing you do this b/c you want to find how many times more than the baseline (which in our base in the median)
  ### gene copy numbers there are
}

#normalAmpMedian <- apply(dfCorrectedCounts[, poolNames], 1, median)
#normalAmpMean <- apply(dfCorrectedCounts[, poolNames], 1, mean)

if (length(poolNames)>1) {
  normalAmpMean <- apply(dfCorrectedCounts[, poolNames], 1, mean);
} else {
  normalAmpMean <- dfCorrectedCounts[[poolNames[1]]];
}



dfCorrectedCounts2 <- dfCorrectedCounts
for (x in allNames) {
  dfCorrectedCounts2[[x]] <- dfCorrectedCounts[[x]]/normalAmpMean
}


dfZscoresCounts <- df
for (x in allNames) {
  dfZscoresCounts[[x]] <- df[[x]]/dfCorrectionRatio[[x]]
}


### did commente below to see average deviation from zero and mean did better
mean((1 - apply(dfCorrectedCounts2[,poolNames], 2, mean)))

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



psamp <- function(samp, ylim=NULL, cex.axis=0.75, ax=TRUE, main=NA, analysisResults=NULL) {
  geneSpacing <- 40;
  chromSpacing <- 0;
  if (is.na(main)) main=samp;
  
  xlim <- c(0, nrow(ampliconInfo) + 1 + geneSpacing*(1 + nrow(geneInfo)) + 23*chromSpacing);
  
  ## plot(y=log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]]), x=(geneSpacing*df$GeneNum + chromSpacing*(df$ChromNum-1) + df$AmpliconIndex), col=df$Color, cex=20.0*sqrt(df$Weights), xaxt="n", main=main, ylab="Log2(CN Ratio)", xlim=NULL, xlab="", ylim=ylim)
  
  ## c(0,nrow(ampliconInfo)+geneSpacing*nrow(geneInfo)+chromSpacing*23), xlab="", ylim=ylim);
  
  plot(y=log2(dfResidRatio[[samp]] / dfCorrectionRatio[[samp]] / sampleMedianGeneEst[[samp]]), x=(geneSpacing*df$GeneNum + chromSpacing*(as.numeric(df$ChromNum)-1) + df$AmpliconIndex), col=df$Color, cex=0.3, xaxt="n", main=main, ylab="Log2(CN Ratio)", xlim=xlim, xlab="", ylim=ylim, xaxs="i", cex.axis=0.9, cex.lab=0.9);
  
  for (gene in geneInfo$Gene) {
    lines(x=geneInfo[gene, c("MinIndex", "MaxIndex")] + geneSpacing*c(-0.5,0.5) + geneSpacing*geneInfo[gene, "GeneNum"] + chromSpacing* as.numeric(geneInfo[gene, "ChromNum"]), y=rep(log2(geneEst[gene, samp]),2));
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



combinedTable <- NULL
# tmpColnames <- colnames(geneEst)
# tmpColnames <- str_remove(tmpColnames, "_MG_X.*")
# tmpColnames <- str_remove(str_remove(tmpColnames, "^X"), "_X.*")
# tmpColnames <- tolower(str_remove_all(str_remove_all(str_remove(tmpColnames, "X.*"), "_"), "\\."))
# tmpColnames  <- str_remove(tmpColnames, "o")
# tcDf <- read.table("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt", sep = "\t",
#                    header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# matchingTc <- tcDf$tc[match(tmpColnames, tcDf$sample)]
# geneEst[,2:ncol(geneEst)] <- sweep(geneEst[,2:ncol(geneEst)], 2, matchingTc[2:length(matchingTc)], "/")

### just do this to regular geneEst[,2:ncol(geneEst)] instead of making geneEst2

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
  setwd("/mnt/DATA5/tmp/kev/testOCAP/")
  pdf(paste("out.", samp, ".pdf", sep=""), width=7 * 1.25, height=5 * 1.25, useDingbats=FALSE);
  if (variantData) pboth(samp, cex.axis=0.5, ylim=c(-2,4), analysisResults=analysisResults) else psamp(samp, cex.axis=0.5, ylim=c(-4,4), analysisResults=analysisResults);
  dev.off();
  write.table(format(analysisResults, scientific=FALSE, digits=4), file=paste("calls.", samp, "_tc.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE);
  combinedTable <- rbind(combinedTable, format(analysisResults, scientific=FALSE, digits=4))
}
write.table(combinedTable, file=paste("combinedCalls.txt"), quote=FALSE, sep="\t", row.names=FALSE);


warnings();






source("/mnt/DATA4/kevhu/scripts/20181205newpsampFunction.R")


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

for (samp in allNames) {
  setwd("/mnt/DATA4/kevhu/tmp2/")
  pdf(file = paste("newsamp.", samp, ".pdf",sep = ""),onefile = TRUE, useDingbats = FALSE,  width=7 * 1.25, height=5 * 1.25)
  newpsamp(samp)
  dev.off()
  rm(testggplot)
}


warnings();










