library(RTCGAToolbox)
library(reshape2)
library(stringr)
library(copynumber)

getFreqData <- function(data){  
  #Check if segments or data:
  if(colnames(data)[1]=="sampleID" || colnames(data)[2]=="arm"){
    #input is segments data frame;
    #could be on a multiseg-format -> convert to uniseg-format:
    #need to convert to appropriate format
    #first find intersection of all breakpoints:
    chr <- unique(data[,2])
    bpts <- matrix(NA,nrow=0,ncol=2)
    for(j in 1:length(chr)){
      subseg <- subsetSegments(data,chrom=chr[j])
      bpts <- rbind(bpts,data.frame(chr[j],sort(unique(c(subseg$start.pos,subseg$end.pos)))))
    } 
    colnames(bpts) <- c("chrom","pos")
    
    #Then interpolate to get pcf-val in all breakpoints
    data = interpolate.pcf(data, bpts) 
  }
  #If not segments, the input data should already be on an appropriate format (chrom, pos, estimates)
  
  return(data) 
} 


reducingFreqBed <- function(df, idx){
  reducedDf <- NULL
  for (i in 1:(length(idx) - 1)) {
    idx1 <- idx[i]
    idx2 <- idx[i+1]
    tmpDf <- df[idx1:idx2,]
    tmpChr <- tmpDf$chr[1]
    tmpStart <- min(tmpDf$pos)
    tmpEnd <- max(tmpDf$pos) - 1
    tmpFreq <- tmpDf[,3][1]
    
    tmpVec <- c("Chr" = tmpChr, "Start" = tmpStart, "End" = tmpEnd, "Freq" = tmpFreq)
    reducedDf <- rbind(reducedDf, tmpVec)
  }
  reducedDf <- data.frame(reducedDf, stringsAsFactors = FALSE)
  return(reducedDf)
}

# used on df with segmentation to get frequency bed
ampsDels <- function(df){
  getFreqOut <- getFreqData(df)
  freqAmp <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] > 0.2)*100
  freqDel <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] < -0.2)*100
  
  freqDf_Amp <- cbind(getFreqOut[,1:2], "amp" = freqAmp)
  freqDf_Del <- cbind(getFreqOut[,1:2], "del" = freqDel)
  
  ampRedIdx <- c(1,1+which(diff(freqDf_Amp$amp)!=0))
  delRedIdx <- c(1,1+which(diff(freqDf_Del$del)!=0))
  
  res <- list(freqDf_Amp, ampRedIdx, freqDf_Del, delRedIdx)
}



# later add functions or something to better get Cn info stratified by genotype

setwd("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/")
getFirehoseData(dataset="OV", runDate="20160128",
                clinical=FALSE, Mutation=TRUE)

getFirehoseData(dataset="OV", runDate="20160128",
                clinical=FALSE, CNASNP=TRUE)

getFirehoseData(dataset="OV", runDate="20160128",
                clinical=FALSE, CNVSNP=TRUE)

ov_mut <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/20160128-OV-Mutations-AllSamples.txt",
                     sep = "\t", header = TRUE, stringsAsFactors = FALSE)

ov_CNA <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/20160128-OV-CNASNPHg19.txt",
                     sep = "\t", header = TRUE, stringsAsFactors = FALSE)

ov_geneLevel <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",
                           sep = "\t", header = TRUE, stringsAsFactors = FALSE)



ov_geneLevel_melt <- melt(ov_geneLevel, id.vars = "Gene.Symbol")
ov_geneLevel_melt$variable <- as.character(ov_geneLevel_melt$variable)
ov_geneLevel_melt$variable <- str_replace_all(ov_geneLevel_melt$variable, "\\.", "\\-")


brca1Losses <- which(ov_geneLevel_melt$Gene.Symbol == "BRCA1" & ov_geneLevel_melt$value <= -1)
brca1Losses2copy <- which(ov_geneLevel_melt$Gene.Symbol == "BRCA1" & ov_geneLevel_melt$value <= -2)

ov_brca1_loss <- ov_geneLevel_melt[brca1Losses,]
ov_brca1_loss_2copy <- ov_geneLevel_melt[brca1Losses2copy,]


brca1MutSamples <- ov_mut$Tumor_Sample_Barcode[which(ov_mut$Hugo_Symbol == "BRCA1")]

# 11/12 of the brca1 muts have a second 1 copy deletion - loh
brca1_mut_1copy <- ov_brca1_loss$variable[grep(paste(brca1MutSamples, collapse = "|"), ov_brca1_loss$variable)]
brca1_2copy <- ov_brca1_loss_2copy$variable

brca1_null_tcga <- c(brca1_mut_1copy, brca1_2copy)



# below should only be for looking at frequnency plots and putting them into
# bed format compatible with synteny plot

ov_CNA_brca1 <- ov_CNA[grep(paste(brca1_null_tcga, collapse = "|"), ov_CNA$Sample), ]
ov_CNA_brca1_noY <- ov_CNA_brca1[-which(ov_CNA_brca1$Chromosome == "24"),]
colnames(ov_CNA_brca1_noY) <- c("sampleID","chrom", "start.pos",
                                "end.pos", "n.probes", "mean")

plotFreq(ov_CNA_brca1_noY, thres.gain = 0.2, thres.loss = -0.2)


ov_brca1_freq_out <- getFreqData(ov_CNA_brca1_noY)
ampDels_brca <- ampsDels(ov_brca1_freq_out)
ov_brca1_amp_bed <- reducingFreqBed(ampDels_brca[[1]],ampDels_brca[[2]])
ov_brca1_del_bed <- reducingFreqBed(ampDels_brca[[3]],ampDels_brca[[4]])


write.table(ov_brca1_amp_bed, "/mnt/DATA5/tmp/kev/misc/20210424tcga_ov_brca1_amp_0.2.bed",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(ov_brca1_del_bed, "/mnt/DATA5/tmp/kev/misc/20210424tcga_ov_brca1_del_0.2.bed",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



ov_CNA_noY <- ov_CNA[-which(ov_CNA$Chromosome == "24"),]
colnames(ov_CNA_noY) <- c("sampleID","chrom", "start.pos",
                          "end.pos", "n.probes", "mean")

plotFreq(ov_CNA_noY, thres.gain = 0.2, thres.loss = -0.2)


ov_freq_out <- getFreqData(ov_CNA_noY)
ampDels_ov <- ampsDels(ov_freq_out)
ov_amp_bed <- reducingFreqBed(ampDels_ov[[1]],ampDels_ov[[2]])
ov_del_bed <- reducingFreqBed(ampDels_ov[[3]],ampDels_ov[[4]])


write.table(ov_amp_bed, "/mnt/DATA5/tmp/kev/misc/20210424tcga_ov_amp_0.2.bed",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(ov_del_bed, "/mnt/DATA5/tmp/kev/misc/20210424tcga_ov_del_0.2.bed",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



# next steps would be to implement these beds into the synteny plots
# not much to do but subset by chrom approriately

# one thing I still might want to do is have a measurement for synteny
# i.e some fisher exact per syntenic block of change?

plotFreq()
