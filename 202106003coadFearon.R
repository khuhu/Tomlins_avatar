### first thing to do is to give standard results i.e heatmap
library(pheatmap)
library(copynumber)
library(GenomicFeatures)
library(dplyr)
library(dbplyr)
library(circlize)
library(stringr)
library(DNAcopy)
library(stringr)
#library(optparse)
library(data.table)

fearonAnno <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20210601Fearon_withData.xlsx")
cnGeneDf <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/cnMatrix_gene.txt",
                       sep = "\t", stringsAsFactors = FALSE, header = TRUE)
cnGeneDf2 <- cnGeneDf[, c(1, grep("EF", colnames(cnGeneDf)))]
newColnames <- colnames(cnGeneDf2)
newColnames <- str_remove(newColnames, "\\_MG.*")
colnames(cnGeneDf2) <- newColnames

fearonAnno$`Sequencing name` <- str_replace(fearonAnno$`Sequencing name`,
                                            "-", "_")


### reordering the samples so they line up with genotype annotation
### and putting genes at the end

cnGeneDf2 <- cnGeneDf2[ ,order(colnames(cnGeneDf2))]
cnGeneDf2 <- cnGeneDf2[,c(ncol(cnGeneDf2), 1:ncol(cnGeneDf2)-1)]
delRows <- cnGeneDf2[grep("Del", cnGeneDf2$Gene), 2:ncol(cnGeneDf2)]
cnGeneDf2 <- cnGeneDf2[-grep("Del", cnGeneDf2$Gene),]
cnGeneNames <- c(cnGeneDf2$Gene, "Trp53Del","ApcDel")
cnGeneDf2 <- cnGeneDf2[,-1]
cnGeneDf2 <- rbind(cnGeneDf2, delRows)
cnGeneDf2[1:nrow(cnGeneDf2),] <- lapply(cnGeneDf2[1:nrow(cnGeneDf2),], function(x) log2(x))
rownames(cnGeneDf2) <- cnGeneNames
fearonCnMat <- t(cnGeneDf2)
fearonCnMat[fearonCnMat < -3] <- -3
fearonCnMat[fearonCnMat > 3] <- 3
fearonCnMat[fearonCnMat > -0.2 & fearonCnMat < 0.2] <- 0

annoTab <- data.frame("Genotype" = fearonAnno$Genotype, stringsAsFactors = FALSE)
rownames(annoTab) <- fearonAnno$`Sample ID`

rownames(fearonCnMat) <- rownames(annoTab)

annoCol <- list("Genotype" = c("CDX2P-CreERT2 Apcfl/+, KrasLSLG12D/+, p53R270H+/ex2-10 fl" = "darkgreen",
                               "CDX2P-CreERT2 Apcfl/+, KrasLSLG12D/+, p53 ex2-10 fl/fl"  = "darkred", 
                               "Wild-type liver" = "white",
                               "wild-type proximal colon mucosa" = "grey"))


heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-3,3,6/1000)


heatmap_graph <- pheatmap(mat = fearonCnMat, cluster_rows = FALSE, cluster_cols = FALSE, color = heatMapCol,
                          breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                          border_color = "black", annotation_row = annoTab, annotation_colors = annoCol)


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210609fearonHeatmap.pdf", useDingbats = TRUE, width = 15)
heatmap_graph
dev.off()


### adjusting the half trp53 fl/ trp53 R270H 
afTcDf <- data.frame("KrasG12D" = fearonAnno$`KRAS (AF)`,
                     "Trp53R270H" = fearonAnno$`TP53 (AF)`,
                     "Tumor content" = round(1 - 2^fearonAnno$`TP53 log2(cnr)`, digits = 2))
afTcDf$Tumor.content[1] <- 0.68
afTcDf$Tumor.content[4] <- 1.0
afTcDf$Tumor.content[7] <- 1.0
tmpAfTcDf <- cbind("samples" = fearonAnno$`Sequencing name`, afTcDf)
write.table(tmpAfTcDf, "/mnt/DATA5/tmp/kev/misc/20210621fearonTc.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

afTcCol <- colorRampPalette(c("#FFFFFF","#CC5500"))(100)
colors.breaks <- seq(-0,1,1/100)


afTc_graph <- pheatmap(mat = afTcDf, cluster_rows = FALSE, cluster_cols = FALSE, color = afTcCol,
                       breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                       border_color = "black")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210609fearonAfTc.pdf", useDingbats = TRUE, width = 15)
afTc_graph
dev.off()

###
###
### functions for getting freq data

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
  # 20210629: for whatever reason cant find the mean of a long logical vector
  # i.e only TRUE AND FALSE
  #freqAmp <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] > 0.2)*100
  #freqDel <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] < -0.2)*100
  freqAmp <- apply(df[,-c(1:2)], 1, function(x) 
    length(which(x > 0.2))/length(x)) * 100
  freqDel <- apply(df[,-c(1:2)], 1, function(x) 
    length(which(x < -0.2))/length(x)) * 100
  
  freqDf_Amp <- cbind(getFreqOut[,1:2], "amp" = freqAmp)
  freqDf_Del <- cbind(getFreqOut[,1:2], "del" = freqDel)
  
  ampRedIdx <- c(1,1+which(diff(freqDf_Amp$amp)!=0))
  delRedIdx <- c(1,1+which(diff(freqDf_Del$del)!=0))
  
  res <- list(freqDf_Amp, ampRedIdx, freqDf_Del, delRedIdx)
}


ampliconData <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/gcCorrectedCounts_matrix.txt",
                           sep = "\t", stringsAsFactors = FALSE, header = TRUE)
bedfile <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/bed.txt", sep = "\t",
                      stringsAsFactors = FALSE, header = FALSE)


varResults <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_anno.txt",
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE)
 
tcga_snp_cn <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/20160128-COAD-CNASNPHg19.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#tcga_seq_cn <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/20160128-COAD-CNAseq.txt",
#                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# 
# tcga_gistic_cn <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/20160128-COAD-all_thresholded.by_genes.txt",
#                           sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# 
# tcga_mut <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/20210608-COAD-mutations.txt",
#                        sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

coad_KAP <- read.table("/mnt/DATA5/tmp/kev/misc/20210628_panCancerCoadKAPIds.txt", sep = "\t",
                       header = TRUE, stringsAsFactors = FALSE)
coad_KAP$samples <- str_replace_all(coad_KAP$samples, "\\.", "\\-")
coad <- read.table("/mnt/DATA5/tmp/kev/misc/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv", sep = "\t",
                   header = TRUE, stringsAsFactors = FALSE)
tcga_arm <- read.table("/mnt/DATA5/tmp/kev/misc/20210628panCancerArmTable_hg19.txt",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
tcga_cna <- read.table("/mnt/DATA5/tmp/kev/misc/20210628panCancerCnaTable_hg19.txt",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

coad_arm <- tcga_arm[which(tcga_arm$Sample %in% coad$Sample.ID),]
coad_cna <- tcga_cna[which(tcga_cna$Sample %in% coad$Sample.ID),]

kap_arm <- tcga_arm[grep(paste(coad_KAP$samples, collapse = "|"), tcga_arm$Sample),]
kap_cna <- tcga_cna[grep(paste(coad_KAP$samples, collapse = "|"), tcga_cna$Sample),]



coad_arm2 <- data.frame("sampleID" = coad_arm$Sample, "chrom" = coad_arm$Chromosome,
                        "start.pos" = coad_arm$Start, "end.pos" = coad_arm$End,
                        "n.probes" = NA, "mean" = coad_arm$newTotalCN, stringsAsFactors = FALSE)

coad_cna2 <- data.frame("sampleID" = coad_cna$Sample, "chrom" = coad_cna$Chromosome,
                        "start.pos" = coad_cna$Start, "end.pos" = coad_cna$End,
                        "n.probes" = NA, "mean" = coad_cna$newTotalCN, stringsAsFactors = FALSE)

kap_arm2 <- data.frame("sampleID" = kap_arm$Sample, "chrom" = kap_arm$Chromosome,
                        "start.pos" = kap_arm$Start, "end.pos" = kap_arm$End,
                        "n.probes" = NA, "mean" = kap_arm$newTotalCN, stringsAsFactors = FALSE)
kap_cna2 <- data.frame("sampleID" = kap_cna$Sample, "chrom" = kap_cna$Chromosome,
                       "start.pos" = kap_cna$Start, "end.pos" = kap_cna$End,
                       "n.probes" = NA, "mean" = kap_cna$newTotalCN, stringsAsFactors = FALSE)


coad_arm_freq_out <- getFreqData(coad_arm2)
ampDels_coad_arm <- ampsDels(coad_arm_freq_out)
coad_arm_amp_bed <- reducingFreqBed(ampDels_coad_arm[[1]], ampDels_coad_arm[[2]])
coad_arm_del_bed <- reducingFreqBed(ampDels_coad_arm[[3]], ampDels_coad_arm[[4]])

coad_cna_freq_out <- getFreqData(coad_cna2)
ampDels_coad_cna <- ampsDels(coad_cna_freq_out)
coad_cna_amp_bed <- reducingFreqBed(ampDels_coad_cna[[1]], ampDels_coad_cna[[2]])
coad_cna_del_bed <- reducingFreqBed(ampDels_coad_cna[[3]], ampDels_coad_cna[[4]])


kap_arm_freq_out <- getFreqData(kap_arm2)
ampDels_kap_arm <- ampsDels(kap_arm_freq_out)
kap_arm_amp_bed <- reducingFreqBed(ampDels_kap_arm[[1]], ampDels_kap_arm[[2]])
kap_arm_del_bed <- reducingFreqBed(ampDels_kap_arm[[3]], ampDels_kap_arm[[4]])

kap_cna_freq_out <- getFreqData(kap_cna2)
ampDels_kap_cna <- ampsDels(kap_cna_freq_out)
kap_cna_amp_bed <- reducingFreqBed(ampDels_kap_cna[[1]], ampDels_kap_cna[[2]])
kap_cna_del_bed <- reducingFreqBed(ampDels_kap_cna[[3]], ampDels_kap_cna[[4]])



write.table(coad_arm_amp_bed, "/mnt/DATA5/tmp/kev/misc/20210629_coadArmAmp.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(coad_arm_del_bed, "/mnt/DATA5/tmp/kev/misc/20210629_coadArmDel.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(coad_cna_amp_bed, "/mnt/DATA5/tmp/kev/misc/20210629_coadCnaAmp.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(coad_cna_del_bed, "/mnt/DATA5/tmp/kev/misc/20210629_coadCnaDel.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


### var analysis & deleted region analysis
###
###



### Vars
mgpVars <- c("hom", "het")
varResults <- varResults[-which(varResults$mm10_mpgpv6_Indels %in% mgpVars),]
varResults$Sample <- str_replace(varResults$Sample, "\\-", "\\_")

segResults <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/segResults.txt",
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t")

segResults <- segResults[grep("EF", segResults$ID),]
sampNames <- unique(segResults$ID)


apcDel <- paste0("AMP_", 3904:3907)
apcDel_df <- ampliconData[which(ampliconData$AmpliconId %in% apcDel), ]
log2(apply(apcDel_df[,2:40], 2, mean))

tp53Del <- paste0("AMP_", 2569:2584)
tp53Del_df <- ampliconData[which(ampliconData$AmpliconId %in% tp53Del),]
log2(apply(tp53Del_df[,2:40], 2, mean))


### thought there was a exon deleted in ATM but the 3 segment region AMPS 2063-65
### 63 and 65 are depleted, but 64 has copy-number ratio of 1

varResults <- varResults[which(varResults$Sample %in% sampNames),]
fdpFilt <- which(varResults$FDP > 20)
faoFilt <- which(varResults$FAO > 5)
freqFilt <- which(varResults$AF > 0.05)
hrunFilt <- which(varResults$HRUN < 4)
strandRatio <- intersect(which(varResults$FSAF/varResults$FSAR > 0.2),
                         which(varResults$FSAF/varResults$FSAR < 5))
qualFilt <- which(varResults$QUAL > 60)
goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, hrunFilt, qualFilt))
varResults<- varResults[goodSamps,]

varResults_exon <- varResults[which(varResults$Func.refGene == "exonic"),]

