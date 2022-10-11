source("/home/kevhu/scripts/20210802syntenyFunctions.R")
library(pheatmap)
library(ggdendro)
library(stringr)
library(gridExtra)

### this 

coadNonzeroTc <- read.table("/mnt/DATA5/tmp/kev/misc/20220521onlyEfdApcTc.txt", sep = "\t",
                            header = TRUE, stringsAsFactors = FALSE)
# coadNonzeroTc <- coadNonzeroTc[which(coadNonzeroTc$tc > 0.2), ]

annotationList <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20211207yingSecondSet.xlsx")
annotationList$sampleStripped <- str_remove(nameStripper(annotationList$`#`), "\\-")
annotationList <- annotationList[-which(annotationList$sampleStripped == "efd36"),]
annotationList$sampleStripped[which(annotationList$sampleStripped %in% paste0("efd", c(1,2,4:7)))] <- paste0("efd0", c(1,2,4:7))


histologyList <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/AKP&CDX2Braftumor.xlsx", skip = 1)
histologyList$sampleStripped <- str_remove(nameStripper(histologyList$Sample), "\\-")
histologyList <- histologyList[1:61,]


fearonVar1 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-76-MG_test1_255_185_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonVar2 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonVar3 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-156-MG_Fearon_20210809_374_382_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonVar4 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)


combinedVars <- rbind(fearonVar1, fearonVar2, fearonVar3, fearonVar4)
mgpVars <- c("hom", "het", "unknown")
combinedVars <- combinedVars[which(combinedVars$mm10_mpgpv6_Indels == ""),]
combinedVars$FSAF <- as.numeric(combinedVars$FSAF)
combinedVars$FSAR <- as.numeric(combinedVars$FSAR)

fdpFilt <- which(combinedVars$FDP > 100)
faoFilt <- which(combinedVars$FAO > 10)
freqFilt <- which(combinedVars$AF > 0.10)
hrunFilt <- which(combinedVars$HRUN < 4)
strandRatio <- intersect(which(combinedVars$FSAF/combinedVars$FSAR > 0.2),
                         which(combinedVars$FSAF/combinedVars$FSAR < 5))
qualFilt <- which(combinedVars$QUAL >= 100)

goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, qualFilt, hrunFilt))
combinedVars_goodsamps <- combinedVars[goodSamps,]
combinedVars_goodsamps_exon <- combinedVars_goodsamps[which(combinedVars_goodsamps$Func.refGene == "exonic"), ]
combinedVars_goodsamps_exon  <- combinedVars_goodsamps_exon[grep("EF", combinedVars_goodsamps_exon$Sample),]

combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[which(combinedVars_goodsamps_exon$name == "."),]                                         
combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[-grep("nonframeshift",
                                                                 combinedVars_goodsamps_exon$ExonicFunc.refGene),]  
combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[-grep("^synonymous SNV",
                                                                 combinedVars_goodsamps_exon$ExonicFunc.refGene),]  

combinedVars_goodsamps_exon$Sample <- str_remove(str_remove(nameStripper(combinedVars_goodsamps_exon$Sample), "\\-"), "x.*")

goodSamps_anno <- annotationList$sampleStripped[-which(annotationList$Qual == "Bad")]

combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[which(combinedVars_goodsamps_exon$Sample %in% goodSamps_anno),]


mutOfInterest <- c("Kras:NM_021284:exon2:c.G35A:p.G12D", "Trp53:NM_001127233:exon8:c.G809A:p.R270H,Trp53:NM_011640:exon8:c.G809A:p.R270H")
mutDf <- combinedVars_goodsamps_exon[grep(paste(mutOfInterest, collapse = "|"), combinedVars_goodsamps_exon$AAChange.refGene),]

combinedVars_somatic <- combinedVars_goodsamps_exon[-grep(paste(mutOfInterest, collapse = "|"), combinedVars_goodsamps_exon$AAChange.refGene),]


mouseNormal <- c("EF_D03_MG_X14", "MG_18X50", "MG_21X53", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")


# fearonCn1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-76-MG_test1_255_185_ApcTrp53/cnMatrix_gene.txt",
#                         sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
# 
# fearonCn2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_ApcTrp53/cnMatrix_gene.txt",
#                         sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
# fearonCn2 <- fearonCn2[,-which(colnames(fearonCn2) %in% mouseNormal)]
# 
# fearonCn3 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-156-MG_Fearon_20210809_374_382_ApcTrp53/cnMatrix_gene.txt",
#                         sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
# fearonCn3 <- fearonCn3[,-which(colnames(fearonCn3) %in% mouseNormal)]
# 
# fearonCn4 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_ApcTrp53/cnMatrix_gene.txt",
#                         sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
# fearonCn4 <- fearonCn4[,-which(colnames(fearonCn4) %in% mouseNormal)]


fearonCn1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-76-MG_test1_255_185_Apc/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

fearonCn2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_Apc/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
fearonCn2 <- fearonCn2[,-which(colnames(fearonCn2) %in% mouseNormal)]

fearonCn3 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-156-MG_Fearon_20210809_374_382_Apc/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
fearonCn3 <- fearonCn3[,-which(colnames(fearonCn3) %in% mouseNormal)]

fearonCn4 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_Apc/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
fearonCn4 <- fearonCn4[,-which(colnames(fearonCn4) %in% mouseNormal)]

fearonCombinedCn <- cbind(fearonCn1, fearonCn2[,2:ncol(fearonCn2)], fearonCn3[,2:ncol(fearonCn3)], fearonCn4[,2:ncol(fearonCn4)])

tmpDf <- melt(fearonCombinedCn)
tmpDf <- tmpDf[grep("Del", tmpDf$Gene),]
# delDf <- dcast(tmpDf, variable ~ Gene)
delDf <- tmpDf[,2:3]
colnames(delDf)[1] <- "Sample"
delDf$Sample <- str_remove(nameStripper(delDf$Sample), "x.*")
delDf$KrasG12D <- 0
delDf$Trp53R270H <- 0

tmpMut <- mutDf[which(mutDf$AAChange.refGene == mutOfInterest[1]),]
delDf$KrasG12D <- tmpMut$AF[match(delDf$Sample, tmpMut$Sample)]
delDf$KrasG12D[which(is.na(delDf$KrasG12D))] <- 0

tmpMut <- mutDf[which(mutDf$AAChange.refGene == mutOfInterest[2]),]
delDf$Trp53R270H<- tmpMut$AF[match(delDf$Sample, tmpMut$Sample)]
delDf$Trp53R270H[which(is.na(delDf$Trp53R270H))] <- 0

colnames(fearonCombinedCn)[2:ncol(fearonCombinedCn)] <- str_remove(nameStripper(colnames(fearonCombinedCn)[2:ncol(fearonCombinedCn)]), "x.*")
cnMat <- fearonCombinedCn[,which(colnames(fearonCombinedCn) %in% annotationList$sampleStripped)]
rownames(cnMat) <- fearonCombinedCn$Gene

### cn heatmap

# combinedCalls1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-76-MG_test1_255_185_ApcTrp53/combinedCalls.txt",
#                              sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
# combinedCalls2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_ApcTrp53/combinedCalls.txt",
#                              sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
# combinedCalls3 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-156-MG_Fearon_20210809_374_382_ApcTrp53/combinedCalls.txt",
#                              sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
# combinedCalls4 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_ApcTrp53/combinedCalls.txt",
#                              sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)


combinedCalls1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-76-MG_test1_255_185_Apc/combinedCalls.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
combinedCalls2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_Apc/combinedCalls.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
combinedCalls3 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-156-MG_Fearon_20210809_374_382_Apc/combinedCalls.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
combinedCalls4 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_Apc/combinedCalls.txt",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)


allcombinedCalls <- rbind(combinedCalls1, combinedCalls2, combinedCalls3, combinedCalls4)
allcombinedCalls$Sample <- str_remove(nameStripper(allcombinedCalls$Sample), "x.*")

geneNames <- rownames(cnMat)

### really poor q-values for event he deletions in this data set .... going to filter just by cnr
xGenes <- c("Bcor", "Enox2", "Ar", "Atrx",
            "Diaph2", "Frmpd4")
delGenes <- c("ApcDel", "Trp53Del")

cnMat2 <- log2(cnMat)
cnMat2[cnMat2 > 2] <- 2
cnMat2[cnMat2 < -2] <- -2
# cnMat2[abs(cnMat2) < 0.2] <- 0
cnMat2 <- t(cnMat2)

cnMat2 <- cnMat2[which(rownames(cnMat2) %in% coadNonzeroTc$sample),]
rownames(cnMat2) <- annotationList$`Sample ID`[match(rownames(cnMat2), annotationList$sampleStripped)]

cnDel <- cnMat2[,which(colnames(cnMat2) %in% delGenes)]
cnMat2 <- cnMat2[,-which(colnames(cnMat2) %in% delGenes)]
cnMat2_noX <- cnMat2[,-which(colnames(cnMat2) %in% xGenes)]

allFearonAnno <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220531allFearonAnno.xlsx")
allFearonAnno$Sample <- nameStripper(allFearonAnno$`Sequencing name`)
allFearonAnno$Sample <- str_remove(allFearonAnno$Sample, "\\-")
allFearonAnno$Sample[1:9] <- paste0("efd0", 1:9)
allFearonAnno$Pathology[1:2] <- "Adenocarcinoma"

annoTab2 <- data.frame("Genotype" = annotationList$Genotype[match(rownames(cnMat2), annotationList$`Sample ID`)])
rownames(annoTab2) <- annotationList$`Sample ID`[match(rownames(cnMat2), annotationList$`Sample ID`)]
annoTab2$Pathology <- allFearonAnno$Pathology[match(rownames(cnMat2), allFearonAnno$`Sample ID`)]

coadNonzeroTc_filt <- coadNonzeroTc[which(coadNonzeroTc$sample %in% annotationList$sampleStripped), ]
coadNonzeroTc_filt$ID <- annotationList$`Sample ID`[match(coadNonzeroTc_filt$sample, annotationList$sampleStripped)]
tcCont <- coadNonzeroTc_filt$tc[match(rownames(annoTab2), coadNonzeroTc_filt$ID)]


tcVar <- NULL
for (i in tcCont) {
  if (i > 0.6) {
    res <- "High (0.6 - 1.0)"
  } else if(i > 0.4 & i < 0.6){
    res <- "Medium (0.4 - 0.6)"
  } else if(i > 0.1 & i < 0.4){
    res <- "Low (0.1 - 0.4)"
  } else if (i < 0.1){
    res <- "None (< 0.1)"
  }
  tcVar <- c(tcVar, res)
}

annoTab2$TumorContent <- tcVar


annoCol <- list("Genotype" = c("CDX2P-CreERT2+  CDX2 fl/fl Braf fl/+" = "darkred", "CreERT2Apcfl/+, KrasLSLG12D/+, p53R270H/ex2-10" = "grey",
                               "CreERT2Apcfl/+, Krasfl/+, p53R270H fl/ex2-10" = "lightblue", "CreERT2Apcfl/+, KrasLSLG12D/+, p532-10 -/-" = "orange",
                               "CreERT2Apcfl/+, KrasLSLG12D/+, p53ex2-10fl/+ Sox9 fl/+" = "darkblue", "CreERT2Apcfl/+, KrasLSLG12D/+" = "darkgreen",
                               "Wild type" = "white", "CDX2P-CreERT2 Apcfl/+, KrasLSLG12D/+, p53R270H+/ex2-10 fl" = "goldenrod4",
                               "CDX2P-CreERT2 Apcfl/+, KrasLSLG12D/+, p53 ex2-10 fl/fl" = "lightgreen"),
                "Pathology" = c("Adenoma" = "gold4", "Adenocarcinoma" = "darkred", "hyperplastic polyp" = "antiquewhite2", "Normal" = "white"),
                "TumorContent" = c("High (0.6 - 1.0)" = "green4", "Medium (0.4 - 0.6)" = "green2", "Low (0.1 - 0.4)" = "slategray4", "None (< 0.1)" = "black"))



heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-2,2,4/1000)

heatmap_graph <- pheatmap(mat = cnMat2_noX, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                          breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                          border_color = "black", annotation_row = annoTab2, annotation_colors = annoCol)

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20220613yingGoodSamplesCnr.pdf", useDingbats = FALSE, height = 8, width = 16)
heatmap_graph
dev.off()

geneNames_noX <- colnames(cnMat2_noX)
yingSigMat <- NULL
for (i in unique(allcombinedCalls$Sample)) {
  tmpDf <- allcombinedCalls[which(allcombinedCalls$Sample == i),]
  tmpDf <- tmpDf[match(geneNames_noX,tmpDf$Gene),]
  tmpVector <- ifelse(10^tmpDf$Log10PValue < 0.05, 1, 0)
  yingSigMat <- rbind(yingSigMat, tmpVector)
}


rownames(yingSigMat) <- unique(allcombinedCalls$Sample)
colnames(yingSigMat) <- geneNames_noX
yingSigMat <- t(yingSigMat)
samps <- annotationList$sampleStripped[match(rownames(cnMat2_noX), annotationList$`Sample ID`)]
yingSigMat <- yingSigMat[, which(colnames(yingSigMat) %in% samps)]
colnames(yingSigMat) <- annotationList$`Sample ID`[match(colnames(yingSigMat), annotationList$sampleStripped)]
yingSigMat <- t(yingSigMat)


cnMat2_noX_qval <- cnMat2_noX
cnMat2_noX_qval[yingSigMat < 1] <- 0

pheatmap(mat = cnMat2_noX_qval, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
         breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
         border_color = "black", annotation_row = annoTab2, annotation_colors = annoCol)




# delDf_good <- delDf[which(delDf$Sample %in% annotationList_good$sampleStripped),]
# delDf_good$Sample <- annotationList_good$`Sample ID`[match(annotationList_good$sampleStripped, delDf_good$Sample)]
# delDf_good[,2:3] <- log2(delDf_good[,2:3])
# 
# delDf_good2 <- delDf_good[,2:5]
# rownames(delDf_good2) <- delDf_good$Sample
# 
# delDf_good2_del <- delDf_good2[,1:2]
# delDf_good2_af <- delDf_good2[,3:4]
# 
# delDf_good3 <- delDf_good
# delDf_good3$anno <- annoTab2$Genotype
# delDf_good_order <- delDf_good3[heatmap_graph$tree_row$order,]
# 
# delDf_good_order$Histology <- histologyList$Histology[match(delDf_good_order$Sample, histologyList$`Tissue Name`)]
# delDf_good_order$sampleStripped <- annotationList_good$sampleStripped[match(delDf_good_order$Sample, annotationList_good$`Sample ID`)]
# delDf_good_order <- delDf_good_order[order(delDf_good_order$anno),]
# delDf_good_order[,2:5] <- signif(delDf_good_order[,2:5], digits = 3)


### making a segmental chart first need to format data to put into getFreqData

mouseNormal <- c("EF_D03_MG_X14", "MG_18X50", "MG_21X53", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")


segRes1 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-76-MG_test1_255_185/segResults.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
segRes2 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/segResults.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
segRes3 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/segResults.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
segRes4 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384/segResults.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

segRes2 <- segRes2[-which(segRes2$ID %in% mouseNormal),]
segRes3 <- segRes3[-which(segRes3$ID %in% mouseNormal),]
segRes4 <- segRes4[-which(segRes4$ID %in% mouseNormal),]

combinedSegRes <- rbind(segRes1, segRes2, segRes3, segRes4)
colnames(combinedSegRes)


# combinedSegRes$seg.mean[which(abs(combinedSegRes$seg.mean) < 0.2)] <- 0
combinedSegRes$seg.mean[which(combinedSegRes$q.val2 > 0.05)] <- 0
combinedSegRes$seg.mean[which(combinedSegRes$seg.mean < -2)] <- -2
combinedSegRes$seg.mean[which(combinedSegRes$seg.mean > 2)] <- 2
combinedSegRes$ID <- nameStripper(combinedSegRes$ID)
combinedSegRes$ID <- str_remove(combinedSegRes$ID, "x.*")
combinedSegRes$length <- combinedSegRes$loc.end - combinedSegRes$loc.start
combinedSegRes$seg.mean[which(combinedSegRes$length < 20e6)] <- 0

selfCalls2 <- cbind("Sample" = rownames(annoTab2), "ID" = allFearonAnno$Sample[match(rownames(annoTab2), allFearonAnno$`Sample ID`)] ,annoTab2)

combinedSegRes2 <- combinedSegRes[which(combinedSegRes$ID %in% selfCalls2$ID),]
combinedSegRes2$chrom <- str_replace(combinedSegRes2$chrom, "23", "20")

combinedSegRes2_form <- combinedSegRes2[c("ID", "chrom", "loc.start", "loc.end",
                                          "num.mark", "seg.mean")]
colnames(combinedSegRes2_form) <- c("sampleID", "chrom", "start.pos",
                                    "end.pos", "n.probes", "mean")
combinedSegRes2_form$n.probes <- NA

combinedSegRes2_freq <- getFreqData(combinedSegRes2_form)
combinedSegRes2$chrom <- str_replace(combinedSegRes2$chrom, "23", "20")

### scale colors according to heatmap colors i.e larger loss darker blue
### 

df_cn <- combinedSegRes2
df_freq <-  combinedSegRes2_freq
chromTextSpec = NULL

plotSegHeatmap <- function(df_cn, df_freq, main = NULL,
                           chromTextSpec = NULL, segsize = 3){
  require(ggplot2)
  
  if(is.null(chromTextSpec)){
    chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "mm10"){
    chromTextdf <- read.table("mnt/DATA5/tmp/kev/misc/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "hg19") {
    chromTextdf <- read.table("/mnt/DATA5/tmp/kev/misc/20210801hg19_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  }
  
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  
  
  ### new way of assigning colors basically searches for where value lies within range of colors
  heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
  
  highVal <- max(abs(df_cn$seg.mean))
  colPalRange <- seq(-highVal, highVal, 2 * highVal/999)
  
  df_cn$loc.start <- df_cn$loc.start/1e6
  df_cn$loc.end <- df_cn$loc.end/1e6
  df_cn$col <- "#000000"
  for (i in seq_along(df_cn$col)) {
    if (df_cn$seg.mean[i] == 0) {
      df_cn$col[i] <- "#FFFFFF"
    } else if(df_cn$seg.mean[i]  > 0) {
      df_cn$col[i] <- heatMapCol[which(colPalRange  == min(colPalRange[which(colPalRange >= df_cn$seg.mean[i])]))]
    } else if (df_cn$seg.mean[i] < 0) {
      df_cn$col[i] <- heatMapCol[which(colPalRange  == max(colPalRange[which(colPalRange <= df_cn$seg.mean[i])]))]
    }
  }
  #df_cn$col[which(df_cn$seg.mean < 0.2 & df_cn$seg.mean > -0.2)] <- "#FFFFFF"
  
  
  
  df_dist <- df_freq[,3:ncol(df_freq)]
  rownames(df_dist) <- paste0("chr", df_freq$chr, ":", df_freq$pos)
  df_dist <- t(df_dist)
  tmpModel <- hclust(dist(df_dist), method = "complete")
  dendo <- as.dendrogram(tmpModel)
  dendo_dat <- dendro_data(dendo, type = "rectangle")
  treeOrder <- dendo_dat$labels$label
  
  ### use tree order to get the dendo order for heatmap
  df_cn2 <- NULL
  for (i in seq_along(treeOrder)) {
    tmpCn <- df_cn[which(df_cn$ID == treeOrder[i]),]
    for (j in unique(tmpCn$chrom)) {
      tmpCn$loc.start[which(tmpCn$chrom == j)] <- tmpCn$loc.start[which(tmpCn$chrom == j)] +
        chromTextdf$graphingStart[which(chromTextdf$chrom == j)]
      tmpCn$loc.end[which(tmpCn$chrom == j)] <- tmpCn$loc.end[which(tmpCn$chrom == j)] +
        chromTextdf$graphingStart[which(chromTextdf$chrom == j)]
    }
    
    tmpCn$ypos <- 1 + i * 0.2
    df_cn2 <- rbind(df_cn2, tmpCn)
  }
  
  df_cn2 <- df_cn2[,c("ID","chrom","loc.start", "loc.end", "ypos","seg.mean", "col")]
  colnames(df_cn2) <- c("id", "chrom", "xstart", "xend", "ypos","cn", "col")
  
  ### assign is temp - need for annoTable
  assign("treeOrder", treeOrder, envir = parent.frame() )
  assign("df_cn2", df_cn2, envir = parent.frame() )
  
  xlabels <- unique(df_cn2$id)
  hline <- unique(df_cn2$ypos) - 0.1
  hline <- c(hline, max(hline) + 0.2)
  chromTextdf$ypos <- max(df_cn2$ypos) + 0.2
  
  
  ### graphing + grob layering
  dendo_graph <- ggplot(ggdendro::segment(dendo_dat)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    coord_flip() + 
    scale_y_reverse(expand = c(0,0)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0), 
                                          add = c(0.5, 1.5))) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank())
  
  
  
  heat_graph <- ggplot(df_cn2) +
    geom_rect(aes(xmin = xstart, xmax = xend,
                  ymin = ypos - 0.1, ymax = ypos + 0.1), fill = df_cn2$col) + 
    geom_vline(xintercept = chromBreak) +
    geom_hline(yintercept = hline, color = "grey") + 
    scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
    scale_y_continuous(breaks=seq(min(df_cn2$ypos), max(df_cn2$ypos), 0.2),
                       labels = xlabels, expand = expansion(mult = c(0, 0.02))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=6)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos, label = chrom), size = 2.5)
  
  dendo_grob <- ggplotGrob(dendo_graph)
  heat_grob <- ggplotGrob(heat_graph)
  
  grobLayout <- rbind(c(1,rep(2,19)),
                      c(1,rep(2,19)),
                      c(1,rep(2,19)),
                      c(1,rep(2,19)))
  grid.arrange(dendo_grob, heat_grob,
               layout_matrix=grobLayout)
}


combinedSegRes2$ID <- annotationList$`Sample ID`[match(combinedSegRes2$ID, annotationList$sampleStripped)]
colnames(combinedSegRes2_freq)[3:ncol(combinedSegRes2_freq)] <- annotationList$`Sample ID`[match(colnames(combinedSegRes2_freq)[3:ncol(combinedSegRes2_freq)], annotationList$sampleStripped)]

plotSegHeatmap(combinedSegRes2, combinedSegRes2_freq, segsize = 3)

combinedSegRes2_noX <- combinedSegRes2[-which(combinedSegRes2$chrom == "20"),]
combinedSegRes2_noX_form <- combinedSegRes2_noX[c("ID", "chrom", "loc.start", "loc.end",
                                                  "num.mark", "seg.mean")]
colnames(combinedSegRes2_noX_form) <- c("sampleID", "chrom", "start.pos",
                                        "end.pos", "n.probes", "mean")
combinedSegRes2_noX_form$n.probes <- NA
combinedSegRes2_noX_freq <- getFreqData(combinedSegRes2_noX_form)

combinedGraph <- plotSegHeatmap(combinedSegRes2_noX, combinedSegRes2_noX_freq, segsize = 1)
selfCalls2 <- selfCalls2[match(treeOrder, selfCalls2$Sample),]

### order samps in selfCalls with tree-order prior to making table of colors and rects

annoTab_coords  <- NULL
for (i in 3:(dim(selfCalls2)[2])) {
  sampVar <- selfCalls2[,i]
  tmpDf <- data.frame("xpos" = 1.2 + 0.2 * (i-2), 
                      "ypos" = seq(1.2, dim(selfCalls2)[1] * .2 + 1, .2),
                      "label"  = sampVar)
  annoTab_coords <- rbind(annoTab_coords, tmpDf)
}

dfAnnoColIdx <- data.frame("labelNames" = c(unique(selfCalls2$Genotype), unique(selfCalls2$Pathology)), 
                           "colors" = c("chartreuse3", "yellow", "white", "lightblue", "darkblue", "darkmagenta","grey",
                                        "darkgreen", "antiquewhite2", "darkred", "white"))

annoTab_coords$col <- "white"
annoTab_coords$col <- dfAnnoColIdx$colors[match(annoTab_coords$label, dfAnnoColIdx$labelNames)]

anno_graph <- ggplot(annoTab_coords) +
  geom_rect(aes(xmin = xpos - 0.1, xmax = xpos + 0.1,
                ymin = ypos - 0.1, ymax = ypos + 0.1), fill = annoTab_coords$col) + 
  geom_hline(yintercept = seq(1.1, 7.3, .2), color = "black") +
  geom_vline(xintercept = c(1.1, 1.3, 1.5), color = "black") + 
  scale_x_continuous(expand = expansion(mult = c(0, 0), 
                                        add = c(0, 0))) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

grobLayout2.1 <- do.call("rbind", replicate(length(unique(df_cn2$ypos)) - 1, c(rep(1,40), 2), simplify = FALSE))
grobLayout2 <- rbind(c(rep(1,40), NA), 
                     grobLayout2.1)


grid.arrange(combinedGraph, anno_graph,
             layout_matrix=grobLayout2)


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20211218yingGoodSamps.pdf", useDingbats = FALSE, height = 5, width = 18)
grid.arrange(combinedGraph, anno_graph,
             layout_matrix=grobLayout2)
dev.off()


### same thing above but trying to do for abs copy number
### idea would be to first get segs and intersect them for each gene

allNoshadSeg <- read.table("/mnt/DATA5/tmp/kev/noshadSnpFiltHclust15AllV2/absoluteRes.txt", sep = "\t",
                           stringsAsFactors = FALSE, header = TRUE, fill = TRUE)

allPloidyCalls <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220517allHgscPloidy.xlsx")
allPloidyCalls <- allPloidyCalls[-which(allPloidyCalls$ploidy_int == "NA"), ]

allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "mg4", "2027lte")
allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "mg15", "13085lt")
allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "mg20", "14399rt")



allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "14150lt", "14150rt")
allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "14154lt", "14154rot")
allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "14656peritnealmt", "14656peritonealmt")
allNoshadSeg$sample <- str_replace_all(allNoshadSeg$sample, "13576rt", "133576rt")





allNoshadSeg_filt <- allNoshadSeg[grep(paste0(paste(allPloidyCalls$sample, round(as.numeric(allPloidyCalls$purityV2), 3),
                                                    round(as.numeric(allPloidyCalls$ploidyV2), 2), sep = "_"), collapse = "|"), allNoshadSeg$sample),]
allNoshadSeg_filt$sC <- as.numeric(allNoshadSeg_filt$sC)
allNoshadSeg_filt[,2:3] <- lapply(allNoshadSeg_filt[,2:3], as.numeric)

### check
paste(allPloidyCalls$sample, round(as.numeric(allPloidyCalls$purityV2), 3),
      round(as.numeric(allPloidyCalls$ploidyV2), 2), sep = "_")[-which(paste(allPloidyCalls$sample, round(as.numeric(allPloidyCalls$purityV2), 3),
                                                                             round(as.numeric(allPloidyCalls$ploidyV2), 2), sep = "_") %in% unique(allNoshadSeg_filt$sample))]
allNoshadSeg_filt$absCn <- NA
for (i in 1:nrow(allPloidyCalls)) {
  tmpString <- paste(allPloidyCalls$sample[i], round(as.numeric(allPloidyCalls$purityV2[i]), 3),
                     round(as.numeric(allPloidyCalls$ploidyV2[i]), 2), sep = "_")
  tmpMatch <- grep(tmpString, allNoshadSeg_filt$sample)
  tmpDf <- allNoshadSeg_filt[tmpMatch,]
  allNoshadSeg_filt$absCn[tmpMatch] <- round(allNoshadSeg_filt$sC[tmpMatch] - as.numeric(allPloidyCalls$ploidy_int[i]))
}


mouseBed <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed",
                       header = FALSE, stringsAsFactors = FALSE, sep = "\t")
mouseBed$V8 <- firstUpper(mouseBed$V8)

ngsGeneNames <- names(table(mouseBed$V8))[which(table(mouseBed$V8) > 2)]
mouseBed <- mouseBed[which(mouseBed$V8 %in% ngsGeneNames),]

ngsMouseGeneInfo <- NULL
i <- unique(mouseBed$V8)[1]
for (i in unique(mouseBed$V8)) {
  tmpDf <- mouseBed[which(mouseBed$V8 %in% i),]
  ngsMouseGeneInfo <- rbind(ngsMouseGeneInfo, c(tmpDf$V1[1], min(tmpDf$V2), max(tmpDf$V3), tmpDf$V8[1]))
}

ngsMouseGeneInfo <- data.frame(ngsMouseGeneInfo, stringsAsFactors = FALSE)
ngsMouseGeneInfo[,2:3] <- lapply(ngsMouseGeneInfo[,2:3], as.numeric)
colnames(ngsMouseGeneInfo) <- c("chrom", "start", "end", "gene")

ngsMouseGeneInfo <- ngsMouseGeneInfo[which(ngsMouseGeneInfo$gene %in% colnames(cnMat2_noX)),]
ngsMouseGeneInfo_grange <- GRanges(seqnames = ngsMouseGeneInfo$chrom,
                                   IRanges(start = ngsMouseGeneInfo$start, end = ngsMouseGeneInfo$end))

coad_abs_geneMat <- cnMat2_noX
coad_abs_geneMat <- coad_abs_geneMat[-which(rownames(coad_abs_geneMat) == "17705 PC1-4 HP"),]
i <- rownames(coad_abs_geneMat)[1]
j <- 1
for (i in rownames(coad_abs_geneMat)) {
  print(i)
  tmpDf <- coad_abs_geneMat[which(rownames(coad_abs_geneMat) == i),]
  tmpEf <- allFearonAnno$Sample[which(allFearonAnno$`Sample ID` == i)]
  tmpNoshad <- allNoshadSeg_filt[grep(tmpEf, allNoshadSeg_filt$sample),]
  tmpAbsGrange <- GRanges(seqnames = tmpNoshad$chr, IRanges(start = tmpNoshad$start, end = tmpNoshad$end))
  
  ### only a problem when a gene is in between two segments, then for loop probably needed
  for (j in seq_along(tmpDf)) {
    tmpOverlap <- subjectHits(findOverlaps(ngsMouseGeneInfo_grange[j], tmpAbsGrange))
    tmpDf[j] <- mean(tmpNoshad$absCn[tmpOverlap])
  }
  coad_abs_geneMat[which(rownames(coad_abs_geneMat) == i),] <- tmpDf
}



colors.breaks.2 <- seq(-3,3,6/1000)
coad_abs_geneMat <- coad_abs_geneMat[-which(rownames(coad_abs_geneMat) %in% c("19904 PC1-4", "8593 kidney", "8593 liver", "8403 liver")),]
heatmap_graph_abs <- pheatmap(mat = coad_abs_geneMat, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                          breaks = colors.breaks.2, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                          border_color = "black", annotation_row = annoTab2, annotation_colors = annoCol)


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20220613yingGoodSamplesAbsCn.pdf", useDingbats = FALSE, height = 12, width = 18)
heatmap_graph_abs
dev.off()

### next doing it by sample
library(RTNsurvival)

#--- semi-supervised clustering

allFearonAnno$bioSampGroup <- unlist(lapply(str_split(allFearonAnno$`Sample ID`, " "), '[[', 1))

vectorBioGroups <- allFearonAnno$bioSampGroup[match(rownames(coad_abs_geneMat), allFearonAnno$`Sample ID`)]
ss_clust <- hclust_semisupervised(data = coad_abs_geneMat,
                                  groups = split(rownames(coad_abs_geneMat), vectorBioGroups))

plot(ss_clust, hang=-1, cex=0.5)

cnMat_dendo <- as.dendrogram(ss_clust)
cnMat_dendo_dat <- dendro_data(cnMat_dendo, type = "rectangle")
ss_treeOrder <- cnMat_dendo_dat$labels


cl_cb <- function(hcl, mat){
  # Recalculate manhattan distances for reorder method
  dists <- dist(mat, method = "euclidean")
  
  # Perform reordering according to OLO method
  hclust_olo <- reorder(hcl, dists)
  return(hclust_olo)
}


# callback_obj <- cl_cb(cnMat_dendo, coad_abs_geneMat)
# pheatmap(mat = coad_abs_geneMat, cluster_rows = FALSE, cluster_cols = FALSE, color = heatMapCol,
#          breaks = .2, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
#          border_color = "black", annotation_row = annoTab2, annotation_colors = annoCol,
#          clustering_callback = callback_obj)

heatmat_graph_abs <- pheatmap(mat = coad_abs_geneMat, cluster_rows = ss_clust, cluster_cols = FALSE, color = heatMapCol,
                              breaks = colors.breaks.2, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                              border_color = "black", annotation_row = annoTab2, annotation_colors = annoCol,
                              clustering_callback = callback_obj)


delDf_efd <- delDf
delDf_efd <- delDf_efd[which(delDf_efd$Sample %in% annotationList$sampleStripped), ]
delDf_efd$Sample <- annotationList$`Sample ID`[match(delDf_efd$Sample, annotationList$sampleStripped)]
delDf_efd <- delDf_efd[which(delDf_efd$Sample %in% rownames(coad_abs_geneMat)), ]

delDf_efd_gene <- delDf_efd[,2:3]
rownames(delDf_efd_gene) <-delDf_efd$Sample
delDf_efd_gene <- log2(delDf_efd_gene)
delDf_efd_mut <- delDf_efd[,4:5]
rownames(delDf_efd_mut) <-delDf_efd$Sample

delDf_efd_gene_order <- delDf_efd_gene[heatmat_graph_abs$tree_row$order,]
delDf_efd_mut_order <- delDf_efd_mut[heatmat_graph_abs$tree_row$order,]


tcContDf <- data.frame("Sample" = coadNonzeroTc_filt$sample[match(rownames(annoTab2), coadNonzeroTc_filt$ID)],
                       "TumorContent" = coadNonzeroTc_filt$tc[match(rownames(annoTab2), coadNonzeroTc_filt$ID)])
tcContDf$Sample <- annotationList$`Sample ID`[match(tcContDf$Sample, annotationList$sampleStripped)]
tcContDf <- tcContDf[which(tcContDf$Sample %in% rownames(delDf_efd_mut_order)), ]
tcContDf$TumorContent[which(tcContDf$TumorContent < 0)] <- 0
tcContDf <- tcContDf[heatmat_graph_abs$tree_row$order, ]


heatmap_del_cn <- pheatmap(mat = delDf_efd_gene_order, cluster_rows = FALSE, cluster_cols = FALSE, color = heatMapCol,
                      breaks = , fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                      border_color = "black")


heatMapColMut <- colorRampPalette(c("#000000","#8B0000"))(100)
colors.breaks.3 <- seq(0, 1 , 1/100)

delDf_efd_mut_order$TummorContent <- tcContDf$TumorContent

heatmap_del_mut <- pheatmap(mat = delDf_efd_mut_order, cluster_rows = FALSE, cluster_cols = FALSE, color = heatMapColMut,
                       breaks = colors.breaks.3, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                       border_color = "black")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20220616fearonCoadefd34to90_heat.pdf", useDingbats = FALSE , height = 12, width = 18)
heatmat_graph_abs
dev.off()

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20220616fearonCoadefd34to90_delcn.pdf", useDingbats = FALSE , height = 12, width = 18)
heatmap_del_cn
dev.off()


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20220616fearonCoadefd34to90_delmut.pdf", useDingbats = FALSE , height = 12, width = 18)
heatmap_del_mut
dev.off()

### above may not be the most accurate way to determine gene because it calculates it from the segments
### alternatrive method is to scale correctly based on the predicted tc and ploidy



cnMat2_noX_tc <- cnMat2_noX_qval
cnMat2_noX_tc <- cnMat2_noX_tc[1:58,]

ploidyTcTable2 <- allPloidyCalls
ploidyTcTable2$sampleId <- annotationList$`Sample ID`[match(ploidyTcTable2$sample, annotationList$sampleStripped)]
ploidyTcTable2_fearon <- ploidyTcTable2[-which(is.na(ploidyTcTable2$sampleId)),]
ploidyTcTable2_fearon$purityV2 <- as.numeric(ploidyTcTable2_fearon$purityV2)
ploidyTcTable2_fearon$ploidy_int <- as.numeric(ploidyTcTable2_fearon$ploidy_int)

cnMat2_noX_tc <- cnMat2_noX_tc[which(rownames(cnMat2_noX_tc) %in% ploidyTcTable2_fearon$sampleId),]
cnMat2_noX_tc <- 2^cnMat2_noX_tc
# cnMat2_noX_tc[cnMat2_noX_tc < 2^0.2 && cnMat2_noX_tc > 2^-0.2] <- 1


i <- 19
for (i in 1:nrow(cnMat2_noX_tc)) {
  tmpVector <- cnMat2_noX_tc[i,]
  tmpIdx <- which(rownames(cnMat2_noX_tc)[i] == ploidyTcTable2_fearon$sampleId)
  tmpTc <- signif(ploidyTcTable2_fearon$purityV2[tmpIdx], 2)
  tmpPloidy <- ploidyTcTable2_fearon$ploidy_int[tmpIdx]
  
  tmpVector[tmpVector > 1] <- tmpVector[tmpVector > 1] / tmpTc
  tmpVector[tmpVector < 1] <- tmpVector[tmpVector < 1] * tmpTc
  
  tmpVector <- tmpVector * tmpPloidy
  tmpVector <- tmpVector - tmpPloidy
  
  # tmpVector[tmpVector > 0] <- tmpVector[tmpVector > 0] / tmpTc
  # tmpVector[tmpVector < 0] <- tmpVector[tmpVector < 0] * tmpTc
  
  tmpVector <- round2(tmpVector)
  
  cnMat2_noX_tc[i, ] <- tmpVector
}

cnMat2_noX_tc <- round2(cnMat2_noX_tc)
quantile(cnMat2_noX_tc, seq(0,1,0.05))


pheatmap(mat = cnMat2_noX_tc, cluster_rows = ss_clust, cluster_cols = FALSE, color = heatMapCol,
         breaks = colors.breaks.2, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
         border_color = "black", annotation_row = annoTab2, annotation_colors = annoCol,
         clustering_callback = callback_obj)

### cluster things phylogenetically or just do unsupervised clustering with features that are only present in at least two samples .... can't be in 0-1 or all samples

combinedVars_goodsamps_fearon <- combinedVars_goodsamps
combinedVars_goodsamps_fearon$Sample <- nameStripper(combinedVars_goodsamps_fearon$Sample)
combinedVars_goodsamps_fearon$Sample <- str_remove_all(combinedVars_goodsamps_fearon$Sample, "x.*")
combinedVars_goodsamps_fearon$Sample <- str_remove_all(combinedVars_goodsamps_fearon$Sample, "\\-")
combinedVars_goodsamps_fearon$Sample <- annotationList$`Sample ID`[match(combinedVars_goodsamps_fearon$Sample, annotationList$sampleStripped)]
combinedVars_goodsamps_fearon <- combinedVars_goodsamps_fearon[-which(is.na(combinedVars_goodsamps_fearon$Sample)),]

acgt_vector <- c("A", "G", "T", "C")
snpTypes <- c("nonsynonymous SNV", "") 

combinedVars_goodsamps_fearon <- combinedVars_goodsamps_fearon[which(combinedVars_goodsamps_fearon$Ref %in% acgt_vector),]
combinedVars_goodsamps_fearon <- combinedVars_goodsamps_fearon[which(combinedVars_goodsamps_fearon$Alt %in% acgt_vector),]
combinedVars_goodsamps_fearon <- combinedVars_goodsamps_fearon[-which(combinedVars_goodsamps_fearon$ExonicFunc.refGene %in% snpTypes),]

### nothing interesting for SNPs so only using features like high level cn gain or deep del + large chromosomal changes

allNoshadSeg_fearon <- allNoshadSeg_filt[grep("efd", allNoshadSeg_filt$sample),]
allNoshadSeg_fearon$sample <- str_remove_all(allNoshadSeg_fearon$sample, "\\_.*")
  
allNoshadSeg_fearon_graph <- allNoshadSeg_fearon[, c("sample", "chr", "start","end", "K", "absCn")]
colnames(allNoshadSeg_fearon_graph) <- c("sampleID", "chrom", "start.pos",
                                 "end.pos", "n.probes", "mean")
allNoshadSeg_fearon_graph$chrom <- str_remove(allNoshadSeg_fearon_graph$chrom, "chr")
allNoshadSeg_fearon_form <- allNoshadSeg_fearon_graph
noshadFilt_fearon_freq <- getFreqData(allNoshadSeg_fearon_form)
allNoshadSeg_fearon_graph$chrom <- str_replace(allNoshadSeg_fearon_graph$chrom, "X", "20")

allNoshadSeg_fearon_graph$sampleID <- annotationList$`Sample ID`[match(allNoshadSeg_fearon_graph$sample, annotationList$sampleStripped)]
colnames(noshadFilt_fearon_freq)[3:ncol(noshadFilt_fearon_freq)] <- annotationList$`Sample ID`[match(colnames(noshadFilt_fearon_freq)[3:ncol(noshadFilt_fearon_freq)] , annotationList$sampleStripped)]

allNoshadSeg_fearon_graph <- allNoshadSeg_fearon_graph[-which(is.na(allNoshadSeg_fearon_graph$sampleID)),]
noshadFilt_fearon_freq <- noshadFilt_fearon_freq[, -which(is.na(colnames(noshadFilt_fearon_freq)))]

colnames(allNoshadSeg_fearon_graph) <- c("ID", "chrom", "loc.start",
                                         "loc.end", "num.mark", "seg.mean")


fearon_cn <- allNoshadSeg_fearon_graph
fearon_freq <-  noshadFilt_fearon_freq
chromTextSpec = NULL

plotSegHeatmap(combinedSegRes2, combinedSegRes2_freq, segsize = 3)

plotSegHeatmap(fearon_cn, fearon_freq, segsize = 3)


### use freq data for features in phylogeny
### might not be much different if I use gene data

fearon_gene_binary <- coad_abs_geneMat
fearon_gene_binary <- round2(fearon_gene_binary)

i <- unique(allFearonAnno$bioSampGroup)[23]
for (i in unique(allFearonAnno$bioSampGroup)) {
  
  tmp <- fearon_gene_binary[grep(i, rownames(fearon_gene_binary)),]
  
  if (is.null(nrow(tmp))) {
    next()
  } else if (dim(tmp)[1] < 3) {
    next()
  }

  
  tmpCount <- tmp
  tmpCount[tmpCount != 0] <- 1
  phyloVars <- colSums(tmpCount)
  phyloVars2 <- which(phyloVars < dim(tmp)[1] & phyloVars > 1)
  tmp2 <- tmp[, which(colnames(tmp) %in% names(phyloVars2))]
  # tmp2 <- apply(tmp2, 2, as.factor)
  distMat <- dist(tmp2)
  tmpTree <- nj(distMat)
  tmpTree$edge.length <- signif(tmpTree$edge.length, digits = 3)
  tmpTree$tip.label <- paste0(tmpTree$tip.label, "(", annoTab2$Pathology[match(tmpTree$tip.label, rownames(annoTab2))],")")
  pdf(file = paste0("/mnt/DATA5/tmp/kev/misc/phylograms/",i , ".pdf"), width = 12, height = 4)
  plot(tmpTree, main = i)
  edgelabels(tmpTree$edge.length, bg="black", col="white", font=2)
  dev.off()
}



