source("/home/kevhu/scripts/20210802syntenyFunctions.R")
library(ggdendro)
library(gridExtra)

annotationList <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20211207yingSecondSet.xlsx")
annotationList$sampleStripped <- str_remove(nameStripper(annotationList$`#`), "\\-")
annotationList <- annotationList[-which(annotationList$sampleStripped == "efd36"),]

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

fdpFilt <- which(combinedVars$FDP > 150)
faoFilt <- which(combinedVars$FAO > 15)
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

### cn 

mouseNormal <- c("EF_D03_MG_X14", "MG_18X50", "MG_21X53", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")


fearonCn1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-76-MG_test1_255_185_ApcTrp53/cnMatrix_gene.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

fearonCn2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_ApcTrp53/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
fearonCn2 <- fearonCn2[,-which(colnames(fearonCn2) %in% mouseNormal)]

fearonCn3 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-156-MG_Fearon_20210809_374_382_ApcTrp53/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
fearonCn3 <- fearonCn3[,-which(colnames(fearonCn3) %in% mouseNormal)]

fearonCn4 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_ApcTrp53/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
fearonCn4 <- fearonCn4[,-which(colnames(fearonCn4) %in% mouseNormal)]


fearonCombinedCn <- cbind(fearonCn1, fearonCn2[,2:ncol(fearonCn2)], fearonCn3[,2:ncol(fearonCn3)], fearonCn4[,2:ncol(fearonCn4)])

tmpDf <- melt(fearonCombinedCn)
tmpDf <- tmpDf[grep("Del", tmpDf$Gene),]
delDf <- dcast(tmpDf, variable ~ Gene)
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


annotationList_good <- annotationList[-which(annotationList$Qual == "Bad"),]
cnMat <- fearonCombinedCn[,which(colnames(fearonCombinedCn) %in% annotationList_good$sampleStripped)]
rownames(cnMat) <- fearonCombinedCn$Gene

### heatmap

combinedCalls1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-76-MG_test1_255_185_ApcTrp53/combinedCalls.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
combinedCalls2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_ApcTrp53/combinedCalls.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
combinedCalls3 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-156-MG_Fearon_20210809_374_382_ApcTrp53/combinedCalls.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
combinedCalls4 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_ApcTrp53/combinedCalls.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

allcombinedCalls <- rbind(combinedCalls1, combinedCalls2, combinedCalls3, combinedCalls4)
allcombinedCalls$Sample <- str_remove(nameStripper(allcombinedCalls$Sample), "x.*")

geneNames <- rownames(cnMat)


### really poor q-values for event he deletions in this data set .... going to filter just by cnr

# yingSigMat <- NULL
# for (i in unique(allcombinedCalls$Sample)) {
#   tmpDf <- allcombinedCalls[which(allcombinedCalls$Sample == i),]
#   tmpDf <- tmpDf[match(geneNames,tmpDf$Gene),]
#   tmpVector <- ifelse(10^tmpDf$Log10PValue < 0.05, 1, 0)
#   yingSigMat <- rbind(yingSigMat, tmpVector)
# }
# 
# rownames(yingSigMat) <- unique(allcombinedCalls$Sample)
# colnames(yingSigMat) <- geneNames
# yingSigMat <- t(yingSigMat)
# 
# sigIdx <- match(colnames(cnMat), colnames(yingSigMat))
# yingSigMat_samps <- yingSigMat[,sigIdx]

xGenes <- c("Bcor", "Enox2", "Ar", "Atrx",
            "Diaph2", "Frmpd4")
delGenes <- c("ApcDel", "Trp53Del")


cnMat2 <- log2(cnMat)
cnMat2[cnMat2 > 3] <- 3
cnMat2[cnMat2 < -3] <- -3
cnMat2[abs(cnMat2) < 0.2] <- 0
cnMat2 <- t(cnMat2)

rownames(cnMat2) <- annotationList_good$`Sample ID`[match(rownames(cnMat2), annotationList_good$sampleStripped)]

cnDel <- cnMat2[,which(colnames(cnMat2) %in% delGenes)]
cnMat2 <- cnMat2[,-which(colnames(cnMat2) %in% delGenes)]
cnMat2_noX <- cnMat2[,-which(colnames(cnMat2) %in% xGenes)]

annoTab2 <- data.frame("Genotype" = annotationList_good$Genotype)
rownames(annoTab2) <- annotationList_good$`Sample ID`
annoCol <- list("Genotype" = c("CDX2P-CreERT2+  CDX2 fl/fl Braf fl/+" = "darkred", "CreERT2Apcfl/+, KrasLSLG12D/+, p53R270H/ex2-10" = "grey",
                               "CreERT2Apcfl/+, Krasfl/+, p53R270H fl/ex2-10" = "lightblue", "CreERT2Apcfl/+, KrasLSLG12D/+, p532-10 -/-" = "orange",
                               "CreERT2Apcfl/+, KrasLSLG12D/+, p53ex2-10fl/+ Sox9 fl/+" = "darkblue", "CreERT2Apcfl/+, KrasLSLG12D/+" = "darkgreen",
                               "Wild type" = "white"))



heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-3,3,6/1000)


heatmap_graph <- pheatmap(mat = cnMat2_noX, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                          breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                          border_color = "black", annotation_row = annoTab2, annotation_colors = annoCol)



delDf_good <- delDf[which(delDf$Sample %in% annotationList_good$sampleStripped),]
delDf_good$Sample <- annotationList_good$`Sample ID`[match(annotationList_good$sampleStripped, delDf_good$Sample)]
delDf_good[,2:3] <- log2(delDf_good[,2:3])

delDf_good2 <- delDf_good[,2:5]
rownames(delDf_good2) <- delDf_good$Sample

delDf_good2_del <- delDf_good2[,1:2]
delDf_good2_af <- delDf_good2[,3:4]

delDf_good3 <- delDf_good
delDf_good3$anno <- annoTab2$Genotype
#delDf_good_order <- delDf_good3[heatmap_graph$tree_row$order,]

delDf_good_order$Histology <- histologyList$Histology[match(delDf_good_order$Sample, histologyList$`Tissue Name`)]
delDf_good_order$sampleStripped <- annotationList_good$sampleStripped[match(delDf_good_order$Sample, annotationList_good$`Sample ID`)]
delDf_good_order <- delDf_good_order[order(delDf_good_order$anno),]
delDf_good_order[,2:5] <- signif(delDf_good_order[,2:5], digits = 3)


write.table(delDf_good_order, "/mnt/DATA5/tmp/kev/misc/20211213yingSampleDelAf.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)


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


combinedSegRes$seg.mean[which(abs(combinedSegRes$seg.mean) < 0.2)] <- 0
combinedSegRes$seg.mean[which(combinedSegRes$q.val2 > 0.05)] <- 0
combinedSegRes$seg.mean[which(combinedSegRes$seg.mean < -2)] <- -2
combinedSegRes$seg.mean[which(combinedSegRes$seg.mean > 2)] <- 2
combinedSegRes$ID <- nameStripper(combinedSegRes$ID)
combinedSegRes$ID <- str_remove(combinedSegRes$ID, "x.*")
combinedSegRes$length <- combinedSegRes$loc.end - combinedSegRes$loc.start
combinedSegRes$seg.mean[which(combinedSegRes$length < 5e6)] <- 0

selfCalls <- read.table("/mnt/DATA5/tmp/kev/misc/20211213yingSampleDelAf.txt", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)

combinedSegRes2 <- combinedSegRes[which(combinedSegRes$ID %in% selfCalls$sampleStripped[-which(selfCalls$call_kh == "no tumor")]),]
#combinedSegRes2 <- combinedSegRes[which(combinedSegRes$ID %in% selfCalls$sampleStripped),]
combinedSegRes2$chrom <- str_replace(combinedSegRes2$chrom, "23", "20")

combinedSegRes2_form <- combinedSegRes2[c("ID", "chrom", "loc.start", "loc.end",
                                          "num.mark", "seg.mean")]
colnames(combinedSegRes2_form) <- c("sampleID", "chrom", "start.pos",
                                    "end.pos", "n.probes", "mean")
combinedSegRes2_form$n.probes <- NA

combinedSegRes2_freq <- getFreqData(combinedSegRes2_form)
combinedSegRes2$chrom <- str_replace(combinedSegRes2$chrom, "23", "20")

selfCalls2 <- selfCalls[-which(selfCalls$call_kh == "no tumor"), c("Sample", "anno", "Histology")]
colnames(selfCalls2) <- c("Sample", "Genotype", "Histology")

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
  colPalRange <- seq(-2, 2, 4/999)
  
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


plotSegHeatmap(combinedSegRes2, combinedSegRes2_freq, segsize = 3)

combinedSegRes2_noX <- combinedSegRes2[-which(combinedSegRes2$chrom == "20"),]
combinedSegRes2_noX_form <- combinedSegRes2_noX[c("ID", "chrom", "loc.start", "loc.end",
                                          "num.mark", "seg.mean")]
colnames(combinedSegRes2_noX_form) <- c("sampleID", "chrom", "start.pos",
                                    "end.pos", "n.probes", "mean")
combinedSegRes2_noX_form$n.probes <- NA
combinedSegRes2_noX_freq <- getFreqData(combinedSegRes2_noX_form)

combinedSegRes2_noX$ID <- selfCalls$Sample[match(combinedSegRes2_noX$ID, selfCalls$sampleStripped)]
colnames(combinedSegRes2_noX_freq)[3:ncol(combinedSegRes2_noX_freq)] <- selfCalls$Sample[match(colnames(combinedSegRes2_noX_freq)[3:ncol(combinedSegRes2_noX_freq)],
                                                                                               selfCalls$sampleStripped)]

combinedGraph <- plotSegHeatmap(combinedSegRes2_noX, combinedSegRes2_noX_freq, segsize = 1)
selfCalls2 <- selfCalls2[match(treeOrder, selfCalls2$Sample),]

### order samps in selfCalls with tree-order prior to making table of colors and rects

annoTab_coords  <- NULL
for (i in 2:(dim(selfCalls2)[2])) {
  sampVar <- selfCalls2[,i]
  tmpDf <- data.frame("xpos" = 1.2 + 0.2 * (i-2), 
                      "ypos" = seq(1.2, dim(selfCalls2)[1] * .2 + 1, .2),
                      "label"  = sampVar)
  annoTab_coords <- rbind(annoTab_coords, tmpDf)
}

dfAnnoColIdx <- data.frame("labelNames" = c(unique(selfCalls2$Genotype), unique(selfCalls2$Histology)), 
                           "colors" = c("chartreuse3", "yellow", "lightblue", "darkblue", "darkmagenta","grey",
                                        "darkgreen", "darkred"))

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

### below was for developing code for heatmap with clustering - and aligning clustering dendogram with heatmap
### problems (1) some small blank (white areas when seg val is zero)
# 
# 
# combinedSegResDistMat <- combinedSegRes2_freq[,3:ncol(combinedSegRes2_freq)]
# rownames(combinedSegResDistMat) <- paste0("chr",combinedSegRes2_freq$chr, ":",combinedSegRes2_freq$pos)
# combinedSegResDistMat <- t(combinedSegResDistMat)
# model <- hclust(dist(combinedSegResDistMat), method = "average")
# dhc <- as.dendrogram(model)
# ddata <- dendro_data(dhc, type = "rectangle")
# p <- ggplot(segment(ddata)) + 
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
#   coord_flip() + 
#   scale_y_reverse(expand = c(0,0)) +
#   scale_x_continuous(expand = expansion(mult = c(0, 0), 
#                                            add = c(0.5, 1.5))) +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         panel.border = element_blank())
# p
# 
# treeOrder <- ddata$labels$label
# treeGrob <- ggplotGrob(p)
# 
# xlabels <- unique(df_cn2$id)
# hline <- unique(df_cn2$ypos) - 0.1
# hline <- c(hline, max(hline) + 0.2)
# chromTextdf$ypos <- max(df_cn2$ypos) + 0.2
# 
# df_cn2
# 
# testHeat <- ggplot(df_cn2) +
#   geom_rect(aes(xmin = xstart, xmax = xend,
#                 ymin = ypos - 0.1, ymax = ypos + 0.1), fill = df_cn2$col) + 
#   geom_vline(xintercept = chromBreak) +
#   geom_hline(yintercept = hline, color = "grey") + 
#   scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
#   scale_y_continuous(breaks=seq(min(df_cn2$ypos), max(df_cn2$ypos), 0.2),
#                      labels = xlabels, expand = expansion(mult = c(0, 0.02))) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()) +
#   geom_text(data = chromTextdf, aes(x = xpos, y = ypos, label = chrom), size = 2.5)
# 
# heatmapGrob <- ggplotGrob(testHeat)
# 
# grobLayout <- rbind(c(1,rep(2,19)),
#                     c(1,rep(2,19)),
#                     c(1,rep(2,19)),
#                     c(1,rep(2,19)))
# grid.arrange(treeGrob, heatmapGrob,
#              layout_matrix=grobLayout)
# 
