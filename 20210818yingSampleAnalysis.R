fearonCn1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonCn2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonVar1 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-156-MG_Fearon_20210809_374_382_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonVar2 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonSeg1 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/segResults.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonSeg2 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/segResults.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonSeg3 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384/segResults.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

annoTable <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20210817samplesForYingAnno.xls")

combinedCalls1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/combinedCalls.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = TRUE)

combinedCalls2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/combinedCalls.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = TRUE)
### cn
samps <- c(paste0("D", 10:33), paste0("D0", c(8,9)))

fearonSeg_combined <- rbind(fearonSeg1, fearonSeg2)
fearonSeg_combined <- fearonSeg_combined[grep(paste0(samps, collapse = "|"),
                                               fearonSeg_combined$ID),]

fearonCombinedCalls <- rbind(combinedCalls1, combinedCalls2)
fearonCombinedCalls <- fearonCombinedCalls[grep(paste0(samps, collapse = "|"),
                                                fearonCombinedCalls$Sample),]
fearonCombinedCalls$Sample <- str_remove(fearonCombinedCalls$Sample, "\\_MG.*")
fearonCombinedCalls$Sample <- str_remove(fearonCombinedCalls$Sample, "\\_X.*")
fearonCombinedCalls <- fearonCombinedCalls[order(fearonCombinedCalls$Sample),]
fearonCombinedCalls <- fearonCombinedCalls[-grep("Del",fearonCombinedCalls$Gene),]

geneNames <- fearonCn1$Gene
geneNames <- geneNames[-grep("Del", geneNames)]
fearonSigMat <- NULL
for (i in unique(fearonCombinedCalls$Sample)) {
  tmpDf <- fearonCombinedCalls[which(fearonCombinedCalls$Sample == i),]
  tmpDf <- tmpDf[match(geneNames,tmpDf$Gene),]
  tmpVector <- ifelse(10^tmpDf$Log10QValue < 0.05, 1, 0)
  fearonSigMat <- rbind(fearonSigMat, tmpVector)
}
rownames(fearonSigMat) <- unique(fearonCombinedCalls$Sample)
colnames(fearonSigMat) <- geneNames

fearonCn_comb <- cbind(fearonCn1, fearonCn2)
fearonCn_comb <- fearonCn_comb[,grep(paste0(samps, collapse = "|"), colnames(fearonCn_comb))]
colnames(fearonCn_comb) <- str_remove(colnames(fearonCn_comb), "\\_MG.*")
colnames(fearonCn_comb) <- str_remove(colnames(fearonCn_comb), "\\_X.*")
fearonCn_comb <- log2(fearonCn_comb)
rownames(fearonCn_comb) <- fearonCn1$Gene
fearonCn_comb <- fearonCn_comb[, order(colnames(fearonCn_comb))]

tc <- 1 - 2^fearonCn_comb[which(rownames(fearonCn_comb) == "Trp53Del"),]
tc_labels <- tc
tc_labels[c(5,6,13)] <- 0

tc[c(5,6,13)] <- 1
tc[which(tc < 0.5)] <- 0.5
fearonCn_comb <- sweep(fearonCn_comb, 2, as.numeric(tc), "/")
fearonCn_comb[fearonCn_comb < -3] <- -3
fearonCn_comb[fearonCn_comb > 3] <- 3
dels <- fearonCn_comb[grep("Del", rownames(fearonCn_comb)),]
fearonCn_comb <- fearonCn_comb[-grep("Del", rownames(fearonCn_comb)),]
fearonCn_comb[fearonCn_comb > log2(2/3) & fearonCn_comb < log(4/3)] <- 0
fearonCn_comb <- t(fearonCn_comb)
rownames(fearonCn_comb) <- annoTable$`New Sample ID`
fearonCn_comb[fearonSigMat == 0] <- 0

xGenes <- c("Frmpd4", "Diaph2", "Atrx", "Med12", "Ar", "Enox2", "Bcor", "Aurkc")
fearonCn_comb <- fearonCn_comb[,-which(colnames(fearonCn_comb) %in% xGenes)]


annoTab2 <- data.frame("Type" = annoTable$Type)
rownames(annoTab2) <- annoTable$`New Sample ID`
annoCol <- list("Genotype" = c("cell-line" = "orange", "normal" = "white", "primary" = "darkred",
                               "lung.met" = "darkblue", "liver.met" = "black"))

heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-3,3,6/1000)

heatmap_graph <- pheatmap(mat = fearonCn_comb, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                          breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                          border_color = "black", annotation_row = annoTab2, annotation_colors = annoCol)

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210917yingHeatmap_genes.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_graph
dev.off()

dels <- t(dels)
heatmap_graph_dels <- pheatmap(dels[heatmap_graph$tree_row$order,], cluster_rows = FALSE, cluster_cols = FALSE,
                               color = heatMapCol,breaks = colors.breaks, fontsize = 5, cellwidth = 5,
                               cellheight = 10,silent = FALSE, border_color = "black")


heatMapCol_tc <- colorRampPalette(c("#FFFFFF","#FFA500"))(100)
colors.breaks_tc <- seq(0,1,1/100)
tc_labels <- t(tc_labels)
heatmap_graph_tc <- pheatmap(tc_labels[heatmap_graph$tree_row$order,], cluster_rows = FALSE, cluster_cols = FALSE,
                               color = heatMapCol_tc, breaks = colors.breaks_tc, fontsize = 5, cellwidth = 5,
                               cellheight = 10,silent = FALSE, border_color = "black")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210917yingHeatmap_del.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_graph_dels
dev.off()

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210917yingHeatmap_tc.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_graph_tc
dev.off()


### Vars

combinedVars <- rbind(fearonVar1, fearonVar2)
mgpVars <- c("hom", "het", "unknown")
combinedVars <- combinedVars[which(combinedVars$mm10_mpgpv6_Indels == ""),]

fdpFilt <- which(combinedVars$FDP > 20)
faoFilt <- which(combinedVars$FAO > 5)
freqFilt <- which(combinedVars$AF > 0.05)
hrunFilt <- which(combinedVars$HRUN < 4)
strandRatio <- intersect(which(combinedVars$FSAF/combinedVars$FSAR > 0.2),
                         which(combinedVars$FSAF/combinedVars$FSAR < 5))
qualFilt <- which(combinedVars$QUAL >= 100)

goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, qualFilt))
combinedVars_goodsamps <- combinedVars[goodSamps,]
combinedVars_goodsamps_exon <- combinedVars_goodsamps[which(combinedVars_goodsamps$Func.refGene == "exonic"), ]

combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[which(combinedVars_goodsamps_exon$name == "."),]                                                  
combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[-grep("nonframeshift",
                                                                 combinedVars_goodsamps_exon$ExonicFunc.refGene),]  
combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[-grep("^synonymous SNV",
                                                                 combinedVars_goodsamps_exon$ExonicFunc.refGene),]  
combinedVars_goodsamps_exon <- combinedVars_goodsamps_exon[grep(paste0(samps, collapse = "|"), combinedVars_goodsamps_exon$Sample),]

