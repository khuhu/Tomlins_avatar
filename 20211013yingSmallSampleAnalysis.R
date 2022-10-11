cnDf1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-76-MG_test1_255_185_3norm/cnMatrix_gene.txt",
                    sep = "\t", stringsAsFactors = FALSE, header = TRUE)


sampleNames <- c("MG_28X60", "MG_29X61", "MG_30X62", "MG_31X63", "MG_32X64",
                 "MG_33X65", "MG_34X66", "MG_35X67")

cnDf2 <- cnDf1[, which(colnames(cnDf1) %in% sampleNames)]
cnDf2 <- cbind("Gene" = cnDf1$Gene,cnDf2)

cnDf1_1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_3norm/cnMatrix_gene.txt",
                    sep = "\t", stringsAsFactors = FALSE, header = TRUE)

sampleNames2 <- c("EF_D91x66", "EF_D92x67", "EF_D93x68")
cnDf2_1 <- cnDf1_1[, which(colnames(cnDf1_1) %in% sampleNames2)]

cnDf1_2 <-  read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_3norm/cnMatrix_gene.txt",
                       sep = "\t", stringsAsFactors = FALSE, header = TRUE)
sampleNames3 <- c("EF_D03_MG_X14")
cnDf2_2 <- cnDf1_2[, which(colnames(cnDf1_2) %in% sampleNames3),]

cnDf3 <- cbind(cnDf2, cnDf2_1, cnDf2_2)

cnDf_mat <- as.matrix(cnDf3[, 2:ncol(cnDf3)])
rownames(cnDf_mat) <- cnDf3$Gene
colnames(cnDf_mat)[12] <- c("EF_D03_MG_X14") 

tcVec <- 1 - cnDf_mat[which(rownames(cnDf_mat) == "ApcDel"),]
tcVec[tcVec < 0.5] <- 0.5
tcVec[9:12]  <- 0

for (i in 1:8) {
  tmpVector <- cnDf_mat[,1]
  tmpVector[tmpVector > 1] <- tmpVector[tmpVector > 1] / tcVec[i]
  tmpVector[tmpVector < 1] <- tmpVector[tmpVector < 1] * tcVec[i]
  cnDf_mat[,i] <- tmpVector
}

cnDf_mat <- log2(cnDf_mat)
cnDf_mat[cnDf_mat > 3] <- 3
cnDf_mat[cnDf_mat < -3] <- -3
cnDf_mat[cnDf_mat > log2(2/3) & cnDf_mat < log2(4/3)] <- 0

### need to get combined call matrices

combinedCalls1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-76-MG_test1_255_185_3norm/combinedCalls.txt",
                    sep = "\t", stringsAsFactors = FALSE, header = TRUE)

combinedCalls2  <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_3norm/combinedCalls.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE)

combinedCalls3  <-  read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_3norm/combinedCalls.txt",
                       sep = "\t", stringsAsFactors = FALSE, header = TRUE)

allcombinedCalls <- rbind(combinedCalls1, combinedCalls2, combinedCalls3)

geneNames <- cnDf1$Gene
#geneNames <- geneNames[-grep("Del", geneNames)]
yingSigMat <- NULL
for (i in unique(allcombinedCalls$Sample)) {
  tmpDf <- allcombinedCalls[which(allcombinedCalls$Sample == i),]
  tmpDf <- tmpDf[match(geneNames,tmpDf$Gene),]
  tmpVector <- ifelse(10^tmpDf$Log10QValue < 0.05, 1, 0)
  yingSigMat <- rbind(yingSigMat, tmpVector)
}

rownames(yingSigMat) <- unique(allcombinedCalls$Sample)
colnames(yingSigMat) <- geneNames
yingSigMat <- t(yingSigMat)

sigIdx <- match(colnames(cnDf_mat), colnames(yingSigMat))

yingSigMat_samps <- yingSigMat[,sigIdx]
yingSigMat_samps <- yingSigMat_samps[-116,]
cnDf_mat_del <- cnDf_mat[116,]
cnDf_mat <- cnDf_mat[-116,]
cnDf_mat[yingSigMat_samps == 0] <- 0


heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-3,3,6/1000)
cnDf_mat <- t(cnDf_mat)

# rownames(cnDf_mat) <- c("23380-PC1", "23388-PC1", "23468-PC1", "23505-PC1",
#                         "24074-PC1", "24094-PC1", "24144-PC1", "24172-PC1",
#                         "EF-D91", "EF-D92", "EF-D93", "EF-D03")

heatmap_genes <- pheatmap(mat = cnDf_mat, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                                 breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                                 border_color = "black")


### variants
samps <- rownames(cnDf_mat)

fearonVar1 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-76-MG_test1_255_185_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

fearonVar2 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_anno.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)

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

### since the copy-number isn't showing too much .. adding vars
cnDf_mat_vars <- cnDf_mat
cnDf_mat_vars <- cbind(cnDf_mat, matrix(0, nrow = nrow(cnDf_mat),ncol = 2))
colnames(cnDf_mat_vars)[129:130] <- c("Mtor:nonsyn", "Ccnd1:stopgain")
cnDf_mat_vars[c(1,2,8), 129] <- c(0.47, 0.52, 0.40)
cnDf_mat_vars[c(5,7), 130] <- c(0.26, 0.33)

#cnDf_mat_vars <- cbind(cnDf_mat_vars, cnDf_mat_del)
#colnames(cnDf_mat_vars)[131] <- c("ApcDel")

heatmap_genes_vars <- pheatmap(mat = cnDf_mat_vars, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                          breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                          border_color = "black")

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/202101013yingHeatmap_genes_vars.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_genes_vars
dev.off()



