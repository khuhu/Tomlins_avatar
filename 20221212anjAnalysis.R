library(stringr)
library(pheatmap)

nameStripper <- function(df){
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
}


### old bed
# ovBed <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_ApcTrp53Del.gc.bed", sep = "\t",
#                     stringsAsFactors = FALSE, header = FALSE)
# ovBed$V8[grep("ApcDel", ovBed$V8)] <- "Apc"
# write.table(ovBed, "/home/kevhu/data/bedFiles/IAD202670_167_Trp53Del.gc.bed", row.names = FALSE, col.names = FALSE,
#             sep = "\t", quote = FALSE)

sunnyAmplicons <-  read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-216-SW_MG_2_469_509/gcCorrectedCounts_matrix.txt",
                              sep = "\t", header = TRUE)

mouseBed <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed",
                       header = FALSE, stringsAsFactors = FALSE)

sunnyAmpliconsNotch1 <- sunnyAmplicons[which(sunnyAmplicons$Gene == "Notch1"),]
sunnyAmpliconsNotch1[, 2:23] <- log2(sunnyAmpliconsNotch1[, 2:23])


sunnyGeneDelTable <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-216-SW_MG_2_469_509/cnMatrix_gene.txt",
                             sep = "\t", header = TRUE)

sunnyCombinedDelTable <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-216-SW_MG_2_469_509/combinedCalls.txt",
                                  sep = "\t", header = TRUE)


sunnySegTable  <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-216-SW_MG_2_469_509/segResults.txt",
                             sep = "\t", header = TRUE)

anjGeneTable <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-264-AD22_MG_517_614/cnMatrix_gene.txt",
                           sep = "\t", header = TRUE)

anjCombinedTable <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-264-AD22_MG_517_614/combinedCalls.txt",
                           sep = "\t", header = TRUE)


anjSegTable <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-264-AD22_MG_517_614/segResults.txt",
                               sep = "\t", header = TRUE)

anjGeneDelTable <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-264-AD22_MG_517_614/cnMatrix_gene.txt",
                           sep = "\t", header = TRUE)

anjCombinedDelTable <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-264-AD22_MG_517_614/combinedCalls.txt",
                               sep = "\t", header = TRUE)

mouseNormal <- c("EF_D03_MG_X14", "MG_18X50", "MG_21X53", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

seqID <- unique(anjCombinedDelTable$Sample)[1:18]
sampleID <- c("15542-V1-W1", "15464-T2-W2", "15464-T3B-W2", "15703-T1-W1", "15704-T7-W2",
              "15804-T1-W1", "18163-KC", "15464-T3B-W2-W", "15464-T2-W2-C",
              "15704-T7-W2-U", "15704-L", "15464-L",
              "BCC-#1-Mouse-1", "BCC-#2-Mouse-1", "Liver-Mouse-1", "BCC-#10-Mouse-2", "BCC-#12-Mouse-2", "BCC-#1-Mouse-3")
nameDf <- data.frame(seqID, sampleID)
cnMat <- anjGeneDelTable[, 2:ncol(anjGeneDelTable)]
rownames(cnMat) <- anjGeneDelTable$Gene
geneNames <- rownames(cnMat)

for (i in 1:nrow(nameDf)) {
  anjCombinedDelTable$Sample[grep(nameDf$seqID[i], anjCombinedDelTable$Sample)] <- nameDf$sampleID[i]
  colnames(cnMat)[which(colnames(cnMat) == nameDf$seqID[i])] <- nameDf$sampleID[i]
}


### really poor q-values for event he deletions in this data set .... going to filter just by cnr

sigMatrix <- reshape2::dcast(anjCombinedDelTable, formula = Gene ~ Sample, value.var = "Log10QValue")
sigMatrix <- sigMatrix[match(rownames(cnMat), sigMatrix$Gene), ]
sigMatrix2 <- sigMatrix[, 2:ncol(sigMatrix)]

### setting to cnr of 1 so it'll get filtered out by the 0.2 filter
cnMat[sigMatrix2 < log10(0.05)] <- 1

xGenes <- c("Bcor", "Enox2", "Ar", "Atrx",
            "Diaph2", "Frmpd4")
delGenes <- c("Trp53Del")
badGenes <- c("Esr1", "Aurkc")

cnMat2 <- log2(cnMat)
cnMat2[cnMat2 > 2] <- 2
cnMat2[cnMat2 < -2] <- -2
cnMat2[abs(cnMat2) < 0.2] <- 0 
cnMat2 <- t(cnMat2)

cnDel <- cnMat2[, -which(colnames(cnMat2) %in% delGenes)]
cnMat2_noX <- cnMat2[,-which(colnames(cnMat2) %in% xGenes)]
cnMat2_noX_del <- cnMat2[-which(rownames(cnMat2) %in% mouseNormal), -which(colnames(cnMat2) %in% c(xGenes, delGenes, badGenes))]

anjNormal <- c("15704-L", "15464-L", "Liver-Mouse-1")
anjBcc <- c("BCC-#1-Mouse-1", "BCC-#2-Mouse-1", "BCC-#10-Mouse-2", "BCC-#12-Mouse-2", "BCC-#1-Mouse-3")
cnMat2_noX_del <- cnMat2_noX_del[-which(rownames(cnMat2_noX_del) %in% anjNormal), ]

cnMat2_noX_del_bcc <- cnMat2_noX_del[which(rownames(cnMat2_noX_del) %in% anjBcc), ]
cnMat2_noX_del_mcc <- cnMat2_noX_del[-which(rownames(cnMat2_noX_del) %in% anjBcc), ]



### need to separate bcc with mcc, combine new bccs with sunny wong samples

cnMat_sw <- sunnyGeneDelTable[, 2:ncol(sunnyGeneDelTable)]
rownames(cnMat_sw) <- sunnyGeneDelTable$Gene

sigMatrix_sw <- reshape2::dcast(sunnyCombinedDelTable, formula = Gene ~ Sample, value.var = "Log10QValue")
sigMatrix_sw <- sigMatrix_sw[match(rownames(cnMat_sw), sigMatrix_sw$Gene), ]
sigMatrix2_sw <- sigMatrix_sw[, 2:ncol(sigMatrix_sw)]


cnMat_sw[sigMatrix2_sw < log10(0.05)] <- 1

cnMat2_sw <- log2(cnMat_sw)
cnMat2_sw[cnMat2_sw > 2] <- 2
cnMat2_sw[cnMat2_sw < -2] <- -2
cnMat2_sw[abs(cnMat2_sw) < 0.2] <- 0 
cnMat2_sw <- t(cnMat2_sw)

cnDel_sw <- cnMat2_sw[, -which(colnames(cnMat2_sw) %in% delGenes)]
cnMat2_noX_sw <- cnMat2_sw[,-which(colnames(cnMat2_sw) %in% xGenes)]
cnMat2_noX_del_sw <- cnMat2_sw[-which(rownames(cnMat2_sw) %in% mouseNormal), -which(colnames(cnMat2_sw) %in% c(xGenes, delGenes, badGenes))]

sw_normal <- paste0("SW_MG_", c("2_X82","4_X84", "6_X86", "8_X88", "10_X90", "13_X93", "15_X95"))
cnMat2_noX_del_sw <- cnMat2_noX_del_sw[-which(rownames(cnMat2_noX_del_sw) %in% sw_normal), ]

cnMat2_noX_del_bcc_combined <- rbind(cnMat2_noX_del_bcc, cnMat2_noX_del_sw)


heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-2, 2, 4/1000)


heatmap_graph_mcc <- pheatmap(mat = cnMat2_noX_del_mcc, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                          breaks = colors.breaks, fontsize = 8, cellwidth = 12, cellheight = 25, silent = FALSE,
                          border_color = "black")

heatmap_graph_bcc <- pheatmap(mat = cnMat2_noX_del_bcc_combined, cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                              breaks = colors.breaks, fontsize = 8, cellwidth = 12, cellheight = 25, silent = FALSE,
                              border_color = "black")

dev.off()
png("/mnt/DATA5/tmp/kev/misc/20230103anjMccHeatmap.png", width = 1800, height = 1100)
heatmap_graph_mcc
dev.off()

dev.off()
png("/mnt/DATA5/tmp/kev/misc/20230103anjBccHeatmap.png", width = 1800, height = 1100)
heatmap_graph_bcc
dev.off()

