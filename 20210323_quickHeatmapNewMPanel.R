library(pheatmap)

cnDf <- read.table("/mnt/DATA5/tmp/kev/newMouse2/cnMatrix_gene.txt",
                   sep = "\t", stringsAsFactors = FALSE, header = TRUE)

combinedCalls <- read.table("/mnt/DATA5/tmp/kev/newMouse2/combinedCalls.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = TRUE)

crossTableIdx <- readxl::read_xlsx("/home/kevhu/data/20201207annotations.xlsx")
crossTableIdx <- crossTableIdx[,1:3]
crossTableIdx$old_name2 <- tolower(crossTableIdx$old_name)
crossTableIdx$old_name2 <- str_remove(crossTableIdx$old_name2 , "[[:punct:]]")

newPanelCalls_names <- as.numeric(str_remove(str_remove(colnames(cnDf)[2:ncol(cnDf)], "MG_"), "X.*"))
newPanelCalls_names2 <- newPanelCalls_names

newPanelCalls_names_cc <- as.numeric(str_remove(str_remove(combinedCalls$Sample, "MG_"), "X.*"))
newPanelCalls_names_cc2 <- newPanelCalls_names_cc

for (i in unique(newPanelCalls_names)) {
  newPanelCalls_names2[which(newPanelCalls_names %in% i)] <- crossTableIdx$old_name2[which(crossTableIdx$mg_id %in% i)]
  newPanelCalls_names_cc2[which(newPanelCalls_names_cc2 %in% i)] <- crossTableIdx$old_name2[which(crossTableIdx$mg_id %in% i)]
}

colnames(cnDf)[2:ncol(cnDf)] <- newPanelCalls_names2
cnDf2 <- cnDf[,-grep("2611n", colnames(cnDf))]

combinedCalls$Sample <- newPanelCalls_names_cc2
combinedCalls2 <- combinedCalls[-grep("2611n", combinedCalls$Sample),]


heatMapCol <- colorRampPalette(c("darkblue","white","darkred"))(1000)
colors.breaks <- seq(-3,3,6/1000)
mouseCNA_mat <- log2(cnDf2[, 2:ncol(cnDf2)])
rownames(mouseCNA_mat) <- cnDf$Gene

combinedCallsDf <- NULL 
geneOrder <- cnDf$Gene
for (i in colnames(mouseCNA_mat)) {
 tmpDf <- combinedCalls2[which(combinedCalls2$Sample %in% i),]
 tmpDf <- tmpDf[match(geneOrder, tmpDf$Gene),]
 tmpDf$var <- 10^(tmpDf$Log10QValue)
 combinedCallsDf <- cbind(combinedCallsDf, tmpDf$var)
}

mouseCNA_mat[combinedCallsDf > 0.05] <- 0

mouseCNA_mat[mouseCNA_mat > 3] <- 3
mouseCNA_mat[mouseCNA_mat < -3] <- -3
#mouseCNA_mat[mouseCNA_mat > 0 & mouseCNA_mat < 0.2] <- 0
#mouseCNA_mat[mouseCNA_mat < 0 & mouseCNA_mat > -0.2] <- 0
mouseCNA_mat <- t(mouseCNA_mat)


pheatmap(mouseCNA_mat,cluster_rows = TRUE,
         cluster_cols = FALSE, color = heatMapCol,
         breaks = colors.breaks, fontsize = 5,
         border_color = "black")

