### making BPRNTA bed file for copy-number data 


# cnRes <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/gcCorrectedCounts_matrix.txt",
#                     sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
# 
# cnRes <- cnRes[which(cnRes$Gene == "Apc" | cnRes$Gene == "Pten" |
#                        cnRes$Gene == "Brca1" | cnRes$Gene == "Trp53" |
#                        cnRes$Gene == "Nf1" | cnRes$Gene == "Rb1"),]
# cnRes <- cnRes[, c(1, 44:47,34, 35, 37, 59:66)]
# 
# apc <- c("18:34310584,34310876") # 4 probes
# pten <- c("19:32799776,32800016") # probes
# 
# 
# apcDel # 3904:3907
# ptenDel # 4063:4066
# brcaDel # 2874:2919
# trp53Del # 2569:2584
# nf1Del # 2684:2687
# rb1Del # 3320:3321
# 
# tmpBed <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.gc.noChr.bed", sep = "\t")
# 
# tmpBed2 <- tmpBed
# tmpBed2$V8[3904:3907] <- "ApcDel"
# tmpBed2$V8[4063:4066] <- "PtenDel"
# tmpBed2$V8[2874:2919] <- "Brca1Del"
# tmpBed2$V8[2569:2584] <- "Trp53Del"
# tmpBed2$V8[2684:2687] <- "Nf1Del"
# tmpBed2$V8[3320:3321] <- "Rb1Del"
# 
# 
# tmpBed2$V1 <- paste0("chr", tmpBed2$V1)
# write.table(tmpBed2, "/home/kevhu/data/bedFiles/IAD202670_167_BPRNTADel.gc.bed", sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# 
# 


firstUpper <- function(gene){
  firstLetter <- toupper(substr(gene, start = 1, stop = 1))
  restOfGene <- tolower(substr(gene, start = 2, stop = nchar(gene)))
  res <- paste0(firstLetter, restOfGene)
  return(res)
}


library(stringr)

tcDf <- read.table("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

### KC-01 to 05 is in mixe run with frearon samples

statsTab1 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-138-MG_cho_20210621_354_343/plugin_out/coverageAnalysis_out.668/Auto_MG_cho_20210621_eros_343.bc_summary.xls",
           sep = "\t", header = TRUE, stringsAsFactors = FALSE)

statsTab2 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/plugin_out/coverageAnalysis_out.681/Auto_MG_cho_202106_3TS_eros_349.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

statsTab3 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-142-MG_cho_20210701_357_353/plugin_out/coverageAnalysis_out.693/Auto_MG_cho_20210701_eros_353.bc_summary.xls",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

allStats <- rbind(statsTab1, statsTab2, statsTab3)
allStats$On.Target <- as.numeric(str_remove(allStats$On.Target, "%"))
allStats$Uniformity <- as.numeric(str_remove(allStats$Uniformity, "%"))

badSamps <- allStats$Sample.Name[which(allStats$On.Target < 80 | allStats$Mean.Depth < 100 | allStats$Uniformity < 80)]


annoTab <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/2021078kathyMouseSelfAnno.xlsx")
annoTab_addtional <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20210717rongAdditionalAnnoSelf.xlsx")

cnGeneDf1 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-138-MG_cho_20210621_354_343/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
cnGeneDf1 <- cnGeneDf1[,1:49]
if(length(which(colnames(cnGeneDf1) %in% badSamps)) > 0){
  cnGeneDf1 <- cnGeneDf1[,-which(colnames(cnGeneDf1) %in% badSamps)]
}

colnames(cnGeneDf1) <- str_remove(colnames(cnGeneDf1), "_.*")
colnames(cnGeneDf1) <- str_remove(colnames(cnGeneDf1), "[[:punct:]]")
colnames(cnGeneDf1) <- str_replace_all(colnames(cnGeneDf1), " ", "")
colnames(cnGeneDf1) <- str_remove(colnames(cnGeneDf1), "O")


cnGeneDf2 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
cnGeneDf2 <- cnGeneDf2[,1:32]
if(length(which(colnames(cnGeneDf2) %in% badSamps)) > 0){
  cnGeneDf2 <- cnGeneDf2[,-which(colnames(cnGeneDf2) %in% badSamps)]
}

colnames(cnGeneDf2) <- str_remove(colnames(cnGeneDf2), "X.*")
colnames(cnGeneDf2) <- str_remove(colnames(cnGeneDf2), "[[:punct:]]")
colnames(cnGeneDf2) <- str_replace_all(colnames(cnGeneDf2), " ", "")
colnames(cnGeneDf2) <- str_remove(colnames(cnGeneDf2), "O")

cnGeneDf3 <- read.table("/mnt/DATA6/mouseData/customCnRuns/Auto_user_AUS5-142-MG_cho_20210701_357_353/cnMatrix_gene.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
cnGeneDf3 <- cnGeneDf3[,1:56]

if(length(which(colnames(cnGeneDf3) %in% badSamps)) > 0){
  cnGeneDf3 <- cnGeneDf3[,-which(colnames(cnGeneDf3) %in% badSamps)]
}

colnames(cnGeneDf3) <- str_remove(colnames(cnGeneDf3), "_X.*")
colnames(cnGeneDf3) <- str_remove(colnames(cnGeneDf3), "[[:punct:]]")
colnames(cnGeneDf3) <- str_replace_all(colnames(cnGeneDf3), " ", "")
colnames(cnGeneDf3) <- str_remove(colnames(cnGeneDf3), "O")

combinedCnDf <- cbind(cnGeneDf1, cnGeneDf2[,2:ncol(cnGeneDf2)], cnGeneDf3[,2:ncol(cnGeneDf3)])

annoNames <- paste0(annoTab$mouse_id, annoTab$side, "T")
annoNames[1] <- "12167met"
annoNames[which(annoNames == "14390RT")] <- c("14390RTE", "14390RTS")
annoNames[which(annoNames == "14433metT")] <- "14433met"

additoinalNames <- str_remove(annoTab_addtional$mouse_id, "_.*")
additoinalNames <- str_remove(additoinalNames, "[[:punct:]]")
additoinalNames <- str_replace_all(additoinalNames, " ", "")
additoinalNames <- str_remove(additoinalNames, "O")

annoNames <- c(annoNames, additoinalNames)

colnames(combinedCnDf)[which(colnames(combinedCnDf) == "12167")] <- "12167met"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "12401RT2")] <- "12401RT"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "14433MT")] <- "14433met"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "14433LTS")] <- "14433LT"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "2027LTS")] <- "2027LT"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "14109RT")] <- "14109LT"
colnames(combinedCnDf)[which(colnames(combinedCnDf) == "2016RT")] <- "2016LT"

allSampNames <- colnames(combinedCnDf)

annoNames[which(annoNames %in% allSampNames)]

combinedCnRong <- combinedCnDf[ , which(colnames(combinedCnDf) %in% annoNames)]
combinedCnRong <- cbind("Gene" = combinedCnDf$Gene, combinedCnRong)

annoNames2 <- annoNames[which(annoNames %in% colnames(combinedCnRong))]
combinedCnRong <- cbind("Genes" = combinedCnRong[,1], combinedCnRong[,match(annoNames2, colnames(combinedCnRong))])

### samples get filtered out by the copy-number pipeline or stats above - 3 removed
### ask Albert about these samples first 12401RT2 is 12401RT1, 14109RT == 14109RT, 14433LTS == 14433LT?
### 2016LT == 2016RT bad?, 2027LTS == 2027LT?, 2013RT bad?



delRows <- combinedCnRong[grep("Del", combinedCnRong$Gene), 2:ncol(combinedCnRong)]
combinedCnRong <- combinedCnRong[-grep("Del", combinedCnRong$Gene),]
combinedCnDf$Gene <- firstUpper(combinedCnDf$Gene)
cnGeneNames <- c(combinedCnRong$Genes, "Trp53Del", "Nf1Del", "Brca1Del",
                 "Rb1Del","ApcDel","PtenDel")
combinedCnRong <- combinedCnRong[,-1]
combinedCnRong <- rbind(combinedCnRong, delRows)
combinedCnRong[1:nrow(combinedCnRong),] <- lapply(combinedCnRong[1:nrow(combinedCnRong),], function(x) log2(x))
rownames(combinedCnRong) <- cnGeneNames
ChoCnMat <- t(combinedCnRong)


ChoCnMat[ChoCnMat < -3] <- -3
ChoCnMat[ChoCnMat > 3] <- 3

#ChoCnMat[ChoCnMat > (log2(1/2) *.8) & ChoCnMat < (log2(3/2) * .8)] <- 0
### use combinedCall sets to filtering based on p-value for individual gene calls - possibly
### use current 1 copy


annoTab2 <- data.frame("Genotype" = annoTab$geno, "Type" = annoTab$type,
                       "Brca1_1" = annoTab$brca1_1, "Brca1_2" = annoTab$brca1_2,
                       "Trp53_1" = annoTab$trp53_1...7, "Trp53_2" = annoTab$trp53_1...8,
                       "Rb1_1" = annoTab$rb1_1, "Rb1_2" = annoTab$rb1_2,
                       "Nf1_1" = annoTab$nf1_1, "Nf1_2" = annoTab$nf1_2,
                       "Pten_1" = annoTab$pten_1, "Pten_2" = annoTab$pten_2,
                       "Apc_1" = annoTab$apc_1, "Apc_2" = annoTab$apc_2, stringsAsFactors = FALSE)

tmpannoTab_addtional <- data.frame("Genotype" = annoTab_addtional$geno, "Type" = annoTab_addtional$type,
                         "Brca1_1" = annoTab_addtional$brca1_1, "Brca1_2" = annoTab_addtional$brca1_2,
                         "Trp53_1" = annoTab_addtional$trp53_1...6, "Trp53_2" = annoTab_addtional$trp53_1...7,
                         "Rb1_1" = annoTab_addtional$rb1_1, "Rb1_2" = annoTab_addtional$rb1_2,
                         "Nf1_1" = annoTab_addtional$nf1_1, "Nf1_2" = annoTab_addtional$nf1_2,
                         "Pten_1" = annoTab_addtional$pten_1, "Pten_2" = annoTab_addtional$pten_2,
                         "Apc_1" = annoTab_addtional$apc_1, "Apc_2" = annoTab_addtional$apc_2, stringsAsFactors = FALSE)

annoTab2 <- rbind(annoTab2, tmpannoTab_addtional)

annoNames[-which(annoNames %in% annoNames2)]
annoTab2 <- annoTab2[which(annoNames %in% annoNames2),]
annoTab2 <- data.frame(lapply(annoTab2, as.character))

annoTab2 <- annoTab2[-which(duplicated(annoNames2)),]
ChoCnMat <- ChoCnMat[-which(duplicated(annoNames2)),]
annoNames2 <- annoNames2[-which(duplicated(annoNames2))]

rownames(annoTab2) <- annoNames2
rownames(ChoCnMat) <- annoNames2


#get annoCol from old paper code
annoCol <- list("Genotype" = c("BPRN" = "red", "PRN" = "darkorange",
                               "BPP" = "darkorange3", "BPN" = "darkred",
                               "AP" = "lightgreen","unk" = "grey"),
                "Type"  = c("hgsc" = "red", "ehgsc" = "red4", "met" = "black",
                            "mmmt" = "dodgerblue", "mmmte" = "lightblue1", "mmmts" = "lightblue3",
                            "endoCa" = "yellow", "pap. endo" = "yellow4"),
                "Brca1_1" = c("1" = "green4", "0" = "white", "2" = "lawngreen"),
                "Brca1_2" = c("1" = "green4", "0" = "white", "2" = "lawngreen"),
                "Trp53_1" = c("1" =  "cyan", "2" = "blue4", "0" = "white"),
                "Trp53_2" = c("1" =  "cyan", "2" = "blue4","0" = "white"),
                "Rb1_1" = c("1" = "purple4", "0" = "white"),
                "Rb1_2" = c("1" = "purple4", "0" = "white"),
                "Nf1_1" = c("1" = "peachpuff", "0" = "white"),
                "Nf1_2" = c("1" = "peachpuff", "0" = "white"),
                "Pten_1" = c("1" = "orange3", "0" = "white"),
                "Pten_2" = c("1" = "orange3", "0" = "white"),
                "Apc_1" = c("1" = "dimgrey", "0" = "white"),
                "Apc_2" = c("1" = "dimgrey", "0" = "white"))


heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
colors.breaks <- seq(-3,3,6/1000)

#ChoCnMat[abs(ChoCnMat) < 0] <- 0

heatmap_graph <- pheatmap(mat = ChoCnMat[,1:(ncol(ChoCnMat) - 6)], cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                          breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                          border_color = "black", annotation_row = annoTab2, annotation_colors = annoCol)




dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210719choHeatmap_genes.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_graph
dev.off()


### create separate 

dels <- ChoCnMat[,grep("Del", colnames(ChoCnMat))]
tc <- apply(dels, 1, function(x) round(1-2^(min(x)), digits = 2))

dels <- dels[,c("Brca1Del", "Trp53Del", "Rb1Del", "Nf1Del", "PtenDel", "ApcDel")]
# this is the order for the clustering
heatmap_graph_dels <- pheatmap(dels[heatmap_graph$tree_row$order,], cluster_rows = FALSE, cluster_cols = FALSE,
                               color = heatMapCol,breaks = colors.breaks, fontsize = 5, cellwidth = 5,
                               cellheight = 10,silent = FALSE, border_color = "black")

tc_df <- data.frame("Tumor content" = tc)
rownames(tc_df) <- rownames(dels)
tcCol <- colorRampPalette(c("#FFFFFF","#CC5500"))(100)
colors.breaks.tc <- seq(-0,1,1/100)

heatmap_graph_tcvals <- pheatmap(tc_df[heatmap_graph$tree_row$order,], cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = tcCol, breaks = colors.breaks.tc, fontsize = 5, cellwidth = 5,
                                 cellheight = 10,silent = FALSE, border_color = "black")


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210719choHeatmap_dels.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_graph_dels
dev.off()


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210711choHeatmap_tc.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_graph_tcvals
dev.off()





### for tc corrected


ChoCnMat_tc <- ChoCnMat
matchingTc <- tcDf$tc[match(tolower(rownames(ChoCnMat_tc)), tcDf$sample)]
which(!is.na(matchingTc))
ChoCnMat_tc[which(!is.na(matchingTc)),] <- sweep(ChoCnMat_tc[which(!is.na(matchingTc)),], 1,
                                                 matchingTc[which(!is.na(matchingTc))], "/")

ChoCnMat_tc[ChoCnMat_tc > (log2(1/2) *.7) & ChoCnMat_tc < (log2(3/2) * .7)] <- 0


heatmap_tc_graph <- pheatmap(mat = ChoCnMat_tc[,1:(ncol(ChoCnMat_tc) - 6)], cluster_rows = TRUE, cluster_cols = FALSE, color = heatMapCol,
                             breaks = colors.breaks, fontsize = 5, cellwidth = 5, cellheight = 10, silent = FALSE,
                             border_color = "black", annotation_row = annoTab2, annotation_colors = annoCol)


dels_tc <- ChoCnMat_tc[,grep("Del", colnames(ChoCnMat_tc))]
dels_tc <- dels_tc[,c("Brca1Del", "Trp53Del", "Rb1Del", "Nf1Del", "PtenDel", "ApcDel")]

heatmap_tc_dels <- pheatmap(dels_tc[heatmap_tc_graph$tree_row$order,], cluster_rows = FALSE, cluster_cols = FALSE,
                            color = heatMapCol,breaks = colors.breaks, fontsize = 5, cellwidth = 5,
                            cellheight = 10,silent = FALSE, border_color = "black")

heatmap_tc_tcvals <- pheatmap(tc_df[heatmap_tc_graph$tree_row$order,], cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = tcCol, breaks = colors.breaks.tc, fontsize = 5, cellwidth = 3,
                                 cellheight = 8,silent = FALSE, border_color = "black")



dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210719choHeatmap_tc_genes.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_tc_graph
dev.off()

dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210719choHeatmap_tc_dels.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_tc_dels
dev.off()


dev.off()
pdf("/mnt/DATA5/tmp/kev/misc/20210719choHeatmap_tc_tcval.pdf", useDingbats = TRUE, width = 15, height = 15)
heatmap_tc_tcvals
dev.off()



### var portion 

annoFile_138 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-138-MG_cho_20210621_354_343_anno.txt",
                           sep = "\t", stringsAsFactors = FALSE, header = TRUE)
annoFile_141 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349_anno.txt",
                           sep = "\t", stringsAsFactors = FALSE, header = TRUE)
annoFile_142 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-142-MG_cho_20210701_357_353_anno.txt",
                           sep = "\t", stringsAsFactors = FALSE, header = TRUE)

combinedAnnoFile <- rbind(annoFile_138, annoFile_141, annoFile_142)
combinedAnnoFile <- combinedAnnoFile[-which(combinedAnnoFile$mm10_mpgpv6_Indels %in% c("hom","het")),]

fdpFilt <- which(combinedAnnoFile$FDP > 20)
faoFilt <- which(combinedAnnoFile$FAO > 5)
freqFilt <- which(combinedAnnoFile$AF > 0.05)
hrunFilt <- which(combinedAnnoFile$HRUN < 4)
strandRatio <- intersect(which(combinedAnnoFile$FSAF/combinedAnnoFile$FSAR > 0.2),
                         which(combinedAnnoFile$FSAF/combinedAnnoFile$FSAR < 5))
qualFilt <- which(combinedAnnoFile$QUAL > 60)
goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, hrunFilt, qualFilt))
combinedAnnoFile<- combinedAnnoFile[goodSamps,]
combinedAnnoFile_exon <- combinedAnnoFile[which(combinedAnnoFile$Func.refGene == "exonic"),]

combinedAnnoFile_exon2 <- combinedAnnoFile_exon[-grep(paste(c("^synonymous SNV", "nonframeshift"), collapse = "|"),
                                                      combinedAnnoFile_exon$ExonicFunc.refGene),]
combinedAnnoFile_exon2 <- combinedAnnoFile_exon2[order(combinedAnnoFile_exon2$ExonicFunc.refGene),]
write.table(combinedAnnoFile_exon2, "/mnt/DATA5/tmp/kev/misc/20210721choMousevariants.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

