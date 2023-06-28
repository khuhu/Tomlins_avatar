source("/home/kevhu/scripts/20210802syntenyFunctions.R")
plotSegHeatmap <- function(df_cn, df_freq, main = NULL,
                           chromTextSpec = NULL, segsize = 3, sizeVar = 6){
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
  dendo_graph <- ggplot(segment(dendo_dat)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    coord_flip() + 
    scale_y_reverse(expand = c(0,0)) +
    #scale_x_continuous(expand = expansion(mult = c(0, 0), 
    #                                      add = c(0.5, 1.5))) +
    
    scale_x_continuous(expand = expansion(mult = c(0, 0), 
                                          add = c(0.5, 1.3 + nrow(df_cn)/2000))) +
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
          axis.text.y = element_text(size=sizeVar)) +
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

### summary files used to remove 

sumFile1 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-76-MG_test1_255_185/plugin_out/coverageAnalysis_out.303/Auto_MG_test1_eros_185.bc_summary.xls",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sumFile2 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/plugin_out/coverageAnalysis_out.577/Auto_MG_EFD4_BBN_eros_304.bc_summary.xls",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sumFile3 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-138-MG_cho_20210621_354_343/plugin_out/coverageAnalysis_out.668/Auto_MG_cho_20210621_eros_343.bc_summary.xls",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sumFile4 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/plugin_out/coverageAnalysis_out.681/Auto_MG_cho_202106_3TS_eros_349.bc_summary.xls",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sumFile5 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-142-MG_cho_20210701_357_353/plugin_out/coverageAnalysis_out.693/Auto_MG_cho_20210701_eros_353.bc_summary.xls",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sumFile6 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/plugin_out/coverageAnalysis_out.757/Auto_MG_Fearon_20210809_eros_382.bc_summary.xls",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sumFile7 <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384/plugin_out/coverageAnalysis_out.761/Auto_MG_Fearon_20210809_2_eros_384.bc_summary.xls",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)

sumFiles <- rbind(sumFile1, sumFile2, sumFile3, 
                  sumFile4, sumFile5, sumFile6,
                  sumFile7)

sumFiles$Uniformity <- as.numeric(str_remove(sumFiles$Uniformity, "\\%"))
sumFiles$Sample.Name <- nameStripper(sumFiles$Sample.Name)
sumFiles$Sample.Name <- str_remove(sumFiles$Sample.Name, "\\-")
sumFiles$Sample.Name <- str_remove(sumFiles$Sample.Name, "x.*")

goodSamps <- sumFiles$Sample.Name[which(sumFiles$Uniformity > 80 & sumFiles$Mean.Depth > 250)]
### EF data

mouseNormal <- c("EF_D03_MG_X14", "MG_18X50", "MG_21X53", "MG_6X38",
                 "MG_8X40","MG_11X43", "MG_13X45")

### I should load all mouse data filt by quality and then make seprate graphs by loading in desiered anno table

segRes1 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-76-MG_test1_255_185/segResults.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
segRes2 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-138-MG_cho_20210621_354_343/segResults.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
segRes3 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-156-MG_Fearon_20210809_374_382/segResults.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
segRes4 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-157-MG_Fearon_20210809_2_375_384/segResults.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
segRes5 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/segResults.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
segRes6 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/segResults.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
segRes7 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-142-MG_cho_20210701_357_353/segResults.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)


segRes2 <- segRes2[-which(segRes2$ID %in% mouseNormal),]
segRes3 <- segRes3[-which(segRes3$ID %in% mouseNormal),]
segRes4 <- segRes4[-which(segRes4$ID %in% mouseNormal),]
segRes5 <- segRes5[-which(segRes5$ID %in% mouseNormal),]
segRes6 <- segRes6[-which(segRes6$ID %in% mouseNormal),]
segRes7 <- segRes7[-which(segRes7$ID %in% mouseNormal),]

combinedSegRes <- rbind(segRes1, segRes2, segRes3, segRes4,
                        segRes5, segRes6, segRes7)

combinedSegRes$ID <- nameStripper(combinedSegRes$ID)
combinedSegRes$ID <- str_remove(combinedSegRes$ID, "x.*")

combinedSegRes <- combinedSegRes[which(combinedSegRes$ID %in% goodSamps), ]


combinedSegRes$seg.mean[which(abs(combinedSegRes$seg.mean) < 0.2)] <- 0
combinedSegRes$seg.mean[which(combinedSegRes$q.val2 > 0.05)] <- 0
combinedSegRes$seg.mean[which(combinedSegRes$seg.mean < -2)] <- -2
combinedSegRes$seg.mean[which(combinedSegRes$seg.mean > 2)] <- 2
combinedSegRes$length <- combinedSegRes$loc.end - combinedSegRes$loc.start
combinedSegRes$seg.mean[which(combinedSegRes$length < 1e6)] <- 0
### cutoff at 35 from ggplot(combinedSegRes, aes(x = num.mark, y = abs(seg.mean))) + geom_point()
combinedSegRes$seg.mean[which(combinedSegRes$num.mark < 35)] <- 0

### below I filter anno

annotationList <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20211207yingSecondSet.xlsx")
annotationList$namestripped <- str_remove(nameStripper(annotationList$`#`), "\\-")
selfCalls <- read.table("/mnt/DATA5/tmp/kev/misc/20220104allEfSampsForSegHeat.txt", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)

badSamps <- selfCalls$sampleStripped[which(selfCalls$call_kh == "no tumor")]

otherAnno <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220107mouseCrcAnno.xlsx")
otherAnno$nameStripped <- str_remove(nameStripper(otherAnno$`#`), "\\-")



combinedSegRes2 <- combinedSegRes[-which(combinedSegRes$ID %in% badSamps),]
combinedSegRes2 <- combinedSegRes2[which(combinedSegRes2$ID %in% otherAnno$nameStripped),]

combinedSegRes2$chrom <- str_replace(combinedSegRes2$chrom, "23", "20")

combinedSegRes2 <- combinedSegRes2[-which(combinedSegRes2$chrom == "20"), ]

combinedSegRes2$ID <- otherAnno$`Sample ID`[match(combinedSegRes2$ID, otherAnno$nameStripped)]

combinedSegRes2_form <- combinedSegRes2[c("ID", "chrom", "loc.start", "loc.end",
                                          "num.mark", "seg.mean")]
colnames(combinedSegRes2_form) <- c("sampleID", "chrom", "start.pos",
                                    "end.pos", "n.probes", "mean")
combinedSegRes2_form$n.probes <- NA

combinedSegRes2_freq <- getFreqData(combinedSegRes2_form)

otherAnno2 <- otherAnno[which(otherAnno$`Sample ID` %in% combinedSegRes2_form$sampleID),
                        c("Sample ID", "Genotype", "Histology", "nameStripped")]

crcSampleNames <- otherAnno2$`Sample ID`[which(otherAnno2$Histology %in% c("adenoC", "cell line"))]
crcSampleNames <- otherAnno$nameStripped[match(crcSampleNames, otherAnno$`Sample ID`)]

adenomaSampleNames <- otherAnno2$`Sample ID`[which(otherAnno2$Histology %in% c("Adenoma"))]
adenomaSampleNames <- otherAnno$nameStripped[match(adenomaSampleNames, otherAnno$`Sample ID`)]


### graph

df_cn <- combinedSegRes2
df_freq <-  combinedSegRes2_freq
chromTextSpec = NULL

combinedGraph <- plotSegHeatmap(combinedSegRes2, combinedSegRes2_freq, segsize = 1)
otherAnno2 <- otherAnno2[match(treeOrder, otherAnno2$`Sample ID`),]
otherAnno2 <- otherAnno2[,1:3]

annoTab_coords  <- NULL
for (i in 2:(dim(otherAnno2)[2])) {
  sampVar <- otherAnno2[,i]
  tmpDf <- data.frame("xpos" = 1.2 + 0.2 * (i-2), 
                      "ypos" = seq(1.2, dim(otherAnno2)[1] * .2 + 1, .2),
                      "label"  = sampVar)
  colnames(tmpDf) <- c("xpos", "ypos", "label")
  annoTab_coords <- rbind(annoTab_coords, tmpDf)
}

### this needs to be changed each time b/c vector of labels if different
dfAnnoColIdx <- data.frame("labelNames" = c(unique(otherAnno2$Genotype), unique(otherAnno2$Histology)), 
                           "colors" = c("pink", "orange", "chartreuse3", "white", "grey",
                                        "lightblue", "darkblue", "darkmagenta",  "grey",
                                        "darkred", "black", "yellow", "white"))

annoTab_coords$col <- "white"
annoTab_coords$col <- dfAnnoColIdx$colors[match(annoTab_coords$label, dfAnnoColIdx$labelNames)]

anno_graph <- ggplot(annoTab_coords) +
  geom_rect(aes(xmin = xpos - 0.1, xmax = xpos + 0.1,
                ymin = ypos - 0.1, ymax = ypos + 0.1), fill = annoTab_coords$col) + 
  geom_hline(yintercept = seq(1.1, 0.2 * dim(annoTab_coords)[1]/(dim(otherAnno2)[2] - 1) + 1.1, .2), color = "black") +
  geom_vline(xintercept = seq(1.1, 0.2 * (dim(annoTab_coords)[2] - 2) + 1.1, 0.2), color = "black") + 
  scale_x_continuous(expand = expansion(mult = c(0, 0), 
                                        add = c(0, 0))) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0),
                                        add = c(0,nrow(otherAnno2)/300 - 0.05))) +
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
pdf("/mnt/DATA5/tmp/kev/misc/20220110newCrcHeatmap.pdf", useDingbats = FALSE, height = 8, width = 18)
grid.arrange(combinedGraph, anno_graph,
             layout_matrix=grobLayout2)
dev.off()

###
###
### annotation for hgsc samples
kathyNameConversion <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20201207annotations.xlsx")
kathyNameConversion$mg_id <- paste0("mg", kathyNameConversion$mg_id)
kathyNameConversion$old_name <- str_remove(nameStripper(kathyNameConversion$old_name), "\\-")
kathyNameConversion$new_name <- paste0(kathyNameConversion$old_name,
                                       "(", kathyNameConversion$mg_id, ")")
kathAnno <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220110allHgscNewPanelAnno.xlsx")
kathAnno$Sample <- str_remove(nameStripper(kathAnno$Sample), "\\-")
kathAnno$Type <- str_replace(kathAnno$Type, "Early HGSC", "eHGSC")
kathAnno$Type <- str_replace(kathAnno$Type, "Tumor", "endoCa")
kathAnno2 <- kathAnno[which(kathAnno$Type %in% c("HGSC","MMMT", "Normal", "eHGSC", "Met", "endoCa")),]
kathyNameConversion$old_name[!is.na(match(kathyNameConversion$old_name, kathAnno2$Sample))]
kathAnno2$Sample[!is.na(match(kathAnno2$Sample, kathyNameConversion$old_name))]

tmpMatch <- kathAnno2$Sample[which(kathAnno2$Sample %in% kathyNameConversion$old_name)]
# check: kathyNameConversion$old_name[match(tmpMatch, kathyNameConversion$old_name)]
kathAnno2$Sample[which(kathAnno2$Sample %in% kathyNameConversion$old_name)] <- kathyNameConversion$mg_id[match(tmpMatch, kathyNameConversion$old_name)]
colnames(kathAnno2) <- c("Sample ID", "Histology","Genotype")

combinedSegRes2 <- combinedSegRes[which(combinedSegRes$ID %in% kathAnno2$`Sample ID`),]
combinedSegRes2$chrom <- str_replace(combinedSegRes2$chrom, "23", "20")

combinedSegRes2 <- combinedSegRes2[-which(combinedSegRes2$chrom == "20"), ]

# combinedSegRes2$ID <- otherAnno$`Sample ID`[match(combinedSegRes2$ID, otherAnno$nameStripped)]

combinedSegRes2_form <- combinedSegRes2[c("ID", "chrom", "loc.start", "loc.end",
                                          "num.mark", "seg.mean")]
colnames(combinedSegRes2_form) <- c("sampleID", "chrom", "start.pos",
                                    "end.pos", "n.probes", "mean")
combinedSegRes2_form$n.probes <- NA

combinedSegRes2_freq <- getFreqData(combinedSegRes2_form)

otherAnno2 <- kathAnno2[which(kathAnno2$`Sample ID` %in% combinedSegRes2_form$sampleID),
                        c("Sample ID", "Genotype", "Histology")]

hgscSampleNames <-otherAnno2$`Sample ID`[which(otherAnno2$Histology %in% c("MMMT", "HGSC"))]

### graph

df_cn <- combinedSegRes2
df_freq <-  combinedSegRes2_freq
chromTextSpec = NULL

combinedGraph <- plotSegHeatmap(combinedSegRes2, combinedSegRes2_freq, segsize = 1, sizeVar = 4)
otherAnno2 <- otherAnno2[match(treeOrder, otherAnno2$`Sample ID`),]
otherAnno2 <- otherAnno2[,1:3]

annoTab_coords  <- NULL
for (i in 2:(dim(otherAnno2)[2])) {
  sampVar <- otherAnno2[,i]
  tmpDf <- data.frame("xpos" = 1.2 + 0.2 * (i-2), 
                      "ypos" = seq(1.2, dim(otherAnno2)[1] * .2 + 1, .2),
                      "label"  = sampVar)
  colnames(tmpDf) <- c("xpos", "ypos", "label")
  annoTab_coords <- rbind(annoTab_coords, tmpDf)
}

### this needs to be changed each time b/c vector of labels if different
dfAnnoColIdx <- data.frame("labelNames" = c(unique(otherAnno2$Genotype), unique(otherAnno2$Histology)), 
                           "colors" = c("green", "white", "green", "green", "purple",
                                        "green", "yellow", "green",  "green","orange",
                                        "green", "red", "green", "green", "green",
                                        "lightblue", "red", "lightblue","pink", "white",
                                        "yellow"))

annoTab_coords$col <- "white"
annoTab_coords$col <- dfAnnoColIdx$colors[match(annoTab_coords$label, dfAnnoColIdx$labelNames)]

anno_graph <- ggplot(annoTab_coords) +
  geom_rect(aes(xmin = xpos - 0.1, xmax = xpos + 0.1,
                ymin = ypos - 0.1, ymax = ypos + 0.1), fill = annoTab_coords$col) + 
  geom_hline(yintercept = seq(1.1, 0.2 * dim(annoTab_coords)[1]/(dim(otherAnno2)[2] - 1) + 1.1, .2), color = "black") +
  geom_vline(xintercept = seq(1.1, 0.2 * (dim(annoTab_coords)[2] - 2) + 1.1, 0.2), color = "black") + 
  scale_x_continuous(expand = expansion(mult = c(0, 0), 
                                        add = c(0, 0))) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0),
                                        add = c(0,nrow(otherAnno2)/300 - 0.05))) +
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
pdf("/mnt/DATA5/tmp/kev/misc/20220110newHgscHeatmap.pdf", useDingbats = FALSE, height = 8, width = 18)
grid.arrange(combinedGraph, anno_graph,
             layout_matrix=grobLayout2)
dev.off()

