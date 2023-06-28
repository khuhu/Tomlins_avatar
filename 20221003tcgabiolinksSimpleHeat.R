plotSegHeatmapV2 <- function(df_cn, df_freq, main = NULL,
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
    chromTextdf <- chromTextdf[1:22,]
  }
  
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  
  
  ### new way of assigning colors basically searches for where value lies within range of colors
  heatMapCol <- colorRampPalette(c("#3852A3","#FFFFFF","#EC1E24"))(1000)
  df_cn$seg.mean[which(df_cn$seg.mean > 3)] <- 3
  df_cn$seg.mean[which(df_cn$seg.mean < -3)] <- -3
  
  highVal <- 3
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
  # df_cn$col[which(df_cn$seg.mean < 0.2 & df_cn$seg.mean > -0.2)] <- "#FFFFFF"
  
  
  
  df_dist <- df_freq[,3:ncol(df_freq)]
  rownames(df_dist) <- paste0("chr", df_freq$chr, ":", df_freq$pos)
  df_dist <- t(df_dist)
  tmpModel <- hclust(dist(df_dist), method = "complete")
  dendo <- as.dendrogram(tmpModel)
  dendo_dat <- dendro_data(dendo, type = "rectangle")
  treeOrder <- dendo_dat$labels$label
  
  ### use tree order to get the dendo order for heatmap
  df_cn2 <- NULL
  i <- 1
  for (i in seq_along(treeOrder)) {
    print(i)
    tmpCn <- df_cn[which(df_cn$ID == treeOrder[i]),]
    for (j in unique(tmpCn$chrom)) {
      print(j)
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
  # dendo_graph <- ggplot(ggdendro::segment(dendo_dat)) + 
  #   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  #   coord_flip() + 
  #   scale_y_reverse(expand = c(0,0)) +
  #   scale_x_continuous(expand = expansion(mult = c(0, 0), 
  #                                         add = c(0.5, 1.5))) +
  #   theme_bw() + 
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(), 
  #         axis.title = element_blank(),
  #         axis.text = element_blank(),
  #         axis.ticks = element_blank(),
  #         panel.border = element_blank())
  # 
  
  
  heat_graph <- ggplot(df_cn2) +
    geom_rect(aes(xmin = xstart, xmax = xend,
                  ymin = ypos - 0.1, ymax = ypos + 0.1),
              fill = df_cn2$col, color = NA) + 
    geom_vline(xintercept = chromBreak) +
    # geom_hline(yintercept = hline, color = "grey") + 
    scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
    scale_y_continuous(breaks=seq(min(df_cn2$ypos), max(df_cn2$ypos), 0.2),
                       labels = xlabels, expand = expansion(mult = c(0, 0.02))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos, label = chrom), size = 2.5)
  
  # dendo_grob <- ggplotGrob(dendo_graph)
  heat_grob <- ggplotGrob(heat_graph)
  
  plot(heat_grob)
  # grobLayout <- rbind(c(1,rep(2,19)),
  #                     c(1,rep(2,19)),
  #                     c(1,rep(2,19)),
  #                     c(1,rep(2,19)))
  # grid.arrange(dendo_grob, heat_grob,
  #              layout_matrix=grobLayout)
}

source("/home/kevhu/scripts/20210802syntenyFunctions.R")

library(TCGAbiolinks)
library(pheatmap)
library(ggdendro)
library(stringr)
library(gridExtra)

### html vignette https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html#Harmonized_data

View(getGDCprojects())
TCGAbiolinks::getDataCategorySummary("TCGA-COAD")

query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Copy Number Variation",
  data.type = "Copy Number Segment",
  legacy = FALSE
)

GDCdownload(query, directory = "/mnt/DATA5/tmp/kev/tmpDbs/tcgaBiolinks/")
data <- GDCprepare(query, directory = "/mnt/DATA5/tmp/kev/tmpDbs/tcgaBiolinks/")


### using data to reconstruct the segment plot
apcMuts <- read.table("/mnt/DATA5/tmp/kev/misc/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv",
                       sep = "\t", stringsAsFactors = FALSE, header = TRUE)
data_reduced <- grep(paste0(apcMuts$Sample.ID, collapse = "|"), data$Sample)
data_reduced_df <- data[data_reduced, ]
data_reduced_df_form <- data_reduced_df[, c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")]
colnames(data_reduced_df_form) <- c("sampleID", "chrom", "start.pos", "end.pos", "n.probes", "mean")
data_reduced_df_form <- data_reduced_df_form[which(data_reduced_df_form$chrom %in% c(1:22)),]
data_reduced_df_form <- as.data.frame(data_reduced_df_form)
data_reduced_freq <- getFreqData(data_reduced_df_form)

df_res <- data.frame(data_reduced_df[, c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")])
colnames(df_res) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
df_res <- df_res[which(df_res$chrom %in% c(1:22)),]

dev.off()
# pdf("/mnt/DATA5/tmp/kev/misc/testTcgaHeat.pdf", height = 12)
pdf("/mnt/DATA5/tmp/kev/misc/testTcgaHeat.pdf")
plotSegHeatmapV2(df_res, data_reduced_freq, segsize = 3, chromTextSpec = "hg19")
dev.off()



### BBN cna data

query <- GDCquery(
  project = "TCGA-BLCA",
  data.category = "Copy Number Variation",
  data.type = "Copy Number Segment",
  legacy = FALSE
)

GDCdownload(query, directory = "/mnt/DATA5/tmp/kev/tmpDbs/tcgaBiolinks/")
data <- GDCprepare(query, directory = "/mnt/DATA5/tmp/kev/tmpDbs/tcgaBiolinks/")

blcaMuts <- read.table("/mnt/DATA5/tmp/kev/misc/blca_withmuts_tcga_clinical_data.tsv",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
data_reduced <- grep(paste0(blcaMuts$Sample.ID, collapse = "|"), data$Sample)
data_reduced_df <- data[data_reduced, ]
data_reduced_df_form <- data_reduced_df[, c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")]
colnames(data_reduced_df_form) <- c("sampleID", "chrom", "start.pos", "end.pos", "n.probes", "mean")
data_reduced_df_form <- data_reduced_df_form[which(data_reduced_df_form$chrom %in% c(1:22)),]
data_reduced_df_form <- as.data.frame(data_reduced_df_form)
data_reduced_freq <- getFreqData(data_reduced_df_form)

df_res <- data.frame(data_reduced_df[, c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")])
colnames(df_res) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
df_res <- df_res[which(df_res$chrom %in% c(1:22)),]

dev.off()
# pdf("/mnt/DATA5/tmp/kev/misc/testTcgaHeat.pdf", height = 12)
pdf("/mnt/DATA5/tmp/kev/misc/blcaHeat.pdf")
plotSegHeatmapV2(df_res, data_reduced_freq, segsize = 3, chromTextSpec = "hg19")
dev.off()


