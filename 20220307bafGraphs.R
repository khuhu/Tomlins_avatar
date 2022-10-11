# original script: /home/kevhu/scripts/20211011tvcSnpAnalysis.R
# original of this script from brahm: /home/kevhu/scripts/20211204tvcSnpAnalysis.R

### function to graph BAF + cn overlap 
###

freqPlot_baf <- function(df, main = "no title", chromTextSpec = NULL){
  require(ggplot2)
  
  if(is.null(chromTextSpec)){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "mm10"){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "hg19") {
    chromTextdf <- read.table("/br_z1/kevin_storage/misc/20210801hg19_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  }
  
  
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  
  ### divide positions by megabase locations
  
  df$CHROM <- as.numeric(str_remove(df$CHROM, "chr"))
  df$point <- df$POS/1e6
  
  for (i in unique(df$CHROM)) {
    df$point[which(df$CHROM == i)] <- df$point[which(df$CHROM == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  
  df <- df[,c("CHROM", "point", "AF")]
  colnames(df) <- c("chrom", "pos", "af")
  
  
  if (nrow(df) == 0) {
    ggplot(df, aes(x = pos, y = af)) + geom_hline(yintercept = c(1,0.75, 0.5, 0.25, 0), color = "#D4D4D4")+
      geom_point(size = 0.25, color = df$color) + geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 1, by = 0.10),
                                                             limits = c(0,1.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + .2, label = chrom), size = 2.5) 
  } else{
    df$color <- "#000000"
    df$color[which(df$af < 0.3)] <- "#FF0000"
    df$color[which(df$af > 0.7)] <- "#FF0000"
    
    ggplot(df, aes(x = pos, y = af)) + geom_hline(yintercept = c(1,0.75, 0.5, 0.25, 0), color = "#D4D4D4")+
      geom_point(size = 0.25, color = df$color) + geom_vline(xintercept=chromBreak) + theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(0, 1, by = 0.10),
                                                             limits = c(0,1.2)) +
      geom_text(data = chromTextdf, aes(x = xpos, y = ypos + .2, label = chrom), size = 2.5) 
  }
}


freqPlot_cn_amp <- function(df_cn, df_amp, main = NULL, chromTextSpec = NULL){
  require(ggplot2)
  
  if(is.null(chromTextSpec)){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "mm10"){
    chromTextdf <- read.table("/home/kevhu/20210801mm10_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  } else if(chromTextSpec == "hg19") {
    chromTextdf <- read.table("/br_z1/kevin_storage/misc/20210801hg19_graphingLimits.txt",
                              sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  }
  
  chromBreak <- c(0, chromTextdf$chromBreaksPos)
  
  df_cn$loc.start <- df_cn$loc.start/1e6
  df_cn$loc.end <- df_cn$loc.end/1e6
  df_cn$col <- "#000000"
  df_cn$col[which(df_cn$seg.mean > 0.2)] <- "#FF0000"
  df_cn$col[which(df_cn$seg.mean < -0.2)] <- "#0000FF"
  
  for (i in unique(df_cn$chrom)) {
    df_cn$loc.start[which(df_cn$chrom == i)] <- df_cn$loc.start[which(df_cn$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
    df_cn$loc.end[which(df_cn$chrom == i)] <- df_cn$loc.end[which(df_cn$chrom == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  df_cn <- df_cn[,c("chrom", "loc.start", "loc.end", "seg.mean", "col")]
  colnames(df_cn) <- c("chrom", "start", "end","cn", "col")
  
  
  df_amp$Pos <- ((df_amp$StartPos + df_amp$EndPos)/2)/1e6
  df_amp$col <- "#64e5d1"
  df_amp$col[grep("SNP", df_amp$Gene)] <- "#006400"
  
  for (i in unique(df_amp$ChromNum)) {
    df_amp$Pos[which(df_amp$ChromNum == i)] <- df_amp$Pos[which(df_amp$ChromNum == i)] +
      chromTextdf$graphingStart[which(chromTextdf$chrom == i)]
  }
  df_amp <- df_amp[,c(1,6,5,7)]
  colnames(df_amp) <- c("chrom", "pos", "cn", "col")
  
  
  ggplot(df_cn) + geom_hline(yintercept = c(-3, -2 , -1 , 0, 1, 2, 3), color = "#D4D4D4") +
    geom_point(data = df_amp, aes(x = pos, y = cn),
               colour = df_amp$col, size = 0.05, inherit.aes = FALSE) +
    geom_segment(aes(x = start, xend = end, y = cn, yend = cn), colour = df_cn$col) + 
    geom_vline(xintercept=chromBreak) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = seq(-3, 3, by = 1),
                                                           limits = c(-3.2, 3.2)) +
    geom_text(data = chromTextdf, aes(x = xpos, y = ypos + 2.3, label = chrom), size = 2.5) + ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5))
}


nameStripper <- function(df){
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
}


tcDf <- read.table("/br_z1/kevin_storage/misc/20210718hgscTcDf.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

non_zero_noPosFilt <- read.table("/br_z1/kevin_storage/misc/20211204non_zero_noPosFilt.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

non_zero_20 <- read.table("/br_z1/kevin_storage/misc/20211203mm10BafsNonzero20.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

allAmps <- read.table("/br_z1/kevin_storage/misc/20211129all_mouse_amps.txt", sep = "\t", stringsAsFactors = FALSE,
                      header = TRUE, check.names = FALSE)

allAmps_tc <- read.table("/br_z1/kevin_storage/misc/20211206all_mouse_amps_tc.txt", sep = "\t", stringsAsFactors = FALSE,
                      header = TRUE, check.names = FALSE)

# allSegs <- read.table("/br_z1/kevin_storage/misc/20211203hgscSegRes.txt", sep = "\t", stringsAsFactors = FALSE,
#                       header = TRUE, check.names = FALSE)

allSegs <- read.table("/br_z1/kevin_storage/misc/20211206all_mouse_seg.txt", sep = "\t", stringsAsFactors = FALSE,
                      header = TRUE, check.names = FALSE)


# allSegs_tc <- read.table("/br_z1/kevin_storage/misc/20211203hgscMm10SegRes_tc.txt", sep = "\t", stringsAsFactors = FALSE,
#                       header = TRUE, check.names = FALSE)

allSegs_tc <- read.table("/br_z1/kevin_storage/misc/20211206all_mouse_seg_tc.txt", sep = "\t", stringsAsFactors = FALSE,
                         header = TRUE, check.names = FALSE)


genoBed <- read.table("/br_z1/kevin_storage/misc/20211120snpsOnIAD202670_167_positions.bed", sep = "\t",
                      stringsAsFactors = FALSE, header = FALSE)

positionString <- paste0(genoBed$V1,":" ,genoBed$V3)

allSegs$seg.mean[which(abs(allSegs$seg.mean) < 0.2)] <- 0
#allSegs$seg.mean[which(allSegs$seg.mean < (log2(5/4)) & allSegs$seg.mean > (log2(3/4)))] <- 0
allSegs$seg.mean[which(allSegs$q.val2 > 0.05)] <- 0
allSegs$seg.mean[which(allSegs$seg.mean < -3)] <- -3
allSegs$seg.mean[which(allSegs$seg.mean > 3)] <- 3
#allSegs$seg.mean[which((allSegs$end.pos - allSegs$start.pos) < (10 * 1e6))] <- 0
allSegs$ID <- nameStripper(allSegs$ID)
allSegs$ID <- str_remove(allSegs$ID, "x.*")


allSegs_tc$seg.mean[which(abs(allSegs_tc$seg.mean) < 0.2)] <- 0
#allSegs_tc$seg.mean[which(allSegs_tc$seg.mean < (log2(5/4)) & allSegs_tc$seg.mean > (log2(3/4)))] <- 0
allSegs_tc$seg.mean[which(allSegs_tc$q.val2 > 0.05)] <- 0
allSegs_tc$seg.mean[which(allSegs_tc$seg.mean < -3)] <- -3
allSegs_tc$seg.mean[which(allSegs_tc$seg.mean > 3)] <- 3
#allSegs_tc$seg.mean[which((allSegs_tc$end.pos - allSegs_tc$start.pos) < (10 * 1e6))] <- 0
allSegs_tc$ID <- str_remove(allSegs_tc$ID, "x.*")


colnames(allAmps)[7:ncol(allAmps)] <- str_remove(colnames(allAmps)[7:ncol(allAmps)] , "x.*")
colnames(allAmps_tc)[7:ncol(allAmps_tc)] <- str_remove(colnames(allAmps)[7:ncol(allAmps_tc)] , "x.*")




### things to add: load in tc so i can add it to the title name - new graphs should just be raw, tc corrected and bafs
### make sure I create a separate mouse_amplicon_tc and z_score for tc - do this in the loop
### make sure only to tc correct things > abs(0.2) so I don't tc correct for noisy data ... seems to happen and doesn't agree with BAFs .. need to check this

i <- "efd64"

for (i in unique(allSegs$ID)) {
  print(i)
  
  
  if (is.null(allAmps[[i]])) {
    next
  }
  
  print(i)
  
  matchingTc <- tcDf$tc[match(tolower(i), tcDf$sample)]
  
  testDf_cn <- allSegs[which(allSegs$ID == i),]
  testDf_tc_cn <- allSegs_tc[which(allSegs_tc$ID == i),]
  testDf_amp <- cbind(allAmps[,c(3:6)], allAmps[[i]])
  b_2 <- freqPlot_cn_amp(testDf_cn, testDf_amp, main = paste(i, "MousePanel, Baf: no filt, Baf w/ filts", "tc:", matchingTc))
  b_1 <- freqPlot_cn_amp(testDf_tc_cn, testDf_amp)

  
  testDf <- non_zero_20[which(non_zero_20$sample == i),]
  a_20 <- freqPlot_baf(testDf)
  
  g20 <- ggplotGrob(a_20)
  gB_2 <- ggplotGrob(b_2)
  gB_1 <- ggplotGrob(b_1)
  
  dev.off()
  # png(filename = paste0("/br_z1/kevin_storage/misc/20211027_tc_zscore_amp", samplesSeg[i], ".png"), width = 1800, height = 500)
  pdf(file = paste0("/br_z1/kevin_storage/misc/20211101noTcAmpAnalysis/20211204_filtCompsBoth__20_", i, ".pdf"), width = 15, height = 7)
  grid::grid.newpage()
  grid::grid.draw(rbind(gB_2, gB_1, g20))
  dev.off()
  
}




