# should be a list of functions needed to be used after
# processing segmentation (input used)
library(GenomicFeatures)
library(dplyr)
library(dbplyr)
library(circlize)
library(stringr)
library(DNAcopy)
library(stringr)
#library(optparse)
library(data.table)

# four variables should be mouse amplicon files, o/e ratio,
# normal vector name for FFPE vs FF, and output path




# data and functions needed for synteny
synteny_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_syntenyDf.txt",
                               sep = "\t", stringsAsFactors = FALSE, header = TRUE)
cyto_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_cytoDf.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = FALSE)
human_gene <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210407humanGeneTable.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)
mouse_gene <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210407mouseGeneTable.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)
geneNameDf <- read.table("/home/kevhu/data/20201021geneNameBiomart.txt", header = TRUE,
                         stringsAsFactors = FALSE, sep = "\t")
mouseBedFile <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.gc.bed", stringsAsFactors = FALSE, sep = "\t",
                           header = FALSE)


tcga_ov_brca1_amp <- read.table("/mnt/DATA5/tmp/kev/misc/20210424tcga_ov_brca1_amp_0.2.bed", sep = "\t",
                                stringsAsFactors = FALSE, header = TRUE)


tcga_ov_brca1_del <- read.table("/mnt/DATA5/tmp/kev/misc/20210424tcga_ov_brca1_del_0.2.bed", sep = "\t",
                                stringsAsFactors = FALSE, header = TRUE)



tcga_ov_amp <- read.table("/mnt/DATA5/tmp/kev/misc/20210424tcga_ov_amp_0.2.bed", sep = "\t",
                                stringsAsFactors = FALSE, header = TRUE)


tcga_ov_del <- read.table("/mnt/DATA5/tmp/kev/misc/20210424tcga_ov_del_0.2.bed", sep = "\t",
                                stringsAsFactors = FALSE, header = TRUE)


# below are the two inputs needed for segmentation
mouseAmplicons <- read.table("/mnt/DATA5/tmp/kev/newMouse2/cnAmplicon_matrix.txt", sep = "\t",
                             stringsAsFactors = FALSE, header = TRUE)

zscore_gc_oe_ratios <- read.table("/mnt/DATA5/tmp/kev/newMouse2/gcCorrectedCounts_matrix.txt", sep = "\t",
                                  stringsAsFactors = FALSE, header = TRUE)


all_probes_grange <- GRanges(seqnames = Rle(zscore_gc_oe_ratios$ChromNum),
                             ranges = IRanges(start = zscore_gc_oe_ratios$StartPos,
                                              end = zscore_gc_oe_ratios$EndPos))

colorVector <- c("#CD5C5C", "#C71585", "#FF8C00", "#F0E68C",
                 "#4B0082", "#32CD32", "#000080", "#800000",
                 "#F5F5DC", "#2F4F4F", "#808000", "#E6E6FA",
                 "#FFD700", "#FFA07A", "#7B68EE", "#008B8B",
                 "#BC8F8F", "#FAEBD7", "#FFE4E1", "#000000",
                 "#00FFFF", "#8FBC8B", "#B22222", "#D0D3D4")

# this needs to be an input vector in the snakefile - get inputs
# mouseNormal <- c("13604n", "14104t", "14154n", "14433n", "2405n", "2519n", "2796n", "3867n","8234n")
# 2796n 13604n - have BRCA1 gene deletions
# need the MG names


mouseNormal <- c("MG_17X49", "MG_18X50", "MG_23X55", "MG_6X38",
                  "MG_8X40","MG_11X43", "MG_13X45")

# used in converted human gene names to mouse
firstUpper <- function(gene){
  firstLetter <- toupper(substr(gene, start = 1, stop = 1))
  restOfGene <- tolower(substr(gene, start = 2, stop = nchar(gene)))
  res <- paste0(firstLetter, restOfGene)
  return(res)
}

maximums <- function(x) which(x - shift(x, 10) > 0  & x - shift(x, 10, type='lead') > 0)
minimums <- function(x) which(x - shift(x, 10) < 0  & x - shift(x, 10, type='lead') < 0)


# function for processing the corrected amplicons into segmentation
# when I find time: parallelize this process

ampSeg <- function(mouseAmplicons, mouseBedFile){
  mouseBedIdx <- as.numeric(str_remove(mouseAmplicons$AmpliconId, "AMP_"))
  mouseAmplicons2 <- cbind(mouseBedFile$V1[mouseBedIdx], mouseBedFile$V2[mouseBedIdx], mouseAmplicons)
  mouseAmplicons2[,6:42] <- log2(mouseAmplicons2[,6:42])
  mouseAmplicons2$`mouseBedFile$V1[mouseBedIdx]` <- as.numeric(mouseAmplicons2$`mouseBedFile$V1[mouseBedIdx]`)
  allMouseCalls <- colnames(mouseAmplicons2)[6:42]
  segResults <- NULL
  for (i in allMouseCalls) {
    tmpCNA_obj <- CNA(cbind(mouseAmplicons2[[i]]),
                      mouseAmplicons2$ChromNum, mouseAmplicons2$`mouseBedFile$V2[mouseBedIdx]`,
                      data.type="logratio",sampleid=i)
    smoothed_tmpCNA <- smooth.CNA(tmpCNA_obj)
    segment_tmpCNA <- segment(smoothed_tmpCNA, verbose = 1, undo.splits = "sdundo", undo.SD = 0.2, alpha = 0.05,
                              p.method = c("hybrid"))
    tmpCalls <- segments.p(segment_tmpCNA)
    tmpCalls_filt <- tmpCalls[which(abs(tmpCalls$seg.mean) > .2),]
    segResults <- rbind(segResults, tmpCalls_filt)
  }
  return(segResults)
}

calcZscore <- function(segRes){
  z_vector <- NULL
  nmean_vector <- NULL
  tmean_vector <- NULL
  for (i in unique(segRes$ID)) {
    tmp_seg_res <- segRes[which(segRes$ID == i), ]
    tmp_seg_ranges <- GRanges(seqnames = Rle(tmp_seg_res$chrom), 
                              ranges = IRanges(start = tmp_seg_res$loc.start,
                                               end = tmp_seg_res$loc.end))
    for (j in seq_along(tmp_seg_ranges)) {
      tmpOverlap <- findOverlaps(all_probes_grange, tmp_seg_ranges[j])
      tmpTest <- zscore_gc_oe_ratios[[i]][queryHits(tmpOverlap)]
      
      tmpNormal <- zscore_gc_oe_ratios[queryHits(tmpOverlap), mouseNormal]
      tmpNormal2 <- unlist(tmpNormal)
      
      z_stat <- (mean(tmpTest) - mean(tmpNormal2))/sd(tmpNormal2)
      normal_seg_mean <- mean(tmpNormal2)
      tumor_seg_mean <- mean(tmpTest)
      z_vector <- c(z_vector, z_stat)
      nmean_vector <- c(nmean_vector, normal_seg_mean)
      tmean_vector <- c(tmean_vector, tumor_seg_mean)
    }
  }
  return(list("z_vector" = z_vector, "nmean_vector"= nmean_vector,
              "tmean_vector" = tmean_vector))
}


### filt could have parameter for size of filtering
segZscoresFilt <- function(segResults, zscores){
  
  # filts for significance, cn-change and length
  segmental_zscores <- cbind(segResults, "z-scores" = zscores[["z_vector"]],
                             "normal_seg_mean" = zscores[["nmean_vector"]],
                             "tumor_seg_mean" = zscores[["tmean_vector"]])
  segmental_zscores2 <- segmental_zscores[which(segmental_zscores$num.mark > 10),]
  segmental_zscores2$length <- segmental_zscores2$loc.end - segmental_zscores2$loc.start
  largeChromCutoff <- 1000000
  segmental_zscores3 <- segmental_zscores2[which(segmental_zscores2$length >= largeChromCutoff),]
  segmental_zscores3 <- segmental_zscores3[which(abs(segmental_zscores3$`z-scores`) > 1.645),]
  segInput <- segmental_zscores3[,c("ID", "chrom", "loc.start", "loc.end", "seg.mean")]
  return(segInput)
}

# input should be segmentation results - should go from mouse to human

syntenyPlotInputs <- function(segInput){
  segInput$chrom <- str_replace_all(segInput$chrom, "23", "X")
  tumor_freq_ranges <- GRanges(seqnames = paste0("chr", segInput$chrom),
                               ranges = IRanges(start = as.numeric(segInput$loc.start),
                                                end = as.numeric(segInput$loc.end)))
  synteny_ranges <- GRanges(seqnames = paste0("chr", synteny_hg38mm10$comp_chr),
                            ranges = IRanges(start = synteny_hg38mm10$comp_start_pos,
                                             end = synteny_hg38mm10$comp_end_pos))
  
  
  res_overlap <- findOverlaps(synteny_ranges, tumor_freq_ranges)
  query_hits <- queryHits(res_overlap)
  subject_hits <- subjectHits(res_overlap)
  
  
  mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg38mm10$comp_chr[query_hits]),
                          "start" = synteny_hg38mm10$comp_start_pos[query_hits],
                          "end" = synteny_hg38mm10$comp_end_pos[query_hits],
                          stringsAsFactors = FALSE)
  
  human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg38mm10$ref_chr[query_hits]),
                          "start" = synteny_hg38mm10$ref_start_pos[query_hits],
                          "end" = synteny_hg38mm10$ref_end_pos[query_hits],
                          stringsAsFactors = FALSE)
  
  
  
  resList <- list("human_bed" = human_bed,
                  "mouse_bed" = mouse_bed)
  
  return(resList)
}

# gene part incomplete b/c odd matching and mouse genes weren't not found ...
# do this later and only for genes on panel ... easiest
# not that important
syntenyGeneInputs <- function(geneList){
  mouseGeneList <- geneList
  humanGeneList <- NULL
  for (i in geneList) {
    humanOrthoName <- geneNameDf$external_gene_name[which(geneNameDf$mmusculus_homolog_associated_gene_name == i)]
    if (length(humanOrthoName)  == 0) {
      humanOrthoName <- toupper(i)
    }
    
    humanGeneList <- c(humanGeneList, humanOrthoName)
  }
  
  tmpH <- human_gene[which(human_gene$external_gene_name %in% humanGeneList), ]
  tmpM <- mouse_gene[which(mouse_gene$external_gene_name %in% geneList), ]
  return(list(tmpH, tmpM))
}



downsampleRegions <- function(df){
  
  finalTbl <- NULL
  for (i in unique(df$chr)) {
    tmpDf <- df[which(df$chr == i),]
    # weird artefact where the last row is nearly entire chromosome
    # for every entry - so removing it
    tmpDf <- tmpDf[1:(nrow(tmpDf) - 1),]
    
    fit <- loess(start ~ cn, degree=1, span = 0.1, data=tmpDf)
    tmpDf$loess <- fit$fitted
    tmpDf_red <- tmpDf[c(maximums(tmpDf$loess),
                      minimums(tmpDf$loess)),]
    tmpDf_red <- tmpDf_red[order(tmpDf_red$start, decreasing = FALSE),]
    redIdx <- c(1,1+which(abs(diff(tmpDf_red$cn)) >= 0.05))
    tmpDf_red2 <- NULL
    for (j in 1:length(redIdx)) {
      if (j == nrow(tmpDf_red)) {
        # if last change point is in last region
        
        tmpIdx <- redIdx[j]
        tmpVec <- c("chr" = i,
                    "start" = tmpDf_red$start[tmpIdx],
                    "end" = tmpDf_red$end[tmpIdx], 
                    "cn" = round(tmpDf_red$cn[tmpIdx],digits = 2))
      } else if (j == length(redIdx)){
        # if last change point is not in last region
        
        tmpIdx <- redIdx[j]
        tmpVec <- c("chr" = i,
                    "start" = tmpDf_red$start[tmpIdx],
                    "end" = tmpDf_red$end[nrow(tmpDf_red)], 
                    "cn" = round(mean(c(tmpDf_red$cn[tmpIdx],
                                        tmpDf_red$cn[nrow(tmpDf_red)])), digits = 2))
      } else {
        tmpIdx <- redIdx[j]
        tmpIdx2 <- redIdx[j + 1]
        tmpVec <- c("chr" = i,
                    "start" = tmpDf_red$start[tmpIdx],
                    "end" = tmpDf_red$start[tmpIdx2] - 1, 
                    "cn" = round(mean(tmpDf_red$cn[tmpIdx:tmpIdx2]),digits = 2))
      }
      tmpDf_red2 <- rbind(tmpDf_red2, tmpVec)
    }
    
    tmpDf_red2 <- data.frame(tmpDf_red2, stringsAsFactors = FALSE)
    tmpDf_red2[,2:4] <- lapply(tmpDf_red2[,2:4], as.numeric)
    
    finalTbl <- rbind(finalTbl, tmpDf_red2)
  }
  print("done")
  return(finalTbl)
}



getFreqBed <- function(amp, del){
  # modified original code for brevity
  amp2 <- amp
  del2 <- del
  
  amp2$Freq <- round(amp2$Freq/100, digits = 2)
  amp2$Chr <- paste0("h_chr", amp2$Chr)
  #amp2 <- amp2[which(amp2$Chr %in% cyto_interest),]
  colnames(amp2) <- c("chr","start","end","cn")
  
  amp3 <- downsampleRegions(amp2)
  amp3$ybot <- 0
  amp3$ytop <- amp3$cn
  
  del2$Freq <- round(del2$Freq/100, digits = 2)
  del2$Chr <- paste0("h_chr", del2$Chr)
  #del2 <- del2[which(del2$Chr %in% cyto_interest),]
  colnames(del2) <- c("chr","start","end","cn")
  del2$cn <- -del2$cn
  
  del3 <- downsampleRegions(del2)
  del3$ytop <- 0
  del3$ybot <- del3$cn

  cn_track_bed <- rbind(amp3, del3)
  
  return(cn_track_bed)
}


# only includes large chromosomal changes 
# the ony change is instead of genes I can use the genomic track instead 
# to represent gains or losses

#segRes <- NULL
#segRes <- ampSeg(mouseAmplicons, mouseBedFile)
#save(segRes, file = "/mnt/DATA5/tmp/kev/misc/20210409mm10SegRes.Robj")
load(file = "/mnt/DATA5/tmp/kev/misc/20210409mm10SegRes.Robj")
segZscores <- calcZscore(segRes)
segZfilt <- segZscoresFilt(segRes, segZscores)


# I AM DUMB, dont need to downsample each time
# just downsample once, and use that file to subset using cyto_of_interest

# downsamping method works well plots for full version using copy-nubmer package
# /mnt/DATA5/tmp/kev/misc/210420tcga_ov_brca1_loss_copynumber.pdf
# /mnt/DATA5/tmp/kev/misc/20210420tcga_ov_freq_copynumber.pdf
tcga_ov_brca1_freq <- getFreqBed(amp = tcga_ov_brca1_amp, del = tcga_ov_brca1_del)
tcga_ov_freq <- getFreqBed(amp = tcga_ov_amp, del = tcga_ov_del)




for(i in unique(segZfilt$ID)){
  print(i)
  i <- "MG_3X35"
  
  segZfilt_sample <- segZfilt[which(segZfilt$ID == i),]
  syntenyBed <- syntenyPlotInputs(segZfilt_sample)
  human_bed <- syntenyBed$human_bed
  mouse_bed <- syntenyBed$mouse_bed
  
  #make sure to floor and ceilign to -1 and 1 for cn later
  
  mouse_cn <- segZfilt_sample[,2:5]
  colnames(mouse_cn) <- c("chr", "start", "end", "cn")
  mouse_cn$chr <- paste0("m_chr", mouse_cn$chr)
  mouse_cn$chr <- str_replace(mouse_cn$chr, "m_chr23", "m_chrX")
  
  mouse_cn$ytop <- ifelse(mouse_cn$cn > 0 , mouse_cn$cn + 0.025, mouse_cn$cn) 
  mouse_cn$ybot <- ifelse(mouse_cn$cn > 0 , mouse_cn$cn, mouse_cn$cn - 0.025) 
  
  cyto_interest <- unique(c(human_bed$chr, mouse_bed$chr))
  cyto_combined_red <- cyto_hg38mm10[which(cyto_hg38mm10$V1 %in% cyto_interest),]
  
  
  
  # get freq bed per sample
  
  cn_track_bed <- rbind(tcga_ov_freq, mouse_cn)
  cn_track_bed <- cn_track_bed[which(cn_track_bed$chr %in% cyto_interest), ]
  # subset cn_track_bed amps and dels by cyto_of_interest
  
  # some weird artefacts in the freq plot for whole-genome
  # looks like an overall arm level measurement, so filter out
  
  assign("tmpVar", cn_track_bed)
  print(dim(cn_track_bed))
  #print(cn_track_bed[which((cn_track_bed$start - cn_track_bed$end) > 0),])
  
  track_color <- rep("#000000", nrow(cn_track_bed))
  track_color <- ifelse(cn_track_bed$cn < 0 , "#00008B", "#8B0000") 
  cn_track_bed$col = track_color
  
  mouse_chrs <- unique(mouse_bed$chr)
  human_chrs <- unique(human_bed$chr)
  
  colorVector2 <- rep("#000000", nrow(mouse_bed))
  for (j in seq_along(mouse_chrs)) {
    colorVector2[which(mouse_bed$chr == mouse_chrs[j])] <- colorVector[j]
  }
  
  
  # plotting each of the graphs
  #pdf(paste0(opt$output,"circos_", i, ".pdf"), width = 6, height = 6, useDingbats = FALSE)
  pdf(paste0("/mnt/DATA5/tmp/kev/testSynteny/","circos_ov_", i, ".pdf"),
      width = 6, height = 6, useDingbats = FALSE)
  circos.par("track.height"= 0.10) + 
    circos.initializeWithIdeogram(cytoband = cyto_combined_red, sort.chr = FALSE,
                                  plotType = c("ideogram", "labels")) + 
    circos.genomicTrack(cn_track_bed, ylim = c(-1.00,1.00),
                        panel.fun = function(region, value, ...) {
                          circos.genomicRect(region, value, col = value$col, 
                                             ytop = value$ytop,
                                             ybottom = value$ybot,
                                             border = NA,...)
                        }) + 
    circos.genomicLink(mouse_bed, human_bed,
                       col = colorVector2,
                       border = NA)
  dev.off()
  readline(prompt="Press [enter] to continue")
}





# below is testing adding various plotting effects to circos plots
# i.e most recent example would be the alteration frequencies



segZfilt_sample <- segZfilt[which(segZfilt$ID == "MG_3X35"),]
syntenyBed <- syntenyPlotInputs(segZfilt_sample)
human_bed <- syntenyBed$human_bed
mouse_bed <- syntenyBed$mouse_bed



cyto_interest <- unique(c(human_bed$chr, mouse_bed$chr))
cyto_combined_red <- cyto_hg38mm10[which(cyto_hg38mm10$V1 %in% cyto_interest),]

# just make each ytop and ybot separately
tcga_ov_brca1_amp2 <- tcga_ov_brca1_amp
tcga_ov_brca1_del2 <- tcga_ov_brca1_del

tcga_ov_brca1_amp2$Freq <- round(tcga_ov_brca1_amp2$Freq/100, digits = 2)
tcga_ov_brca1_amp2$Chr <- paste0("h_chr", tcga_ov_brca1_amp2$Chr)
tcga_ov_brca1_amp2 <- tcga_ov_brca1_amp2[which(tcga_ov_brca1_amp2$Chr %in% cyto_interest),]
colnames(tcga_ov_brca1_amp2) <- colnames(mouse_bed)
tcga_ov_brca1_amp2$ybot <- 0
tcga_ov_brca1_amp2$ytop <- tcga_ov_brca1_amp2$cn

tcga_ov_brca1_del2$Freq <- round(tcga_ov_brca1_del2$Freq/100, digits = 2)
tcga_ov_brca1_del2$Chr <- paste0("h_chr", tcga_ov_brca1_del2$Chr)
tcga_ov_brca1_del2 <- tcga_ov_brca1_del2[which(tcga_ov_brca1_del2$Chr %in% cyto_interest),]
colnames(tcga_ov_brca1_del2) <- colnames(mouse_bed)
tcga_ov_brca1_del2$cn <- -tcga_ov_brca1_del2$cn
tcga_ov_brca1_del2$ytop <- 0
tcga_ov_brca1_del2$ybot <- tcga_ov_brca1_del2$cn

mouse_bed2 <- mouse_bed
mouse_bed2$ytop <- ifelse(mouse_bed2$cn > 0 , mouse_bed2$cn + 0.05, mouse_bed2$cn) 
mouse_bed2$ybot <- ifelse(mouse_bed2$cn > 0 , mouse_bed2$cn, mouse_bed2$cn -0.05) 

cn_track_bed <- rbind(tcga_ov_brca1_amp2, tcga_ov_brca1_del2, mouse_bed2)


#color_cn = colorRamp2(c(-1.0, 0, 1.0), c("#00008B", "#000000", "#8B0000"))
#track_color <- color_cn(cn_track_bed$cn)
track_color <- rep("#000000", nrow(cn_track_bed))
track_color <- ifelse(cn_track_bed$cn < 0 , "#00008B", "#8B0000") 

cn_track_bed$color <- track_color

mouse_chrs <- unique(mouse_bed$chr)
human_chrs <- unique(human_bed$chr)

colorVector2 <- rep("#000000", nrow(mouse_bed))
for (j in seq_along(mouse_chrs)) {
  colorVector2[which(mouse_bed$chr == mouse_chrs[j])] <- colorVector[j]
}



# so the way to have both frequency on human and segments on the same track
# is to use rect

# the big question now is how to have to rectangles with colored losses


circos.par("track.height"= 0.10) + 
  circos.initializeWithIdeogram(cytoband = cyto_combined_red, sort.chr = FALSE,
                                plotType = c("ideogram", "labels")) + 
  circos.genomicTrack(cn_track_bed, ylim = c(-1.00,1.00),
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, col = value$col,
                                    ytop = value$ytop, ybottom = value$ybot, border = NA,...)
                      }) + 
  circos.genomicLink(mouse_bed, human_bed,
                     col = colorVector2,
                     border = NA)



# need to find way to first create a line plot per chromosome
# then find inflection points - subset to a total of less than 5k per gains or losses
# note cn is not copy-number just mislabeling for the freq represened as height
# i.e freq of loss is -1 * freq

# idea is to make bed into line plot the looks for inflection points

library(data.table)

maximums <- function(x) which(x - shift(x, 10) > 0  & x - shift(x, 10, type='lead') > 0)
minimums <- function(x) which(x - shift(x, 10) < 0  & x - shift(x, 10, type='lead') < 0)

# I was an idiot and did they the smoothing + inflection point together
# they should be done separately

# maybe loess smoothing, inflection point finding, then make bed from close

#gain
tcga_ov_brca1_amp3 <- tcga_ov_brca1_amp2[which(tcga_ov_brca1_amp2$chr == "h_chr7"),]
tcga_ov_brca1_amp3 <- tcga_ov_brca1_amp3[1:(nrow(tcga_ov_brca1_amp3)-1),]
fit <- loess(start ~ cn, degree=1, span = 0.1, data=tcga_ov_brca1_amp3)
tcga_ov_brca1_amp3$loess <- fit$fitted
tcga_ov_brca1_amp3_red <- tcga_ov_brca1_amp3[c(maximums(tcga_ov_brca1_amp3$loess),
                                               minimums(tcga_ov_brca1_amp3$loess)),]
tcga_ov_brca1_amp3_red <- tcga_ov_brca1_amp3_red[order(tcga_ov_brca1_amp3_red$start,
                                                       decreasing = FALSE),]

red2Idx <- c(1,1+which(abs(diff(tcga_ov_brca1_amp3_red$cn)) >= 0.05))
h_chr1_amp_red2 <- NULL
for (i in 1:length(red2Idx)) {
  print(i)
  if(i == nrow(tcga_ov_brca1_amp3_red)){
    # if changepoint is on last region
    tmpIdx <- red2Idx_2[i]
    tmpVec <- c("chr" = tcga_ov_brca1_amp3_red$chr[tmpIdx],
                "start" = tcga_ov_brca1_amp3_red$start[tmpIdx],
                "end" = tcga_ov_brca1_amp3_red$end[tmpIdx], 
                "cn" = round(mean(tcga_ov_brca1_amp3_red$cn[tmpIdx]),digits = 2))
    h_chr1_amp_red2 <- rbind(h_chr1_amp_red2, tmpVec)
  } else if(i == length(red2Idx)){
    # last change point that is not the last region
    tmpIdx <- red2Idx[i]
    tmpVec <- c("chr" = tcga_ov_brca1_amp3_red$chr[tmpIdx],
                "start" = tcga_ov_brca1_amp3_red$start[tmpIdx],
                "end" = max(tcga_ov_brca1_amp3_red$end), 
                "cn" = round(tcga_ov_brca1_amp3_red$cn[tmpIdx],digits = 2))
    h_chr1_amp_red2 <- rbind(h_chr1_amp_red2, tmpVec)
  } else {
    # every other case
    tmpIdx <- red2Idx[i]
    tmpIdx2 <- red2Idx[i + 1]
    tmpVec <- c("chr" = tcga_ov_brca1_amp3_red$chr[tmpIdx],
                "start" = tcga_ov_brca1_amp3_red$start[tmpIdx],
                "end" = tcga_ov_brca1_amp3_red$start[tmpIdx2] - 1, 
                "cn" = round(mean(tcga_ov_brca1_amp3_red$cn[tmpIdx:tmpIdx2]),digits = 2))
    h_chr1_amp_red2 <- rbind(h_chr1_amp_red2, tmpVec)
  }
}

h_chr1_amp_red2 <- data.frame(h_chr1_amp_red2, stringsAsFactors = FALSE)
h_chr1_amp_red2[,2:4] <- lapply(h_chr1_amp_red2[,2:4], as.numeric)

#order check
h_chr1_amp_red2$diff <- h_chr1_amp_red2$start - h_chr1_amp_red2$end

# loss
tcga_ov_brca1_del3 <- tcga_ov_brca1_del2[which(tcga_ov_brca1_del2$chr == "h_chr7"),]
fit2 <- loess(start ~ cn, degree=1, span = 0.1, data=tcga_ov_brca1_del3)
tcga_ov_brca1_del3$loess <- fit2$fitted
tcga_ov_brca1_del3_red <- tcga_ov_brca1_del3[c(maximums(tcga_ov_brca1_del3$loess),
                                               minimums(tcga_ov_brca1_del3$loess)),]
tcga_ov_brca1_del3_red <- tcga_ov_brca1_del3_red[order(tcga_ov_brca1_del3_red$start, decreasing = FALSE),]
red2Idx_2 <- c(1,1+which(abs(diff(tcga_ov_brca1_del3_red$cn)) >= 0.05))
h_chr1_del_red2 <- NULL
for (i in 1:(length(red2Idx_2))) {
  print(i)
  if(i == nrow(tcga_ov_brca1_del3_red)){
    # if changepoint is on last region
    tmpIdx <- red2Idx_2[i]
    tmpVec <- c("chr" = tcga_ov_brca1_del3_red$chr[tmpIdx],
                "start" = tcga_ov_brca1_del3_red$start[tmpIdx],
                "end" = tcga_ov_brca1_del3_red$end[tmpIdx], 
                "cn" = round(mean(tcga_ov_brca1_del3_red$cn[tmpIdx]),digits = 2))
    h_chr1_del_red2 <- rbind(h_chr1_del_red2, tmpVec)
  } else if(i == length(red2Idx_2)){
    # last change point that is not the last region
    tmpIdx <- red2Idx_2[i]
    tmpVec <- c("chr" = tcga_ov_brca1_del3_red$chr[tmpIdx],
                "start" = tcga_ov_brca1_del3_red$start[tmpIdx],
                "end" = max(tcga_ov_brca1_del3_red$end), 
                "cn" = round(tcga_ov_brca1_del3_red$cn[tmpIdx],digits = 2))
    h_chr1_del_red2 <- rbind(h_chr1_del_red2, tmpVec)
  } else {
    # every other case
    tmpIdx <- red2Idx_2[i]
    tmpIdx2 <- red2Idx_2[i + 1]
    tmpVec <- c("chr" = tcga_ov_brca1_del3_red$chr[tmpIdx],
                "start" = tcga_ov_brca1_del3_red$start[tmpIdx],
                "end" = tcga_ov_brca1_del3_red$start[tmpIdx2] - 1, 
                "cn" = round(mean(tcga_ov_brca1_del3_red$cn[tmpIdx:tmpIdx2]),digits = 2))
    h_chr1_del_red2 <- rbind(h_chr1_del_red2, tmpVec)
  }
}

h_chr1_del_red2 <- data.frame(h_chr1_del_red2, stringsAsFactors = FALSE)
h_chr1_del_red2[,2:4] <- lapply(h_chr1_del_red2[,2:4], as.numeric)

h_chr1_del_red2$diff <- h_chr1_del_red2$start - h_chr1_del_red2$end 
# graphs showing conservation of pattern

arm1_amp <- ggplot(tcga_ov_brca1_amp2_red[1:10000,]) + geom_area(aes(x = start, y = cn))
arm2_amp <- ggplot(tcga_ov_brca1_amp2_red[10000:19697,]) + geom_area(aes(x = start, y = cn)) 
reduced_amp <- ggplot(h_chr1_amp_red2) + geom_area(aes(x = start, y = cn))

grid.arrange(arm1_amp, arm2_amp, reduced_amp,  ncol = 1)


arm1_del <- ggplot(tcga_ov_brca1_del2_red[1:15000,]) + geom_area(aes(x = start, y = cn))
arm2_del <- ggplot(tcga_ov_brca1_del2_red[15000:27028,]) + geom_area(aes(x = start, y = cn)) 
reduced_del <- ggplot(h_chr1_del_red2) + geom_area(aes(x = start, y = cn))

grid.arrange(arm1_del, arm2_del, reduced_del,  ncol = 1)
