# should be a list of functions needed to be used after
# processing segmentation (input used)
library(GenomicFeatures)
library(dplyr)
library(dbplyr)
library(circlize)
library(stringr)
library(DNAcopy)
library(stringr)
library(optparse)

# four variables should be mouse amplicon files, o/e ratio,
# normal vector name for FFPE vs FF, and output path

option_list = list(
  make_option(c("-a", "--amplicon"), type="character", default=NULL, 
              help="amplicon output", metavar="character"),
  make_option(c("-r", "--ratios"), type="character", default=NULL, 
              help="ratio output", metavar = ),
  make_option(c("-n", "--normals"), type="character", default=NULL, 
              help="should be either FFPE or FF", metavar="character"),
  make_option(c("-o", "--out"), type="character", default = NULL, 
              help="output path", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



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


# below are the two inputs needed for segmentation
# mouseAmplicons <- read.table("/mnt/DATA5/tmp/kev/newMouse2/cnAmplicon_matrix.txt", sep = "\t", stringsAsFactors = FALSE,
#                             header = TRUE)

#zscore_gc_oe_ratios <- read.table("/mnt/DATA5/tmp/kev/newMouse2/gcCorrectedCounts_matrix.txt", sep = "\t", stringsAsFactors = FALSE,
#                                  header = TRUE)

mouseAmplicons <- read.table(opt$amplicon, sep = "\t", stringsAsFactors = FALSE,
                             header = TRUE)

zscore_gc_oe_ratios <- read.table(opt$ratios,
                                  sep = "\t", stringsAsFactors = FALSE, header = TRUE)

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

mouseNormal <- NULL
if (opt$normal == "FFPE") {
  mouseNormal <- c("MG_17X49", "MG_18X50", "MG_23X55", "MG_6X38",
                   "MG_8X40","MG_11X43", "MG_13X45")
} else if(opt$normal == "FF"){
  # plaeceholder currently, no FF samples
  # still need a method to automatically 
  # differentiate in the pipeline
}


# used in converted human gene names to mouse
firstUpper <- function(gene){
  firstLetter <- toupper(substr(gene, start = 1, stop = 1))
  restOfGene <- tolower(substr(gene, start = 2, stop = nchar(gene)))
  res <- paste0(firstLetter, restOfGene)
  return(res)
}


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
  
  # when comparing to human, change to human specific cn
  
  mouse_bed <- data.frame("chr" = paste0("m_chr", synteny_hg38mm10$comp_chr[query_hits]),
                          "start" = synteny_hg38mm10$comp_start_pos[query_hits],
                          "end" = synteny_hg38mm10$comp_end_pos[query_hits],
                          "cn" = segInput$seg.mean[subject_hits],
                          stringsAsFactors = FALSE)
  
  human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg38mm10$ref_chr[query_hits]),
                          "start" = synteny_hg38mm10$ref_start_pos[query_hits],
                          "end" = synteny_hg38mm10$ref_end_pos[query_hits],
                          "cn" = segInput$seg.mean[subject_hits],
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

# only includes large chromosomal changes 
# the ony change is instead of genes I can use the genomic track instead 
# to represent gains or losses

segRes <- NULL
segRes <- ampSeg(mouseAmplicons, mouseBedFile)
#save(segRes, file = "/mnt/DATA5/tmp/kev/misc/20210409mm10SegRes.Robj")
#load(file = "/mnt/DATA5/tmp/kev/misc/20210409mm10SegRes.Robj")
segZscores <- calcZscore(segRes)
segZfilt <- segZscoresFilt(segRes, segZscores)


for(i in unique(segZfilt$ID)){
  segZfilt_sample <- segZfilt[which(segZfilt$ID == i),]
  syntenyBed <- syntenyPlotInputs(segZfilt_sample)
  human_bed <- syntenyBed$human_bed
  mouse_bed <- syntenyBed$mouse_bed
  
  cn_track_bed <- rbind(mouse_bed, human_bed)
  cn_track_bed2 <- cn_track_bed[, 1:3]
  
  cyto_interest <- unique(c(human_bed$chr, mouse_bed$chr))
  cyto_combined_red <- cyto_hg38mm10[which(cyto_hg38mm10$V1 %in% cyto_interest),]
  
  color_cn = colorRamp2(c(-1.0, 0, 1.0), c("#00008B", "#000000", "#8B0000"))
  track_color <- color_cn(cn_track_bed$cn)
  cn_track_bed$col = track_color
  
  mouse_chrs <- unique(mouse_bed$chr)
  human_chrs <- unique(human_bed$chr)
  
  colorVector2 <- rep("#000000", nrow(mouse_bed))
  for (j in seq_along(mouse_chrs)) {
    colorVector2[which(mouse_bed$chr == mouse_chrs[j])] <- colorVector[j]
  }
  
  # plotting each of the graphs
  pdf(paste0(opt$output,"circos_", i, ".pdf"), width = 6, height = 6, useDingbats = FALSE)
  #pdf(paste0("/mnt/DATA5/tmp/kev/newMouse2/","circos_", i, ".pdf"), width = 6, height = 6, useDingbats = FALSE)
  circos.par("track.height"= 0.05) + 
    circos.initializeWithIdeogram(cytoband = cyto_combined_red, sort.chr = FALSE) + 
    circos.genomicTrack(cn_track_bed, ylim = c(-1,1),
                        panel.fun = function(region, value, ...) {
                          circos.genomicLines(region, value, col = value$col, type = "segment",...)
                        }) + 
    circos.genomicLink(mouse_bed, human_bed,
                       col = colorVector2,
                       border = NA)
  dev.off()
}

# should also create the output for the seg output wiht z scores
# the segRes should be the output file snakemake uses

