library(GenomicFeatures)
library(dplyr)
library(dbplyr)
library(circlize)
library(stringr)

firstUpper <- function(gene){
  firstLetter <- toupper(substr(gene, start = 1, stop = 1))
  restOfGene <- tolower(substr(gene, start = 2, stop = nchar(gene)))
  res <- paste0(firstLetter, restOfGene)
  return(res)
}

### getting syntenic blocks 
syntenic_blocks <- DBI::dbConnect(RSQLite::SQLite(), "/mnt/DATA5/tmp/kev/programs/syntenybrowser/synteny.db")
src_dbi(syntenic_blocks)

syntenic_blocks_data <- tbl(syntenic_blocks, "syntenic_block")
head(syntenic_blocks_data)

humanStart <- syntenic_blocks_data %>%
  filter(ref_taxonid == 9606) %>%
  select(ref_taxonid, ref_chr, ref_start_pos, ref_end_pos, comp_taxonid,
         comp_chr, comp_start_pos, comp_end_pos, symbol)

humanStart <- as.data.frame(humanStart)

hg38 <- read.cytoband(species = "hg38")
hg38_cyto <- hg38$df
hg38_cyto$tmpString <- paste0(hg38_cyto$V1, hg38_cyto$V4)

hg38toMm10Chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/hg38ToMm10.over.chain")

hg38biomartTable <- read.table("/home/kevhu/data/20201030hg38KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

mm10biomartTable <- read.table("/home/kevhu/data/20201030Mm10KnownCanbiomartQuery.txt", header = TRUE,
                               stringsAsFactors = FALSE, sep = "\t")

ov_cna_freq <- read.table("/mnt/DATA5/tmp/kev/misc/20210320tcgaOvAltGenesFreq.txt",
                          sep =  "\t", stringsAsFactors = FALSE, skip = 1)
ov_cna_freq <- ov_cna_freq[,c(1,3:ncol(ov_cna_freq))]
colnames(ov_cna_freq) <- c("gene", "cytoband", "CNA", "Profiled", "SampleNumber",
                           "Freq", "IsCancerGene")
ov_cna_freq$Freq <- as.numeric(str_remove_all(ov_cna_freq$Freq, "\\%"))
ov_cna_freq_2 <- ov_cna_freq[which(ov_cna_freq$IsCancerGene == "Yes"),]
ov_cna_freq_2 <- ov_cna_freq_2[which(ov_cna_freq_2$Freq > 10),]

ov_cyto <- ov_cna_freq_2$cytoband[-grep("-", ov_cna_freq_2$cytoband)]
ov_cyto2 <- ov_cna_freq_2$cytoband[grep("-", ov_cna_freq_2$cytoband)]
ov_cyto <- c(ov_cyto, "3q28", "3q29", "3q27.3", "3q28", "19p13.12", "19p13.11")
ov_cyto <- ov_cyto[-which(duplicated(ov_cyto))]

ov_freq_df <- data.frame("chr" = substr(ov_cyto, 1, 1), 
                         "cyto" = substr(ov_cyto, 2, nchar(ov_cyto)))

ov_freq_df$tmpString <- paste0("chr" ,ov_freq_df$chr, ov_freq_df$cyto)
cyto_coords <- which(hg38_cyto$tmpString %in% ov_freq_df$tmpString)
ov_freq_df$start <- hg38_cyto$V2[cyto_coords]
ov_freq_df$end <- hg38_cyto$V3[cyto_coords]

### use grange to find intersection of cytoband locations and syteny blocks
### just need to find which blocks I need to map for the circos plot
### bed dataframes are created for links on circos plot

ov_freq_ranges <- GRanges(seqnames = paste0("chr", ov_freq_df$chr),
                          ranges = IRanges(start = ov_freq_df$start,
                                           end = ov_freq_df$end))

synteny_ranges <- GRanges(seqnames = paste0("chr", humanStart$ref_chr),
                          ranges = IRanges(start = humanStart$ref_start_pos,
                                           end = humanStart$ref_end_pos))


res_overlap <- findOverlaps(synteny_ranges, ov_freq_ranges)
query_nondup <- queryHits(res_overlap)
query_nondup <- query_nondup[-which(duplicated(query_nondup))]

human_bed <- data.frame("chr" = humanStart$ref_chr[query_nondup],
                        "start" = humanStart$ref_start_pos[query_nondup],
                        "end" = humanStart$ref_end_pos[query_nondup],
                        stringsAsFactors = FALSE)

human_bed$chr <- paste0("h_chr", human_bed$chr)

mouse_bed <- data.frame("chr" = humanStart$comp_chr[query_nondup],
                        "start" = humanStart$comp_start_pos[query_nondup],
                        "end" = humanStart$comp_end_pos[query_nondup],
                        stringsAsFactors = FALSE)

mouse_bed$chr <- paste0("m_chr", mouse_bed$chr)

### below for obtaining gene locations for human and mice
### this will be a link table in the shiny app, just making sure
### getting genes as rectangles looks as intended

ov_cna_freq_gene <- ov_cna_freq[which(ov_cna_freq$cytoband %in% c(ov_cyto, ov_cyto2)),]
ov_cna_freq_gene <- ov_cna_freq_gene[which(ov_cna_freq_gene$IsCancerGene == "Yes"),]
genesToMap <- ov_cna_freq_gene$gene

hgGeneBed <- NULL
for (i in genesToMap) {
  tmpDf <- hg38biomartTable[which(hg38biomartTable$external_gene_name == i),]
  tmpStart <- min(tmpDf$exon_chrom_start)
  tmpEnd <- max(tmpDf$exon_chrom_end)
  tmpChrom <- tmpDf$chromosome_name[1]
  
  hgGeneBed <- rbind(hgGeneBed, c(i, tmpChrom, tmpStart, tmpEnd))
}

hgGeneBed <- data.frame(hgGeneBed, stringsAsFactors = FALSE)
hgGeneBed <- hgGeneBed[-which(is.na(hgGeneBed$X2)),]
hgGeneBed$X2 <- paste0("h_chr", hgGeneBed$X2)

mmGeneBed <- NULL
for (i in hgGeneBed$X1) {
  tmpDf <- mm10biomartTable[which(mm10biomartTable$external_gene_name == firstUpper(i)),]
  tmpStart <- min(tmpDf$exon_chrom_start)
  tmpEnd <- max(tmpDf$exon_chrom_end)
  tmpChrom <- tmpDf$chromosome_name[1]
  
  mmGeneBed <- rbind(mmGeneBed, c(i, tmpChrom, tmpStart, tmpEnd))
}

mmGeneBed <- data.frame(mmGeneBed, stringsAsFactors = FALSE)
mmGeneBed <- mmGeneBed[-which(is.na(mmGeneBed$X2)),]
mmGeneBed$X2 <- paste0("m_chr", mmGeneBed$X2)


geneCombinedBed <- rbind(mmGeneBed, hgGeneBed)
geneCombinedBed <- geneCombinedBed[,2:4]
colnames(geneCombinedBed) <- c("chr", "start", "end")
#geneCombinedBed$value <- runif(n = nrow(geneCombinedBed),min = -1,max = 1)
geneCombinedBed$start <- as.numeric(geneCombinedBed$start)
geneCombinedBed$end <- as.numeric(geneCombinedBed$end)



mm10 <- read.cytoband(species = "mm10")
mm10$df$V1 <-  paste0("m_", mm10$df$V1)
hg38 <- read.cytoband(species = "hg38")
hg38$df$V1 <- paste0("h_", hg38$df$V1)

cyto_combined <- rbind(mm10$df, hg38$df)
length(grep("h_", unique(cyto_combined$V1)))
length(grep("m_", unique(cyto_combined$V1)))



cyto_interest <- c(unique(human_bed$chr), unique(mouse_bed$chr))
cyto_combined_red <- cyto_combined[which(cyto_combined$V1 %in% cyto_interest),]


# shiny app steps

### automate reduction of to chromosomes used
### make sure to automate the creation of the color vector
### also automate the creation of the gaps
### last thing is to test out gene markers 
### 

dev.off()
circos.clear()
#circos.par("gap.degree" = c(rep(2, 23), 10, rep(2,20), 10))
circos.par("track.height"= 0.05)
circos.initializeWithIdeogram(cytoband = cyto_combined_red)
circos.genomicTrack(geneCombinedBed, ylim = c(0,0.5),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop = 0, ybottom = 0.5)
                    })
circos.genomicLink(human_bed, mouse_bed,
                   col = c(rep("#A84A32",6), rep("#32A869",11), rep("#A8327F",10)),
                   border = NA)


### going to try and do this by chromosome now
chr8_idx <- which(humanStart$ref_chr == "8")
human_bed_8 <- data.frame("chr" = humanStart$ref_chr[chr8_idx],
                          "start" = humanStart$ref_start_pos[chr8_idx],
                          "end" = humanStart$ref_end_pos[chr8_idx])
human_bed_8$chr <- paste0("h_chr", human_bed_8$chr)
mouse_bed_8 <- data.frame("chr" = humanStart$comp_chr[chr8_idx],
                          "start" = humanStart$comp_start_pos[chr8_idx],
                          "end" = humanStart$comp_end_pos[chr8_idx])
mouse_bed_8$chr <- paste0("m_chr", mouse_bed_8$chr)


dev.off()
circos.clear()
circos.par("gap.degree" = c(rep(2, 23), 10, rep(2,20), 10))
circos.initializeWithIdeogram(cytoband = cyto_combined)
circos.genomicLink(human_bed_8, mouse_bed_8, col = rand_color(nrow(human_bed_8), transparency = 0.5),
                   border = NA)


