### functions when needed to perform synteny analysis
library(data.table)
library(foreach)
library(doParallel)

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

segZscoresFilt_zeroOut <- function(segResults, zscores){
  
  # filts for significance, cn-change and length
  segmental_zscores <- cbind(segResults, "z-scores" = zscores[["z_vector"]],
                             "normal_seg_mean" = zscores[["nmean_vector"]],
                             "tumor_seg_mean" = zscores[["tmean_vector"]])
  segmental_zscores2 <- segmental_zscores
  segmental_zscores2$seg.mean[which(segmental_zscores2$num.mark < 10)] <- 0
  segmental_zscores2$length <- segmental_zscores2$loc.end - segmental_zscores2$loc.start
  largeChromCutoff <- 1000000
  segmental_zscores2$seg.mean[which(segmental_zscores2$length <= largeChromCutoff)] <- 0
  segmental_zscores2$seg.mean[which(abs(segmental_zscores2$`z-scores`) < 1.645)] <- 0
  segInput <- segmental_zscores2[,c("ID", "chrom", "loc.start", "loc.end", "seg.mean")]
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


syntenyPlotInputsFreq <- function(segInput){
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
                          "freq" = segInput$seg.mean[subject_hits],
                          stringsAsFactors = FALSE)
  
  human_bed <- data.frame("chr" = paste0("h_chr", synteny_hg38mm10$ref_chr[query_hits]),
                          "start" = synteny_hg38mm10$ref_start_pos[query_hits],
                          "end" = synteny_hg38mm10$ref_end_pos[query_hits],
                          stringsAsFactors = FALSE)
  
  
  
  resList <- list("human_bed" = human_bed,
                  "mouse_bed" = mouse_bed)
  
  return(resList)
}

getFreqData <- function(data){  
  #Check if segments or data:
  if(colnames(data)[1]=="sampleID" || colnames(data)[2]=="arm"){
    #input is segments data frame;
    #could be on a multiseg-format -> convert to uniseg-format:
    #need to convert to appropriate format
    #first find intersection of all breakpoints:
    chr <- unique(data[,2])
    bpts <- matrix(NA,nrow=0,ncol=2)
    for(j in 1:length(chr)){
      subseg <- subsetSegments(data,chrom=chr[j])
      bpts <- rbind(bpts,data.frame(chr[j],sort(unique(c(subseg$start.pos,subseg$end.pos)))))
    } 
    colnames(bpts) <- c("chrom","pos")
    
    #Then interpolate to get pcf-val in all breakpoints
    data = interpolate.pcf(data, bpts) 
  }
  #If not segments, the input data should already be on an appropriate format (chrom, pos, estimates)
  
  return(data) 
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


reducingFreqBed <- function(df, idx){
  reducedDf <- NULL
  for (i in 1:(length(idx) - 1)) {
    idx1 <- idx[i]
    idx2 <- idx[i+1]
    tmpDf <- df[idx1:idx2,]
    tmpChr <- tmpDf$chr[1]
    tmpStart <- min(tmpDf$pos)
    tmpEnd <- max(tmpDf$pos) - 1
    tmpFreq <- tmpDf[,3][1]
    
    tmpVec <- c("Chr" = tmpChr, "Start" = tmpStart, "End" = tmpEnd, "Freq" = tmpFreq)
    reducedDf <- rbind(reducedDf, tmpVec)
  }
  reducedDf <- data.frame(reducedDf, stringsAsFactors = FALSE)
  return(reducedDf)
}

# used on df with segmentation to get frequency bed
ampsDels <- function(df){
  getFreqOut <- getFreqData(df)
  freqAmp <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] > 0.2)*100
  freqDel <- rowMeans(getFreqOut[,-c(1:2),drop=FALSE] < -0.2)*100
  
  freqDf_Amp <- cbind(getFreqOut[,1:2], "amp" = freqAmp)
  freqDf_Del <- cbind(getFreqOut[,1:2], "del" = freqDel)
  
  ampRedIdx <- c(1,1+which(diff(freqDf_Amp$amp)!=0))
  delRedIdx <- c(1,1+which(diff(freqDf_Del$del)!=0))
  
  res <- list(freqDf_Amp, ampRedIdx, freqDf_Del, delRedIdx)
}


round2 = function(x, n =  0) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5 + sqrt(.Machine$double.eps)
  z = trunc(z)
  z = z/10^n
  z*posneg
}

separateSegments_tcga <- function(df, ploidyDf){
  # separates the segments prior to generating frequencies - diff for tcga data and amplicon
  # whole chromosome/arms (covers ~80% of max chr coverage from panel); changes > 1MB
  res <- NULL
  res2 <- NULL
  
  # this was genius and got rid of the loop - literally infinitely faster *pats self on back*
  df$newTotalCN <- round2(df$Copynumber - ploidyDf$ploidy[match(df$Sample, ploidyDf$array)], 0)
  
  for (i in unique(df$Chromosome)) {
    tmpChrDf <- df[which(df$Chromosome == i),]
    tmpChrDf$newTotalCN[which(tmpChrDf$length < 1000000)] <- 0
    armTable <- tmpChrDf 
    cnaTable <- tmpChrDf
    tmpCyto <- human_arm[which(human_arm$chromosome == i),]

    minLength <- min(tmpCyto$length80)
    
    armTable$newTotalCN[which(tmpChrDf$length <= minLength)] <- 0
    cnaTable$newTotalCN[which(tmpChrDf$length > minLength)] <- 0
    
    res <- rbind(res, armTable)
    res2 <- rbind(res2, cnaTable)
  }
  return(list(res, res2))
}

separateSegments_m <- function(df){
  # separates the segments prior to generating frequencies - diff for tcga data and amplicon
  # whole chromosome/arms (covers ~80% of max chr coverage from panel); changes > 1MB
  
}


synteny_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_syntenyDf.txt",
                               sep = "\t", stringsAsFactors = FALSE, header = TRUE)
cyto_hg38mm10 <- read.table("/mnt/DATA6/kevin_recovery/apps/circos/20210401_mm10_hg38_cytoDf.txt",
                            sep = "\t", stringsAsFactors = FALSE, header = FALSE)
#hg38tohg19_chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/hg38ToHg19.over.chain")

tcga_seg <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/TCGA_mastercalls.abs_segtabs.fixed.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE)

tcga_ploidy <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
                          sep = "\t", stringsAsFactors = FALSE, header = TRUE)

human_cyto <- cyto_hg38mm10[grep("h_chr",cyto_hg38mm10$V1),]
human_arm <- NULL
for (i in unique(human_cyto$V1)) {
  tmpDf <- human_cyto[which(human_cyto$V1 == i),]
  centIdx <- which(tmpDf$V5 == "acen")
  arm1 <- c("chromosome" = tmpDf$V1[1], "start" = tmpDf$V2[1],
            "end" = tmpDf$V2[min(centIdx) - 1])
  arm2 <- c("chromosome" = tmpDf$V1[1], "start" = tmpDf$V2[max(centIdx) + 1],
            "end" = tmpDf$V2[nrow(tmpDf)])
  human_arm <- rbind(human_arm, arm1, arm2)
}

human_arm <- data.frame(human_arm, stringsAsFactors = FALSE)
human_arm$start <- as.numeric(human_arm$start)
human_arm$end <- as.numeric(human_arm$end)
human_arm$chromosome <- str_remove(human_arm$chromosome, "h_chr")
human_arm$length <- human_arm$end  - human_arm$start
human_arm$length80 <- human_arm$length * 0.8


hg38ChromLocs <- NULL
for (i in unique(human_cyto$V1)) {
  tmpDf <- human_cyto[which(human_cyto$V1 == i),]
  hg38ChromLocs <- rbind(hg38ChromLocs, c("chrom" = str_remove(tmpDf$V1[1], "h_chr"),
                                             "start" = min(tmpDf$V2), "end" = max(tmpDf$V3)))
}

hg38ChromLocs <- data.frame(hg38ChromLocs, stringsAsFactors = FALSE)
hg38ChromLocs$start <- as.numeric(hg38ChromLocs$start)
hg38ChromLocs$end <- as.numeric(hg38ChromLocs$end)

write.table(hg38ChromLocs, "/mnt/DATA5/tmp/kev/misc/20210628chromList_hg38.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

### honestly most of the large blocks don't change - get exact coordinates later
loRes <- NULL
indexToSkip <- NULL
for (i in 1:nrow(synteny_hg38mm10)) {
  tmpGrange <- GRanges(seqnames = paste0("chr",synteny_hg38mm10$ref_chr[i]),
                       IRanges(start = synteny_hg38mm10$ref_start_pos[i],
                               end = synteny_hg38mm10$ref_end_pos[i]))
  tmpLo <- data.frame(liftOver(tmpGrange, hg38tohg19_chain))
  if (nrow(tmpLo) == 0) {
    indexToSkip <- c(indexToSkip, i)
    next()
  }
  tmpLo$hg38start <- synteny_hg38mm10$ref_start_pos[i]
  tmpLo$hg38end <- synteny_hg38mm10$ref_end_pos[i]
  loRes <- rbind(loRes, tmpLo)
}


cl <- makeCluster(20)
registerDoParallel(cl)
finalTable <- NULL
finalTable <- foreach(i = unique(tcga_seg$Sample), .packages = "data.table", .combine = "rbind") %dopar% {
  res <- NULL
  sampleDf <- tcga_seg[which(tcga_seg$Sample == i),]
  for (j in unique(sampleDf$Chromosome)) {
    sampleDf_chr <- sampleDf[which(sampleDf$Chromosome == j),]
    indices <- c(1,1+which(abs(diff(sampleDf_chr$Modal_Total_CN)) >= 1))
    tmpDf <- data.frame("Sample" = rep(i,length(indices)), "Chromosome" = rep(j,length(indices)),
                        "Start" = sampleDf_chr$Start[indices],
                        "End" = sampleDf_chr$End[shift(indices,1, type = "lead")-1],
                        "Copynumber" = sampleDf_chr$Modal_Total_CN[indices])
    # because of how i use shift and lead, without doubling the last entry, I get NA
    tmpDf$End[nrow(tmpDf)] <- max(sampleDf_chr$End)
    res <- rbind(res, tmpDf)
  }
  res
}
stopCluster(cl)


finalTable$Sample <- as.character(finalTable$Sample)
write.table(finalTable, "/mnt/DATA5/tmp/kev/misc/20210628panCancerAllSegReduced_hg19.txt", quote = FALSE, col.names = TRUE, row.names = FALSE,
            sep = "\t")

finalTable <-  read.table("/mnt/DATA5/tmp/kev/misc/20210628panCancerAllSegReduced_hg19.txt", sep = "\t",
                          stringsAsFactors = FALSE, header = TRUE)

finalTable$length <- finalTable$End - finalTable$Start
finalTable2 <- finalTable[which(finalTable$Chromosome %in% c(1:22)),]
armRes <- separateSegments_tcga(finalTable2, tcga_ploidy)
armTable <- armRes[[1]]
cnaTable <- armRes[[2]]

write.table(armTable, "/mnt/DATA5/tmp/kev/misc/20210628panCancerArmTable_hg19.txt", sep = "\t", col.names = TRUE,
           row.names = FALSE, quote = FALSE)
write.table(cnaTable, "/mnt/DATA5/tmp/kev/misc/20210628panCancerCnaTable_hg19.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)

### overall test I should have is Cochran-Mazel Test (Chi-squared) or stratified
### for stratification just make a table for amplified and deleted genes per pathway
### 2x2 table amp & deleted by mouse and human. also note this is hg19, not hg38
