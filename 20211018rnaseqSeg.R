library(copynumber)
library(stringr)
library(gg3D)

### code to get z-scores from RNA


calcZscore <- function(segRes){
  sampZ1 <- NULL
  sampZ2 <- NULL
  z1_vector <- NULL
  z2_vector <- NULL
  p_vector1 <- NULL
  p_vector2 <- NULL
  nmean_vector <- NULL
  tmean_vector <- NULL
  q_vector1 <- NULL
  q_vector2 <- NULL
  
  uniqIdx <- unique(match(str_replace_all(segRes$ID, "[[:punct:]]", ""), 
                          str_replace_all(str_replace_all(colnames(zscore_gc_oe_ratios), "[[:punct:]]", ""), " ", "")))
  zScoreColnames <- colnames(zscore_gc_oe_ratios)[uniqIdx]
  for (i in seq_along(unique(segRes$ID))) {
    tmp_seg_res <- segRes[which(segRes$ID == unique(segRes$ID)[i]), ]
    tmp_seg_ranges <- GRanges(seqnames = Rle(tmp_seg_res$chrom), 
                              ranges = IRanges(start = tmp_seg_res$loc.start,
                                               end = tmp_seg_res$loc.end))
    print(unique(segRes$ID)[i])
    print(zScoreColnames[i])
    
    sampZ1 <- NULL
    sampZ2 <- NULL
    
    for (j in seq_along(tmp_seg_ranges)) {
      tmpOverlap <- findOverlaps(all_probes_grange, tmp_seg_ranges[j])
      tmpTest <- zscore_gc_oe_ratios[[zScoreColnames[i]]][queryHits(tmpOverlap)]
      
      if (!is.numeric(tmpTest)) {
        sampZ1 <- c(sampZ1, NA)
        sampZ2 <- c(sampZ2, NA)
        nmean_vector <- c(nmean_vector, NA)
        tmean_vector <- c(tmean_vector, NA)
        next()
      } else{
        tmpNormal <- zscore_gc_oe_ratios[queryHits(tmpOverlap), mouseNormal]
        tmpNormal2 <- unlist(apply(tmpNormal, 2, mean))
        
        z1_stat <- (mean(tmpTest) - mean(tmpNormal2))/sd(tmpNormal2)
        z2_stat <- (mean(tmpTest) - mean(tmpNormal2))/(sd(tmpNormal2) * mean(tmpTest))
        pVal_z1 <- 2.0*pnorm(-abs(z1_stat))
        pVal_z2 <- 2.0*pnorm(-abs(z2_stat))
        
        
        normal_seg_mean <- mean(tmpNormal2)
        tumor_seg_mean <- mean(tmpTest)
        sampZ1 <- c(sampZ1, pVal_z1)
        sampZ2 <- c(sampZ2, pVal_z2)
        z1_vector <- c(z1_vector, z1_stat)
        z2_vector <- c(z2_vector, z2_stat)
        p_vector1 <- c(p_vector1, pVal_z1)
        p_vector2 <- c(p_vector2, pVal_z2)
        nmean_vector <- c(nmean_vector, normal_seg_mean)
        tmean_vector <- c(tmean_vector, tumor_seg_mean)
      }
    }
    qVal_z1 <- (p.adjust(sampZ1, method = "BH"))
    qVal_z2 <- (p.adjust(sampZ2, method = "BH"))
    q_vector1 <- c(q_vector1, qVal_z1)
    q_vector2 <- c(q_vector2, qVal_z2)
    
  }
  return(list("z1_vector" = z1_vector, "z2_vector" = z2_vector,
              "p1_vector" = p_vector1, "p2_vector" = p_vector2,
              "q1_vector" = q_vector1, "q2_vector" = q_vector2,
              "nmean_vector"= nmean_vector, "tmean_vector" = tmean_vector))
}

segZscoresFilt <- function(segResults, zscores){
  
  # filts for significance, cn-change and length
  segmental_zscores <- cbind(segResults, "z-scores1" = zscores[["z1_vector"]], "z-scores2" = zscores[["z2_vector"]],
                             "p-val1" = zscores[["p1_vector"]], "p-val2" = zscores[["p2_vector"]],
                             "q-val1" = zscores[["q1_vector"]], "q-val2" = zscores[["q2_vector"]],
                             "normal_seg_mean" = zscores[["nmean_vector"]], "tumor_seg_mean" = zscores[["tmean_vector"]])
  return(segmental_zscores)
}



tcDf <- read.table("/mnt/DATA5/tmp/kev/misc/20210718hgscTcDf.txt", sep = "\t", stringsAsFactors = FALSE,
                   header = TRUE)

setwd("/mnt/DATA5/tmp/kev/programs/igordot_copynumber/")
files.sources = list.files()
files.sources <- files.sources[-which(files.sources == "sysdata.rda")]
load("sysdata.rda")
sapply(files.sources, source)


# tmp <- read.table("/mnt/DATA6/mouseData/mouseBedsIdx.txt",sep = "\t", header = TRUE)
# 
# allSum <- NULL
# allSum2 <- NULL
# for (i in tmp$summaryFile) {
#   tmpFile <- read.table(i, sep = "\t", header = TRUE)
#   if (ncol(tmpFile) < 4) {
#     next()
#   } else if (ncol(tmpFile) == 6) {
#     allSum <- rbind(allSum, tmpFile)
#   } else if (ncol(tmpFile) == 7) {
#     allSum2 <- rbind(allSum2, tmpFile)
#   }
#   #allSum <- rbind(allSum, tmpFile)
# }
# 

### need to double check the tumor sequenced match - i.e just says tumor and not left and right
# rnaSeqSamples <- c("2285LOTX50", "2519LOTX52", "2611LOTX53", 
#                    "2611ROTX54", "2942ROTX57", "10879LT_X58")
rsemDir <- '/mnt/DATA5/tmp/kev/tmpDbs/starIndexes/'
setwd(rsemDir)
rsem_files <- system('ls *.rsem.genes.results', intern = TRUE)
rsem_files <- paste0(rsemDir, rsem_files)

rsemValTable <- NULL
listOfSRRs <- NULL
for (i in seq_along(rsem_files)) {
  sampleName <- str_remove(str_remove(rsem_files[i], "/mnt/DATA5/tmp/kev/tmpDbs/starIndexes/"), "\\.rsem.*")
  tmpTable <- read.table(rsem_files[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  ### col 5 is expected counts, 6, is TPM, 7 is FPKM
  tmpTable <- tmpTable[,c(1,2,7)]
  if (i == 1) {
    rsemValTable <- rbind(rsemValTable, tmpTable)
  } else {
    rsemValTable <- cbind(rsemValTable, tmpTable[,3])
  }
  
  listOfSRRs <- c(listOfSRRs, sampleName)
}

colnames(rsemValTable)[3:ncol(rsemValTable)] <- listOfSRRs

listOfGroups <- list(c("SRR9933670", "SRR9933671", "SRR9933672"),
                     c("SRR9933673", "SRR9933674", "SRR9933675"),
                     c("SRR9933676", "SRR9933677", "SRR9933678"),
                     c("SRR9933679", "SRR9933680", "SRR9933681"),
                     c("SRR9933697", "SRR9933698", "SRR9933699"),
                     c("SRR9933706", "SRR9933707", "SRR9933708"),
                     c("SRR9933751", "SRR9933752", "SRR9933753"),
                     c("SRR9933754", "SRR9933755", "SRR9933756"),
                     c("SRR9933772", "SRR9933773", "SRR9933774"))

listOfSampleNames <- c("Ov1", "Ov2", "Ov3", "Ov4",
                "2611lt","2942rt", "2285lt", "2519lt",
                "10879lt") 

newRsemTable <- rsemValTable[,1:2]
for (i in seq_along(listOfGroups)) {
  newCol <- apply(rsemValTable[,listOfGroups[[i]]], 1, sum)
  newRsemTable <- cbind(newRsemTable, newCol)
}


colnames(newRsemTable)[3:ncol(newRsemTable)] <- listOfSampleNames
rsemValTable <- newRsemTable

### order is done so i can calculate non-expression, make 1 min value for any expressed genes




rsemValTable[,3:ncol(rsemValTable)] <- log2(rsemValTable[,3:ncol(rsemValTable)])
rsemMat <- as.matrix(rsemValTable[,3:ncol(rsemValTable)])
rsemMat[is.infinite(rsemMat)] <- 0
rsemMat_TF <- rsemMat >= 1
expressedRatios <- apply(rsemMat_TF, 1, mean)
rsemMat[rsemMat < 1] <- 1
rsemValTable[,3:ncol(rsemValTable)] <- rsemMat
rsemValTable2 <- rsemValTable[-which(expressedRatios < 0.2),]

medNormVals <- apply(rsemValTable2[,c(3:6)], 1, median)

for (i in 3:ncol(rsemValTable2)) {
  rsemValTable2[,i] <- rsemValTable2[,i] - medNormVals
}


geneVariability <- apply(rsemValTable2[,3:ncol(rsemValTable2)], 1, var)
quantile( geneVariability, seq(0,1, 0.05))

rsemValTable2 <- rsemValTable2[-which(geneVariability >= 4.54),]


# added tc correction
matchingTc <- tcDf$tc[match(tolower(colnames(rsemValTable2)[3:11]), tcDf$sample)]
rsemValTable2[3:11] <- 2^rsemValTable2[3:11]
for (i in which(!is.na(matchingTc))) {
  ### columns are off for rsemTable, add 2 
  tmpVec <- rsemValTable[,i + 2]
  tmpVec[tmpVec > 1] <- tmpVec[tmpVec > 1] / matchingTc[i]
  tmpVec[tmpVec < 1] <- tmpVec[tmpVec < 1] * matchingTc[i]
  rsemValTable[,i + 2] <- tmpVec
}
rsemValTable2[3:11] <- log2(rsemValTable2[3:11])


# end of tc-correction if wanted to change back

### query unique geneids with biomart
library(biomaRt)
uniqGeneIds <- unique(rsemValTable2$gene_id)

ensemblMm <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# filtsMm <- listFilters(ensemblMm)
# attrMm <- listAttributes(ensemblMm)

# attr_tmp <- c("external_gene_name", "chromosome_name", "start_position", "end_position")
# 
# tmpMm <- getBM(mart = ensemblMm, attributes = attr_tmp, filters = "wikigene_name",
#                values = uniqGeneIds)

attr_tmp2 <- c("wikigene_name", "chromosome_name", "start_position", "end_position")
tmpMm2 <- getBM(mart = ensemblMm, attributes = attr_tmp2, filters = "wikigene_name",
               values = uniqGeneIds)
tmpMm2_1 <- tmpMm2[-which(duplicated(tmpMm2$wikigene_name)),]

rsemValTable3 <- rsemValTable2[which(rsemValTable2$gene_id %in% tmpMm2_1$wikigene_name), ]
rsemValTable3 <- cbind(tmpMm2_1[match(rsemValTable3$gene_id,tmpMm2_1$wikigene_name),], 
                       rsemValTable3)
chromNames <- c(1:19, "X")
rsemValTable3 <- rsemValTable3[which(rsemValTable3$chromosome_name %in% chromNames), ]

rsemValTable3$pos <- apply(rsemValTable3[,3:4], 1, mean)
rsemValTable4 <- rsemValTable3[-which(duplicated(paste(rsemValTable3$chromosome_name, rsemValTable3$pos))), ]


allPcf <- NULL
for (i in 7:(ncol(rsemValTable4)-1)) {
  tmp <- rsemValTable3[,c(2,ncol(rsemValTable4),i)]
  tmp <- tmp[order(tmp$chromosome_name, tmp$pos), ]
  colnames(tmp)[1] <- "chrom"
  tmp.wins <- winsorize(data= tmp,verbose=FALSE, assembly = "mm10", gamma = 18, k = 30)
  tmp.seg <- pcf(data= tmp.wins,verbose=FALSE, gamma = 18,
                 assembly = "mm10", fast = FALSE)
  allPcf <- rbind(allPcf, tmp.seg)
}

colnames(allPcf)[c(1,2,4,5,7)] <- c("ID", "chrom", "loc.start", "loc.end", "seg.mean")


mouseNormal <- c("Ov1", "Ov2", "Ov3", "Ov4")
all_probes_grange <- GRanges(seqnames = Rle(rsemValTable4$chromosome_name),
                             ranges = IRanges(start = rsemValTable4$start_position,
                                              end = rsemValTable4$end_position))
zscore_gc_oe_ratios <- rsemValTable4
allRnaPcf_zscore <- calcZscore(allPcf)
allRnaPcf_zscoreFilt <- segZscoresFilt(allPcf, allRnaPcf_zscore)



write.table(allRnaPcf_zscoreFilt, "/mnt/DATA5/tmp/kev/misc/20211025allPcfRnaseq_zscore_tc.txt", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)



### trying cbs
library(DNAcopy)

allRnaSeg <- NULL
for (i in 7:(ncol(rsemValTable4)-1)) {
  tmpCNA_obj <- CNA(cbind(rsemValTable4[,i]),
                    rsemValTable4$chromosome_name, rsemValTable4$pos,
                    data.type="logratio",sampleid=colnames(rsemValTable4)[i])
  smoothed_tmpCNA <- smooth.CNA(tmpCNA_obj)
  segment_tmpCNA <- segment(smoothed_tmpCNA, verbose = 1, undo.splits = "sdundo",
                            p.method = c("hybrid"), undo.SD = 2, min.width = 2)
  
  
  tmpCalls <- segments.p(segment_tmpCNA)
  allRnaSeg <- rbind(allRnaSeg, tmpCalls)
}

allRnaSeg$ID <- str_remove(allRnaSeg$ID, "^X")


mouseNormal <- c("Ov1", "Ov2", "Ov3", "Ov4")
all_probes_grange <- GRanges(seqnames = Rle(rsemValTable4$chromosome_name),
                             ranges = IRanges(start = rsemValTable4$start_position,
                                              end = rsemValTable4$end_position))
zscore_gc_oe_ratios <- rsemValTable4
allRnaSeg_zscore <- calcZscore(allRnaSeg)
allRnaSeg_zscoreFilt <- segZscoresFilt(allRnaSeg, allRnaSeg_zscore)


write.table(allRnaSeg_zscoreFilt, "/mnt/DATA5/tmp/kev/misc/20211109allCbsRnaseq_zscore_tc.txt", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)

