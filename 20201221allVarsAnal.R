source("/mnt/DATA6/kevin_recovery/scripts/fastReadFile.R")

library(GenomicRanges)
crossTableIdx <- readxl::read_xlsx("/home/kevhu/data/20201207annotations.xlsx")
crossTableIdx <- crossTableIdx[,1:3]
crossTableIdx$old_name2 <- tolower(crossTableIdx$old_name)
crossTableIdx$old_name2 <- str_remove(crossTableIdx$old_name2 , "[[:punct:]]")

convertedCosmic <- faster.readfile("/mnt/DATA5/tmp/kev/misc/20201103convertedCosmic_par.txt", 0)
convertedCosmic <- convertedCosmic[-which(convertedCosmic$hgGene == "hgGene"),]
convertedCosmic <- convertedCosmic[-which(convertedCosmic$hgGene == "CEL"),]

reversed <- which((as.numeric(convertedCosmic$mm10end) - as.numeric(convertedCosmic$mm10start)) < 0)
reverseStrandsEnds <- convertedCosmic$mm10end[reversed]
convertedCosmic$mm10end[reversed] <- convertedCosmic$mm10start[reversed]
convertedCosmic$mm10start[reversed] <- reverseStrandsEnds
convertedCosmicGrange <- GRanges(seqnames=Rle(convertedCosmic$mm10chr),
                                 ranges=IRanges(as.numeric(convertedCosmic$mm10start),
                                                as.numeric(convertedCosmic$mm10end)))



covertedCosmic_hgChr <- paste0("chr",unlist(lapply(strsplit(convertedCosmic$Gr38_position, ":"), "[[", 1)))
covertedCosmic_hgPos <- unlist(lapply(strsplit(convertedCosmic$Gr38_position, ":"), "[[", 2))
covertedCosmic_hgStart <- unlist(lapply(strsplit(covertedCosmic_hgPos, "-"), "[[", 1))
covertedCosmic_hgEnd <- unlist(lapply(strsplit(covertedCosmic_hgPos, "-"), "[[", 2))

convertedCosmicGrange_hg <- GRanges(seqnames=Rle(covertedCosmic_hgChr),
                                 ranges=IRanges(as.numeric(covertedCosmic_hgStart),
                                                as.numeric(covertedCosmic_hgEnd)))

### I think for this to work better - I just map each of the converted mutations to a protein domain location 
### and use that instead of just 3 or 5 surrouding base pair. more justification than only doing arbitrary amounts of amino acids
prot2geneHg <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/prot2Hg/download",
                          stringsAsFactors = FALSE, header = TRUE, sep = ";")

mappedChrs <- c(paste0("chr", 1:22), "chrX", "chrY")
prot2geneHg  <- prot2geneHg [which(prot2geneHg$X38_chr %in% mappedChrs),]

prot2geneHgGrange <- GRanges(seqnames = Rle(prot2geneHg$X38_chr),
                             ranges = IRanges(as.numeric(prot2geneHg$X38_chr_start),
                                              as.numeric(prot2geneHg$X38_chr_end)))


subjectRes_cosmic <- subjectHits(findOverlaps(prot2geneHgGrange, convertedCosmicGrange_hg))
queryRes_prot <- queryHits(findOverlaps(prot2geneHgGrange, convertedCosmicGrange_hg))

convertedCosmic$protAnno <- "none"
convertedCosmic$protAnno[subjectRes_cosmic] <- prot2geneHg$feature_name[queryRes_prot]


### quickly creating bed file for the genotyping SNPs - use intersect output bed 
###

#oneBpBedAmpliseqInput <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/20200805mousePanelDraftNoAll_1bp.bed",
#                                    sep = "\t", stringsAsFactors = FALSE, header = FALSE)
#inputSnpLocations <- oneBpBedAmpliseqInput[which(oneBpBedAmpliseqInput$V4 == "SNP_Genotyping"),]
#write.table(inputSnpLocations, "/mnt/DATA6/kevin_recovery/ComprehensiveMouse/20210118SnpsOnMouse.bed", row.names = FALSE, col.names = FALSE,
#            quote = FALSE, sep = "\t")

### did some bed intersect with the final SNPs used and got 490 hits - below is result of that
### doing further work to create a hotspot file specific to torrentsuite - need anchor base which is base before allele of question
###

#genoTypeSnps <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/20210118genoSnpsOnPanel.bed", header = FALSE,
#                           stringsAsFactors = FALSE)
#genoTypeSnps2 <- genoTypeSnps
#genoTypeSnps2$V2 <- genoTypeSnps$V2 - 1
#write.table(genoTypeSnps2, "/mnt/DATA6/kevin_recovery/ComprehensiveMouse/20210118genoGetFasta.bed", col.names = FALSE,
#            row.names = FALSE, quote = FALSE, sep = "\t")

###next I will make the bed file without the header - add manually after
###

#vcf490 <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/20210118genoHotspot.vcf", sep = "\t", header = FALSE,
#                     stringsAsFactors = FALSE)
#bedAnchor <- unlist(read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/20210118genoWithAnchor.bed",
#                                  sep = " ", stringsAsFactors = FALSE))
#bedAmplicon <- read.table("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/20210118bedIntersectWithAms.bed",
#                          sep = "\t", header = FALSE, stringsAsFactors = FALSE)

#bedAnchor_nucs <- bedAnchor[seq(2,980, 2)]
#tmpCol <- str_remove_all(vcf490$V7, pattern = ".*_")
#tmpCol_REF <- paste0("REF=",sapply(str_split(tmpCol, "/"), "[[", 1))
#tmpCol_OBS<- paste0("OBS=",sapply(str_split(tmpCol, "/"), "[[", 2))
#tmpCol_ANC <- paste0("ANCHOR=",toupper(substr(bedAnchor_nucs, 1, 1)))
#tmpCol_full3 <- paste(tmpCol_REF, tmpCol_OBS, tmpCol_ANC, sep = ";")
#noHeaderHotspot <- data.frame(vcf490[1:3], paste0(vcf490$V4, 1:490), tmpCol_full3, bedAmplicon$V8)
#write.table(noHeaderHotspot, "/mnt/DATA6/kevin_recovery/ComprehensiveMouse/20210118hotspot_noHeader.bed",
#            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

### hom/het definition are easy to understand; the unkown definition is maps transcripts mapt o coding regions - no complete ORFs though
###

albertVarsAll <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/MG_test1 variants (Albert).xlsx", sheet = 1)
albertVarsFilt <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/MG_test1 variants (Albert).xlsx", sheet = 2)


mgpVars <- c("hom", "het", "unknown")

allAnnoNew <- read.table("/mnt/DATA5/tmp/kev/misc/20201221newPanelVar.txt",sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tmp <- allAnnoNew[which(allAnnoNew$mm10_mpgpv6_Indels %in% c("het")),]
allAnnoNew <- allAnnoNew[-which(allAnnoNew$mm10_mpgpv6_Indels %in% mgpVars),]

newPanelCalls_names <- as.numeric(str_remove(str_remove(allAnnoNew$Sample, "MG_"), "X.*"))
newPanelCalls_names2 <- newPanelCalls_names



for (i in unique(newPanelCalls_names)) {
  newPanelCalls_names2[which(newPanelCalls_names %in% i)] <- crossTableIdx$old_name2[which(crossTableIdx$mg_id %in% i)]
}
allAnnoNew$Sample <- newPanelCalls_names2

### looking at cross-species annotation
mouseVarGrange <- GRanges(seqnames=Rle(allAnnoNew$Chr),
                         ranges=IRanges(allAnnoNew$Start,
                                       allAnnoNew$End))

#convertedCosmicHits <- subjectHits(findOverlaps(mouseVarGrange, convertedCosmicGrange))
#allAnnoNew_Hits <- queryHits(findOverlaps(mouseVarGrange, convertedCosmicGrange))



### the overlapping hits - might be useful later for looking at hotspot regions - not directly orthologous
### similar to above idea - dan and scott did within 3 or 5 basepairs changes in a mutation ... could be somethingb to look at here
### similar to above part 2: domain specific mutations between orthologous genes

subject_mouseToHuman <- subjectHits(findOverlaps(mouseVarGrange, convertedCosmicGrange))
query_mouseToHuman <- queryHits(findOverlaps(mouseVarGrange, convertedCosmicGrange))

#allAnnoNew_exactMatches <- which(paste0(allAnnoNew$Chr,allAnnoNew$Start, allAnnoNew$End) %in% paste0(convertedCosmic$mm10chr,convertedCosmic$mm10start,convertedCosmic$mm10end))
#convertedCosmic_exactMatches <-  which(paste0(convertedCosmic$mm10chr,convertedCosmic$mm10start,convertedCosmic$mm10end) %in% paste0(allAnnoNew$Chr,allAnnoNew$Start, allAnnoNew$End))

#match_exact <- match(paste0(allAnnoNew$Chr,allAnnoNew$Start, allAnnoNew$End), paste0(convertedCosmic$mm10chr,convertedCosmic$mm10start,convertedCosmic$mm10end))
#match_exact <- match_exact[-which(is.na(match_exact))]

convertedDf <- data.frame("hgGene" = rep("none", nrow(allAnnoNew)), "ENST_ID_VER" = rep("none", nrow(allAnnoNew)),
                          "hgCdsMutation" = rep("none", nrow(allAnnoNew)),
                          "hgAaMutation" = rep("none", nrow(allAnnoNew)),
                          "CosmicID" = rep("none", nrow(allAnnoNew)), 
                          "EXAC_AF" = rep("none", nrow(allAnnoNew)),
                          "GNOMAD_AF" = rep("none", nrow(allAnnoNew)),
                          "Gr38_position" = rep("none", nrow(allAnnoNew)),
                          "mm10chr" = rep("none", nrow(allAnnoNew)),
                          "mm10_start" = rep(0, nrow(allAnnoNew)),
                          "mm10end" = rep(0, nrow(allAnnoNew)),
                          "domainInfo" = rep("none", nrow(allAnnoNew)),stringsAsFactors = FALSE)


#convertedDf[allAnnoNew_exactMatches,] <- convertedCosmic[match_exact,]
convertedDf[query_mouseToHuman,] <- convertedCosmic[subject_mouseToHuman,]
allAnnoNew <- cbind(allAnnoNew, convertedDf)
allAnnoNew <- allAnnoNew[-which(allAnnoNew$Sample == "2611n"),]

bedErrors <- paste0(allAnnoNew$Chr, allAnnoNew$Start, allAnnoNew$End)
bedErrors_freq <- table(bedErrors)/length(unique(allAnnoNew$Sample))
badPositions <- names(bedErrors_freq)[which(bedErrors_freq > 0.05)]
allAnnoNew <- allAnnoNew[-which(bedErrors %in% badPositions),]


### take out indels to separate them into unqiue entries for filtering
allAnnoNew_indels <- allAnnoNew[grep(allAnnoNew$FAO, pattern = ","),]
allAnnoNew <- allAnnoNew[-grep(allAnnoNew$FAO, pattern = ","),]


newIndels <- NULL
colsOfInterest <- which(colnames(allAnnoNew_indels) %in% c("AF", "FAO", "FSAF", "FSAR"))
for (i in 1:nrow(allAnnoNew_indels)) {
  tmpAf <- unlist(str_split(allAnnoNew_indels$AF[i], ","))
  tmpFao <- unlist(str_split(allAnnoNew_indels$FAO[i], ","))
  tmpFsaf <- unlist(str_split(allAnnoNew_indels$FSAF[i], ","))
  tmpFsar <- unlist(str_split(allAnnoNew_indels$FSAR[i], ","))
  oneEntry <- allAnnoNew_indels[i, -colsOfInterest]
  newEntry <- oneEntry[rep(1, each = length(tmpAf)), ]
  newEntry2 <- data.frame(newEntry, "AF" = tmpAf, "FAO" = tmpFao,
                          "FSAF" = tmpFsaf, "FSAR" = tmpFsar, stringsAsFactors = FALSE)
  newIndels <- rbind(newIndels, newEntry2)
}


newIndels$FDP <- as.numeric(newIndels$FDP)
newIndels$FAO <- as.numeric(newIndels$FAO)
newIndels$GQ <- as.numeric(newIndels$GQ)
newIndels$AF <- as.numeric(newIndels$AF)
newIndels$FSAF <- as.numeric(newIndels$FSAF)
newIndels$FSAR <- as.numeric(newIndels$FSAR)

indel_fdpFilt <- which(newIndels$FDP > 20)
indel_faoFilt <- which(newIndels$FAO > 5)
indel_freqFilt <- which(newIndels$AF > 0.05)
indel_strandRatio <- intersect(which(newIndels$FSAF/newIndels$FSAR > 0.2),
                         which(newIndels$FSAF/newIndels$FSAR < 5))
indel_goodSamps <- Reduce(intersect, list(indel_fdpFilt, indel_faoFilt, indel_freqFilt, indel_strandRatio))
newIndels_goodsamps <- newIndels[indel_goodSamps,]
newIndels_goodsamps2 <- newIndels_goodsamps[which(newIndels_goodsamps$Func.refGene == "exonic"),]


### last filter will get rid of substituion errors that are probably wrong  - debating whether to get rid of all substituions or just long ones
indelSub <- newIndels_goodsamps2[grep("substitution", newIndels_goodsamps2$ExonicFunc.refGene),]
deletions_goodsamps <- newIndels_goodsamps2[-grep("substitution", newIndels_goodsamps2$ExonicFunc.refGene),]


allAnnoNew$FDP <- as.numeric(allAnnoNew$FDP)
allAnnoNew$FAO <- as.numeric(allAnnoNew$FAO)
allAnnoNew$GQ <- as.numeric(allAnnoNew$GQ)
allAnnoNew$AF <- as.numeric(allAnnoNew$AF)
allAnnoNew$FSAF <- as.numeric(allAnnoNew$FSAF)
allAnnoNew$FSAR <- as.numeric(allAnnoNew$FSAR)

fdpFilt <- which(allAnnoNew$FDP > 20)
faoFilt <- which(allAnnoNew$FAO > 5)
freqFilt <- which(allAnnoNew$AF > 0.05)
strandRatio <- intersect(which(allAnnoNew$FSAF/allAnnoNew$FSAR > 0.2),
                 which(allAnnoNew$FSAF/allAnnoNew$FSAR < 5))
goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio))
allAnnoNew_goodsamps <- allAnnoNew[goodSamps,]

allAnnoNew_goodsampsExonic <- allAnnoNew_goodsamps[-which(allAnnoNew_goodsamps$ExonicFunc.refGene == ""),]
additionalSubs <- allAnnoNew_goodsampsExonic[grep("substitution", allAnnoNew_goodsampsExonic$ExonicFunc.refGene),]
allAnnoNew_goodsampsExonic2 <- allAnnoNew_goodsampsExonic[-grep("substitution", allAnnoNew_goodsampsExonic$ExonicFunc.refGene),]


finalTable_noSubs <- rbind(allAnnoNew_goodsampsExonic2, deletions_goodsamps)
finalSubs <- rbind(additionalSubs, indelSub)

write.table(finalTable_noSubs, "/mnt/DATA5/tmp/kev/misc/20210215_mgpVarList_noSub.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(finalSubs, "/mnt/DATA5/tmp/kev/misc/20210215_mgpVarList_Subs.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


### workflow from old samps
###

allAnnoOld <- read.table("/mnt/DATA5/tmp/kev/misc/20201221oldPanelVar.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
allAnnoOld <- allAnnoOld[-which(allAnnoOld$Sample == "Sample"),]
allAnnoOld <- allAnnoOld[-which(allAnnoOld$mm10_mpgpv6_Indels %in% mgpVars),]
allAnnoOld$Sample <- str_remove(allAnnoOld$Sample, "X")
allAnnoOld$Sample <- tolower(allAnnoOld$Sample)
allAnnoOld$Sample <- str_replace_all(allAnnoOld$Sample, "[[:punct:]]", "")
allAnnoOld_indels <- allAnnoOld[grep(allAnnoOld$FAO, pattern = ","),]
allAnnoOld <- allAnnoOld[-grep(allAnnoOld$FAO, pattern = ","),]


allAnnoOld$string <- paste0(allAnnoOld$Gene.refGene, allAnnoOld$Start, allAnnoOld$End)
bedBadString2 <- names(table(allAnnoOld$string)[which(table(allAnnoOld$string) > round(length(unique(allAnnoOld$Sample)) * 0.01))])
bedErrorIdx2 <- which(allAnnoOld$string %in% bedBadString2)
allAnnoOld <- allAnnoOld[-bedErrorIdx2,]

allAnnoOld$FDP <- as.numeric(allAnnoOld$FDP)
allAnnoOld$FAO <- as.numeric(allAnnoOld$FAO)
allAnnoOld$GQ <- as.numeric(allAnnoOld$GQ)
allAnnoOld$AF <- as.numeric(allAnnoOld$AF)
allAnnoOld$FSAF <- as.numeric(allAnnoOld$FSAF)
allAnnoOld$FSAR <- as.numeric(allAnnoOld$FSAR)

fdpFilt_old <- which(allAnnoOld$FDP > 20)
faoFilt_old <- which(allAnnoOld$FAO > 5)
freqFilt_old <- which(allAnnoOld$AF > 0.05)
strandRatio_old <- intersect(which(allAnnoOld$FSAF/allAnnoOld$FSAR > 0.2),
                         which(allAnnoOld$FSAF/allAnnoOld$FSAR < 5))
goodSamps_old <- Reduce(intersect, list(fdpFilt_old, faoFilt_old, freqFilt_old, strandRatio_old))
allAnnoOld_goodsamps <- allAnnoOld[goodSamps_old,]


allAnnoNew_goodsamps_subset <- allAnnoNew_goodsamps[which(allAnnoNew_goodsamps$Sample %in% allAnnoOld_goodsamps$Sample),]
allAnnoOld_goodsamps_subset <- allAnnoOld_goodsamps[which(allAnnoOld_goodsamps$Sample %in% allAnnoNew_goodsamps_subset$Sample),]

normSampList <- c("2796n", "3867n", "8234n", "13604n", "14104t", "14154n", "14433n", "2405n", "2519n", "2611n")

allAnnoNew_goodsamps_subset_n <- allAnnoNew_goodsamps_subset[which(allAnnoNew_goodsamps_subset$Sample %in% normSampList),]
allAnnoNew_goodsamps_subset_t  <- allAnnoNew_goodsamps_subset[-which(allAnnoNew_goodsamps_subset$Sample %in% normSampList),]

allAnnoOld_goodsamps_subset_n <- allAnnoOld_goodsamps_subset[which(allAnnoOld_goodsamps_subset$Sample %in% normSampList),]
allAnnoOld_goodsamps_subset_t <- allAnnoOld_goodsamps_subset[-which(allAnnoOld_goodsamps_subset$Sample %in% normSampList),]

#allAnnoOld_goodsamps_subset_t2 <- allAnnoOld_goodsamps_subset_t[which(allAnnoOld_goodsamps_subset_t$Gene.refGene %in% allAnnoNew_goodsamps_subset_t$Gene.refGene),]
#allAnnoNew_goodsamps_subset_t2 <- allAnnoNew_goodsamps_subset_t[which(allAnnoNew_goodsamps_subset_t$Gene.refGene %in% allAnnoOld_goodsamps_subset_t2$Gene.refGene),]


### testing differences in affected pathways


outputOld <- run_pathfindR()

tmpBed <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.bed", sep = "\t", header = FALSE, skip = 1, stringsAsFactors = FALSE)
tmpBed$V1 <- stringr::str_remove(tmpBed$V1, "chr")
write.table(tmpBed, "/home/kevhu/data/bedFiles/IAD202670_167_Designed_noCHR.bed", sep = "\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE)
