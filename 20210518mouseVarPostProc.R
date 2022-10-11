source("/mnt/DATA6/kevin_recovery/scripts/fastReadFile.R")

library(GenomicRanges)
library(rtracklayer)

cosmic <- faster.readfile("/mnt/DATA5/tmp/kev/misc/cmc_export.v92.edited.tsv", sep = "\t", 0)
cosmic_red <- cosmic[which(cosmic$COSMIC_SAMPLE_MUTATED > 30),]
mm10Tohg38Chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/mm10ToHg38.over.chain")
mgpVars <- c("hom", "het", "unknown")

allAnnoNew <- read.table("/home/kevhu/scripts/newMousePanelPipeline/reportAnno/Auto_user_AUS5-76-MG_test1_255_185_anno.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)

allAnnoNew2 <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-76-MG_test1_255_185_anno.txt",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#hetSites <- allAnnoNew[which(allAnnoNew$mm10_mpgpv6_Indels %in% c("het")),]
allAnnoNew <- allAnnoNew[-which(allAnnoNew$mm10_mpgpv6_Indels %in% mgpVars),]

fdpFilt <- which(allAnnoNew$FDP > 20)
faoFilt <- which(allAnnoNew$FAO > 5)
freqFilt <- which(allAnnoNew$AF > 0.05)
hrunFilt <- which(allAnnoNew$HRUN > 4)
strandRatio <- intersect(which(allAnnoNew$FSAF/allAnnoNew$FSAR > 0.2),
                         which(allAnnoNew$FSAF/allAnnoNew$FSAR < 5))
goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio))
allAnnoNew_goodsamps <- allAnnoNew[goodSamps,]
allAnnoNew_goodsamps_exon <- allAnnoNew_goodsamps[which(allAnnoNew_goodsamps$Func.refGene == "exonic"), ]

# this looking at a simple check to see if there is also one in humans
# I think I can do the reverse lookup for conserved nonsynomous mutations - this makes more sense to me

mouseVarRange <- GRanges(seqnames = allAnnoNew_goodsamps_exon$Chr[2], 
                         IRanges(start = allAnnoNew_goodsamps_exon$Start[2],
                                 end = allAnnoNew_goodsamps_exon$End[2]))

loOut <- liftOver(mouseVarRange, mm10Tohg38Chain)

# make the converted vars into the mutation coord format
loOutDf <- data.frame(mouseVarRange, stringsAsFactors = FALSE)
loOutDf$seqnames <- as.character(loOutDf$seqnames)
loOutDf$matchString <- paste0(str_remove(loOutDf$seqnames, "chr"), ":", loOutDf$start, "-", loOutDf$end)

which(loOutDf$matchString %in% cosmic$`Mutation genome position GRCh38`)
