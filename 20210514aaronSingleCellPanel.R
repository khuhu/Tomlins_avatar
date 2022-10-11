source("/home/kevhu/scripts/20210512functionListMousePanel.R")

listOfChr <- c(as.character(1:22), "X", "Y")

geneCdsHs <- read.table("/home/kevhu/data/20210604hg19KnownCanbiomartQuery.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE)
geneCdsHs <-  geneCdsHs[which(geneCdsHs$chromosome_name %in% listOfChr),]

inputTable <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/aaronSingleCell.xlsx",
                                sheet = 1)

allTblsList <- parseInputTable(inputTable)

tmp <- getAmpliconTargets(geneList = allTblsList[[3]], ampSize = 225, numAmps = 6)

# the NA's are from empty 1amp processing
tmp <- tmp[-which(is.na(tmp$chromosome)),]
table(tmp$gene)

write.table(tmp, "/mnt/DATA5/tmp/kev/misc/20210604aaronSingleCell6amps.bed", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


# getting cbioportal based CDS's
tilingGeneCbio <- cBioportalIds$Tid[which(cBioportalIds$Gene %in% allTblsList[[4]])]
tilingGeneDat <- geneCdsHs[which(geneCdsHs$ensembl_transcript_id %in% tilingGeneCbio),]
tmp2 <- cdsConversion(tilingGeneDat)
tilingGeneDat2 <- tmp2[,c(1:11,14,15)]
colnames(tilingGeneDat2)[12:13] <- c("cds_genomic_start", "cds_genomic_end")


write.table(tilingGeneDat2, "/mnt/DATA5/tmp/kev/misc/20210604aaronSingleCellCds.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")



# new cds check, YHWAZ and TERT are neg. strand, WWOX is pos. strand
tmpYWHAZ <- geneCdsHs[which(geneCdsHs$ensembl_transcript_id == "ENST00000353245"),]
tmpYWHAZ_cds <- cdsConversion(tmpYWHAZ)
tmpYWHAZ_cds$tmp <- abs(tmpYWHAZ_cds$strand_cds_start_conv - tmpYWHAZ_cds$strand_cds_end_conv)
tmpYWHAZ_cds$cds_end - tmpYWHAZ_cds$cds_start

tmpWWOX <- geneCdsHs[which(geneCdsHs$ensembl_transcript_id == "ENST00000566780"),]
tmpWWOX_cds <- cdsConversion(tmpWWOX)
tmpWWOX_cds$tmp <- abs(tmpWWOX_cds$strand_cds_start_conv - tmpWWOX_cds$strand_cds_end_conv)
tmpWWOX_cds$cds_end - tmpWWOX_cds$cds_start

tmpTERT <- geneCdsHs[which(geneCdsHs$ensembl_transcript_id == "ENST00000310581"),]
tmpTERT_cds <- cdsConversion(tmpTERT)
tmpTERT_cds$tmp <- abs(tmpTERT_cds$strand_cds_start_conv - tmpTERT_cds$strand_cds_end_conv)
tmpTERT_cds$cds_end - tmpTERT_cds$cds_start


