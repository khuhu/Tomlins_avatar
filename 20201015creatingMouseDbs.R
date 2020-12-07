### new script to query biomart for gene names, cdna coordinates etc

library(biomaRt)

#listOfEnsIdsHg <- read.table("/mnt/DATA5/tmp/kev/programs/ensemblApi/src/ensembl/misc-scripts/human_canonical_transcripts_noid_v99.fa",
#                           header = FALSE, sep = "\t")
#listOfEnsIdsMm <- read.table("/mnt/DATA5/tmp/kev/programs/ensemblApi/src/ensembl/misc-scripts/mouse_canonical_transcripts_noid_v99.fa",
#                             header = FALSE, sep = "\t")

knownCanHg38 <- read.table("/home/kevhu/data/20201015hg38UcscXrefKnownCan.txt", stringsAsFactors = FALSE,
                           sep = "\t", header = TRUE)
knownCanMm10 <- read.table("/home/kevhu/data/20201020mm10UcscXrefKnownCan.txt", stringsAsFactors = FALSE,
                           sep = "\t", header = TRUE)


ucsc_txn_ids <- sub(x = knownCanHg38$hg38.knownCanonical.transcript, pattern = "\\..*", replacement = "")
ucsc_txn_mm10_ids <- sub(x = knownCanMm10$mm10.knownCanonical.transcript, pattern = "\\..*",
                         replacement = "")

ensemblHg <- useMart("ensembl", dataset= "hsapiens_gene_ensembl")
ensemblMm <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") 

filtsHg <- listFilters(ensemblHg)
filtsMm <- listFilters(ensemblMm)
attrHg <- listAttributes(ensemblHg)
attrMm <- listAttributes(ensemblMm)

#attrHg_prot <- c("external_gene_name", "ensembl_transcript_id",
#                 "cdna_coding_start", "cdna_coding_end", "chromosome_name",
#                 "start_position", "end_position", "strand", "rank")


#attrHg_prot <- c("external_gene_name", "ensembl_transcript_id",
#                 "cdna_coding_start", "cdna_coding_end", "chromosome_name",
#                 "start_position", "end_position", "strand", "rank")

attrHg_prot <- c("external_gene_name", "ensembl_transcript_id",
                 "cdna_coding_start", "cdna_coding_end", "chromosome_name",
                 "cds_start", "cds_end", "strand", "rank",
                 "exon_chrom_start", "exon_chrom_end")


attrHg_geneName <- c("external_gene_name", "ensembl_transcript_id",
                 "mmusculus_homolog_associated_gene_name")

attr_peptide <- c("external_gene_name", "ensembl_transcript_id",
                     "peptide")



### test with a TP53 gene txn id before running batch
chromFilters <- c(as.character(1:22), "X", "Y")

#cdnaLevel <- getBM(mart = ensemblHg, attributes = attrHg_prot, filters = "ensembl_transcript_id",
#      values = listOfEnsIdsHg$V1)

#cdnaLevelHgFilt <- cdnaLevel[-which(is.na(cdnaLevel$cdna_coding_start)),]
#cdnaLevelHgFilt <- cdnaLevelHgFilt[which(cdnaLevelHgFilt$chromosome_name %in% chromFilters),]


#write.table(cdnaLevelHgFilt, "/home/kevhu/data/20201016hg38biomartQuery.txt", sep = "\t",
#            quote = FALSE, col.names = TRUE, row.names = FALSE)

cdnaLevelMm <- getBM(mart = ensemblMm, attributes = attrHg_prot, filters = "external_gene_name",
                   values = knownCanMm10$mm10.kgXref.geneSymbol)

cdnaLevelMmFilt <- cdnaLevelMm[-which(is.na(cdnaLevelMm$cdna_coding_start)),]
cdnaLevelMmFilt <- cdnaLevelMmFilt[which(cdnaLevelMmFilt$chromosome_name %in% chromFilters),]


write.table(cdnaLevelMmFilt, "/home/kevhu/data/20201030Mm10KnownCanbiomartQuery.txt", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)


cdnaLevelHgUcsc <- getBM(mart = ensemblHg, attributes = attrHg_prot, filters = "external_gene_name",
                         values = knownCanHg38$hg38.kgXref.geneSymbol)

cdnaLevelHgUcscFilt <- cdnaLevelHgUcsc[-which(is.na(cdnaLevelHgUcsc$cdna_coding_start)),]
cdnaLevelHgUcscFilt <- cdnaLevelHgUcscFilt[which(cdnaLevelHgUcscFilt$chromosome_name %in% chromFilters),]

write.table(cdnaLevelHgUcscFilt, "/home/kevhu/data/20201030hg38KnownCanbiomartQuery.txt", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)





### need to get cdna to check if 

cdnaSeqHgUcsc <- getBM(mart = ensemblHg, attributes = attrHg_cdna, filters = "ensembl_transcript_id",
                         values = ucsc_txn_ids)

cdnaSeqMmUcsc <- getBM(mart = ensemblMm, attributes = attrHg_cdna, filters = "ensembl_transcript_id",
                       values = ucsc_txn_mm10_ids)


write.table(cdnaSeqHgUcsc, "/home/kevhu/data/20201021HgUcscCdna.txt",
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)

write.table(cdnaSeqMmUcsc, "/home/kevhu/data/20201021MmUcscCdna.txt",
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)


geneName <- getBM(ensemblHg, attributes = attrHg_geneName, filters = "ensembl_transcript_id",
                  values = ucsc_txn_ids)

write.table(geneName, "/home/kevhu/data/20201021geneNameBiomart.txt",
             quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)


proteinHgSeq <- getBM(ensemblHg, attributes = attr_peptide, filters = "external_gene_name",
                    values = knownCanHg38$hg38.kgXref.geneSymbol)

write.table(proteinHgSeq, "/home/kevhu/data/20201030proteinHg.txt",
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE) 


proteinMmSeq <- getBM(ensemblMm, attributes = attr_peptide, filters = "external_gene_name",
                      values = knownCanMm10$mm10.kgXref.geneSymbol)


write.table(proteinMmSeq, "/home/kevhu/data/20201030proteinMm.txt",
            quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE) 



