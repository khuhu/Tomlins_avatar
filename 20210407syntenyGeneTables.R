# quick and dirty way to get df of genes for synteny mapping
library(biomaRt)

cbioportal <- read.table("/mnt/DATA5/tmp/kev/misc/20210404cbioBiorMart_genes_tid.txt",
                         sep = "\t", stringsAsFactors = FALSE, header = TRUE)
cbioportal2 <- cbioportal[which(cbioportal$Gene %in% geneDf$Gene),]

ensemblHg <- useMart("ensembl", dataset= "hsapiens_gene_ensembl")
ensemblMm <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") 

filtsHg <- listFilters(ensemblHg)
filtsMm <- listFilters(ensemblMm)
attrHg <- listAttributes(ensemblHg)
attrMm <- listAttributes(ensemblMm)

attrHg_hom <- c("external_gene_name", "ensembl_transcript_id",
                "mmusculus_homolog_ensembl_gene",
                "mmusculus_homolog_ensembl_peptide")

orthoRes <- getBM(mart = ensemblHg, attributes = attrHg_hom, filters = "ensembl_transcript_id",
                  values = cbioportal$Tid)

orthoRes2 <- orthoRes[-which(orthoRes$mmusculus_homolog_ensembl_gene == ""),]
orthoRes2 <- orthoRes2[-which(duplicated(orthoRes2$external_gene_name)),]


attrMm_gene <- c("external_gene_name", "chromosome_name",
                "start_position", "end_position")

attrHg_gene <- c("external_gene_name", "chromosome_name",
                 "start_position", "end_position")



geneResMm <- getBM(mart = ensemblMm, attributes = attrMm_gene,
                   filters = "ensembl_gene_id", values = orthoRes$mmusculus_homolog_ensembl_gene)
which(duplicated(geneResMm$external_gene_name))
geneResMm2 <- geneResMm[-which(duplicated(geneResMm$external_gene_name)),]
  
geneResHg <- getBM(mart = ensemblHg, attributes = attrHg_gene,
                   filters = "ensembl_transcript_id", values = orthoRes2$ensembl_transcript_id)
which(duplicated(geneResHg$external_gene_name))

write.table(geneResMm2, "/mnt/DATA6/kevin_recovery/apps/circos/20210407mouseGeneTable.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(geneResHg, "/mnt/DATA6/kevin_recovery/apps/circos/20210407humanGeneTable.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

