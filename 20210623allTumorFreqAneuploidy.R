# I think I should formally redo creating the freq matrices here too 
# with the new annotations from cbioportal, along with with the aneuploidy comparisons
# should be boxplots comparing aneuploidy events between genotype
# I can separate human non-cell data into specific genotype if I want


library(readxl)

# remake freq table and for all subsequent samples
cbio_anno_coad <- read.table("/mnt/DATA5/tmp/kev/misc/PATIENT_DATA_oncoprint_KAP.tsv",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)

cbio_anno_coad2 <- cbio_anno_coad[4:9, 3:ncol(cbio_anno_coad)]
rownames(cbio_anno_coad2) <- c("APC_cn", "TP53_cn", "KRAS_cn", 
                               "APC_mut", "TP53_mut", "KRAS_mut")
cbio_anno_coad2 <- t(cbio_anno_coad2)
cbio_anno_coad2 <- cbind("samples" = rownames(cbio_anno_coad2), cbio_anno_coad2)
rownames(cbio_anno_coad2) <- NULL
cbio_anno_coad2 <- data.frame(cbio_anno_coad2, stringsAsFactors = FALSE)

cbio_anno_bprn <- read.table("/mnt/DATA5/tmp/kev/misc/PATIENT_DATA_oncoprint_BPRN.tsv",
                             sep = "\t", stringsAsFactors = FALSE, header = TRUE)

cbio_anno_bprn2 <- cbio_anno_bprn[4:11, 3:ncol(cbio_anno_bprn)]
rownames(cbio_anno_bprn2) <- c("Brca1_cn", "Tp53_cn", "Nf1_cn", "Rb1_cn",
                               "Brca1_mut", "Tp53_mut", "Nf1_mut", "Rb1_mut")
cbio_anno_bprn2 <- t(cbio_anno_bprn2)
cbio_anno_bprn2 <- cbind("samples" = rownames(cbio_anno_bprn2), cbio_anno_bprn2)
rownames(cbio_anno_bprn2) <- NULL
cbio_anno_bprn2 <- data.frame(cbio_anno_bprn2, stringsAsFactors = FALSE) 

coad_KAP <- which((cbio_anno_coad2$APC_cn == "homdel_rec" | grepl("driver" , cbio_anno_coad2$APC_mut)) &
        (cbio_anno_coad2$TP53_cn == "homdel_rec" | grepl("driver" , cbio_anno_coad2$TP53_mut)) &
          grepl("driver", cbio_anno_coad2$KRAS_mut))
cbio_anno_coad3 <- cbio_anno_coad2[coad_KAP,c(1:3,5:7)]

write.table(cbio_anno_coad3, "/mnt/DATA5/tmp/kev/misc/20210628_panCancerCoadKAPIds.txt", sep = "\t", col.names = TRUE,
            quote = FALSE, row.names = FALSE)
# so these four changes never co-occur in the tcga dataset ....
ov_BPRN <- which((cbio_anno_bprn2$Brca1_cn == "homdel_rec" | grepl("driver", cbio_anno_bprn2$Brca1_mut)) &
                     (cbio_anno_bprn2$Tp53_cn == "homdel_rec" | grepl("driver", cbio_anno_bprn2$Tp53_mut)) &
                     (cbio_anno_bprn2$Nf1_cn == "homdel_rec" | grepl("driver", cbio_anno_bprn2$Nf1_mut)) &
                     (cbio_anno_bprn2$Rb1_cn == "homdel_rec" | grepl("driver", cbio_anno_bprn2$Rb1_mut)))



# will be interesting to see the phenotypes of the cho lab samples that are more similar
# in pathogenesis (compared to human)


# mainly because we have two types of models with coad and ov - coad which alters genes which
# mimics the pathongenesis in human, whereas BPRN never co-occurs in hgscs. hgsc model 
# does show similarity to human tumor i.e McCool paper. 

# maybe if I look at frequency data between all and certain genotypes ... and there is no
# difference ..... shows how it's a secondary feature of tumor development




# CCLE aneuploidy numbers - make comparisons to human TCGA data
# tangent, but along with this I should also have something for the
# drawing protein done by then too

ccle_cn <- read_xlsx("/mnt/DATA5/tmp/kev/tmpDbs/CCLE/CCLE_ABSOLUTE_combined_20181227.xlsx")
ccle_anno <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/CCLE/Cell_lines_annotations_20181226.txt",
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
ccle_ov <- ccle_cn[grep("OVARY", ccle_cn$sample),]
ccle_coad <- ccle_cn[grep("INTESTINE", ccle_cn$sample),]
