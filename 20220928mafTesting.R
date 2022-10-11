library(vcfR)
library(maftools)
#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')

laml = read.maf(maf = laml.maf,
                verbose = FALSE)


#By default the function plots top20 mutated genes
oncoplot(maf = laml, draw_titv = TRUE)


### can only run individual maf files
### so combining them first

setwd("/mnt/DATA5/tmp/kev/testMm10Vep/")
filenames <- system("ls *.maf", intern = TRUE)

i <- filenames[1]
allMaf <- NULL
for (i in filenames) {
  tmpFile <- read.table(i, skip = 1, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  allMaf <- rbind(allMaf, tmpFile)
}

cat("#version 2.4 \n",file="combinedMaf.maf")
write.table(allMaf, "./combinedMaf.maf", sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

genes <- c("Trp53", "Kmt2d", "Nf1", "Kmt2c", "Apc", "Crebbp", "Mtor", "Notch2", "Ptch1", "Setd2",
           "Atm", "Brca1", "Brca2", "Cdh1", "Cdk12", "Erbb3", "Pik3ca")

testFile <- "/mnt/DATA5/tmp/kev/testMm10Vep/combinedMaf.maf"
testMaf <- read.maf(testFile, verbose = TRUE)
oncoplot(maf = testMaf, genes = genes)

### bbn maf - downloaded set of files from firehose


setwd("/mnt/DATA5/tmp/kev/tmpDbs/genomicDataCommons/gdac.broadinstitute.org_BLCA.Mutation_Packager_Calls.Level_3.2016012800.0.0/")
filenames <- system("ls *maf*", intern = TRUE)
filenames <- filenames[-which(filenames %in% c("TCGA-BT-A0YX-01.maf.txt", "TCGA-BT-A2LD-01.maf.txt",
                                               "TCGA-CU-A0YN-01.maf.txt", "TCGA-G2-A3IB-01.maf.txt",
                                               "TCGA-GC-A3RB-01.maf.txt", "combinedMaf.maf"))]
humanMaf <- NULL
for (i in filenames) {
  print(i)
  tmpFile <- read.table(i, sep = "\t", header = TRUE, 
                        stringsAsFactors = FALSE, fill = TRUE)
  tmpFile <- tmpFile[, 1:50]
  humanMaf <- rbind(humanMaf, tmpFile)
}

write.table(humanMaf, "./combinedMaf.maf", sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

hGenes <- toupper(c("Tp53", "Mll3", "Nf1", "Mll2", "Apc", "Crebbp", "Mtor", "Notch2", "Ptch1", "Setd2",
                    "Atm", "Brca1", "Brca2", "Cdh1", "Cdk12", "Erbb3", "Pik3ca"))

humanFile <- "./combinedMaf.maf"
humanMaf2 <- read.maf(humanFile, verbose = TRUE)
oncoplot(maf = humanMaf2, genes = hGenes)

### aaron bbn samples

### can only run individual maf files
### so combining them first

setwd("/mnt/DATA5/tmp/kev/testMm10Vep/aaronFilt/")
filenames <- system("ls *.maf", intern = TRUE)

i <- filenames[1]
allMaf <- NULL
for (i in filenames) {
  tmpFile <- read.table(i, skip = 1, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  allMaf <- rbind(allMaf, tmpFile)
}

tmpVarList <- read.table("../20221004aaronVarList.txt", sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE)

tmpVarList$string2 <- paste(tmpVarList$Sample, tmpVarList$Chr, tmpVarList$Start, tmpVarList$End)
allMaf2 <- allMaf
allMaf2$samplename <- str_remove(allMaf2$Tumor_Sample_Barcode, "\\_.*")
allMaf2$string <- paste(allMaf2$samplename, paste0("chr", allMaf2$Chromosome), allMaf2$Start_Position, allMaf2$End_Position)

allMaf  <- allMaf[which(allMaf2$string %in% tmpVarList$string2), ]
tmpVarList$string2[-which(tmpVarList$string2 %in% allMaf2$string)]


cat("#version 2.4 \n",file="combinedMaf.maf")
write.table(allMaf, "./combinedMaf.maf", sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

genes <- c("Trp53", "Kmt2d", "Nf1", "Kmt2c", "Crebbp", "Notch1","Setd2", "Apc",
           "Atm", "Mtor", "Brca1", "Notch2", "Ptch1", "Arid1a",
           "Brca2", "Cdh1", "Cdk12", "Erbb3", "Pik3ca")

testFile <- "/mnt/DATA5/tmp/kev/testMm10Vep/aaronFilt/combinedMaf.maf"
testMaf <- read.maf(testFile, verbose = TRUE)
oncoplot(maf = testMaf, barcode_mar = 10, 
         showTumorSampleBarcodes = TRUE)

oncoplot(maf = testMaf, barcode_mar = 10, genes = genes,  
         showTumorSampleBarcodes = TRUE)


tmpVarList$string2[-which(tmpVarList$string2 %in% allMaf2$string)]
tmpVarList$Gene.refGene[-which(tmpVarList$string2 %in% allMaf2$string)]


### looking at variants from two bbn runs
###
###

bbn239 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-239-BBN_mouse_bladder_MG_493_562/combinedCalls.txt",
                     sep = "\t", stringsAsFactors = FALSE, header = TRUE)

bbn260 <- read.table("/mnt/DATA6/mouseData/copynumber/Auto_user_AUS5-260-BBN_mouse_bladder_MG_2_514_605/combinedCalls.txt",
                     sep = "\t", stringsAsFactors = FALSE, header = TRUE)


bbn239$logCn <- log2(bbn239$CopyNumberRatio)
bbn239$logCn <- bbn239$logCn / 0.5


bbn260$logCn <- log2(bbn260$CopyNumberRatio)
bbn260$logCn <- bbn260$logCn / 0.5


bbn239_1copy <- bbn239[which(abs(2^bbn239$logCn *2 - 2) > 0.999),]
bbn260_1copy <- bbn260[which(abs(2^bbn260$logCn *2 - 2) > 0.999),]
bbn239_1copy$sample2 <- str_remove(bbn239_1copy$Sample, "\\_D.*")
bbn260_1copy$sample2 <- str_remove(bbn260_1copy$Sample, "\\_D.*")

bbn239_1copy <- bbn239_1copy[which(bbn239_1copy$Log10QValue < -1.30103),]
bbn260_1copy <- bbn260_1copy[which(bbn260_1copy$Log10QValue < -1.30103),]

bbn239_1copy$string <- paste(bbn239_1copy$sample2, bbn239_1copy$Gene)
bbn260_1copy$string <- paste(bbn260_1copy$sample2, bbn260_1copy$Gene)

bbn260_1copy$string[which(bbn260_1copy$string %in% bbn260_1copy$string)]


concordantCnChanges <- bbn260_1copy[which(bbn260_1copy$string %in% bbn260_1copy$string),]

mouseBed <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed",
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE)

xGenes <- unique(mouseBed$V8[which(mouseBed$V1 == "chrX")])

concordantCnChanges <- concordantCnChanges[-which(concordantCnChanges$Gene %in% xGenes),]



write.table(concordantCnChanges, "/mnt/DATA5/tmp/kev/misc/20221005bbnConcordantCns.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

