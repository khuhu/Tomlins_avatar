### tvc was ran using tvc docker file
### files transfered from brahm to avatar for a quick comparison

nameStripper <- function(df){
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
}

library(reshape2)
library(stringr)

### don't know why but the fam file is off with vcfs from admixtureFilesV2
refDir <- "/mnt/DATA6/mouseData/admixtureFiles/"
# refDir <- "/mnt/DATA6/mouseData/admixtureFilesV2/"
setwd(refDir)
# refFiltFiles <- paste0(refDir, "", system("ls *.vcf.gz", intern = TRUE))
# refFiltFiles <- paste0(refDir, "", system("ls *485snp.withchr.vcf.gz", intern = TRUE))
refFiltFiles <- paste0(refDir, "", system("ls *.485snpsCallOut.mm10.vcf.gz", intern = TRUE))
# refFiltFiles <- refFiltFiles[grep("filt", refFiltFiles)]
# refFiltFiles  <- refFiltFiles[-which(refFiltFiles == "/mnt/DATA6/mouseData/admixtureFiles/KC-05_MG_X64.filt.vcf.gz")]


setwd("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/")
system("cp /mnt/DATA6/mouseData/admixtureFiles/admixtureSnpNameList.txt ./")
listOfFies <- system("ls *geno.vcf", intern = TRUE)

filtFiles <- NULL
for (i in listOfFies) {
  strippedName <- str_remove(i, ".geno.vcf")
  cmd <- paste("bcftools view --exclude-types indels,mnps,other -i'ID=@admixtureSnpNameList.txt'", i, "-Oz -o", paste0(strippedName, ".filt.vcf.gz"))
  system(cmd)
  filtFiles <- c(filtFiles, paste0(strippedName, ".filt.vcf.gz"))
}


cmd2 <- "for file in *.filt.vcf.gz ; do   tabix -f -p vcf ${file}  ; done"
system(cmd2)


filtFiles2 <- c(filtFiles, refFiltFiles)
writeLines(filtFiles2, "listOfFiltFiles.txt")

cmd6 <- "bcftools merge --force-samples --file-list listOfFiltFiles.txt -Oz -o admixturetest.vcf.gz"
system(cmd6)

system("tabix -f -p vcf admixturetest.vcf.gz")
# for reruns
setwd("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/")
# system("plink --vcf admixturetest.vcf.gz --make-bed --out admixtureBed --const-fid 0 --allow-no-sex --no-parents --geno 0.2 --mind 0.9")
system("plink --vcf admixturetest.vcf.gz --make-bed --out admixtureBed --const-fid 0 --allow-no-sex --no-parents --geno 0.5 --mind 0.9")

tmpTable <- read.table("admixtureBed.fam", sep = " ", 
                       stringsAsFactors = FALSE, header = FALSE)
popFile <- tmpTable[,1:2]
popFile$V1 <- "-"
# popFile$V1[grep("129", popFile$V2)] <- "129"
# popFile$V1[grep("BALB", popFile$V2)] <- "BALB"
# popFile$V1[grep("C57B", popFile$V2)] <- "C57BL"
# popFile$V1[grep("C3H", popFile$V2)] <- "C3H"
# popFile$V1[grep("FVB", popFile$V2)] <- "FVB"

popFile$V1[grep("ERR9878960", popFile$V2)] <- "129"
popFile$V1[grep("ERR9880211", popFile$V2)] <- "BALB"
popFile$V1[grep("ERR9880473", popFile$V2)] <- "C3H"
popFile$V1[grep("ERR9880493", popFile$V2)] <- "C57BL"
popFile$V1[grep("ERR9880654", popFile$V2)] <- "FVB"

write.table(popFile[,1], "admixtureBed.pop", col.names = FALSE,
            row.names = FALSE, quote = FALSE)


system("admixture admixtureBed.bed 5 --supervised")


### same but for next directory
### should be automated later on
### 

setwd("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-142-MG_cho_20210701_357_353/")
system("cp /mnt/DATA6/mouseData/admixtureFiles/admixtureSnpNameList.txt ./")
listOfFies <- system("ls *geno.vcf", intern = TRUE)


filtFiles <- NULL
for (i in listOfFies) {
  strippedName <- str_remove(i, ".geno.vcf")
  cmd <- paste("bcftools view --exclude-types indels,mnps,other -i'ID=@admixtureSnpNameList.txt'", i, "-Oz -o", paste0(strippedName, ".filt.vcf.gz"))
  system(cmd)
  filtFiles <- c(filtFiles, paste0(strippedName, ".filt.vcf.gz"))
}


cmd2 <- "for file in *.filt.vcf.gz ; do   tabix -f -p vcf ${file}  ; done"
system(cmd2)


filtFiles2 <- c(filtFiles, refFiltFiles)
writeLines(filtFiles2, "listOfFiltFiles.txt")

cmd6 <- "bcftools merge --force-samples --file-list listOfFiltFiles.txt -Oz -o admixturetest.vcf.gz"
system(cmd6)

system("tabix -f -p vcf admixturetest.vcf.gz")
setwd("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-142-MG_cho_20210701_357_353/")
system("plink --vcf admixturetest.vcf.gz --make-bed --out admixtureBed --const-fid 0 --allow-no-sex --no-parents --geno 0.5  --mind 0.9")


tmpTable <- read.table("admixtureBed.fam", sep = " ", 
                       stringsAsFactors = FALSE, header = FALSE)
popFile <- tmpTable[,1:2]
popFile$V1 <- "-"
# popFile$V1[grep("129", popFile$V2)] <- "129"
# popFile$V1[grep("BALB", popFile$V2)] <- "BALB"
# popFile$V1[grep("C57B", popFile$V2)] <- "C57BL"
# popFile$V1[grep("C3H", popFile$V2)] <- "C3H"
# popFile$V1[grep("FVB", popFile$V2)] <- "FVB"

popFile$V1[grep("ERR9878960", popFile$V2)] <- "129"
popFile$V1[grep("ERR9880211", popFile$V2)] <- "BALB"
popFile$V1[grep("ERR9880473", popFile$V2)] <- "C3H"
popFile$V1[grep("ERR9880493", popFile$V2)] <- "C57BL"
popFile$V1[grep("ERR9880654", popFile$V2)] <- "FVB"


write.table(popFile[,1], "admixtureBed.pop", col.names = FALSE,
            row.names = FALSE, quote = FALSE)


system("admixture admixtureBed.bed 5 --supervised")




setwd("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-138-MG_cho_20210621_354_343/")
system("cp /mnt/DATA6/mouseData/admixtureFiles/admixtureSnpNameList.txt ./")
listOfFies <- system("ls *geno.vcf", intern = TRUE)


filtFiles <- NULL
for (i in listOfFies) {
  strippedName <- str_remove(i, ".geno.vcf")
  cmd <- paste("bcftools view --exclude-types indels,mnps,other -i'ID=@admixtureSnpNameList.txt'", i, "-Oz -o", paste0(strippedName, ".filt.vcf.gz"))
  system(cmd)
  filtFiles <- c(filtFiles, paste0(strippedName, ".filt.vcf.gz"))
}


cmd2 <- "for file in *.filt.vcf.gz ; do   tabix -f -p vcf ${file}  ; done"
system(cmd2)


filtFiles2 <- c(filtFiles, refFiltFiles)
writeLines(filtFiles2, "listOfFiltFiles.txt")

cmd6 <- "bcftools merge --force-samples --file-list listOfFiltFiles.txt -Oz -o admixturetest.vcf.gz"
system(cmd6)

system("tabix -f -p vcf admixturetest.vcf.gz")
setwd("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-138-MG_cho_20210621_354_343/")
system("plink --vcf admixturetest.vcf.gz --make-bed --out admixtureBed --const-fid 0 --allow-no-sex --no-parents --geno 0.5  --mind 0.9")


tmpTable <- read.table("admixtureBed.fam", sep = " ", 
                       stringsAsFactors = FALSE, header = FALSE)
popFile <- tmpTable[,1:2]
popFile$V1 <- "-"
# popFile$V1[grep("129", popFile$V2)] <- "129"
# popFile$V1[grep("BALB", popFile$V2)] <- "BALB"
# popFile$V1[grep("C57B", popFile$V2)] <- "C57BL"
# popFile$V1[grep("C3H", popFile$V2)] <- "C3H"
# popFile$V1[grep("FVB", popFile$V2)] <- "FVB"

popFile$V1[grep("ERR9878960", popFile$V2)] <- "129"
popFile$V1[grep("ERR9880211", popFile$V2)] <- "BALB"
popFile$V1[grep("ERR9880473", popFile$V2)] <- "C3H"
popFile$V1[grep("ERR9880493", popFile$V2)] <- "C57BL"
popFile$V1[grep("ERR9880654", popFile$V2)] <- "FVB"


write.table(popFile[,1], "admixtureBed.pop", col.names = FALSE,
            row.names = FALSE, quote = FALSE)


system("admixture admixtureBed.bed 5 --supervised")




### reading res and comparing to haploqa


fam141 <- read.table("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/admixtureBed.fam",
                     sep = " ", stringsAsFactors = FALSE, header = FALSE)
qtable141 <- read.table("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/admixtureBed.5.Q",
                        sep = " ", stringsAsFactors = FALSE, header = FALSE)
rownames(qtable141) <- fam141$V2
rownames(qtable141)[(nrow(qtable141)-4):nrow(qtable141)] <- c("129", "BALB", "C3H", "C57Bl/6j", "FVB")

fam142 <- read.table("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-142-MG_cho_20210701_357_353/admixtureBed.fam",
                     sep = " ", stringsAsFactors = FALSE, header = FALSE)
qtable142 <- read.table("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-142-MG_cho_20210701_357_353/admixtureBed.5.Q",
                        sep = " ", stringsAsFactors = FALSE, header = FALSE)
rownames(qtable142) <- fam142$V2
rownames(qtable142)[(nrow(qtable142)-4):nrow(qtable142)] <- c("129", "BALB", "C3H", "C57Bl/6j", "FVB")


fam138 <- read.table("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-138-MG_cho_20210621_354_343/admixtureBed.fam",
                     sep = " ", stringsAsFactors = FALSE, header = FALSE)
qtable138 <- read.table("/mnt/DATA6/mouseData/testAdmixtureVcfs/Auto_user_AUS5-138-MG_cho_20210621_354_343/admixtureBed.5.Q",
                        sep = " ", stringsAsFactors = FALSE, header = FALSE)
rownames(qtable138) <- fam138$V2
rownames(qtable138)[(nrow(qtable138)-4):nrow(qtable138)] <- c("129", "BALB", "C3H", "C57Bl/6j", "FVB")

combinedQTable <- rbind(qtable141, qtable142, qtable138)
colnames(combinedQTable)[1:5] <- c("129", "BALB", "C3H", "C57Bl/6j", "FVB")
combinedQTable <- combinedQTable * 100


haploQaList <- paste0("/mnt/DATA5/tmp/kev/misc/20221001haploQA", "/", system("ls /mnt/DATA5/tmp/kev/misc/20221001haploQA", intern = TRUE))

i <- haploQaList[1]
fullHaplo <- NULL
for (i in haploQaList) {
  tmpFile <- read.table(i, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  ### splitting the two haplotypes and adding them
  tmpFile$percent_of_genome2 <- tmpFile$percent_of_genome/2
  tmpFile2 <- data.frame(rbind(cbind(tmpFile$haplotype_1, tmpFile$percent_of_genome2),
                               cbind(tmpFile$haplotype_2, tmpFile$percent_of_genome2)))
  tmpFile2$X2 <- as.numeric(tmpFile2$X2)
  res <- NULL
  for (j in unique(tmpFile2$X1)) {
    res <- rbind(res, c(tmpFile$original_sample_id[1], j, sum(tmpFile2$X2[which(tmpFile2$X1 == j)])))
  }
  res <- data.frame(res)
  fullHaplo <- rbind(fullHaplo, res)
}

fullHaplo$string <- paste(fullHaplo$X1, fullHaplo$X2, fullHaplo$X3)
fullHaplo <- fullHaplo[-which(duplicated(fullHaplo$string)),]

combinedQTable$sample <- rownames(combinedQTable)

refQTable <- combinedQTable[grep("129|FVB|BALB|C3H|SRA", combinedQTable$sample),]
combinedQTable2 <- combinedQTable[-grep("129|FVB|BALB|C3H|SRA", combinedQTable$sample), ]
combinedQTable2$sample <- str_remove(nameStripper(combinedQTable2$sample), "\\-")

fullHaplo$X3 <- as.numeric(fullHaplo$X3)
colnames(fullHaplo)[1:3] <- c("sample", "strain", "percent") 
fullHaplo2 <- dcast(data = fullHaplo, formula = sample ~ factor(strain), value.var = "percent",
                    fun.aggregate = sum)
fullHaplo2$sample <- nameStripper(fullHaplo2$sample)


combinedQTable3 <- combinedQTable2[which(combinedQTable2$sample %in% fullHaplo2$sample), ]
fullHaplo3 <- fullHaplo2[which(fullHaplo2$sample %in% combinedQTable3$sample), ]

compTable <- cbind(combinedQTable3, fullHaplo3[match(combinedQTable3$sample, fullHaplo3$sample),])
compTable2 <- compTable[, c(6, 1:5, 8:12)]
rownames(compTable2) <- NULL
colnames(compTable2)[2:6] <- paste0("ngs_", colnames(compTable2)[2:6])
colnames(compTable2)[7:11] <- paste0("snp_", colnames(compTable2)[7:11])
compTable2[,7:11] <- compTable2[,7:11] * -1
compTable3 <- melt(compTable2)

compTable3 <- compTable3[order(compTable3$value, decreasing = TRUE),]
compTable3$value <- signif(compTable3$value, digits = 2)
colnames(compTable3) <- c("SampleName", "label","Fraction")
compTable3$SampleName <- factor(compTable3$SampleName, levels = unique(compTable3$SampleName))
compTable3$label <- as.character(compTable3$label)
compTable3$label <- str_remove(compTable3$label, "X")

compTable3$color <- "#808080"
compTable3$color[grep("129", compTable3$label)] <- "#009E73"
compTable3$color[grep("BALB", compTable3$label)] <- "#D55E00"
compTable3$color[grep("C3H", compTable3$label)] <- "#001a9e"
compTable3$color[grep("C57", compTable3$label)] <- "#9e007e"
compTable3$color[grep("FVB", compTable3$label)] <- "#9e9600"

# colPalette <- c("#009E73","#D55E00", "#001a9e", "#9e007e", "#9e9600")

ggplot(compTable3, aes(x = SampleName, y = Fraction, fill=color)) +
  geom_bar(position="fill", stat="identity", color = "black") + 
  theme_bw() + 
  # scale_fill_manual(values=colPalette) + theme_bw() + 
  theme(axis.text.x=element_text(angle=90,hjust=1))



### simplest way to show how we miss certain areas of the haplotype map is maybe the biplot showing
### all the mutations and how they group to certain genotypes
