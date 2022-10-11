### for the snakemake as a whole, probably best to just do it by report and have the target file be list of files to process i.e the norm.vcf.gz
### then the output file be some type of combined vcf
### to do this though, I would need to to run bcftools from an Rscript which is a bit odd ..... 

### if I don't want to do it by report can just have one directory and the input is same, output target would be the finished filtered file
### this would be way harder downstream to figure out which combinedVcf for admixture belong to which 


### thing to note is I need to have separate module for admixture since I'd need to create docker for bcftools and docker
library(stringr)

setwd("/mnt/DATA6/mouseData/admixtureVcfs/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/")
listOfFies <- system("ls *norm.vcf.gz", intern = TRUE)


filtFiles <- NULL
for (i in listOfFies) {
  strippedName <- str_remove(i, ".norm.vcf.gz")
  cmd <- paste("bcftools view --exclude-types indels,mnps,other -i'ID=@admixtureSnpNameList.txt'", i, "-Oz -o", paste0(strippedName, ".filt.vcf.gz"))
  system(cmd)
  filtFiles <- c(filtFiles, paste0(strippedName, ".filt.vcf.gz"))
}


cmd2 <- "for file in *.filt.vcf.gz ; do   tabix -f -p vcf ${file}  ; done"
system(cmd2)



### below is commented out because I only used it once to filter vcf files for SNPs
# refFileDir <- "/mnt/DATA6/mouseData/admixtureFiles/"
# setwd(refFileDir)
# cmd3 <- "ls *_sorted.vcf.gz"
# refFiles <- system(cmd3, intern = TRUE)
# refFiles <- paste0(refFileDir, refFiles)


### don't need this, found vcfs with chr
# withChrFiles <- NULL
# for (i in refFiles) {
#   strippedName <- str_remove(i, "_sorted.vcf.gz")
#   cmd <- paste("bcftools annotate --rename-chrs chromosomes.txt",  i, "-Oz -o", paste0(strippedName, ".with_chr.vcf.gz"))
#   system(cmd)
#   withChrFiles <- c(withChrFiles, paste0(strippedName, ".with_chr.vcf.gz"))
# }
# cmd4 <- "for file in *.with_chr.vcf.gz ; do   tabix -f -p vcf ${file}  ; done"
# system(cmd4)
# 
# refFiltFiles <- NULL
# for (i in withChrFiles) {
#   strippedName <- str_remove(i, ".with_chr.vcf.gz")
#   cmd <- paste("bcftools view -R", "KC-05_MG_X64.filt.vcf.gz",i, "-Oz -o", paste0(strippedName, ".filt.vcf.gz"))
#   system(cmd)
#   refFiltFiles <- c(refFiltFiles, paste0(strippedName, ".filt.vcf.gz"))
# }
# 
# cmd5 <- "for file in *.filt.vcf.gz ; do   tabix -f -p vcf ${file}  ; done"
# system(cmd5)

setwd("/mnt/DATA6/mouseData/admixtureVcfs/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/")
filtFiles2 <- c(filtFiles, refFiltFiles)
writeLines(filtFiles2, "listOfFiltFiles.txt")

cmd6 <- "bcftools merge --force-samples --file-list listOfFiltFiles.txt -Oz -o admixturetest.vcf.gz"
system(cmd6)

system("tabix -f -p vcf admixturetest.vcf.gz")

### need to run twice because I can't filter by geno since it's not ordered. the new file is ordered, so 
### rerun command on new file
system("plink --vcf admixturetest.vcf.gz --make-bed --out admixtureBed --const-fid 0 --allow-no-sex --no-parents --geno 0.1")


tmpTable <- read.table("admixtureBed.fam", sep = " ", 
                       stringsAsFactors = FALSE, header = FALSE)
popFile <- tmpTable[,1:2]
popFile$V1 <- "-"
popFile$V1[grep("129", popFile$V2)] <- "129"
popFile$V1[grep("BALB", popFile$V2)] <- "BALB"
popFile$V1[grep("C57B", popFile$V2)] <- "C57BL"
popFile$V1[grep("C3H", popFile$V2)] <- "C3H"
popFile$V1[grep("FVB", popFile$V2)] <- "FVB"


write.table(popFile[,1], "admixtureBed.pop", col.names = FALSE,
            row.names = FALSE, quote = FALSE)


system("admixture admixtureBed.bed 5 --supervised")



### top portion was for beginning to end including making necessary files
### input needs to be report or dir cause everything else seems static
### changed slightly since sometimes SNPs aren't annotated

setwd("/mnt/DATA6/mouseData/admixtureVcfs/Auto_user_AUS5-141-MG_cho_202106_3TS_356_349/")
system("cp /mnt/DATA6/mouseData/admixtureFiles/admixtureSnpNameList.txt ./")
system("cp /mnt/DATA6/mouseData/admixtureFiles/KC-05_MG_X64.filt.vcf.gz ./snpPost.vcf.gz")
system("cp /mnt/DATA6/mouseData/admixtureFiles/KC-05_MG_X64.filt.vcf.gz.tbi ./snpPost.vcf.gz.tbi")
listOfFies <- system("ls *norm.vcf.gz", intern = TRUE)


filtFiles <- NULL
for (i in listOfFies) {
  strippedName <- str_remove(i, ".norm.vcf.gz")
  cmd <- paste("bcftools view -R snpPost.vcf.gz --exclude-types indels,mnps,other", i, "-Oz -o", paste0(strippedName, ".filt.vcf.gz"))
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

### need to run twice because I can't filter by geno since it's not ordered. the new file is ordered, so 
### rerun command on new file
system("plink --vcf admixturetest.vcf.gz --make-bed --out admixtureBed --const-fid 0 --allow-no-sex --no-parents --geno 0.1")


tmpTable <- read.table("admixtureBed.fam", sep = " ", 
                       stringsAsFactors = FALSE, header = FALSE)
popFile <- tmpTable[,1:2]
popFile$V1 <- "-"
popFile$V1[grep("129", popFile$V2)] <- "129"
popFile$V1[grep("BALB", popFile$V2)] <- "BALB"
popFile$V1[grep("C57B", popFile$V2)] <- "C57BL"
popFile$V1[grep("C3H", popFile$V2)] <- "C3H"
popFile$V1[grep("FVB", popFile$V2)] <- "FVB"


write.table(popFile[,1], "admixtureBed.pop", col.names = FALSE,
            row.names = FALSE, quote = FALSE)


system("admixture admixtureBed.bed 5 --supervised")

