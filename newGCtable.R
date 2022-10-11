###create gc bed 
gcbed <- read.table("/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-166-Cho_mouse_1_353_410/plugin_out/coverageAnalysis_out.721/local_beds/IAD124056_167_Designed.gc.bed", 
           skip = 1)



gcbed <- read.table("/home/hovelson/lists/bed/OCP.20150630.designed.noTrack.GC.bed", 
                    skip = 1)

### with pool number
col4 <- sub("*;","",sub("GENE_ID=","", gcbed$V5))
col4 <- gsub("[[:punct:]]","",col4)

### without pool number
col4 <- sub("*;......","",sub("GENE_ID=","", gcbed$V5))

totalBase <- gcbed$V3 - gcbed$V2
percentage <- gcbed$V6/totalBase

gcbed <- cbind(gcbed, col4, totalBase, percentage)

newGcBed <- gcbed[,c(1,2,3,4,6,8,9,7)]

write.table(newGcBed, "/home/kevhu/data/IAD124056_167_Designed.test.gc.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



write.table(newGcBed, "/home/kevhu/data/IAD124056_167_Designed.test.gc.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


write.table(newGcBed, "/home/kevhu/data/IAD124056_167_Designed.nopool.gc.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(newGcBed, "/home/kevhu/data/IAD124056_167_Designed.pool.gc.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")




###testing for OCPV3
###first OCPv4 was from Hayes runs path: /mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-130-Hayes_CTC_AC05_Run1_OCPv3_295_328/plugin_out/coverageAnalysis_out.621/local_beds/OCP.20150630.designed.gc.bed
gcbed <- read.table("/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-108-SQBladder_DNA_RNA_270_266/plugin_out/coverageAnalysis_out.476/local_beds/OCP3.20140506.designed.gc.bed",
                    skip = 1)

### with pool number
col4 <- gsub("*GENE_ID=","", gcbed$V5)
col4 <- strsplit(col4, ";")
finCol4 <- NULL
for(i in seq_along(col4)){
  dummy <- paste(col4[[i]][1],col4[[i]][2],sep = "")
  finCol4 <- c(finCol4, dummy)
}



totalBase <- gcbed$V3 - gcbed$V2
percentage <- gcbed$V6/totalBase

newgcbed <- cbind(gcbed[,c(1,2,3,4,6)], totalBase, percentage, finCol4)
newgcbed$finCol4 <- gsub("[[:punct:]]", "",newgcbed$finCol4)




write.table(newgcbed, "/home/kevhu/data/OCPv3.test.gc.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


write.table(newGcBed, "/home/kevhu/data/IAD124056_167_Designed.nopool.gc.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


write.table(newGcBed, "/mnt/DATA4/kevhu/testData/OCP1c.nopool.gc.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")




###OCP1v


gcbed <- read.table("/home/hovelson/lists/bed/OCP.20150630.designed.noTrack.GC.bed", 
                    skip = 1)


col4 <- gsub("*GENE_ID=","", gcbed$V8)
col4 <- strsplit(col4, ";")
finCol4 <- NULL
for(i in seq_along(col4)){
  dummy <- paste(col4[[i]][1],col4[[i]][2],sep = "")
  finCol4 <- c(finCol4, dummy)
}



#totalBase <- gcbed$V3 - gcbed$V2
#percentage <- gcbed$V6/totalBase

newgcbed <- cbind(gcbed[,c(1:7)], finCol4)
newgcbed$finCol4 <- gsub("[[:punct:]]*", "",newgcbed$finCol4)
newgcbed$finCol4 <- gsub("Pool[0-9]", "",newgcbed$finCol4)


write.table(newgcbed, "/mnt/DATA4/kevhu/testData/OCP1c.nopool.gc.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

###this is for the new PAN GU bed file - slgihtly different format so I got lazy and used the the gc file from IonTorrent
gcbed <- read.table("/home/kevhu/data/WG_IAD127899.20170720.designed.gc.bed",
                    skip = 1, stringsAsFactors = FALSE)

### without pool number
totalBase <- gcbed$V3 - gcbed$V2
percentage <- gcbed$V6/totalBase
#gene <- gsub("*GENE_ID=","", gcbed$V5)
#gene <- gsub("[[:punct:]]*", "",gene)
#gene <- gsub("Pool[0-9]", "",gene)

###going to make own gene table b/c the current one is off... 

genecov <- read.table("/home/kevhu/data/PanGUpanel_coverage.csv", sep = ",", stringsAsFactors = FALSE)
genecov <- genecov[2:nrow(genecov),1:12]
#genecov$newchr <- gsub("chr","",genecov$V1)
#genecov$newchr <- as.numeric(gsub("X","23", genecov$newchr))
#genecov.reord <- genecov[order(genecov$newchr,genecov$V2),]

gcbed$newGene <- NULL
for(i in 1:nrow(gcbed)){
  gcbed$newGene[i] <- genecov$V10[which(genecov$V4 == gcbed$V4[i])]
  print(genecov$V10[which(genecov$V4 == gcbed$V4[i])])
}

gene <- gcbed$newGene
gcbed$newGene <- NULL

gcbed <- cbind(gcbed, gene, totalBase, percentage)

newGcBed <- gcbed[,c(1,2,3,4,6,8,9,7)]

write.table(newGcBed, "/home/kevhu/data/WG_IAD127899.20170720.designed.forscript.GC.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")





