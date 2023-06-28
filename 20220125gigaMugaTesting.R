library("MUGAExampleData")
data(geno)
data("FinalReport1")
data("FinalReport2")


### for tomorrow at least parse the file reports and create the QA metrics mentioned in seminal paper 
### this function should only work for these .rda files. the real text files i just need to read.table and skip pertinent lines

headerName <- unlist(strsplit(FinalReport1[10], "\t"))




finalReport1_table <- do.call(rbind, sapply(FinalReport1[11:length(FinalReport1)], function(x) parsedLine <- strsplit(x, "\t")))
finalReport2_table <- do.call(rbind, sapply(FinalReport2[11:length(FinalReport2)], function(x) parsedLine <- strsplit(x, "\t")))

rownames(finalReport1_table) <- NULL
rownames(finalReport2_table) <- NULL

colnames(finalReport1_table) <- headerName
colnames(finalReport2_table) <- headerName

finalReport1_table <- data.frame(finalReport1_table, stringsAsFactors = FALSE)
finalReport2_table <- data.frame(finalReport2_table, stringsAsFactors = FALSE)

finalReport1_table[,5:11] <- lapply(finalReport1_table[,5:11], as.numeric)
finalReport2_table[,5:11] <- lapply(finalReport2_table[,5:11], as.numeric)

finalReport1_table$R2 <- sqrt(finalReport1_table$X + finalReport1_table$Y)
finalReport2_table$R2 <- sqrt(finalReport2_table$X + finalReport2_table$Y)

testDist <- rnorm(10000, mean = 0.97, sd = 0.42)

sampleStats <- NULL
for (i in unique(finalReport1_table$Sample.ID)) {
  tmp <- finalReport1_table[which(finalReport1_table$Sample.ID == i),]
  ksRes <- ks.test(tmp$R2, testDist)
  res <- c(i, mean(tmp$R2), sd(tmp$R2), ksRes$statistic, ksRes$p.value)
  sampleStats <- rbind(sampleStats, res)
}

sampleStats2 <- NULL
for (i in unique(finalReport1_table$Sample.ID)) {
  tmp <- finalReport1_table[which(finalReport1_table$Sample.ID == i),]
  ksRes <- ks.test(tmp$R2, testDist)
  res <- c(i, mean(tmp$R2), sd(tmp$R2), ksRes$statistic, ksRes$p.value)
  sampleStats2 <- rbind(sampleStats2, res)
}





library(argyle)
data(ex)
summary(ex)

## see the marker map
map <- markers(ex)
head(map)

## see sample metadata
mice <- samples(ex)
head(mice)

intens <- intensity(ex)
lapply(intens, dim)


load("/mnt/DATA5/tmp/kev/tmpDbs/uncGigaMuga/clusters.gigamuga.Rdata")
load("/mnt/DATA5/tmp/kev/tmpDbs/uncGigaMuga/snps.gigamuga.Rdata")

