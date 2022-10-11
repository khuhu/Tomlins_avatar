id.tab <- read.table("/mnt/DATA4/kevhu/tmp/combined.allSamples.CNcalls.txt", stringsAsFactors = FALSE,
                     header = TRUE)
id.tab2 <- id.tab[,c(1,2,4)]
listGenes <- unique(id.tab$Gene)

library(reshape2)


dummy <- dcast(id.tab2, formula = Gene ~ Sample, value.var = "CopyNumberRatio")
write.table(dummy, "/mnt/DATA4/kevhu/tmp/geneMatrix.IE.txt", col.names = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE)


id.ie.tab <- read.table("/mnt/DATA4/kevhu/tmp3/combined.allSamples.CNcalls.txt", stringsAsFactors = FALSE,
                     header = TRUE)

dummy2 <- dcast(id.ie.tab, formula = Gene ~ Sample, value.var = "CopyNumberRatio")
write.table(dummy2, "/mnt/DATA4/kevhu/tmp3/geneMatrix.ID.txt", col.names = TRUE, quote = FALSE, sep = "\t",
            row.names = FALSE)
