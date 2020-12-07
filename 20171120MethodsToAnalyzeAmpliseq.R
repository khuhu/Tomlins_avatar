library(ggplot2)
library(gridExtra)

nonUniqAmp <- read.table("/mnt/DATA4/kevhu/LorenaData/20171121CombinedAmpUT.nonUniq.txt", header = TRUE, stringsAsFactors = FALSE,
                         sep = "\t")
nonUniqAmp$tot_cov <- apply(nonUniqAmp[,c("fwd_cov","rev_cov")], 1, sum)
nonUniqAmp$logtot_cov <- log10(nonUniqAmp$tot_cov +1)
nonUniqAmp$logtot_reads <- log10(nonUniqAmp$total_reads + 1)
plot1 <- ggplot(data = nonUniqAmp, aes(logtot_reads,logtot_cov)) + geom_point(alpha=0.1)


uniqAmp <- read.table("/mnt/DATA4/kevhu/LorenaData/20171127CombinedAmp.Uniq.txt", header = TRUE, stringsAsFactors = FALSE,
                      sep = "\t")
uniqAmp$tot_cov <- apply(uniqAmp[,c("fwd_cov","rev_cov")], 1, sum)
uniqAmp$logtot_cov <- log10(uniqAmp$tot_cov +1)
uniqAmp$logtot_reads <- log10(uniqAmp$total_reads + 1)
uniqAmp$tot_covReads <- uniqAmp$fwd_cov + uniqAmp$rev_cov
uniqAmp$logtot_covReads <- log10(uniqAmp$tot_covReads + 1)
plot2 <- ggplot(data = uniqAmp, aes(logtot_reads,logtot_cov)) + geom_point(alpha=0.1)


nonUniqCov <- read.table("/mnt/DATA4/kevhu/LorenaData/20171121CombinedCovUT.txt", header = TRUE, stringsAsFactors = FALSE,
                         sep = "\t")
nonUniqCov$tot_e2e <- apply(nonUniqCov[,c("fwd_e2e","rev_e2e")], 1, sum)
nonUniqCov$logtot_e2e <- log10(nonUniqCov$tot_e2e +1)
nonUniqCov$logtot_reads <- log10(nonUniqCov$total_reads + 1)
plot3 <- ggplot(data = nonUniqCov, aes(logtot_reads,logtot_e2e)) + geom_point(alpha=0.1)




grid.arrange(plot1,plot2,plot3, ncol = 2)

dummyDf <- data.frame("attritubtes" = nonUniqAmp$attributes,"log10tot_e2e" = nonUniqCov$logtot_e2e, "log10tot_uniqCov" = uniqAmp$logtot_cov,
                      "log10tot_nonUniqTot" = nonUniqAmp$logtot_reads, "log10totCov_nonUniqAmp" = nonUniqAmp$logtot_cov)
dummyDf$gene <- sub("GENE_ID=*","",dummyDf$attritubtes)
dummyDf$gene <- sub(";(.*)","",dummyDf$gene)
totVsE2e <- ggplot(data = dummyDf, aes(log10tot_nonUniqTot,log10tot_e2e)) + geom_point(alpha=0.1) + 
  geom_text(data = dummyDf[which(grepl("AMY1", dummyDf$attritubtes)),],aes(log10tot_nonUniqTot,log10tot_e2e,label = gene),inherit.aes = FALSE)
totVsUniq <- ggplot(data = dummyDf, aes(log10tot_nonUniqTot,log10tot_uniqCov)) + geom_point(alpha=0.1) +
  geom_text(data = dummyDf[which(grepl("AMY1", dummyDf$attritubtes)),],aes(log10tot_nonUniqTot,log10tot_uniqCov,label = gene),inherit.aes = FALSE)
e2eVsUniq <- ggplot(data = dummyDf, aes(log10tot_e2e,log10tot_uniqCov)) + geom_point(alpha=0.1) +
  geom_text(data = dummyDf[which(grepl("AMY1", dummyDf$attritubtes)),],aes(log10tot_e2e,log10tot_uniqCov,label = gene),inherit.aes = FALSE)
nonUniqCovVsUniqCov <- ggplot(data = dummyDf, aes(log10totCov_nonUniqAmp,log10tot_uniqCov)) + geom_point(alpha=0.1) +
  geom_text(data = dummyDf[which(grepl("AMY1", dummyDf$attritubtes)),],aes(log10tot_e2e,log10tot_uniqCov,label = gene),inherit.aes = FALSE)

grid.arrange(totVsE2e, totVsUniq, nonUniqCovVsUniqCov,e2eVsUniq, ncol = 2)



which(grepl("AMY1", dummyDf$attritubtes))



### Analysis for Kelly's Data
kellyAmpUniq <- read.table("/mnt/DATA4/kevhu/KellyData/20171121CombinedAmp.Uniq.txt", header = TRUE, stringsAsFactors = FALSE,
                      sep = "\t")

kellyAmpUniq$tot_cov <- apply(kellyAmpUniq[,c("fwd_cov","rev_cov")], 1, sum)
kellyAmpUniq$logtot_cov <- log10(kellyAmpUniq$tot_cov +1)
kellyAmpUniq$logtot_reads <- log10(kellyAmpUniq$total_reads + 1)

kellyCovNonUniq <- read.table("/mnt/DATA4/kevhu/KellyData/20171121CombinedCov.nonUniq.txt", header = TRUE, stringsAsFactors = FALSE,
                      sep = "\t")

kellyCovNonUniq$tot_e2e <- apply(kellyCovNonUniq[,c("fwd_e2e","rev_e2e")], 1, sum)
kellyCovNonUniq$logtot_e2e <- log10(kellyCovNonUniq$tot_e2e +1)
kellyCovNonUniq$logtot_reads <- log10(kellyCovNonUniq$total_reads + 1)

kellyAmpNonUniq <- read.table("/mnt/DATA4/kevhu/KellyData/20171121CombinedAmp.nonUniq.txt", header = TRUE, stringsAsFactors = FALSE,
                           sep = "\t")

kellyAmpNonUniq$tot_cov <- apply(kellyAmpNonUniq[,c("fwd_cov","rev_cov")], 1, sum)
kellyAmpNonUniq$logtot_cov <- log10(kellyAmpNonUniq$tot_cov +1)
kellyAmpNonUniq$logtot_reads <- log10(kellyAmpNonUniq$total_reads + 1)

dummyDf <- data.frame("attritubtes" = kellyAmpNonUniq$attributes,"log10tot_e2e" = kellyCovNonUniq$logtot_e2e, "log10tot_uniqCov" = kellyAmpUniq$logtot_cov,
                      "log10tot_nonUniqTot" =kellyAmpNonUniq$logtot_reads)
dummyDf$gene <- sub("GENE_ID=*","",dummyDf$attritubtes)
dummyDf$gene <- sub(";(.*)","",dummyDf$gene)

totVsE2e <- ggplot(data = dummyDf, aes(log10tot_nonUniqTot,log10tot_e2e)) + geom_point(alpha=0.1)
totVsUniq <- ggplot(data = dummyDf, aes(log10tot_nonUniqTot,log10tot_uniqCov)) + geom_point(alpha=0.1)
e2eVsUniq <- ggplot(data = dummyDf, aes(log10tot_e2e,log10tot_uniqCov)) + geom_point(alpha=0.1)


grid.arrange(totVsE2e, totVsUniq, e2eVsUniq, ncol = 2)





###analysis for Lorena single cohort 
uniqCov.1coh <- read.table("/mnt/DATA4/kevhu/LorenaData/20171128CombinedCov.1cohort.Uniq.txt", header = TRUE, stringsAsFactors = FALSE,
                            sep = "\t")
uniqCov.1coh$tote2e <- apply(uniqCov.1coh[,c("fwd_e2e","rev_e2e")],1,sum)
uniqE2eVsTot <- ggplot(data = uniqCov.1coh, aes(log10(total_reads+1),log10(tote2e +1)))
uniqE2eVsTot <- uniqE2eVsTot + geom_point(alpha=0.1)


uniqAmp.1coh <- read.table("/mnt/DATA4/kevhu/LorenaData/20171128CombinedAmp.1cohort.Uniq.txt", header = TRUE, stringsAsFactors = FALSE,
                           sep = "\t")
uniqAmp.1coh$totcov <- apply(uniqAmp.1coh[,c("fwd_cov","rev_cov")],1,sum)



nonUniqCov.1coh <- read.table("/mnt/DATA4/kevhu/LorenaData/20171128CombinedCov.1cohort.nonUniq.txt", header = TRUE, stringsAsFactors = FALSE,
                           sep = "\t")
nonUniqCov.1coh$tote2e <- apply(nonUniqCov.1coh[,c("fwd_e2e","rev_e2e")],1,sum)

dummyDf <- data.frame("uniqE2e" = uniqCov.1coh$tote2e, "nonUniqE2e" = nonUniqCov.1coh$tote2e, "uniqAmpCov" = uniqAmp.1coh$totcov)

plot1 <- ggplot(data = dummyDf, aes(log10(nonUniqE2e+1), log10(uniqE2e+1))) + geom_point(alpha = 0.1) 
plot2 <- ggplot(data = dummyDf, aes(log10(uniqAmpCov+1),log10(uniqE2e +1))) + geom_point(alpha = 0.1)

grid.arrange(uniqE2eVsTot, plot1, plot2,ncol = 2)
