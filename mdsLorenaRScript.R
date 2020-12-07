.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))

print(.libPaths())
library(edgeR)

tableOfDat <- read.table(file = "/mnt/DATA4/kevhu/LorenaData/20171130UT_e2e_Combined.Uniq.txt", sep = "\t", header = TRUE)
tableOfDat <- tableOfDat[,c(1:6,which(grepl("total_cov",colnames(tableOfDat))))]
tableOfDat$Gene <- gsub("\\;.*","", tableOfDat$attributes)
tableOfDat$Gene <- gsub("GENE_ID=","", tableOfDat$Gene)

edgeRmat <- tableOfDat[,grepl("total_cov",colnames(tableOfDat))]
rownames(edgeRmat) <- tableOfDat$Gene

dummyMat <- NULL
edgeRmat.cpm <- cpm(edgeRmat)
for(i in 1:nrow(edgeRmat)){
  dummyRow <- NULL
  for(j in 1:(ncol(edgeRmat)/2)){
    if((sum(edgeRmat.cpm[i,c(2*j-1,2*j)] > 5) == 2 | sum(edgeRmat.cpm[i,c(2*j-1,2*j)] > 5) == 0) == TRUE){
      k <- 1
    }
    else{
      k <- 0
    }
    dummyRow <- c(dummyRow, k)
  }
  dummyMat <- rbind(dummyMat, dummyRow)
}


keepFirstFil <- which(rowSums(dummyMat) >= quantile(rowSums(dummyMat), .05))
edgeRmat2 <- edgeRmat[keepFirstFil,]
d <- dist(edgeRmat2)
fit <- cmdscale(d,eig=TRUE, k=2)

saveRDS(fit, "/mnt/DATA4/kevhu/mdsPlotObj")


