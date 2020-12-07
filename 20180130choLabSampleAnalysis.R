
library(gplots)
library(RColorBrewer)
dataSetRawCalls <- rawGeneEst
dataSetRawCalls.2 <- as.matrix(dataSetRawCalls[,2:ncol(dataSetRawCalls)], nrow= nrow(dataSetRawCalls))

heatmap((dataSetRawCalls.2), Colv = hclust((dataSetRawCalls.2)))




heatmap.2(dataSetRawCalls.2, Rowv = FALSE, trace = "none", symm=F,symkey=F,symbreaks=T, scale="none", dendrogram = c("column"),col=my_palette, 
          breaks=colors, density.info="none")


geneEst.2 <- as.matrix(geneEst[,2:ncol(geneEst)], nrow = nrow(geneEst))


heatmap.2(dataSetRawCalls.2, col=my_palette, 
          density.info="none", trace="none", 
          dendrogram=c("column"), symm=F,symkey=F,symbreaks=T, scale="none")

hmcol<-brewer.pal(8,"YlOrRd")
heatmap(data,col=hmcol)

geneEst.3 <- t(log2(geneEst.2))
geneEst.3 <- (geneEst.3[,c("Brca1","Trp53","Nf1","Rb1")])


dataSetRawCalls.3 <- t(dataSetRawCalls.2)
dataSetRawCalls.3 <- dataSetRawCalls.3[,c("Brca1","Trp53","Nf1","Rb1")]

heatmap.2(dataSetRawCalls.3, col = hmcol, Colv = FALSE,
          density.info="none", trace="none", scale = "column",
          dendrogram=c("row"), cexRow = 0.5)



