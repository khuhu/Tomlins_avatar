library(edgeR)

tableOfDat <- read.table(file = "/mnt/DATA4/kevhu/KellyData/LNCap_Transcriptome_Compiled.countMat.csv", sep = ",", header = TRUE)
tableOfDat <- tableOfDat[,c(1:17)]

tableOfDat$Gene <- gsub("\\;.*","", tableOfDat$attributes)
tableOfDat$Gene <- gsub("GENE_ID=","", tableOfDat$Gene)

tableOfDat <- tableOfDat[-c(20813,20814),]
edgeRmat <- tableOfDat[,6:(ncol(tableOfDat)-1)]
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

rownames(dummyMat) <- rownames(edgeRmat)
keepFirstFil <- which(rowSums(dummyMat) >= quantile(rowSums(dummyMat), .05))


#edgeRmat.pooled <- NULL
#new.names <- NULL
###pooling technical replicates, because Marioni et at 2008 states that techincial variance i.e library replicates follows a Poisson distribution
#for(i in seq(2,ncol(edgeRmat), 2)){
#  pooled.name <- colnames(edgeRmat)[i-1]
#  new.names <- c(new.names, pooled.name)
#  pooled <- rowSums(edgeRmat[,c(i-1, i)])
#  edgeRmat.pooled <- cbind(edgeRmat.pooled, pooled)
#}

#colnames(edgeRmat.pooled) <- new.names

dgeObj <- DGEList(counts = edgeRmat)
dgeObj$samples
bioExpressFilt <- c(unlist(which(rowSums(cpm(dgeObj) > 5) >= 2)))

#keep <- intersect(bioExpressFilt, corrFilter)
keep <- intersect(bioExpressFilt, keepFirstFil)

dgeObj <- dgeObj[keep, , keep.lib.sizes = FALSE]

dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples

groups <- factor(c(rep("scr",3),rep("scrHRas",3), rep("sh.p53",3), rep("p53Hars",3)))
#groups <- factor(c(rep("scr",3),rep("scrHRas",3), rep("none",6)))

design <- model.matrix(~0 + groups, data=dgeObj$samples)
design
dgeObj <- estimateDisp(dgeObj, design)
fit <- glmFit(dgeObj, design)
lrt <- glmLRT(fit, contrast=c(0,-1,1,0))
#lrt <- glmLRT(fit, contrast=c(0,-1,1))
topTags(lrt)
lrt.table <- lrt$table
lrtPvals <- lrt$table$PValue
lrtPVals.ad <- p.adjust(lrtPvals, method = c("BH"))
lrt$table$QVal <- lrtPVals.ad
lrt.tab <- lrt$table
lrt.tab$Gene <- rownames(lrt.tab)
lrt.tab$color <- NULL

for(i in seq_along(lrt.tab$logFC)){
  if(abs(lrt.tab$logFC[i]) > 1 & lrt.tab$QVal[i] < 0.05){
    lrt.tab$color[i]  <- "red"
  }
  else{
    lrt.tab$color[i] <- "blue"
  }
}


plot2 <- ggplot(data = lrt.tab,aes(logFC, -log10(QVal)))
plot2 <- plot2 +  geom_point(pch= 20, color = lrt.tab$color, size = 1.25) + geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.25)
#plot2 <-  plot2 + geom_text(data = lrt.tab[which(lrt.tab$color == "red"),], aes(label = Gene, vjust =-1.0), size = 3, alpha = 1.0) + ggtitle("Precursor vs Cancer")
plot2 <-  plot2 + ggtitle("scr scrHras")


plot2 <- plot2 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot2


sig.genes <- lrt.tab[which(lrt.tab$color == "red"),]
write.table(sig.genes,"/mnt/DATA4/kevhu/KellyData/edgeRsigGenesScrHras.txt", quote = FALSE, row.names = FALSE, sep = "\t")




design <- model.matrix(~0 + groups, data=dgeObj$samples)
dgeObj <- estimateDisp(dgeObj,design)
design
fit <- glmFit(dgeObj, design)
lrt <- glmLRT(fit, contrast=c(1,0,0,-1))
topTags(lrt)
lrt.table <- lrt$table
lrtPvals <- lrt$table$PValue
lrtPVals.ad <- p.adjust(lrtPvals, method = c("BH"))
lrt$table$QVal <- lrtPVals.ad
lrt.tab <- lrt$table
lrt.tab$Gene <- rownames(lrt.tab)

lrt.tab$color <- NULL
for(i in seq_along(lrt.tab$logFC)){
  if(abs(lrt.tab$logFC[i]) > 1 & lrt.tab$QVal[i] < 0.05){
    lrt.tab$color[i]  <- "red"
  }
  else{
    lrt.tab$color[i] <- "blue"
  }
}


plot3 <- ggplot(data = lrt.tab,aes(logFC, -log10(QVal)))
plot3 <- plot3 +  geom_point(pch= 20, color = lrt.tab$color, size = 1.25) + geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.25)
#plot3 <-  plot3 + geom_text(data = lrt.tab[which(lrt.tab$color == "red"),], aes(label = Gene, vjust =-1.5), size = 3.0) + ggtitle("CTNBB1 Wt vs Mut")
plot3 <-  plot3 + ggtitle("CTNBB1 Wt vs Mut")

plot3  


sig.genes <- lrt.tab[which(lrt.tab$color == "red"),]
write.table(sig.genes,"/mnt/DATA4/kevhu/KellyData/edgeRsigGenesP53Hras.txt", quote = FALSE, row.names = FALSE, sep = "\t")



