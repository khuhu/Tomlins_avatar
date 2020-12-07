#tableOfDat <- read.table(file = "/mnt/DATA4/kevhu/LorenaData/20171103LorenaRNAseq.csv", sep = ",", header = TRUE, skip = 7)
library(edgeR)
library(readxl)
library(pheatmap)
tableOfDat <- read.table(file = "/mnt/DATA4/kevhu/LorenaData/20171130UT_e2e_Combined.Uniq.txt", sep = "\t", header = TRUE)

tableOfDat <- tableOfDat[,c(1:6,which(grepl("total_cov",colnames(tableOfDat))))]
tableOfDat$Gene <- gsub("\\;.*","", tableOfDat$attributes)
tableOfDat$Gene <- gsub("GENE_ID=","", tableOfDat$Gene)
#edgeRmat <- tableOfDat[,6:29]
which(grepl("total_cov",colnames(tableOfDat)))

edgeRmat <- tableOfDat[,grepl("total_cov",colnames(tableOfDat))]
rownames(edgeRmat) <- tableOfDat$Gene

#keep <- rowSums(cpm(edgeRmat) > 5) >= 10
#edgeRmat <- edgeRmat[keep, ]

###edit this for so that it keeps 2 and 0. i. concordance between the two technical replicates - bit too conservative
#for(i in 1:(ncol(edgeRmat)/2)){
#  rem <- rowSums(cpm(edgeRmat[,c(2*i-1,2*i)]) > 3) == 2 | rowSums(cpm(edgeRmat[,c(2*i-1,2*i)]) > 3) == 0
#  print(length(rem))
#  edgeRmat <- edgeRmat[rem,]
#}


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




hist(rowSums(dummyMat))
quantile(rowSums(dummyMat), 0.05)
summary(rowSums(dummyMat))

###filter using - how I separate the techinical replicates is a makeshift method and probably need to think harder about it
#listOfCors <- NULL
#set1 <- edgeRmat[,grep("L2", colnames(edgeRmat))]
#set2 <- edgeRmat[,-grep("L2", colnames(edgeRmat))]
#for(i in 1:nrow(edgeRmat)){
#  a <- cor(c(unlist(set1[i,])), c(unlist(set2[i,])), method = c("pearson"))
#  listOfCors <- c(listOfCors, a)
#}
#listOfCors.na <- which(is.na(listOfCors))
#listOfCors2 <- listOfCors[-which(is.na(listOfCors))]

#looking at cors by samp

listOfCors <- NULL
set1 <- edgeRmat[,grep("L2", colnames(edgeRmat))]
set2 <- edgeRmat[,-grep("L2", colnames(edgeRmat))]
for(i in 1:ncol(set1)){
  a <- cor(c(unlist(set1[,i])), c(unlist(set2[,i])), method = c("pearson"))
  listOfCors <- c(listOfCors, a)
}

summary(listOfCors)




#quantile(listOfCors2, 0.25)
#keepFirstFil <- which(listOfCors > quantile(listOfCors2, 0.25))

#edgeRmat <- edgeRmat[which(listOfCors > quantile(listOfCors, 0.25))



###cutoff is 1 sd away and emprically testing when sampling errors for things with cpm < 5
#order of cutoffs is pearson, spearman, and kendall
#edgeRmat <- edgeRmat[which(listOfCors > 0.6697699),]
#edgeRmat <- edgeRmat[which(listOfCors > 0.6457861),]
#edgeRmat <- edgeRmat[which(listOfCors > 0.5123475),]
#edgeRmat <- edgeRmat[which(listOfCors > quantile(listOfCors, 0.1)),]



###testing new filter
#testDf <- edgeRmat
#for(i in 1:(ncol(testDf)/2)){
#  rem <- rowSums(cpm(testDf[,c(2*i-1,2*i)]) > 5) == 2 | rowSums(cpm(testDf[,c(2*i-1,2*i)]) > 5) == 0
#  print(length(rem))
#  testDf <- testDf[rem,]
#}


###Test to see how filtering works
#test <- edgeRmat[10:15,]
#test.keep <- cpm(test) > 5
#rowSums(test.keep) >=2



edgeRmat.pooled <- NULL
new.names <- NULL
###pooling technical replicates, because Marioni et at 2008 states that techincial variance i.e library replicates follows a Poisson distribution
for(i in seq(2,ncol(edgeRmat), 2)){
  pooled.name <- colnames(edgeRmat)[i-1]
  new.names <- c(new.names, pooled.name)
  pooled <- rowSums(edgeRmat[,c(i-1, i)])
  edgeRmat.pooled <- cbind(edgeRmat.pooled, pooled)
}

colnames(edgeRmat.pooled) <- new.names
#rownames(edgeRmat.pooled) <- tableOfDat$Gene

#library(edgeR)

dgeObj <- DGEList(counts = edgeRmat.pooled)
dgeObj$samples
bioExpressFilt <- c(unlist(which(rowSums(cpm(dgeObj) > 5) >= 2)))

#keep <- intersect(bioExpressFilt, corrFilter)
keep <- intersect(bioExpressFilt, keepFirstFil)

dgeObj <- dgeObj[keep, , keep.lib.sizes = FALSE]

##checking what the actual count cutoff is with a cpm of 5
#which(dgeObj$samples$lib.size == min(dgeObj$samples$lib.size))



#can conclude the count cutoff for the smallest library is 50 counts with cpm 5

#keep <- rowSums(cpm(dgeObj) > 5) >= 10
#dgeObj <- dgeObj[keep, , keep.lib.sizes = FALSE]
dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples

### diff groups
preCan <- factor(c("Pre","Can","Can","Can","Pre","Can","Can","Can","Can","Can","Pre","Pre"))
CTNBB1 <- factor(c("WT","Mut","WT","WT","Mut","Mut","WT","Mut","Mut","Mut","WT","WT"))
AdenoSquam <- factor(c("None","Squa","Adeno","Adeno","None","Adeno","Adeno","Adeno","Adeno","Squa","None","None"))
Squam <- factor(c("None","Squa","None","None","None","None","None","None","None","Squa","None","None"))

###preCan

design <- model.matrix(~0 + preCan, data=dgeObj$samples)
design
dgeObj <- estimateDisp(dgeObj, design)
fit <- glmFit(dgeObj, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
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
plot2 <-  plot2 + ggtitle("Precursor vs Cancer")


plot2 <- plot2 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot2


setwd("/mnt/DATA4/kevhu/LorenaData/")
pdf(file = "20180108.PreCancerConcordance.pdf",onefile = TRUE, width = 3)
plot2
dev.off()

setwd("/mnt/DATA4/kevhu/LorenaData/")
png(file = "20180108.PreCancerConcordance.png", width = 3.5, height = 5,units="in", res = 200)
plot2
dev.off()




###write.table
sig.genes <- lrt.tab[which(lrt.tab$color == "red"),]
write.table(sig.genes,"20171222.PreCancerWSquamConcordance.siggenes.txt", quote = FALSE, row.names = FALSE, sep = "\t")


###CTBNN1
dgeObj$samples$group <- CTNBB1
design <- model.matrix(~0+group, data=dgeObj$samples)
colnames(design) <- levels(dgeObj$samples$group)
dgeObj <- estimateDisp(dgeObj,design)
design
fit <- glmFit(dgeObj, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
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

plot3 <- plot3 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot3



pdf(file = "20180108CTNBB1.e2e.ConcordanceFilter.pdf",onefile = TRUE)
plot3
dev.off()


png(file = "20180108.CTNBB1Concordance.png", width = 3.5, height = 5,units="in", res = 200)
plot3
dev.off()

sig.genes <- lrt.tab[which(lrt.tab$color == "red"),]
write.table(sig.genes,"20180108CTNBB1.e2e.ConcordanceFilter.siggenes.txt", quote = FALSE, row.names = FALSE, sep = "\t")




###rewrite for adeno
dgeObj$samples$group <- AdenoSquam
design <- model.matrix(~0+group + CTNBB1, data=dgeObj$samples)
colnames(design) <- levels(dgeObj$samples$group)
design
dgeObj <- estimateDisp(dgeObj,design)
fit <- glmFit(dgeObj, design)
lrt <- glmLRT(fit, contrast=c(1,0,-1,0))
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


plot4 <- ggplot(data = lrt.tab,aes(logFC, -log10(QVal)))
plot4 <- plot4 +  geom_point(pch= 20, color = lrt.tab$color, size = 1.25) + geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.25)
#plot4 <-  plot4 + geom_text(data = lrt.tab[which(lrt.tab$color == "red"),], aes(label = Gene, vjust =-1.5), size = 3.0) + ggtitle("Adeno vs Squam")
plot4 <-  plot4 + ggtitle("Adeno vs Squam")

plot4 <- plot4 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot4




pdf(file = "20171216AdenoSquam.e2e.SpearmanFilter.pdf",onefile = TRUE)
plot4
dev.off()

png(file = "20180108.AdenoSquamConcordance.png", width = 3.5, height = 5,units="in", res = 200)
plot4
dev.off()


sig.genes <- lrt.tab[which(lrt.tab$color == "red"),]
write.table(sig.genes,"20171216AdenoSquamSpearman.siggenes.txt", quote = FALSE, row.names = FALSE, sep = "\t")


###rewrite for adeno + CTNBB1
#dgeObj$samples$group <- AdenoSquam

design <- model.matrix(~0 + CTNBB1 + Squam, data=dgeObj$samples)
design
dgeObj <- estimateDisp(dgeObj,design)
fit <- glmFit(dgeObj, design)
lrt <- glmLRT(fit, contrast=c(-1,1,0))
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


plot4 <- ggplot(data = lrt.tab,aes(logFC, -log10(QVal)))
plot4 <- plot4 +  geom_point(pch= 20, color = lrt.tab$color, size = 1.25) + geom_vline(xintercept = c(-1,1), linetype = "dashed", size = 0.25) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.25)
plot4 <-  plot4 + geom_text(data = lrt.tab[which(lrt.tab$color == "red"),], aes(label = Gene, vjust =-1.5), size = 3.0) + ggtitle("CTBNN1")
plot4 <- plot4 + theme_bw() + theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
plot4




pdf(file = "/mnt/DATA4/kevhu/LorenaData/20180106CTBNN1wSquamConcordance.pdf",onefile = TRUE)
plot4
dev.off()

sig.genes <- lrt.tab[which(lrt.tab$color == "red"),]
setwd("/mnt/DATA4/kevhu/LorenaData/")
write.table(sig.genes,"20180106AdenoSquamWCTNNB1Concordance.lowerThres.siggenes.txt", quote = FALSE, row.names = FALSE, sep = "\t")

###
###
###
###20180917 - making the graphs for revisions/edits
### try median centering, but with all genes, instead of reduced set ..... 

geneList.heatmap <- lrt.tab[which(lrt.tab$color == "red"),]$Gene
geneList.heatmap <- c(geneList.heatmap, "CXCL14","GAD1","NOTUM","NKD1","SLC16A1","DKK4","ZNRF3")
heatmapDf <- cpm(dgeObj)
colnames(heatmapDf) <- str_replace(colnames(heatmapDf), "L.*","")
heatmapDf <- log2(heatmapDf + 1)
heatmap.medians <- apply(heatmapDf, 2, median)
heatmapDf <- heatmapDf - heatmap.medians
#heatmapDf <- apply(heatmapDf, 1, scale)
heatmapDf.2 <- heatmapDf[which(rownames(heatmapDf) %in% geneList.heatmap),]


#rownames(heatmapDf.2) <- colnames(heatmapDf)
#heatmapDf.2 <- t(heatmapDf.2)
quantile.range <- quantile(heatmapDf.2, probs = seq(0, 1, 0.01))
colors.breaks <- seq(-8,8,16/1000)

caseAnno <- c("02","02","03","03","03","06","06","08","08","10","10","10")
lesionAnno <- c("CAH","ECsq","EC","EC","CAH","EC","EC","EC","EC","ECsq","CAH","CAH")
mutStat <- c("WT","Mut","WT","WT","Mut","Mut","WT","Mut","Mut","Mut","WT","WT")
anno.col <- list("Lesion"  = c("ECsq" = "#228B22","EC" = "#104E8B","CAH" = "#FF3030"),
                 "Case" = c("02" = "#8B4513","06" = "#32CD32","10" = "#8B4789","03" = "#8B6914","08" = "#FFEBCD"),
                 "CTNNB1_status" = c("WT" = "#CD2990","Mut" = "#4876FF"))

anno.Dat <- data.frame("Case" = factor(c("02","02","03","03","03","06","06","08","08","10","10","10")),
                       "Lesion" = factor(c("CAH","ECsq","EC","EC","CAH","EC","EC","EC","EC","ECsq","CAH","CAH")),
                       "CTNNB1" = factor(c("WT","Mut","WT","WT","Mut","Mut","WT","Mut","Mut","Mut","WT","WT")))

rownames(anno.Dat) <- colnames(heatmapDf.2)
heatmapDf.2 <- heatmapDf.2[,order(anno.Dat$CTNNB1, anno.Dat$Lesion, anno.Dat$Case, decreasing = TRUE)]
heatMapCol <- colorRampPalette(c("blue","white","red"))(1000)

pdf(file = "/mnt/DATA4/kevhu/LorenaData/20180918CTNNB1medCent.pdf", width = 7, height = 8, useDingbats = TRUE)
pheatmap(heatmapDf.2, color = heatMapCol,
         cluster_cols = FALSE, cluster_rows = TRUE,fontsize = 10,
         breaks = colors.breaks, annotation_colors = anno.col,
         annotation_col = anno.Dat, border_color = "black",
         cellwidth = 25, cellheight = 25, treeheight_row = 0,
         treeheight_col = 0)
dev.off()

### majority of genes about 46/49 found to be in the dataset are reduced to ~6-7 after filtering
### i.e: which(rownames(edgeRmat.pooled) %in% TCGA.diffGenes)
TCGA.diff <- read_xlsx(path = "/mnt/DATA4/kevhu/LorenaData/Final_LUSQvsLUADTCGA_SAT.xlsx")
TCGA.diffGenes <- TCGA.diff$LU_ADSYMBOL
heatmapDf.TCGA <- heatmapDf[which(rownames(heatmapDf) %in% TCGA.diffGenes),]
#a <- read_xlsx(path = "/mnt/DATA4/kevhu/LorenaData/Luadvssq.xlsx", sheet = "Sheet1")
#b <- a[,c(1,2,997:1001)]
#c <- b[9:525,]
#heatmapDf.TCGA <- heatmapDf[which(rownames(heatmapDf) %in% c$LU_ADSYMBOL),]
heatmapDf.TCGA <- heatmapDf.TCGA[,order(anno.Dat$CTNNB1, anno.Dat$Lesion, anno.Dat$Case, decreasing = TRUE)]

pdf(file = "/mnt/DATA4/kevhu/LorenaData/20180918TcgaComp.pdf", width = 7, height = 7, useDingbats = TRUE)
pheatmap(heatmapDf.TCGA, color = heatMapCol,
         cluster_cols = FALSE, cluster_rows = TRUE,fontsize = 10,
         breaks = colors.breaks, annotation_colors = anno.col,
         annotation_col = anno.Dat, border_color = "black",
         cellwidth = 25, cellheight = 25, treeheight_row = 0,
         treeheight_col = 0)
dev.off()

### Fischer on the two sets of DE genes - use adenoSquam comp for data: fisher comp is in both (ours + LU comp) vs only in LU
### ref in setting this up is ... https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
### p-value < 2.2e-16 with 95CI (0.0016 - 0.0059), OR: 0.0033

b.1 <- b[-which(is.na(b$FoldChangeUpSquam)),]
b.1$qVal <- p.adjust(b.1$TtestpvalueUpSquam, method = "BH")
b.DE <- b.1[which(abs(b.1$FoldChangeUpSquam) > 1 & b.1$qVal < 0.05),]
b.nonDE <- b.1[which(abs(b.1$FoldChangeUpSquam) < 1 | b.1$qVal > 0.05),]
length(which(lrt.tab$Gene[which(lrt.tab$color == "red")] %in% b.DE$LU_ADSYMBOL))
length(which(lrt.tab$Gene[which(lrt.tab$color == "blue")] %in% b.nonDE$LU_ADSYMBOL))

preVsCan.fisher <- matrix(c(11, 6177-11, 6439, 11954 - 6439), nrow = 2)
fisher.test(preVsCan.fisher, alternative = "two.sided")


### doing the clustering for genes found in clinical Sig B-catenin EEC paper - they scale by gene
#which(rownames(edgeRmat.pooled) %in% paperGenes)
paperGenes <- c("XBP1","ESR1","MYB","PGR","CDH2","PDGFA","WNT5A","WNT5B",
                "FOXM1","VEGFA","CCNB1","CDC20","GIMAP5","GIMAP7","LCK","STAT1")
#heatmapDf.paper <- cpm(dgeObj)
#colnames(heatmapDf.paper) <- str_replace(colnames(heatmapDf.paper), "L.*","")
#heatmapDf.paper <- apply(heatmapDf.paper, 1, scale)
#heatmapDf.paper <- t(heatmapDf.paper)
#colnames(heatmapDf.paper) <- colnames(heatmapDf)
#heatmapDf.paper2 <- heatmapDf.paper[match(paperGenes ,rownames(heatmapDf.paper)),]

heatmapDf.paper2 <- heatmapDf[match(paperGenes,rownames(heatmapDf)),]
heatmapDf.paper2 <- heatmapDf.paper2[-which(is.na(heatmapDf.paper2[,1])),]

#colors.breaks2 <- seq(-2,2,4/1000)
colors.breaks2 <- seq(-7,7,14/1000)

#anno.Dat.paper <- anno.Dat[order(anno.Dat$CTNNB1, anno.Dat$Lesion, anno.Dat$Case, decreasing = TRUE),]
heatmapDf.paper2 <- heatmapDf.paper2[,order(anno.Dat$CTNNB1, anno.Dat$Lesion, anno.Dat$Case, decreasing = TRUE)]

pdf(file = "/mnt/DATA4/kevhu/LorenaData/20180918cluster2Comp.pdf", width = 7, height = 7, useDingbats = TRUE)
pheatmap(heatmapDf.paper2, color = heatMapCol,
         cluster_cols = FALSE, cluster_rows = FALSE,fontsize = 10,
         clustering_distance_rows = "correlation",
         breaks = colors.breaks2, annotation_colors = anno.col,
         annotation_col = anno.Dat, border_color = "black",
         cellwidth = 25, cellheight = 25, treeheight_row = 0,
         treeheight_col = 0)
dev.off()

### grabbed data from CTNNB1 comp
#logFC   logCPM        LR      PValue      QVal  Gene color
#CDH2  -1.8587579 3.079579 8.9883674 0.002717036 0.3145403  CDH2  blue
#PDGFA -0.4593408 5.515410 0.8910582 0.345190704 0.9333833 PDGFA  blue
#WNT5A -0.8700021 5.496136 2.9006686 0.088542818 0.8083246 WNT5A  blue



####
#### Making the GSEA file - make rank file where gene is first column and second column is some metric
#### one metric that incorporate both FC (kind of) and p-value is the following: log10(P-value)/sign(logFC)
#### kind of disregards the magnitude of FC but we'll see ..... 

#gseaTxtTable <- data.frame(rownames(gseaTxtTable), rep("na", nrow(gseaTxtTable)), cpm(dgeObj),
#                           stringsAsFactors = FALSE)
#colnames(gseaTxtTable)[1:2] <- c("NAME","DESCRIPTION")
#gseaTxtTable <- gseaTxtTable[which(gseaTxtTable$NAME %in% lrt.tab$Gene[which(lrt.tab$color == "red")]),]
gseaTxtTable <- data.frame(lrt.tab$Gene, log10(lrt.tab$QVal)/sign(lrt.tab$logFC), stringsAsFactors = FALSE)

write.table(gseaTxtTable, file = "/mnt/DATA4/kevhu/LorenaData/20180919gseaAdenoSquam.rnk",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
