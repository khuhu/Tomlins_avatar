
###making the necessary format for edgeR i.e count matrix
tableOfDat <- read.table(file = "/mnt/DATA4/kevhu/20171103LorenaRNAseq.csv", sep = ",", header = TRUE, skip = 7)
tableOfDat$Gene <- gsub("\\;.*","", tableOfDat$attributes)
tableOfDat$Gene <- gsub("GENE_ID=","", tableOfDat$Gene)
length(unique(tableOfDat$Gene))

#What I previously thought was correct, each amplicon is sequencing a differnt transript


edgeRmat <- tableOfDat[,6:29]
colnames(edgeRmat) <- colnames(tableOfDat)[6:29]
rownames(edgeRmat) <- tableOfDat$Gene

library(edgeR)

dgeObj <- DGEList(counts = edgeRmat)
dgeObj$samples

### cpm is counts/lib size * 1 milli . Need to find a cutoff for cpm which equals 5 counts for the smallest library
which(dgeObj$samples$lib.size == min(dgeObj$samples$lib.size))
rownames(dgeObj$samples)[1]

###grabbing cpm value of 
### maybe show distribution of of before and after each filter -> how it shows overral RNA sampled
which(edgeRmat$UT19 == 5)[1]
dge.cpm <- cpm(edgeRmat)
counts5filter <- cpm(edgeRmat)[30,1]
keep <- rowSums(cpm(dgeObj) > counts5filter) >= 2
dgeObj <- dgeObj[keep, , keep.lib.sizes = FALSE]

dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples$group <- preCan


preCan <- factor(c("Pre","Pre","Can","Can","Can","Can","Can","Can","Pre","Pre","Can","Can","Can","Can","Can","Can",
            "Can","Can","Can","Can","Can","Can","Pre","Pre"))

CTNBB1 <- factor(c("WT","WT","Mut","Mut","WT","WT","WT","WT","Mut","Mut","Mut","Mut","WT","WT","Mut","Mut","Mut",
            "Mut","Mut","Mut","WT","WT","WT","WT"))




###FisherExact

###Cancer vs preCan conparison

design <- model.matrix(~preCan)
dgeObj <- estimateDisp(dgeObj,design)
et <- exactTest(dgeObj, pair = c("Pre","Can"))
etPVals <- et$table$PValue
etPVals.ad <- p.adjust(etPVals, method = c("BH"))
et$table$QVal <- etPVals.ad
which(et$table$QVal < 0.05)
et.tab <- et$table
et.tab$Gene <- rownames(et.tab)

et.tab$color <- NULL
for(i in seq_along(et.tab$logFC)){
  if(abs(et.tab$logFC[i]) > 1 & et.tab$QVal[i] < 0.05){
    et.tab$color[i]  <- "red"
  }
  else{
    et.tab$color[i] <- "blue"
  }
}


plot1 <- ggplot(data = et.tab,aes(logFC, -log10(QVal)))
plot1 <- plot1 +  geom_point(pch= 20, color = et.tab$color) + geom_vline(xintercept = c(-1,1), linetype = "dashed") + geom_hline(yintercept = -log10(0.05), linetype = "dashed")
plot1 <-  plot1 + geom_text(data = et.tab[which(et.tab$color == "red"),], aes(label = Gene, vjust =-2), size = 2)

###GLM approach for pre vs Can
preCan <- factor(c("Pre","Pre","Can","Can","Can","Can","Can","Can","Pre","Pre","Can","Can","Can","Can","Can","Can",
                   "Can","Can","Can","Can","Can","Can","Pre","Pre"))
dgeObj$samples$group <- preCan
design <- model.matrix(~0+group, data=dgeObj$samples)
colnames(design) <- levels(dgeObj$samples$group)
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
plot2 <- plot2 +  geom_point(pch= 20, color = lrt.tab$color) + geom_vline(xintercept = c(-1,1), linetype = "dashed") + geom_hline(yintercept = -log10(0.05), linetype = "dashed")
plot2 <-  plot2 + geom_text(data = lrt.tab[which(lrt.tab$color == "red"),], aes(label = Gene, vjust =-1.5), size = 2) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Precursor vs Cancer")
plot2


#testing CTNBB1 with glm method

CTNBB1 <- factor(c("WT","WT","Mut","Mut","WT","WT","WT","WT","Mut","Mut","Mut","Mut","WT","WT","Mut","Mut","Mut",
            "Mut","Mut","Mut","WT","WT","WT","WT"))

dgeObj$samples$group <- CTNBB1
design <- model.matrix(~0+group, data=dgeObj$samples)
colnames(design) <- levels(dgeObj$samples$group)
dgeObj <- estimateDisp(dgeObj,design)

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
plot3 <- plot3 +  geom_point(pch= 20, color = lrt.tab$color) + geom_vline(xintercept = c(-1,1), linetype = "dashed") + geom_hline(yintercept = -log10(0.05), linetype = "dashed")
plot3 <-  plot3 + geom_text(data = lrt.tab[which(lrt.tab$color == "red"),], aes(label = Gene, vjust =-1.5), size = 2) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("CTNBB1 Wt vs Mut")
plot3


#glm for adeno vs squamos
AdenoSquam <- factor(c("Nein","Nein","Squa","Squa","Nein","Nein","Nein","Nein","Nein","Nein","Adeno","Adeno","Nein","Nein","Adeno","Adeno",
                       "Adeno","Adeno","Squa","Squa","Nein","Nein","Nein","Nein"))

dgeObj$samples$group <- AdenoSquam
design <- model.matrix(~0+group, data=dgeObj$samples)
colnames(design) <- levels(dgeObj$samples$group)
dgeObj <- estimateDisp(dgeObj,design)

fit <- glmFit(dgeObj, design)
lrt <- glmLRT(fit, contrast=c(1,0,-1))
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
plot4 <- plot4 +  geom_point(pch= 20, color = lrt.tab$color) + geom_vline(xintercept = c(-1,1), linetype = "dashed") + geom_hline(yintercept = -log10(0.05), linetype = "dashed")
plot4 <-  plot4 + geom_text(data = lrt.tab[which(lrt.tab$color == "red"),], aes(label = Gene, vjust =-1.5), size = 1.5) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Adeno vs Squuaaaaaa")
plot4


pdf(file = "AdenoVsSquaa.pdf",onefile = TRUE)
plot4
dev.off()


pdf(file = "CTNBB1.pdf",onefile = TRUE)
plot3
dev.off()

pdf(file = "PrecursorVCancer.pdf",onefile = TRUE)
plot2
dev.off()



testList <- lrt.tab[which(lrt.tab$color == "red"),]
rownames(testList) <- NULL

write.table(testList, file = "/mnt/DATA2/share/Users/Lorena/testList.csv",quote = FALSE, row.names = FALSE,sep = ",")
