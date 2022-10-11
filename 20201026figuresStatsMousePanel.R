library(readxl)

newMouseBed <- read.table("/home/kevhu/data/bedFiles/IAD202296_167_Designed.bed", sep = "\t",
                          header = FALSE, stringsAsFactors = FALSE, skip = 1)

geneSheet <- read_xlsx("/home/kevhu/data/20200803mouseInput.xlsx", sheet = 1)

head(newMouseBed)


### getting distribution of targets
geneNames <- sub(x = newMouseBed$V6, pattern = ".*SUBMITTED_REGION=", replacement = "")
geneNames <- sub(x = geneNames, pattern = "exon.*", replacement = "")
geneNames <- sub(x = geneNames, pattern = "Exon.*", replacement = "")
geneNames <- sub(x = geneNames, pattern = "\\;.*", replacement = "")
geneNames <- sub(x = geneNames, pattern = ".*GENE_ID=", replacement = "")
head(geneNames)

table(geneNames)
unique(geneNames)

tsgs <- c("Apc", "Arid1a", "Atm", "Atrx", "Brca1", "Brca2","Cdh1", "Cdk12","Cdkn2a", "Crebbp",
          "Erbb2", "Errb3", "Fat1", "Fbxw7", "Foxa1", "Gata3", "Kmt2c", "Kmt2d",
          "Mtor", "Nf1", "Nf2", "Notch1", "Notch2", "Pbrm1","Pik3ca", "Pik3r1", "Ppp6c", "Ptch1", 
          "Pten", "Ptpn14", "Rb1", "Setd2", "Smad4", "Smo", "Trp53", "Vhl")


tsgGeneNames <- geneNames[which(geneNames %in% tsgs)]
geneNames2 <- geneNames[-which(geneNames %in% tsgs)]
cnGeneNames <- names(which(table(geneNames2) > 9))
geneNames3 <- geneNames2[which(geneNames2 %in% cnGeneNames)]
other <- geneNames2[-which(geneNames2 %in% cnGeneNames)]

table(tsgGeneNames)
table(geneNames3)
table(other)

slices <- c(2659, 1014, 490, 89)
lbls <- c("36 Genes fully tiled (TSGs)", "80 Oncogenes/Gene targets (> 10 Amps)",
          "490 SNP Genotyping", "Other")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels

colorScheme <- c("#F16A70", "#B1D877", "#8CDCDA", "#4D4D4D")
pie(slices,labels = lbls, col=colorScheme,
    main="Distribution of Amplicons on Mouse Panel (4250 Amplicons)")

### solid pie chart for makeup of panel
dev.off()
pdf("/home/kevhu/data/20201026pueChartMouse.pdf",useDingbats = TRUE)
pie(slices,labels = lbls, col=colorScheme,
    main="Distribution of Amplicons on Mouse Panel")
dev.off()
