library(ggplot2)

mutationalLoad <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/mutation-load-updated.txt", sep = "\t",
                             header = TRUE)
mutationalLoad$totalLoad <- mutationalLoad$Silent.per.Mb + mutationalLoad$Non.silent.per.Mb


mouseBed <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed", sep = "\t",
                       stringsAsFactors = FALSE)

bbnmut <- read.table("/mnt/DATA5/tmp/kev/sigProfileExtractor/20221004concordMutMat_AllMut.txt",
                  sep = "\t", stringsAsFactors = FALSE, header = TRUE)
bbnMutLoad <- apply(bbnmut[,2:ncol(bbnmut)], 2, sum)/(sum(mouseBed$V3 - mouseBed$V2)/1e6)

tmpBbnMutLoadDf <- data.frame("Cohort" = rep("BBN", ncol(bbnmut) - 1), "Patient_ID" = colnames(bbnmut)[2:ncol(bbnmut)], "Tumor_Sample_ID" = colnames(bbnmut)[2:ncol(bbnmut)], 
                              "Silent.per.Mb" = rep(0, ncol(bbnmut) - 1), "Non.silent.per.Mb" =  rep(0, ncol(bbnmut) - 1), "totalLoad" = bbnMutLoad)

mutationalLoad <- rbind(mutationalLoad, tmpBbnMutLoadDf)

mutationalLoad$totalLoad[which(mutationalLoad$totalLoad > 100)] <- 100
### add mutation load of BBN mouse samples
### defintions of tmb high from foundationCdxOne study is <= 10 muts/mb
### from https://pubmed.ncbi.nlm.nih.gov/32919526/

customColor <- rep("black", length(unique(mutationalLoad$Cohort)))
customColor[2] <- "darkblue"

ggplot(mutationalLoad, aes(x = Cohort, totalLoad, color = Cohort)) + geom_boxplot() + ggtitle("TMB for 10000 TCGA samples") + 
  scale_y_continuous(breaks=seq(0, 100, 10)) + geom_hline(yintercept = 10, linetype = "dashed", color  = "red") + 
  scale_color_manual(values = customColor) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5),
        legend.position="none")


### going to analyze these sample by separating each into low, intermediate and high - subset the maf files for these
### then run through sigprofiler for each set
### note in the paper response was even higher in 

mutationalLoad_high <- mutationalLoad[which(mutationalLoad$totalLoad > 10), ]
mutationalLoad_low <- mutationalLoad[which(mutationalLoad$totalLoad < 10), ]

mutationalLoad_high_blca <- mutationalLoad_high[which(mutationalLoad_high$Cohort == "BLCA"), ]
mutationalLoad_high_LUAD <- mutationalLoad_high[which(mutationalLoad_high$Cohort == "LUAD"), ]
mutationalLoad_high_LUSC <- mutationalLoad_high[which(mutationalLoad_high$Cohort == "LUSC"), ]
mutationalLoad_high_SKCM <- mutationalLoad_high[which(mutationalLoad_high$Cohort == "SKCM"), ]
mutationalLoad_high_UCEC <- mutationalLoad_high[which(mutationalLoad_high$Cohort == "UCEC"), ]



### read in maf file and then separate it

mafFile <- read.table("/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/mc3.v0.2.8.PUBLIC.maf",
                      sep = "\t", stringsAsFactors = FALSE, header = TRUE)

mafFile_blca <- mafFile[grep(paste0(mutationalLoad_high_blca$Tumor_Sample_ID,collapse = "|"),
                             mafFile$Tumor_Sample_Barcode), ]
write.table(mafFile_blca, "/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/20230116blcaHighTmb.maf",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


mafFile_luad <- mafFile[grep(paste0(mutationalLoad_high_LUAD$Tumor_Sample_ID,collapse = "|"),
                             mafFile$Tumor_Sample_Barcode), ]
write.table(mafFile_luad, "/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/20230116luadHighTmb.maf",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


mafFile_lusc <- mafFile[grep(paste0(mutationalLoad_high_LUSC$Tumor_Sample_ID,collapse = "|"),
                             mafFile$Tumor_Sample_Barcode), ]
write.table(mafFile_lusc, "/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/20230116luscHighTmb.maf",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


mafFile_skcm <- mafFile[grep(paste0(mutationalLoad_high_SKCM$Tumor_Sample_ID,collapse = "|"),
                             mafFile$Tumor_Sample_Barcode), ]
write.table(mafFile_skcm, "/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/20230116skcmHighTmb.maf",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


mafFile_ucec <- mafFile[grep(paste0(mutationalLoad_high_UCEC$Tumor_Sample_ID,collapse = "|"),
                             mafFile$Tumor_Sample_Barcode), ]
write.table(mafFile_ucec, "/mnt/DATA5/tmp/kev/tmpDbs/cbioportal/20230116ucecHighTmb.maf",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


### use signalProfilerMatix generator for maf -> matrix -> run signalProfilerExtractor
### reference genome is hg19 i.e Grch37, loading tables, then picking highest signatures to show for all samples

cbioBlcaAct <- read.table("/mnt/DATA5/tmp/kev/sigProfilerAssignment/cbio_blca/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
                          sep = "\t", stringsAsFactors = FALSE, header = TRUE)
apply(cbioBlcaAct[, 2:ncol(cbioBlcaAct)], 2, sum)[order(apply(cbioBlcaAct[, 2:ncol(cbioBlcaAct)], 2, sum), decreasing = TRUE)]
cbioBlcaSigs <- c("SBS2", "SBS13", "SBS5", "SBS10b", "SBS1", "SBS10a", "SBS7a", "SBS15")



cbioLuadAct <- read.table("/mnt/DATA5/tmp/kev/sigProfilerAssignment/cbio_luad/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
                          sep = "\t", stringsAsFactors = FALSE, header = TRUE)
apply(cbioLuadAct[, 2:ncol(cbioLuadAct)], 2, sum)[order(apply(cbioLuadAct[, 2:ncol(cbioLuadAct)], 2, sum), decreasing = TRUE)]
cbioLuadSigs <- c("SBS4", "SBS5", "SBS13", "SBS2", "SBS1", "SBS24", "SBS39", "SBS45", "SBS6")



cbioLuscAct <- read.table("/mnt/DATA5/tmp/kev/sigProfilerAssignment/cbio_lusc/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
                          sep = "\t", stringsAsFactors = FALSE, header = TRUE)
apply(cbioLuscAct[, 2:ncol(cbioLuscAct)], 2, sum)[order(apply(cbioLuscAct[, 2:ncol(cbioLuscAct)], 2, sum), decreasing = TRUE)]
cbioLuscSigs <- c("SBS4", "SBS5", "SBS13", "SBS2", "SBS7b", "SBS1", "SBS39", "SBS15", "SBS45", "SBS7a", "SBS24", "SBS87", "SBS30")


cbioSkcmAct <- read.table("/mnt/DATA5/tmp/kev/sigProfilerAssignment/cbio_skcm/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
                          sep = "\t", stringsAsFactors = FALSE, header = TRUE)
apply(cbioSkcmAct[, 2:ncol(cbioSkcmAct)], 2, sum)[order(apply(cbioSkcmAct[, 2:ncol(cbioSkcmAct)], 2, sum), decreasing = TRUE)]
cbioSkcmSigs <- c("SBS7b", "SBS7a", "SBS45", "SBS10b", "SBS1", "SBS31", "SBS49", "SBS38", "SBS5")


cbioUcecAct <- read.table("/mnt/DATA5/tmp/kev/sigProfilerAssignment/cbio_ucec/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
                          sep = "\t", stringsAsFactors = FALSE, header = TRUE)
apply(cbioUcecAct[, 2:ncol(cbioUcecAct)], 2, sum)[order(apply(cbioUcecAct[, 2:ncol(cbioUcecAct)], 2, sum), decreasing = TRUE)]
cbioUcecSigs <- c("SBS10b", "SBS10a", "SBS6", "SBS15", "SBS1", "SBS5", "SBS14", "SBS20", "SBS54", "SBS28", "SBS21", "SBS26", "SBS42", "SBS7b", "SBS23", "SBS46", "SBS19")


allNoteWorthSigs <- unique(c(cbioBlcaSigs, cbioLuadSigs, cbioLuscSigs, cbioSkcmSigs, cbioUcecSigs))

### testing dotplot
library(ggdendro)
library(cowplot)
library(tidyverse)
library(ggtree)
library(patchwork) 
gene_cluster <- read.table("/mnt/DATA5/tmp/kev/misc/scRNA_dotplot_data.tsv.gz", sep = "\t",
                        stringsAsFactors = FALSE, header = TRUE)

dotplot <- gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')


tmp <-  gene_cluster %>% mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100)


plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')


