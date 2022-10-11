### these functions will be useful later on if I ever needed them, but for now using Sri's samples


### use function to remap names correctly

remapBCs <- function(x,y, ....){
  for (i in seq_along(y$CorrectBarcode)) {
    if(nchar(y$CorrectBarcode[i]) == 2){
      y$CorrectBarcode[i] <- paste0("IonXpress_0", y$CorrectBarcode[i])
    }
    else if(nchar(y$CorrectBarcode[i]) == 1){
      y$CorrectBarcode[i] <- paste0("IonXpress_00", y$CorrectBarcode[i])
    }
  }
  firstRow <- x[1, 4:ncol(x)]
  x[1, 4:ncol(x)] <- y$SampleName[match(firstRow, y$CorrectBarcode)]
  x[,which(is.na(x[1, 4:ncol(x)]))] <- firstRow[,which(is.na(x[1, 4:ncol(x)]))]
  return(x)
}

### 173_029. Doesn't make sense .... there are 39 samples listed in excel sheet + not listed what is swapped. i.e is barcode 17 25 and vice versa

run_170_023_remap <- remapBCs(run_170_023, trinhAnno_170_023)
run_173_029_remap <- remapBCs(run_173_029, trinhAnno_173_029)
run_174_033_remap <- remapBCs(run_174_033, trinhAnno_174_033)
run_175_035_remap <- remapBCs(run_175_035, trinhAnno_175_035)
run_176_037_remap <- remapBCs(run_176_037, trinhAnno_176_037)
run_177_039_remap <- remapBCs(run_177_039, trinhAnno_177_039)



### andi's samples - on target % > 60, 50% mapped reads being end-to-end, > 300,000 maped reads

rnaMatrixProcessing <- function(x){
  colnames(x) <- x[1,]
  totalReads <- as.numeric(x[2,4:ncol(x)])
  x <- x[-c(1,2),]
  x[,4:ncol(x)] <- lapply(x[,4:ncol(x)], as.numeric)
  listToRet <- list(x, totalReads)
  notEnoughReads <- which(totalReads < 300,000)
  badE2eRatio <- which(apply(x[,4:ncol(x)], 2, sum)/totalReads < 0.50)
  badSamples <- unique(c(notEnoughReads, badE2eRatio))
  print(paste("not enough reads:" , colnames(x)[(notEnoughReads + 4)]))
  print(paste("bad e2e ratio:", colnames(x)[(badE2eRatio + 4)]))
  if (length(badSamples) > 0) {
    x  <- x[,-(badSamples + 4)]
  }
  return(x)
}



### reading in Sri's table of 125 samples

urineData <- read.table("/mnt/DATA4/kevhu/urineRNA/20200225urine_counts_shiny_Sri.csv", sep = ",", header = FALSE, stringsAsFactors = FALSE)
urineData <- urineData[-c(3,4),]
namesForCol <- urineData[1, ]
namesForCol <- str_remove(namesForCol, "X")
urineData[1, ] <- namesForCol
urineData_processed <- rnaMatrixProcessing(urineData)
e2eRpmFactor <- apply(urineData_processed[, 4:ncol(urineData_processed)], 2, sum)/1e6
urineData_processed_libNorm <- cbind(urineData_processed[,1:2],
                                     urineData_processed[,4:ncol(urineData_processed)]/ e2eRpmFactor )



urineData <- read_xlsx("/mnt/DATA4/kevhu/urineRNA/All urine samps hdat.xlsx", sheet = 1)
colnames(urineData)[1] <- "Gene"
urineData_processed <- urineData[,1:180]
e2eRpmFactor <- apply(urineData_processed[, 3:ncol(urineData_processed)], 2, sum)/1e6
urineData_processed_libNorm <- cbind(urineData_processed[,1:2],
                                     urineData_processed[,3:ncol(urineData_processed)]/ e2eRpmFactor )

grep("102", colnames(urineData_processed_libNorm))
a <- ggplot(data = urineData_processed_libNorm) + 
  geom_point(aes(x = urineData_processed_libNorm[,3], y = urineData_processed_libNorm[,63])) + 
  xlab("102") + ylab("102.1")
grep("103", colnames(urineData_processed_libNorm))
b <- ggplot(data = urineData_processed_libNorm) + 
  geom_point(aes(x = urineData_processed_libNorm[,4], y = urineData_processed_libNorm[,64])) + 
  xlab("103") + ylab("103.1")
grep("105", colnames(urineData_processed_libNorm))
c <- ggplot(data = urineData_processed_libNorm) + 
  geom_point(aes(x = urineData_processed_libNorm[,5], y = urineData_processed_libNorm[,66])) + 
  xlab("105") + ylab("105.1")
grep("113", colnames(urineData_processed_libNorm))
d <- ggplot(data = urineData_processed_libNorm) + 
  geom_point(aes(x = urineData_processed_libNorm[,74], y = urineData_processed_libNorm[,121])) + 
  xlab("113") + ylab("113.1")
grep("122", colnames(urineData_processed_libNorm))
e <- ggplot(data = urineData_processed_libNorm) + 
  geom_point(aes(x = urineData_processed_libNorm[,9], y = urineData_processed_libNorm[,89])) + 
  xlab("122") + ylab("122.1")
grep("124", colnames(urineData_processed_libNorm))
f <- ggplot(data = urineData_processed_libNorm) + 
  geom_point(aes(x = urineData_processed_libNorm[,7], y = urineData_processed_libNorm[,91])) + 
  xlab("124") + ylab("124.1")

gridExtra::grid.arrange(a,b,c,d,e,f, ncol = 3)

cor(urineData_processed_libNorm[,3], urineData_processed_libNorm[,63])
cor(log2(urineData_processed_libNorm[,3] + 1), log2(urineData_processed_libNorm[,63] + 1))
cor(urineData_processed_libNorm[,5], urineData_processed_libNorm[,66])
cor(log2(urineData_processed_libNorm[,5] + 1), log2(urineData_processed_libNorm[,66] + 1))
cor(urineData_processed_libNorm[,74], urineData_processed_libNorm[,121])
cor(log2(urineData_processed_libNorm[,74] + 1), log2(urineData_processed_libNorm[,121] + 1))
cor(urineData_processed_libNorm[,9], urineData_processed_libNorm[,89])
cor(log2(urineData_processed_libNorm[,9] + 1), log2(urineData_processed_libNorm[,89] + 1))
cor(urineData_processed_libNorm[,7], urineData_processed_libNorm[,91])
cor(log2(urineData_processed_libNorm[,7] + 1), log2(urineData_processed_libNorm[,91] +1 ))


urineData_libNormed_melted <- melt(urineData_processed_libNorm[,2:ncol(urineData_processed_libNorm)])
urineData_libNormed_melted$value <- log2(urineData_libNormed_melted$value + 1)
pdf("/mnt/DATA4/kevhu/urineRNA/20200225boxplotSamples.pdf", useDingbats = FALSE, width = 10)
ggplot(data = urineData_libNormed_melted, aes(x = Gene, y = value)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) + ylab("log2(RPM + 1)") 
dev.off()


tmprss2fusions <- c("TMPRSS2-ERG.T1E2.COSF123", "TMPRSS2-ERG.T1E4.COSF125", "TMPRSS2-ERG.T2E4.COSF128", 
                    "TMPRSS2-ERG.T1E5.COSF26","TMPRSS2-ERG.T3E4.COSF130")

urineData_tmprss2 <- urineData_processed_libNorm[c(which(urineData_processed_libNorm$Gene %in% tmprss2fusions)),
                                                 3:ncol(urineData_processed_libNorm)]
rownames(urineData_tmprss2) <- tmprss2fusions
urineData_tmprss2 <- t(urineData_tmprss2)
urineData_tmprss2 <- apply(urineData_tmprss2, 2, as.numeric)

urineData_tmprss2 <- cbind(c(1:nrow(urineData_tmprss2)), urineData_tmprss2)
urineData_tmprss2 <- data.frame(urineData_tmprss2, stringsAsFactors = FALSE)

t1e2 <- urineData_tmprss2[order(urineData_tmprss2$TMPRSS2.ERG.T1E2.COSF123, decreasing = FALSE),]
t1e2$V1 <- factor(t1e2$V1, levels = t1e2$V1)
t1e4 <- urineData_tmprss2[order(urineData_tmprss2$TMPRSS2.ERG.T1E4.COSF125, decreasing = FALSE),]
t1e4$V1 <- factor(t1e4$V1, levels = t1e4$V1)
t2e4 <- urineData_tmprss2[order(urineData_tmprss2$TMPRSS2.ERG.T2E4.COSF128, decreasing = FALSE),]
t2e4$V1 <- factor(t2e4$V1, levels = t2e4$V1)
t1e5 <- urineData_tmprss2[order(urineData_tmprss2$TMPRSS2.ERG.T1E5.COSF26, decreasing = FALSE),]
t1e5$V1 <- factor(t1e4$V1, levels = t1e4$V1)
t3e4 <- urineData_tmprss2[order(urineData_tmprss2$TMPRSS2.ERG.T3E4.COSF130, decreasing = FALSE),]
t3e4$V1 <- factor(t3e4$V1, levels = t3e4$V1)

t1e2_plot <- ggplot(data = t1e2) + geom_point(aes(x = V1, y = log2(t1e2$TMPRSS2.ERG.T1E2.COSF123 + 1)), color = "cornflowerblue", alpha = 0.5) + 
  geom_hline(yintercept = 5) + 
  theme_bw() + theme(axis.line = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(), axis.title.y = element_text(size = 8),
                     plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("log2(TMPRSS2-T1E2 RPM + 1)") + ggtitle("T1E2") 

t1e4_plot <- ggplot(data = t1e4) + geom_point(data = t1e4, aes(x = V1, y = log2(t1e4$TMPRSS2.ERG.T1E4 + 1)), color = "red4", alpha = 0.2) + 
  geom_hline(yintercept = 5) +
  theme_bw() + theme(axis.line = element_blank(),axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.y = element_text(size = 8),
                     plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("log2(TMPRSS2-T1E4 RPM + 1)") + ggtitle("T1E4") 

t2e4_plot <- ggplot(data = t2e4) + geom_point(data = t2e4, aes(x = V1, y = log2(t2e4$TMPRSS2.ERG.T2E4.COSF128+ 1)), color = "gold4", alpha = 0.2) + 
  geom_hline(yintercept = 5) +
  theme_bw() + theme(axis.line = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(), axis.title.y = element_text(size = 8),
                     plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("log2(TMPRSS2-T2E4 RPM + 1)") + ggtitle("T2E4") 

t1e5_plot <- ggplot(data = t1e5) + geom_point(data = t1e5, aes(x = V1, y = log2(t1e5$TMPRSS2.ERG.T1E5.COSF26 + 1)), color = "mediumpurple2", alpha = 0.2) + 
  geom_hline(yintercept = 5) +
  theme_bw() + theme(axis.line = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(), axis.title.y = element_text(size = 8),
                     plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("log2(TMPRSS2-T1E5 RPM + 1)") + ggtitle("T1E5") 

t3e4_plot <- ggplot(data = t3e4) + geom_point(data = t3e4, aes(x = V1, y = log2(t3e4$TMPRSS2.ERG.T3E4.COSF130 + 1)), color = "olivedrab4", alpha = 0.2) + 
  geom_hline(yintercept = 5) +
  theme_bw() + theme(axis.line = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(), axis.title.y = element_text(size = 8),
                     plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("log2(TMPRSS2-T3E4 RPM + 1)") + ggtitle("T3E4") 


pdf("/mnt/DATA4/kevhu/urineRNA/20200303t2ergFusions.pdf", useDingbats = FALSE, width = 10)
grid.arrange(t1e2_plot, t1e4_plot, t2e4_plot,
             t1e5_plot, t3e4_plot, nrow = 2)
dev.off()



