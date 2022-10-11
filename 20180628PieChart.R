library(xlsx)
library(stringr)
library(ggplot2)
library(ggthemes)
library(scales)


###switch between the two
#df <- read.xlsx("/mnt/DATA2/share/Users/Kevin/LorenaStuff/ID_IE_Final Zygosity Table for Kevin.xlsx", sheetIndex = 1, sheetName = "ID_order1",stringsAsFactors=FALSE)
df <- read.xlsx("/mnt/DATA2/share/Users/Kevin/LorenaStuff/ID_IE_Final Zygosity Table for Kevin.xlsx", sheetIndex = 3, sheetName = "IE",stringsAsFactors=FALSE)



colnames(df)[1] <- "IDs"
df$TP53[is.na(df$TP53)] <- "None"
df$CDKN2A[is.na(df$CDKN2A)] <- "None"

typeLesion <- unique(df$Type)


dummyTable <- NULL
dummyTable <- rbind(dummyTable, c("TP53Homo","TP53Het","CDKN2AHomo","CDKN2AHet"))


for(i in 1:nrow(df)){
  P53.het.counter <- 0
  P53.homo.counter <- 0
  CDKN2A.het.counter <- 0
  CDKN2A.homo.counter <- 0
  
  P53.het.counter <- sum(grepl("Het",unlist(strsplit(df$TP53[i],","))))
  P53.homo.counter <- sum(grepl("Homo",unlist(strsplit(df$TP53[i],","))))
  CDKN2A.het.counter <- sum(grepl("Het",unlist(strsplit(df$CDKN2A[i],","))))
  CDKN2A.homo.counter <- sum(grepl("Homo",unlist(strsplit(df$CDKN2A[i],","))))
  
  finalRow <- c(P53.homo.counter, P53.het.counter, CDKN2A.homo.counter, CDKN2A.het.counter)
  dummyTable <- rbind(dummyTable, finalRow)
}

colnames(dummyTable) <- dummyTable[1,] 
dummyTable <- dummyTable[-1,] 
rownames(dummyTable) <- df$IDs


typeLesion <- unique(df$Type)
typeLesion <- typeLesion[order(typeLesion)]
AK <- which(df$Type == typeLesion[1])
IS <- which(df$Type == typeLesion[2])
In <- which(df$Type == typeLesion[3])

dummyTable2 <- apply(dummyTable, 2, as.numeric)

AK.tabe <- dummyTable2[AK,]
IS.tabe <- dummyTable2[IS,]
In.tabe <- dummyTable2[In,]

listOfTab <- c("AK.tabe", "IS.tabe", "In.tabe")

finalMatrix <- NULL
maxRows <- max(c(nrow(AK.tabe),nrow(IS.tabe), nrow(In.tabe)))
for(i in seq_along(listOfTab)){
  newRows <- maxRows - nrow(eval(as.name(listOfTab[i])))
  newMat <- matrix(rep(NA, newRows * 4), ncol = 4)
  
  finalMatrix <- cbind(finalMatrix,rbind(eval(as.name(listOfTab[i])), newMat))
  #rbind these NAs to the respective matrix
  #can combne the matrices afterwards
}

finalMatrix <- data.frame(finalMatrix)

### need to make grouping for melted data
grouping <- c(rep("Actinic Keratosis",maxRows),rep("In Situ",maxRows),rep("Invasive",maxRows))
grouping <- rep(grouping,4)



listOfTypes <- c("Actinic Keratosis CDKN2A Het","In Situ CDKN2A Het","Invasive CDKN2A Het",
                 "Actinic Keratosis CDKN2A Homo","In Situ CDKN2A Homo","Invasive CDKN2A Homo",
                 "Actinic Keratosis TP53 Het","In Situ TP53 Het","Invasive TP53 Het",
                 "Actinic Keratosis TP53 Homo","In Situ TP53 Homo","Invasive TP53 Homo")

dummyCol2 <- c("0 Mutations","1 Mutation ","2 Mutations","3 Mutations")


finalMatrix.ordered <- finalMatrix[,order(colnames(finalMatrix))]
colnames(finalMatrix.ordered) <- listOfTypes

finalMatrix.ordered.melted <- melt(finalMatrix.ordered)
finalMatrix.ordered.melted <- cbind(finalMatrix.ordered.melted,grouping)





finalMatrix.pie <- NULL
for(i in 1:ncol(finalMatrix)){
  dummyTable <- table(finalMatrix.ordered[,i])
  
  zeroPer <- dummyTable[names(dummyTable) == "0"]/sum(dummyTable)
  onePer <- dummyTable[names(dummyTable) == "1"]/sum(dummyTable)
  twoPer <- dummyTable[names(dummyTable) == "2"]/sum(dummyTable)
  threePer <- dummyTable[names(dummyTable) == "3"]/sum(dummyTable)
  
  if(length(zeroPer) == 0){
    zeroPer <- 0
  }
  if(length(onePer) == 0){
    onePer<- 0
  }
  if(length(twoPer) == 0){
    twoPer <- 0
  }
  if(length(threePer) == 0){
    threePer <- 0
  }
  
  dummyCol1 <- c(zeroPer, onePer, twoPer, threePer)
  dummyCol3 <- c(rep(listOfTypes[i],4))
  dummyCombinedCol <- cbind(dummyCol1, dummyCol2, dummyCol3)
  finalMatrix.pie <- rbind(finalMatrix.pie, dummyCombinedCol)
}

colnames(finalMatrix.pie) <- c("Freq","NumMut", "Type")
finalMatrix.pie <- data.frame(finalMatrix.pie, stringsAsFactors = FALSE)
finalMatrix.pie$Freq <- sapply(finalMatrix.pie$Freq, as.numeric)

finalMatrix.pie$NumMut <- factor(finalMatrix.pie$NumMut)
finalMatrix.pie$Type <- factor(finalMatrix.pie$Type)
finalMatrix.pie$Freq <- round(finalMatrix.pie$Freq, 3)

finalMatrix.pie$Freq.labels <- finalMatrix.pie$Freq * 100
finalMatrix.pie$Freq.labels <- as.character(finalMatrix.pie$Freq.labels)

for(i in seq_along(finalMatrix.pie$Freq.labels)){
  if(finalMatrix.pie$Freq.labels[i] == "0"){
    finalMatrix.pie$Freq.labels[i] <- ""
  }
  else{
    finalMatrix.pie$Freq.labels[i] <- paste(finalMatrix.pie$Freq.labels[i], "%")
  }
}

pdf(file = "/mnt/DATA4/kevhu/LorenaData/20180702IEZygosityPieChart.pdf", useDingbats = FALSE)
ggplot(data = finalMatrix.pie, aes(x = factor(1) , y = Freq,fill = NumMut)) +
  geom_bar(width = 1, stat = "identity") + facet_wrap(~Type)+ coord_polar(theta = "y")+ theme_bw() + ggtitle("IE Zygosity Plot") +
  theme(axis.text = element_blank(), panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5, size = 30),
        strip.background = element_blank(), axis.title.x = element_text(size = 30), strip.text = element_blank()) + 
  labs(x = "", y="") + geom_text(aes(label = Freq.labels), position = position_stack(vjust = 0.5)) +
  scale_x_discrete(name = "Invasive   In Situ     Actinic Keratosis") 
dev.off()



