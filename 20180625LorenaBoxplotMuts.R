library(xlsx)
library(stringr)
library(ggplot2)
library(ggthemes)
library(scales)

#df <- read.xlsx("/mnt/DATA2/share/Users/Kevin/LorenaStuff/ID_IE_Final Zygosity Table for Kevin.xlsx", sheetIndex = 2, sheetName = "IE", stringsAsFactors=FALSE)
df <- read.xlsx("/mnt/DATA2/share/Users/Kevin/LorenaStuff/ID_IE_Final Zygosity Table for Kevin.xlsx", sheetIndex = 1, sheetName = "ID_order1",stringsAsFactors=FALSE)



colnames(df)[1] <- "IDs"
df$TP53[is.na(df$TP53)] <- "None"
df$CDKN2A[is.na(df$CDKN2A)] <- "None"


boxplot.table <- NULL
boxplot.table <- rbind(boxplot.table, c("TP53Homo","TP53Het","CDKN2AHomo","CDKN2AHet"))

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
  boxplot.table <- rbind(boxplot.table, finalRow)
}
                                                                                                                                                                                                                                                                                                                                                 
colnames(boxplot.table) <- boxplot.table[1,] 
boxplot.table <- boxplot.table[-1,] 
rownames(boxplot.table) <- df$IDs


typeLesion <- unique(df$Type)
AK <- which(df$Type == typeLesion[1])
IS <- which(df$Type == typeLesion[2])
In <- which(df$Type == typeLesion[3])

boxplot.table2 <- apply(boxplot.table, 2, as.numeric)

AK.tabe <- boxplot.table2[AK,]
IS.tabe <- boxplot.table2[IS,]
In.tabe <- boxplot.table2[In,]

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

finalMatrix.TP53 <- finalMatrix[,which(grepl("TP53",colnames(finalMatrix)))]
finalMatrix.TP53 <- finalMatrix.TP53[,order(colnames(finalMatrix.TP53))]
finalMatrix.TP53.melted <- melt(finalMatrix.TP53)
finalMatrix.TP53.melted <- cbind(finalMatrix.TP53.melted, grouping)

finalMatrix.CDKN2A <- finalMatrix[,which(grepl("CDKN2A",colnames(finalMatrix)))]
finalMatrix.CDKN2A <- finalMatrix.CDKN2A[,order(colnames(finalMatrix.CDKN2A))]


finalMatrix.ordered <- finalMatrix[,order(colnames(finalMatrix))]
finalMatrix.ordered.melted <- melt(finalMatrix.ordered)
finalMatrix.ordered.melted <- cbind(finalMatrix.ordered.melted,grouping)

###outlier shape set to NA to avoid plotting it twice
q <- ggplot(data = finalMatrix.ordered.melted, aes(x = variable, y = value)) + geom_boxplot(aes(color = grouping), outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.25, height = 0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color=guide_legend("Tumor Types")) + xlab("Mutation Types") +
  ylab("Number of Mutations") + scale_x_discrete(labels = c("","CDK2NA Het","","","CDK2NA Homo","","","TP53 Het","","","TP53 Homo","")) + 
  scale_y_continuous(limits = c(0,4), breaks = 1:4) +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.ticks.x = element_blank(), axis.title = element_text(size = 17, face = "bold"),legend.key.size = unit(1.5,"cm")) +
  geom_vline(xintercept = xCoords, color = "grey80", size = 1)

q

q_build <- ggplot_build(q)
#next part is for getting coordinate to build vertical grey lines

xCoords <- NULL
for(i in c(3,6,9)){
    xmax <- q_build[1]$data[[1]]$xmax[i]
    xCoords <- c(xCoords,xmax)
}

### added scalar so it's evenly between the two sets
xCoords <- xCoords + 0.15



write.table(finalMatrix.ordered, "/mnt/DATA4/kevhu/LorenaData/ID.FinalMatrix.txt",sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)






####### Pie charts

df <- read.xlsx("/mnt/DATA2/share/Users/Kevin/LorenaStuff/ID_IE_Final Zygosity Table for Kevin.xlsx", sheetIndex = 1, sheetName = "ID_order1",stringsAsFactors=FALSE)

colnames(df)[1] <- "IDs"
df$TP53[is.na(df$TP53)] <- "None"
df$CDKN2A[is.na(df$CDKN2A)] <- "None"

typeLesion <- unique(df$Type)

AK <- which(df$Type == typeLesion[1])
IS <- which(df$Type == typeLesion[2])
In <- which(df$Type == typeLesion[3])

AK.tabe <- boxplot.table2[AK,]
IS.tabe <- boxplot.table2[IS,]
In.tabe <- boxplot.table2[In,]

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
  
  
  dummyCol1 <- c(zeroPer,onePer,twoPer,threePer)
  dummyCol3 <- c(rep(listOfTypes[i],4))
  dummyCombinedCol <- cbind(dummyCol1, dummyCol2, dummyCol3)
  finalMatrix.pie <- rbind(finalMatrix.pie, dummyCombinedCol)
}

colnames(finalMatrix.pie) <- c("Freq","NumMut", "Type")
finalMatrix.pie <- data.frame(finalMatrix.pie, stringsAsFactors = FALSE)
finalMatrix.pie$Freq <- sapply(finalMatrix.pie$Freq, as.numeric)

finalMatrix.pie$NumMut <- factor(finalMatrix.pie$NumMut)
finalMatrix.pie$Type <- factor(finalMatrix.pie$Type)


ggplot(data = finalMatrix.pie, aes(x = factor(1) , y = Freq,fill = factor(NumMut))) +
  geom_bar(width = 1, stat = "identity") + facet_wrap(~Type) + coord_polar(theta = "y") + theme_bw() +
  theme(axis.text = element_blank()) 
  #theme(axis.text.x = element_text(size = 10), axis.ticks.x = element_blank(), legend.key.size = unit(1.5,"cm"))





