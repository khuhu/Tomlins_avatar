#df <- read.table("/mnt/DATA4/kevhu/LorenaData/ZygosityTableWithClusteringOrder.csv", sep = ",", stringsAsFactors = FALSE)
library(xlsx)
library(stringr)
df <- read.xlsx("/mnt/DATA2/share/Users/Kevin/LorenaStuff/ID_IE_Final Zygosity Table for Kevin.xlsx", sheetIndex = 2, sheetName = "IE", stringsAsFactors=FALSE)
#df <- read.xlsx("/mnt/DATA2/share/Users/Kevin/LorenaStuff/ID_IE_Final Zygosity Table for Kevin.xlsx", sheetIndex = 1, sheetName = "ID_order1",stringsAsFactors=FALSE)

colnames(df)[1] <- "IDs"
#df[1,1] <- "IDs"
#df[1,1] <- "IDs"
#df$V5 <- NULL
#colnames(df) <- df[1,]
#df <- df[-1,]

### required for downstream graphs
df$TP53.1 <- df$TP53
df$CDKN2A.1 <- df$CDKN2A

df$TP53 <- gsub("Het", 1, df$TP53)
df$TP53 <- gsub("Homo", 2, df$TP53)
df$TP53 <- gsub(",","+",df$TP53)
df$CDKN2A <- gsub("Het", 1, df$CDKN2A)
df$CDKN2A <- gsub("Homo", 2, df$CDKN2A)
df$CDKN2A <- gsub(",","+",df$CDKN2A)

df$IDs <- str_replace_all(df$IDs, "[[:punct:]]", "")
df$IDs <- trimws(df$IDs)

df$TP53 <- sapply(df$TP53, function(x) sum(eval(parse(text = x))))
df$CDKN2A <- sapply(df$CDKN2A, function(x) sum(eval(parse(text = x))))

df <- df[order(df$Type),]

###IMPORTANT LINE OF CODE IF YOU WANT TO IMPOSE ORDERING OF X AXIS LABELS
df$IDs.ordered <- factor(df$IDs, levels = df$IDs)

#df2 <- df[order(df$Type),]
##trying base plot

#par(mar=c(3,3,0.5,3),cex=1, mfrow=c(2, 1))

#bplt <- barplot(df[[2]],col="red", ylab = "TP53", ylim = c(0,round(max(df[[2]]),max(df[[3]])) + 1))
#axis(1, bplt, df[[1]], lty=0, las=2, cex.axis=1.0)
#barplot(df[[3]], col="blue", ylab = "CDKN2A", axes=FALSE)
#axis(2, -seq(0, max(df[[2]], 1)))


###trying ggplot method
library(reshape2)
library(ggplot2)

df$TP53[is.na(df$TP53)] <- 0
df$CDKN2A[is.na(df$CDKN2A)] <- 0
#newDf <- NULL
#for(i in 1:nrow(df)){
#  for(j in 2:3){
#    a <- rep(df$IDs[i], df[i,j])
#    Genes <- rep(colnames(df)[j], df[i,j])
#    c <- cbind(a,Genes)
#    newDf <- rbind(newDf, c)
#  }
#}

#newDf <- data.frame(newDf, stringsAsFactors = FALSE)

#ggplot(data = newDf, aes(x = a, fill = Genes)) + 
#  geom_bar(data = subset(newDf, Genes == "TP53"), width = 0.6) + 
#  geom_bar(data = subset(newDf, Genes == "CDKN2A"), 
#           mapping = aes(y = - ..count.. ),
#           position = "identity", width = 0.6) +
#  scale_y_continuous(labels = abs) + 
#  xlab("Samples") + 
#  ylab("Mutational Burden") + 
#  theme(axis.text.x = element_text(angle = 90, hjust = -1)) + 
#  theme(axis.line = element_line(colour = "black"),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank()) +
#  geom_hline(yintercept = 0)


#View(mtcars)


###ggplot(df, aes(IDs,y= TP53)) + geom_col()
###
###
###

#customPal <- c("darkblue", "purple")

q <- ggplot(data = df) + 
  geom_col(data = df,aes(IDs.ordered, y=TP53), fill = "darkblue", width = 0.6) + 
  geom_col(data=df,aes(IDs.ordered, y=-CDKN2A), fill = "darkblue", width = 0.6) +
  #scale_colour_manual(values = c("TP53" = "darkblue", "CDKN2A" = "purple")) +
  scale_y_continuous(labels = abs) +
  xlab("Samples") + 
  ylab("Mutational Burden") + 
  #theme(axis.text.x = element_text(angle = 90, hjust = -1)) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle=90,margin = margin(t = 0, r = 0, b = 5, l = 0))) +
  geom_hline(yintercept = 0)
q


graphProps <- ggplot_build(q)[[1]][1:2]
xCoords <- unlist(graphProps[[1]]['x'])
xCoords2 <- (graphProps[[1]])
#graphProps.CDK <- (graphProps[[2]])


TP53.coords <- NULL
for(i in seq_along(xCoords)){
  if(df$TP53[i] < 2){
    next()
  }
  genoCounts <- length(unlist(strsplit(df$TP53.1[i],",")))
  iters <- genoCounts - 1
  if(iters < 1){
    iters <- 1
  }
  len <- floor(df$TP53[i]/genoCounts)
  for(j in 1:iters){
    print(paste("i:",i,"j:",j,"iters:",iters,"len",len))
    coord <- NULL
    coord <- c(xCoords2$xmin[i],xCoords2$xmax[i], len * j)
    TP53.coords <- rbind(TP53.coords, coord)
  }
}


CDKN2A.coords <- NULL
for(i in seq_along(xCoords)){
  if(df$CDKN2A[i] < 2){
    next()
  }
  genoCounts <- length(unlist(strsplit(df$CDKN2A.1[i],",")))
  iters <- genoCounts - 1
  if(iters < 1){
    iters <- 1
  }
  len <- floor(df$CDKN2A[i]/genoCounts)
  for(j in 1:iters){
    print(paste("i:",i,"j:",j,"iters:",iters,"len",len))
    coord <- NULL
    coord <- c(xCoords2$xmin[i],xCoords2$xmax[i], len * j)
    CDKN2A.coords <- rbind(CDKN2A.coords, coord)
  }
}

colnames(TP53.coords) <- c("xmin","xmax","y")
TP53.coords <- data.frame(TP53.coords, stringsAsFactors = FALSE)
colnames(CDKN2A.coords) <- c("xmin","xmax","y")
CDKN2A.coords <- data.frame(CDKN2A.coords, stringsAsFactors = FALSE)


pdf("/mnt/DATA4/kevhu/LorenaData/20180621IEMutPlot", width=7 * 1.25, height=5 * 1.25, useDingbats=FALSE)

q + 
  geom_segment(data = TP53.coords, aes(x = xmin, y = y, xend = xmax, yend= y), col = "red") +
  geom_segment(data = CDKN2A.coords, aes(x = xmin, y = -y, xend = xmax, yend= -y), col = "red") + 
  ggtitle("IE samples") +
  theme(plot.title = element_text(hjust=0.5))

dev.off()
