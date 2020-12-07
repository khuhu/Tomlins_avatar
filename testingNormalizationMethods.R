###testing out the correcting on subset of high amplicon count -> ATM1

ATM1.test <- dfExpected[grep("ATM", dfExpected$Gene),]
for(i in seq_along(ATM1.test$Gene)){
  if(grepl("Pool1", ATM1.test$Gene[i])){
    ATM1.test$PoolNumber[i] <- "1" 
  }
  if(grepl("Pool2",  ATM1.test$Gene[i])){
    ATM1.test$PoolNumber[i] <- "2" 
  }
}


ATM1.test.Pool1 <- ATM1.test[which(ATM1.test$PoolNumber == 1),]
ATM1.test.Pool2 <- ATM1.test[which(ATM1.test$PoolNumber == 2),]

library(ggplot2)

a <- ggplot(data=ATM1.test,
       aes(ATM1.test$PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ,group=as.factor(ATM1.test$PoolNumber), color=as.factor(ATM1.test$PoolNumber))) + geom_histogram() 

b <- ggplot(data=ATM1.test,
       aes(ATM1.test$PR74.28prime.bc2.bat5.OCPv1i2_UM_28_repeat__MK5AH,group=as.factor(ATM1.test$PoolNumber), color=as.factor(ATM1.test$PoolNumber))) + geom_histogram() 

c <- ggplot(data=ATM1.test,
       aes(ATM1.test$PR66.39.bc55.bat2.OCPv1i2_UM_39__1SJZM,group=as.factor(ATM1.test$PoolNumber), color=as.factor(ATM1.test$PoolNumber))) + geom_histogram() 

d<- ggplot(data=ATM1.test,
       aes(ATM1.test$PR66.39prime.bc4.bat5.OCPv1i2_UM_39_repeat__MYBBI,group=as.factor(ATM1.test$PoolNumber), color=as.factor(ATM1.test$PoolNumber))) + geom_histogram() 

ggplot(data=ATM1.test,
       aes(ATM1.test$LU23.45.bc61.bat2.OCPv1i2_UM_45__ZL816,group=as.factor(ATM1.test$PoolNumber), color=as.factor(ATM1.test$PoolNumber))) + geom_density() 


ggplot(data=ATM1.test,
       aes(ATM1.test$LU23.45prime.bc8.bat5.OCPv1i2_UM_45_repeat2__6WZS0,group=as.factor(ATM1.test$PoolNumber), color=as.factor(ATM1.test$PoolNumber))) + geom_density() 


library(gridExtra)

grid.arrange(a,b)
grid.arrange(c,d)


for(x in allNames){
  dumVar <- ATM1.test[[x]]
  matrix(dumVar[which(),])
}


minMaxNorm <- function(x){
  (x - min(x))/(max(x) - min(x))
}

ATM1.test[["PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ"]] 
scale(ATM1.test.Pool1[["PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ"]], scale = FALSE)




med1 <- median(ATM1.test.Pool1[["PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ"]])
med2 <- median(ATM1.test.Pool2[["PR74.28prime.bc2.bat5.OCPv1i2_UM_28_repeat__MK5AH"]])
diffmed <- abs(med1 - med2)
mean1 <- mean(ATM1.test.Pool1[["PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ"]])
mean2 <- mean(ATM1.test.Pool2[["PR74.28prime.bc2.bat5.OCPv1i2_UM_28_repeat__MK5AH"]])
diffmean <- abs(mean1 - mean2)


a <- ggplot(data=ATM1.test.Pool1,
            aes(ATM1.test.Pool1$PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ)) + geom_histogram()
b <- ggplot(data=ATM1.test.Pool2,
            aes(ATM1.test.Pool2$PR74.28prime.bc2.bat5.OCPv1i2_UM_28_repeat__MK5AH)) + geom_histogram()

grid.arrange(a,b)

#higher one gets subtracted, while lower one gets added

ATM1.test.Pool1[["PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ"]]  <- ATM1.test.Pool1[["PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ"]] - diffmean/2 - diffmed/2
ATM1.test.Pool2[["PR74.28prime.bc2.bat5.OCPv1i2_UM_28_repeat__MK5AH"]] <- ATM1.test.Pool2[["PR74.28prime.bc2.bat5.OCPv1i2_UM_28_repeat__MK5AH"]] + diffmean/2 + diffmed/2


ATM1.test.Pool1 <- ATM1.test[which(ATM1.test$PoolNumber == 1),]
ATM1.test.Pool2 <- ATM1.test[which(ATM1.test$PoolNumber == 2),]
ATM1.test.Pool1[["PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ"]]  <- sqrt(ATM1.test.Pool1[["PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ"]])
ATM1.test.Pool2[["PR74.28prime.bc2.bat5.OCPv1i2_UM_28_repeat__MK5AH"]] <- sqrt(ATM1.test.Pool2[["PR74.28prime.bc2.bat5.OCPv1i2_UM_28_repeat__MK5AH"]])
ATM1.test.Pool1[["PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ"]]  <- log(ATM1.test.Pool1[["PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ"]])
ATM1.test.Pool2[["PR74.28prime.bc2.bat5.OCPv1i2_UM_28_repeat__MK5AH"]] <- log(ATM1.test.Pool2[["PR74.28prime.bc2.bat5.OCPv1i2_UM_28_repeat__MK5AH"]])


minMaxNorm(ATM1.test[["PR74.28.bc44.bat2.OCPv1i2_UM_28__MLQZQ"]])



###how about this - b/c it is a per gene difference in pools. I would calculate the difference of means between the two pools and add or subtract that difference from either the lower or higher amplicon count respectively




####Below is the start of the code I plan to implement

ATM1.test.Pool1 <- as.matrix(ATM1.test[which(ATM1.test$PoolNumber == 1),][,4:31])
ATM1.test.Pool2 <- as.matrix(ATM1.test[which(ATM1.test$PoolNumber == 2),][,4:31])


for(i in seq_along(ATM1.test.Pool1[,1])){
  assign(rownames(ATM1.test.Pool1)[i], ATM1.test.Pool1[,1][i])
}

dim(ATM1.test.Pool1)
dim(ATM1.test.Pool2)

dim(ATM1.test.Pool1)[1] - dim(ATM1.test.Pool2)[1]
dim(ATM1.test.Pool2)[1] - dim(ATM1.test.Pool1)[1]
dummyMat <- matrix(data = rep(NA, abs(dim(ATM1.test.Pool2)[1] - dim(ATM1.test.Pool1)[1]) * dim(ATM1.test.Pool1)[2]),
                   nrow = dim(ATM1.test.Pool2)[1] - dim(ATM1.test.Pool1)[1], ncol = dim(ATM1.test.Pool1)[2])

###
biocLite('preprocessCore')
library(preprocessCore)
normalize.quantiles()

x1length <- nrow(ATM1.test.Pool1)
x2length <- nrow(ATM1.test.Pool2)

if(x1length > x2length){
  ATM1.test.Pool2 <- rbind(ATM1.test.Pool2, dummyMat)
}
if(x2length > x1length){
  ATM1.test.Pool1 <- rbind(ATM1.test.Pool1, dummyMat)
}


combMatPool1.names <- rownames(ATM1.test.Pool1)
combMatPool2.names <- rownames(ATM1.test.Pool2)
combMat <- cbind(ATM1.test.Pool1[,1], ATM1.test.Pool2[,2])
combMat.normed <- normalize.quantiles(combMat)

pool1.normed <- data.frame(ampliconNames = combMatPool1.names, counts = as.numeric(combMat.normed[,1]))
pool2.normed <- data.frame(ampliconNames = combMatPool2.names, counts = as.numeric(combMat.normed[,2]))

pool1.normed <- na.omit(pool1.normed)
pool2.normed <- na.omit(pool2.normed)

library(reshape2)
library(stringi)
pool.combined <- rbind(pool1.normed, pool2.normed)
pool.combined <- pool.combined[stri_order(pool.combined$ampliconNames),]

ATM1.test[,4] <- pool.combined[,2]




###Actual implementation of above - same logic

for (x in allNames) {
  dfExpected[[x]] <- df$Weights * sum(df[[x]]);
}

###simple thing to note is that we just use a gc bed file with pool names attached to gene names
###after this normalization step we can remove it (the pool # - gsub)
for(i in seq_along(dfExpected$Gene)){
  if(grepl("Pool1", dfExpected$Gene[i])){
    dfExpected$PoolNumber[i] <- "1" 
  }
  if(grepl("Pool2",  dfExpected$Gene[i])){
    dfExpected$PoolNumber[i] <- "2" 
  }
  if(grepl("Pool3",  dfExpected$Gene[i])){
    dfExpected$PoolNumber[i] <- "3" 
  }
  if(grepl("Pool4",  dfExpected$Gene[i])){
    dfExpected$PoolNumber[i] <- "4" 
  }
}

geneList <- unique(dfExpected$Gene)
geneList <- unique(gsub("Pool[0-9]", "", geneList))

library(stringi)
library(preprocessCore)

for(x in allNames){
  for(i in seq_along(geneList)){
    geneMatches <- grep(geneList[i], dfExpected$Gene)
    genesubset <- dfExpected[geneMatches,]
    Pool1 <- genesubset[which(genesubset$PoolNumber == 1),]
    Pool2 <- genesubset[which(genesubset$PoolNumber == 2),]
    dummyMat <- matrix(data = rep(NA, abs(dim(Pool2)[1] - dim(Pool1)[1]) * dim(Pool1)[2]),
                       nrow = abs(dim(Pool2)[1] - dim(Pool1)[1]), ncol = dim(ATM1.test.Pool1)[2])
    Pool1length <- nrow(Pool1)
    Pool2length <- nrow(Pool2)
    if(x1length > x2length){
      Pool2 <- rbind(Pool2, dummyMat)
    }
    else if(x2length > x1length){
      Pool1 <- rbind(Pool1, dummyMat)
    }
    combMatPool1.names <- rownames(Pool1)
    combMatPool2.names <- rownames(Pool2)
    combMat <- cbind(Pool1[,1], Pool2[,2])
    combMat.normed <- normalize.quantiles(combMat)
    pool1.normed <- data.frame(ampliconNames = combMatPool1.names, counts = as.numeric(combMat.normed[,1]))
    pool2.normed <- data.frame(ampliconNames = combMatPool2.names, counts = as.numeric(combMat.normed[,2]))
    pool1.normed <- na.omit(pool1.normed)
    pool2.normed <- na.omit(pool2.normed)
    pool.combined <- rbind(pool1.normed, pool2.normed)
    pool.combined <- pool.combined[stri_order(pool.combined$ampliconNames),]
    dfExpected[[x]] <- pool.combined[,2]
  }
}



###Notes from lab meetin
###Load in common snps using IonReporter to force in SNP calls and then evaluate the other SNPs 
###for PanGU
