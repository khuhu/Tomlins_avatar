###gene list created a bit prior to quantile normalization, but i will implement it right after the normalization
geneList
dfCopynumber <- df;
for (x in allNames) {
  dfCopynumber[[x]] <- log2(dfResidRatio[[x]] / dfCorrectionRatio[[x]] / sampleMedianGeneEst[[x]])
}


i <- 100;
cooksOutlier <- NULL
###it's funny i can probably just add this into the quantile normalization function
for(i in seq_along(geneList)){
  geneMatches <- grep(geneList[i], df$Gene)
  genesubset.cooks <-  df[geneMatches,allNames]
  genesubset.cooks <- genesubset.cooks[,-which(colnames(genesubset.cooks) %in% c("PR171","PR176","PR340","PR342","PR-310","PR-342"))]
  #if(nrow(genesubset.cooks == 0)){
  #  next()
  #}
  print(i)
  for(k in 1:nrow(genesubset.cooks)){
    idx5 <- NULL
    if(k == 1){
      strt <- 1
      end <- 5
      idx5 <- c(strt:end)
    }
    if(k == 2){
      strt <- 1
      end <- 5
      idx5 <- c(strt:end)
      cooksD <- cooks.distance(lm(form, data = nearest5))
    }
    if(k == nrow(genesubset.cooks)){
      strt <- nrow(genesubset.cooks) - 4
      end <- nrow(genesubset.cooks)
      idx5 <- c(strt:end)
    }
    if(k == nrow(genesubset.cooks) - 1){
      strt <- nrow(genesubset.cooks) - 4
      end <- nrow(genesubset.cooks)
      idx5 <- c(strt:end)
    }
    if(is.null(idx5)){
      strt <- k - 2
      end <- k + 2
      idx5 <- c(strt:end)
    }
    nearest5 <- genesubset.cooks[idx5,]
    nearest5 <- data.frame(t(nearest5), stringsAsFactors = FALSE)
    divList <- sample(1:ncol(nearest5), ncol(nearest5))
    var1 <- round(length(divList)/3)
    var2 <- 2* var1
    var3 <- length(divList)
    #below will be the right number wen I subset the idx
    var1.2 <- apply(nearest5[,divList[1:var1]], 1, mean)
    var2.2 <- apply(nearest5[,divList[(var1 + 1):var2]],1,mean)
    var3.2 <- apply(nearest5[,divList[(var2 +1):var3]],1,mean)
    nearest5$var1 <- var1.2
    nearest5$var2 <- var2.2
    nearest5$var3 <- var3.2
    #nearest5 <- data.frame((nearest5), stringsAsFactors = FALSE)
    #nearest5$med <- apply(nearest5, 1, median)
    
    ###this next part on the forumla I need to think of ... do I make another variable as the response i.e a median of the 3 agglomerated ones?
    nearest5$meanVar <- apply(nearest5[,1:56],1,mean)
    form <- as.formula(paste("meanVar","~", paste(c("var1","var2","var3"), collapse = "+")))
    #plot(cooks.distance(lm(form, data = nearest5)))
    cooksD <- cooks.distance(lm(form, data = nearest5))
    cutoff <- 4/(length(nearest5) - 4)
    print(which(cooksD > cutoff))
  }
}





###testing crap
nearest5 <- genesubset.cooks[1:10,]
nearest5 <- data.frame((nearest5), stringsAsFactors = FALSE)
#nearest5 <- data.frame((nearest5), stringsAsFactors = FALSE)
nearest5$med <- apply(nearest5, 1, median)
form <- as.formula(paste(colnames(nearest5)[ncol(nearest5)],"~", paste(colnames(nearest5)[1:8], collapse = "+")))
form
fit <- lm(form, data = nearest5)
cooksD <- cooks.distance(lm(form, data = nearest5))
summary(fit)
plot(cooksD)
cooksD



### using a modified z-score method i.e instead of mean/sd you use median/(MAD * conversion factor) the conversion factor approximates the scale factor b/w sd and MAD
###pruning bad samples 

###OCPv3
df2 <- df[,-which(colnames(df) %in% c("PR171","PR176","PR340","PR342","PR-310","PR-342"))]

###OCP1c
df2 <- df[,-which(colnames(df) %in% c("UT-23","UT-24","UT-57","UT-92","None","MD1510293","SU155330","AS-5"))]

geneList <- unique(df$Gene)
geneList


#function for calculating modified z score
###Below is current implementation for this method
### regular modz function wrong b/c it doesn't include the scale factor for MAD to SD
modifiedZ <- function(x, ...){
  medians <- apply(x, 1, median)
  mads <- apply(x, 1, mad)
  x2 <- sweep(x, 1, medians, "-")
  x3 <- sweep(x2, 1, mads, "/")
}


modifiedZ2 <- function(x, ...){
  medians <- median(x)
  mads <- mad(x)
  x2 <- (x-medians)/(1.486 * mads)
  return(x2)
}

avg <- function(x,...){
  means <- apply(x, 1, mean)
}


alteredZ <- function(x,y){
  z <- (x - mean(y))/sd(y)
  return(z)
}

modzscoreList <- NULL
modzscore <- NULL
outlierList.modz<- NULL

zscoresList.alt <- NULL
outlierList.alt <-NULL
alt.scores <- NULL


zscoresList <- NULL
outlierList <- NULL

for(i in seq_along(geneList)){
  if(length(grep(geneList[i], df$Gene)) < 5){
    next()
  }
  #print(i)
  geneMatches <- grep(geneList[i], df$Gene)
  genesubset.cooks <-  dfCopynumber[geneMatches,allNames]
  genesubset.cooks <- genesubset.cooks[,-which(colnames(genesubset.cooks) %in% c("UT-23","UT-24","UT-57","UT-92","None","MD1510293","SU155330","AS-5"))]
  #make below into a function. would save space - i.e the sampling from IQR
  geneEst2 <- geneEst[,-which(colnames(geneEst) %in% c("UT-23","UT-24","UT-57","UT-92","None","MD1510293","SU155330","AS-5"))]
  geneEstIQR <- unlist(geneEst2[grep(geneList[i], geneEst2$Gene), 2:ncol(geneEst2)])
  #geneEstIQR2 <- subset(geneEstIQR, geneEstIQR >= quantile(geneEstIQR, 0.25) & geneEstIQR <= quantile(geneEstIQR, 0.75))
  #geneEstIQR3 <- sample(geneEstIQR2, length(geneEstIQR2))
  geneEstIQR3 <- sample(geneEstIQR, 60)
  print(i)
  for(k in 1:nrow(genesubset.cooks)){
    print(k)
    idx5 <- NULL
    if(k == 1){
      strt <- 1
      end <- 3
      idx5 <- c(strt:end)
    }
    if(k == 2){
      strt <- 1
      end <- 3
      idx5 <- c(strt:end)
    }
    if(k == nrow(genesubset.cooks)){
      strt <- nrow(genesubset.cooks) - 2
      end <- nrow(genesubset.cooks)
      idx5 <- c(strt:end)
    }
    if(k == nrow(genesubset.cooks) - 1){
      strt <- nrow(genesubset.cooks) - 3
      end <- nrow(genesubset.cooks) - 1
      idx5 <- c(strt:end)
    }
    if(is.null(idx5)){
      strt <- k - 1
      end <- k + 1
      idx5 <- c(strt:end)
    }
    nearest5 <- genesubset.cooks[idx5,]
    nearest5 <- data.frame(nearest5, stringsAsFactors = FALSE)
    nearest5.means <- apply(nearest5, 1, mean)
    nearest5.meds <- apply(nearest5, 1, median)
    dummyDist <- c(geneEstIQR3, nearest5.means)
    dummyDist2 <- c(geneEstIQR3, nearest5.meds)
    modzscore <- modifiedZ2(dummyDist2)
    alt.scores <- alteredZ(nearest5.means, geneEstIQR3)
    zscores <- scale(dummyDist)
    zscoresList <- c(zscoresList, zscores)
    zscoresList.alt <- c(zscoresList.alt, alt.scores)
    modzscoreList <- c(modzscoreList, modzscore)
    outlierList <- c(outlierList,rownames(zscores)[which(abs(zscores) > 3)])
    outlierList.alt <- c(outlierList.alt,names(alt.scores)[which(abs(alt.scores) > 3)])
    #outlierListMeds <- c(outlierListMeds, rownames(zscores)[which(abs(modzscore) > 3)])
    outlierList.modz <- c(outlierList.modz, names(modzscore)[which(abs(modzscore) > 3)])
  }
}

for(i in seq_along(geneList)){
  print(length(grep(geneList[i], df$Gene)))
}

outlierList9 <- outlierList
outlierListAll <- outlierList

###notes on above, everything becomes an outlier if I calculate z score where mean and sd is solely based on the geneEst





###looking at just t-test 
###idea would be not using a window, but sample ~15 from each (might increase when I add more samples). calculate t-test and then correct for multiplicity 

tot.pval <- NULL
rownames.tot <- NULL
for(i in seq_along(geneList)){
  geneMatches <- grep(geneList[i], df$Gene)
  if(length(geneMatches) == 0){
    next()
  }
  
  genesubset.cooks <-  dfCopynumber[geneMatches,allNames]
  genesubset.cooks <- genesubset.cooks[,-which(colnames(genesubset.cooks) %in% c("UT-23","UT-24","UT-57","UT-92","None","MD1510293","SU155330","AS-5"))]
  geneEst2 <- geneEst[,-which(colnames(geneEst) %in% c("UT-23","UT-24","UT-57","UT-92","None","MD1510293","SU155330","AS-5"))]
  geneEst2.1 <- unlist(geneEst2[grep(geneList[i], geneEst2$Gene), 2:ncol(geneEst2)])
  print(i)
  for(j in 1:nrow(genesubset.cooks)){
    dummyTMat <- NULL
    testCase <- sample(genesubset.cooks[j,], 80)
    refCase <- sample(geneEst2.1, 80)
    dummyTMat <- cbind(refCase, testCase)
    pval <- t.test(refCase,testCase, mu = 0, alternative = "two.sided", conf.level = 0.9973)$p.value
    ampId <- rownames(genesubset.cooks[j,])
    tot.pval <- c(tot.pval, pval)
    rownames.tot <- c(rownames.tot, ampId)
  }
}


wilcox.test(refCase, t(testCase), paired = FALSE)




names(tot.pval) <- rownames.tot
ad.tot.pval <- p.adjust(tot.pval, method = "BH")
length(which(ad.tot.pval < 0.0027))

hist(tot.pval)
hist(ad.tot.pval)


totMeans <- NULL
for(i in 1:10000){
 a <- mean(t(sample(genesubset.cooks[2,], 80)))
 totMeans <- c(totMeans, a)
}
hist(totMeans)

