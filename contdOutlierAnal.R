###analysis of distributions of gene Est and amplicons

###below just looks to see whether or not the individually sampled amplicons (w/in a gene) follow the same distribution
geneList
gMatch <- grep(geneList[86], df$Gene)
testGene1 <- dfCopynumber[gMatch,allNames]
testGene1 <- t(testGene1)
par(mar=c(1,1,1,1))
par(mfrow=c(5,5))
for(i in 1:25){
  qqplot(testGene1[,sample(1:ncol(testGene1),1)],testGene1[,sample(1:ncol(testGene1),1)])
  abline(a=0,b=1)
}


###okay now I will see whether or not the distribution of geneEsts for that gene are concordant for amplicons
###conclusion: no they are not distributed in the same manner ... so this is why all the t-tests give siginificant p-values
geneEst2 <- geneEst[,-which(colnames(geneEst) %in% c("UT-23","UT-24","UT-57","UT-92","None","MD1510293","SU155330","AS-5"))]
geneEstTestSet <- unlist(geneEst2[grep(geneList[86], geneEst2$Gene), 2:ncol(geneEst2)])
par(mar=c(1,1,1,1))
par(mfrow=c(5,5))
for(i in 1:25){
  qqplot(geneEstTestSet,testGene1[,sample(1:ncol(testGene1),1)])
  abline(a=0,b=1)
}




### next i will need to reprogram what's below in order to see which comparison i want to make ..... I can just do one overall t-test like before.....
### or try windowed method .........


tot.pval <- NULL
rownames.tot <- NULL
for(i in seq_along(geneList)){
  geneMatches <- grep(geneList[i], df$Gene)
  if(length(geneMatches) == 0){
    next()
  }
  genesubset.cooks <-  dfCopynumber[geneMatches,allNames]
  genesubset.cooks <- genesubset.cooks[,-which(colnames(genesubset.cooks) %in% c("UT-23","UT-24","UT-57","UT-92","None","MD1510293","SU155330","AS-5"))]
  print(i)
  for(j in 1:nrow(genesubset.cooks)){
    dummyTMat <- NULL
    testCase <- genesubset.cooks[j,]
    refCase <- genesubset.cooks[-j,]
    dummyTMat <- cbind(refCase, testCase)
    pval <- t.test(testCase, mu = 0, alternative = "two.sided")$p.value
    ampId <- rownames(genesubset.cooks[j,])
    tot.pval <- c(tot.pval, pval)
    rownames.tot <- c(rownames.tot, ampId)
  }
}


names(tot.pval) <- rownames.tot
ad.tot.pval <- p.adjust(tot.pval, method = "BH")
length(which(ad.tot.pval < 0.0027))

hist(tot.pval)
hist(ad.tot.pval)



### will just do a modified-z for b/w amplicons - changed the way the function works

modifiedZ2 <- function(x, y, ...){
  x2 <- (x- median(y))/(1.486 * mad(y))
  return(x2)
}


modzscoreList <- NULL
modzscore <- NULL
outlierList.modz<- NULL

for(i in seq_along(geneList)){
  if(length(grep(geneList[i], df$Gene)) < 3){
    next()
  }
  geneMatches <- grep(geneList[i], df$Gene)
  genesubset.cooks <-  dfCopynumber[geneMatches,allNames]
  genesubset.cooks <- genesubset.cooks[,-which(colnames(genesubset.cooks) %in% c("UT-23","UT-24","UT-57","UT-92","None","MD1510293","SU155330","AS-5"))]
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

