###Investigating dispersion estimates

y <- fit$counts
mu <- fit$fitted.values
phi <- fit$dispersion
v <- mu*(1+phi*mu)
resid.pearson <- (y-mu) / sqrt(v)
raw.resid <- y - mu

boxplot(resid.pearson)
plot(density(resid.pearson, na.rm = TRUE))
plot(density(raw.resid))



###comparing with residuals from arabidopsis data
library(NBPSeq)
data(arab)
head(arab)

Treat <- factor(substring(colnames(arab),1,4))
Treat <- relevel(Treat, ref="mock")
Time <- factor(substring(colnames(arab),5,5))

arabData <- DGEList(counts=arab, group=Treat)
keep2 <- rowSums(cpm(y)>2) >= 3
arabData <- arabData[keep2,,keep.lib.sizes = FALSE]
arabData <- calcNormFactors(arabData)
arabData$samples

design <- model.matrix(~0+Time+Treat)
rownames(design) <- colnames(arabData)
arabData <- estimateDisp(arabData, design, robust=TRUE)

arabData$common.dispersion
plotBCV(arabData)

arabFit <- glmQLFit(arabData, design, robust=TRUE)

arabObs <- arabFit$counts
arabExp <- arabFit$fitted.values
arabDisp  <- arabFit$dispersion
arabTrueVar <- arabExp * (1 + arabDisp*arabExp)
arabResid.pearson <- (arabObs - arabExp)/sqrt(arabTrueVar)
plot(density(arabResid.pearson))
arabRaw.Resid <- arabObs - arabExp


###getting another orthogonal set of data from Pickrell
library(tweeDEseqCountData)
data(pickrell1)
Counts <- exprs(pickrell1.eset)
Gender <- pickrell1.eset$gender
rm(pickrell1.eset)
data(annotEnsembl63)
annot <- annotEnsembl63[,c("Symbol","Chr")]
rm(annotEnsembl63)
pickrellData <- DGEList(counts=Counts, genes=annot[rownames(Counts),])
isexpr <- rowSums(cpm(pickrellData)>1) >= 20
hasannot <- rowSums(is.na(pickrellData$genes))==0
pickrellData <- pickrellData[isexpr & hasannot, , keep.lib.sizes=FALSE]
pickrellData <- calcNormFactors(pickrellData)
design <- model.matrix(~0+Gender)
pickrellData <- estimateDisp(pickrellData, design, robust=TRUE)
pickrellFit <- glmQLFit(pickrellData, design, robust=TRUE)
qlf <- glmQLFTest(pickrellFit)

pickrellObs <- pickrellFit$counts
pickrellExp <- pickrellFit$fitted.values
pickrellDisp <- pickrellFit$dispersion
pickrellSampVar <- pickrellExp * (1 + pickrellDisp*pickrellExp)
pickrelPearson <- (pickrellObs - pickrellExp)/sqrt(pickrellSampVar)

lower.quant <- quantile(pickrelPearson, 0.25)
upper.quant <- quantile(pickrelPearson, 0.75)

###looking at the UPenn stats page, anything with an aboslute value larger than 2 or 3 for pearson residual is a bad fit or values 
#https://onlinecourses.science.psu.edu/stat504/node/86

par(mfrow=c(2,2))
hist(resid.pearson)
hist(arabResid.pearson)
hist(raw.resid)
hist(arabRaw.Resid)

par(mfrow=c(2,2))
boxplot(resid.pearson)
boxplot(arabResid.pearson)
boxplot(pickrelPearson)
boxplot(pickrelPearson, ylim = c(lower.quant-1, upper.quant+1))

###below looking at percentage of fits which are badly fitted with the model
length(which(abs(resid.pearson) > 2))/length(resid.pearson)
length(which(abs(resid.pearson) > 3))/length(resid.pearson)
length(which(abs(arabResid.pearson) > 2))/length(arabResid.pearson)
length(which(abs(arabResid.pearson) > 3))/length(arabResid.pearson)
length(which(abs(pickrelPearson) > 2))/length(pickrelPearson)
length(which(abs(pickrelPearson) > 3))/length(pickrelPearson)

chisq.test(arabRaw.Resid, resid.pearson)
ks.test(arabRaw.Resid, resid.pearson)




