### 20230516: changed directories after server crash/wipe


library(readxl)
library(ggplot2)
library(pROC)
library(reshape2)
library(gridExtra)
library(caret)
library(stringr)
library(Rcpp)

### loading pretrained models
load("/mnt/DATA6/kevin_recovery/urineRNA/20180808urineLogGeneModel.Robj")
load("/mnt/DATA6/kevin_recovery/urineRNA/20180808urineMipsModel.Robj")
load("/mnt/DATA6/kevin_recovery/urineRNA/20180808urinePsaModel.Robj")

regModel.coeffs <- c(-4.466250739,0.357331671,0.577094327,0.194125163)
highRisk.coeffs <- c(-5.85877915,0.590375797,0.553164487,0.143707409)

### set number of genes and loaded data with set groups
load("/mnt/DATA6/kevin_recovery/urineRNA/20190620urineG1Subsets.Rdata")
load("/mnt/DATA6/kevin_recovery/urineRNA/20190620urineG2Subsets.Rdata")
load("/mnt/DATA6/kevin_recovery/urineRNA/20180807randoForestGenes.Robj")

### making a dataframe for Andi - probability for it being high-grade cancer



a <- read_xlsx("/mnt/DATA6/kevin_recovery/urineRNA/AC18.1_byKLK3_SAT.JStoKevin.xlsx")
a <- a[1:(nrow(a) -1),]
a.1 <- log2(a[,genes30]+1)
twoClass <- as.character(a$`Pathology Intermediate/Low vs. VeryLow/Benign`)
twoClass.1 <- str_replace_all(twoClass, "1", "Cancer")
twoClass.1 <- factor(str_replace_all(twoClass.1,"0","NoCancer"))

### above is for bening + very low vs low and inter; below is first 3 vs inter

#twoClass <- as.character(a$`Pathology Risk Group(1=Benign, 2=Very , 3=, 4=iate)`)
#twoClass.1 <- str_replace_all(twoClass, "4", "Cancer")
#twoClass.1 <- str_replace_all(twoClass.1,"3","NoCancer")
#twoClass.1 <- str_replace_all(twoClass.1,"2","NoCancer")
#twoClass.1 <- str_replace_all(twoClass.1,"1","NoCancer")
#twoClass.1 <- factor(twoClass.1)

logMod.prob <- predict(fiveFoldModel.log, newdata = a.1)
#confusionMatrix(data = logMod.prob, twoClass.1, positive = "Cancer")
probROC.logMod <- predict(fiveFoldModel.log, type = c("prob"),newdata = a.1)
rocGraph.logMod <- roc(twoClass.1 ~ probROC.logMod$`Cancer`)
plot(rocGraph.logMod)
rocGraph.logMod

ci.auc(rocGraph.logMod, conf.level = 0.95, method = "delong")

tableOfCIs <- data.frame("Model" = "AS log model","AUC" = 66,"Lower bound" = 50, "Upper bound" = 83, stringsAsFactors = FALSE)
probDfTestingAS <- data.frame(probROC.logMod$Cancer, stringsAsFactors = FALSE)
rownames(probDfTestingAS) <- a$SampleID



PSA <- a$`Serum PSA`
a.2 <- data.frame(twoClass.1,PSA)
a.2$ones <- rep(1,nrow(a.2))

model_weights <- ifelse(a.2$twoClass.1 == "Cancer",
                        (1/table(a.2$twoClass.1)[1]) * 0.5,
                        (1/table(a.2$twoClass.1)[2]) * 0.5)
train_control.combined <- trainControl(method = "cv", number = 5, savePredictions = TRUE,
                                       classProbs = TRUE, summaryFunction = twoClassSummary)

set.seed(2222)
fiveFoldModel.psa <- train(twoClass.1 ~., data = a.2, trControl=train_control.combined, method="glmnet",
                           metric = "ROC", family = "binomial", weights = model_weights)



fiveFoldModel.psa
psaPredAs <- predict(fiveFoldModel.psa, newdata = a.2)

probROC.psa <- predict(fiveFoldModel.psa, type = c("prob"))
rocGraph.psa <- roc(twoClass.1 ~ probROC.psa$`Cancer`)
plot(rocGraph.psa)
rocGraph.psa$auc


#confusionMatrix(data = psaPredAs, twoClass.1, positive = "Cancer")
ci.auc(rocGraph.psa, conf.level = 0.95, method = "delong")

tableOfCIs <- rbind(tableOfCIs, c("AS PSA model",53, 35, 71))
probDfTestingAS <- cbind(probDfTestingAS, probROC.psa$Cancer)





a.3 <- data.frame(twoClass.1, PSA, log2(a[,c("TMPRSS2-ERG.T1E4.COSF125","PCA3.E2E3")] + 1))
colnames(a.3)[3:4] <- c("TMPRSS2:ERG","PCA3")

probROC.Mips <- predict(fiveFoldModel.mips, type = c("prob"),newdata = a.3)
rocGraph.Mips <- roc(twoClass.1 ~ probROC.Mips$`Cancer`)
plot(rocGraph.Mips)
rocGraph.Mips

probDfTestingAS <- cbind(probDfTestingAS, probROC.Mips$Cancer)

predAsRetrained <- predict(fiveFoldModel.mips, newdata = a.3)
#confusionMatrix(predAsRetrained, twoClass.1)
ci.auc(rocGraph.Mips, conf.level = 0.95, method = "delong")
tableOfCIs <- rbind(tableOfCIs, c("AS Mips model",43, 27, 65))


### should be PSA, PCA3, and T2ERG in that order thats why I switched it when calculating the original scores

### looking at normal mips model
a.4 <- a.3
a.4 <- a.4[,c(1,2,4,3)]
a.4$twoClass.1 <- relevel(a.4$twoClass.1, "NoCancer")
tmaMips.model <- glm(formula = twoClass.1 ~ .,data = a.4, family = binomial(link = "logit"))
tmaMips.model.reg <- tmaMips.model
tmaMips.model.highRisk <- tmaMips.model
tmaMips.model.reg$coefficients <- regModel.coeffs
tmaMips.model.highRisk$coefficients <- highRisk.coeffs

fullHigh.predict <-  predict.glm(object = tmaMips.model.highRisk, type=c("response"), newdata = a.4)
c <- roc(a.4$twoClass.1 ~ fullHigh.predict)
plot(c)
fullReg.predict <-  predict.glm(object = tmaMips.model.reg, type=c("response"), newdata = a.4)
d <- roc(a.4$twoClass.1 ~ fullReg.predict)
plot(d)


probDfTestingAS <- cbind(probDfTestingAS, fullHigh.predict)
probDfTestingAS <- cbind(rownames(probDfTestingAS), probDfTestingAS)
colnames(probDfTestingAS) <- c("Samples","29_Transcript_Model", "PSA", "Retrained_Mips", "High-grade_Mips")
rownames(probDfTestingAS) <- NULL


### I think I need to do self cutoffs for these things
#mipsHighRisk <- predict.glm(object = tmaMips.model.highRisk, type = "response",newdata = a.4)
#mipsHighRisk[which(mipsHighRisk > 0.5)] <- "Cancer"
#mipsHighRisk <- as.factor(mipsHighRisk)
#confusionMatrix(mipsHighRisk, twoClass.1)
ci.auc(c, conf.level = 0.95, method = "delong")
tableOfCIs <- rbind(tableOfCIs, c("AS high-risk model", 56, 37, 75))

#write.table(probDfTestingAS, "/mnt/DATA6/kevin_recovery/urineRNA/20190801asProbTable.txt", row.names = FALSE, col.names = TRUE,
#            quote = FALSE, sep = "\t")


#### making a basic graph for the AS data
#pdf(file = "/mnt/DATA6/kevin_recovery/urineRNA/20190816AScombinedGraphs3GroupSplit.pdf", onefile = TRUE, useDingbats = TRUE, pointsize = 18, width = 10, height = 10)
plot(rocGraph.psa, col = alpha("blue", alpha = 0.5), main = "ROCs: Active Surveillance n = 45", ylab = "", xlab="", lty = 1, asp = FALSE)
par(new = TRUE)
plot.roc(c, col = alpha("orange", alpha = 0.5), asp = FALSE)
par(new = TRUE)
plot.roc(rocGraph.Mips, main=NULL, asp = FALSE, col = alpha("brown4",alpha = 0.5))
par(new=TRUE)
plot.roc(rocGraph.logMod, col = alpha("red", alpha = 0.5), asp = FALSE)
par(new = TRUE)
legend("bottomright", legend = c("PSA - AUC: 0.47",
                                 "High grade risk Mips - AUC: 0.57",
                                 "Retrained Seq Mips - AUC: 0.40",
                                 "29 transcript model - AUC: 0.65"),
       col = c("blue","orange","brown4","red"), lty=c(rep(1,3)), cex = .70)
dev.off()






### this is for the original training data and testing data
###
###

### loading data, list of genes and then makign the partitions for data
load("/mnt/DATA6/kevin_recovery/urineRNA/20180808urinePsaModel.Robj")

a <- read_xlsx("/mnt/DATA6/kevin_recovery/urineRNA/AC17.2ExtDes2_mastermatrix_allsamples_byKLK3.xlsx")
b <- read_xlsx("/mnt/DATA6/kevin_recovery/urineRNA/20189_07_26_AC17.2ExtDes2_mastermatrix.xlsx")
a <- a[-which(a$`Grade group` == 2),]
a.1 <- a[1:109,3:86]
rownames(a.1) <- as.vector(a$Contig_ID[1:109])
ids <- a$`Grade group bin`[1:109]
b <- b[match(rownames(a.1), b$Contig_ID),]
PSA <- b$LAB_RESNUM

combined.subset <- c(g1.subset,g2.subset)
held.out <- a.1[combined.subset,]
a.2 <- a.1[-combined.subset,]
ids.2 <- factor(ids[-combined.subset])
ids.3 <- str_replace_all(ids.2, "GG2-5", "Cancer")
ids.3 <- str_replace_all(ids.3, "Benign/GG1", "NoCancer")
a.2.log <- log2(a.1[-combined.subset,]+1)
a.2.log.30 <- cbind(ids.3,a.2.log[,genes30])


###logged data
### note the initial ordering does seem to matter, but may need ordering for case weights

#a.7 <- a.2.log.30
#a.7 <- a.7[order(a.7$ids.3),]
ids.4 <- ids.3
ids.3 <- ids.3[order(ids.3)]


### subsetting for the training data
ids.train <- as.factor(ids.3)
train.dat <- a.2.log.30

### subsetting for the testing data
ids.held.out <- ids[combined.subset]
ids.held.out <- str_replace_all(ids.held.out, "GG2-5", "Cancer")
ids.held.out <- str_replace_all(ids.held.out, "Benign/GG1", "NoCancer")
ids.held.out <- as.factor(ids.held.out)
held.out.log <- log2(held.out[,genes30] + 1)
prob.log <- predict(fiveFoldModel.log, newdata = held.out.log)
probROC.log <- predict(fiveFoldModel.log, type = c("prob"),newdata = held.out.log)
rocGraph.log <- roc(ids.held.out ~ probROC.log$`Cancer`)

plot(rocGraph.log, main = "PSA with model ROC")
rocGraph.log$auc

testingPred.log <- fiveFoldModel.log$pred
testingPred.log <- testingPred.log[which(testingPred.log$lambda > 0.04765),]
testingPred.log <- testingPred.log[which(testingPred.log$alpha == 0.55),]
plot.roc(pROC::roc(testingPred.log$obs ~ testingPred.log$Cancer))
trainingLogRoc <- roc(testingPred.log$obs ~ testingPred.log$Cancer)

### getting CI intervals for paper
#confusionMatrix(data = factor(prob.log), ids.held.out, positive = "Cancer")
ci.auc(rocGraph.log, conf.level = 0.95, method = "delong")
tableOfCIs <- rbind(tableOfCIs, c("Testing NGS-Mips", 82, 65, 98))


ci.auc(trainingLogRoc, conf.level = 0.95, method = "delong")
#confusionMatrix(data = factor(prob.train.log), ids.train, positive = "Cancer")
tableOfCIs <- rbind(tableOfCIs, c("Training NGS-Mips", 90, 83, 97))



### PSA

PSA.training <- PSA[-combined.subset]
PSA.training <- PSA.training[order(ids.4)]
a.5 <- data.frame(ids.3,PSA.training)
a.5$ones <- rep(1,nrow(a.5))

held.out.30.psa <- data.frame(PSA[combined.subset])
held.out.30.psa <- cbind(ids.held.out,held.out.30.psa)
held.out.30.psa$ones <- rep(1, nrow(held.out.30.psa))
colnames(held.out.30.psa) <- colnames(a.5)

prob.psa <- predict(fiveFoldModel.psa, newdata = held.out.30.psa)
confusionMatrix(data = factor(prob.psa), ids.held.out, positive = "Cancer")
probROC.psa <- predict(fiveFoldModel.psa, type = c("prob"),newdata = held.out.30.psa)
rocGraph.psa <- roc(ids.held.out ~ probROC.psa$`Cancer`)
plot(rocGraph.psa)
rocGraph.psa$auc

testingPred.psa <- fiveFoldModel.psa$pred
testingPred.psa <- testingPred.psa[which(testingPred.psa$alpha == 0.1),]
testingPred.psa <- testingPred.psa[which(testingPred.psa$lambda >= 0.031),]
plot.roc(pROC::roc(testingPred.psa$obs ~ testingPred.psa$Cancer))

trainingPsaRoc <- roc(testingPred.psa$obs ~ testingPred.psa$Cancer)

### PSA CIs
#confusionMatrix(data = factor(prob.psa), ids.held.out, positive = "Cancer")
ci.auc(rocGraph.psa, conf.level = 0.95, method = "delong")
tableOfCIs <- rbind(tableOfCIs, c("Testing PSA", 69, 51, 87))
ci.auc(trainingPsaRoc, conf.level = 0.95, method = "delong")
#confusionMatrix(data = factor(trainingPSAOrig), factor(ids.3), positive = "Cancer")
tableOfCIs <- rbind(tableOfCIs, c("Training PSA", 65, 51, 78))



###



### mips
###
###


tmprss2.fusion.training <- log2(a.2[,c("TMPRSS2-ERG.T1E4.COSF125")] + 1)[order(ids.4),]
PCA3 <- as.vector(log2(a.2[,which(grepl("PCA3", colnames(a.2)))] + 1))[order(ids.4),]
mips.genes <- data.frame(cbind(PSA.training, PCA3, tmprss2.fusion.training),stringsAsFactors = FALSE)
colnames(mips.genes) <- c("PSA","PCA3","TMPRSS2:ERG")
mips.genes <- cbind(ids.3, mips.genes) 



### making auc for testing
testingPred.mips <- fiveFoldModel.mips$pred
testingPred.mips <- testingPred.mips[which(testingPred.mips$lambda > 0.00047),]
testingPred.mips <- testingPred.mips[which(testingPred.mips$alpha == 0.55),]
plot.roc(pROC::roc(testingPred.mips$obs ~ testingPred.mips$Cancer))

trainingMipsRoc <- roc(testingPred.mips$obs ~ testingPred.mips$Cancer)

###

held.out.tmrpss <- log2(held.out[,c("TMPRSS2-ERG.T1E4.COSF125")] + 1)
held.out.pca <- log2(held.out[,which(grepl("PCA3", colnames(held.out)))] + 1)
held.out.mips <- cbind(ids.held.out,
                       held.out.30.psa[,2],
                       held.out.pca,
                       held.out.tmrpss)
colnames(held.out.mips) <- colnames(mips.genes)
prob.mips <- predict(fiveFoldModel.mips, newdata = held.out.mips)
confusionMatrix(data = factor(prob.mips), ids.held.out, positive = "Cancer")
probROC.mips <- predict(fiveFoldModel.mips, type = c("prob"),newdata = held.out.mips)
rocGraph.mips <- roc(ids.held.out ~ probROC.mips$`Cancer`)
plot(rocGraph.mips)
rocGraph.mips$auc


### Retrained Mips CIs
ci.auc(rocGraph.mips, conf.level = 0.95, method = "delong")
#confusionMatrix(data = factor(prob.mips), ids.held.out, positive = "Cancer")
tableOfCIs <- rbind(tableOfCIs, c("Retrained Mips Testing", 74, 55, 92))
ci.auc(trainingMipsRoc, conf.level = 0.95, method = "delong")
#confusionMatrix(data = factor(trainingRetMips), factor(ids.4), positive = "Cancer")
tableOfCIs <- rbind(tableOfCIs, c("Retrained Mips Training", 80, 74, 87))


### high risk mips
fullHighMipHeldTrain.predict <-  predict.glm(object = tmaMips.model.highRisk, type=c("response"), newdata = mips.genes)
y <- roc(ids.3 ~ fullHighMipHeldTrain.predict)
plot(y)
fullHighMipHeld.predict <-  predict.glm(object = tmaMips.model.highRisk, type=c("response"), newdata = held.out.mips)
z <- roc(ids.held.out ~ fullHighMipHeld.predict)
plot(z)


## CI intervals on original testing and training for high risk mips


ci.auc(y, conf.level = 0.95, method = "delong")
#confusionMatrix(factor(traingOrigHighGMips, levels = c("Cancer", "NoCancer")), factor(ids.3))
tableOfCIs <- rbind(tableOfCIs, c("Orig Training High risk mips", 72, 60, 84))

ci.auc(z, conf.level = 0.95, method = "delong")
#confusionMatrix(factor(testingOrigHighGMips, levels = c("Cancer", "NoCancer")), factor(ids.held.out))
tableOfCIs <- rbind(tableOfCIs, c("Orig Testing High risk mips", 69, 50, 87))

#write.table(tableOfCIs, "/mnt/DATA6/kevin_recovery/urineRNA/20190903ciTables.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#### training graph
#pdf(file = "/mnt/DATA6/kevin_recovery/urineRNA/20190915urineModelsTraining.pdf", onefile = TRUE, useDingbats = TRUE, pointsize = 18, width = 10, height = 10)
plot.roc(pROC::roc(testingPred.psa$obs ~ testingPred.psa$Cancer), col = alpha("blue", alpha = 0.5), asp = FALSE,
         main = "ROCs: Training n=73", ylab = "", xlab="")
par(new = TRUE)
plot.roc(pROC::roc(testingPred.log$obs ~ testingPred.log$Cancer), col = alpha("red", alpha = 0.5), asp = FALSE)
par(new = TRUE)
plot.roc(pROC::roc(testingPred.mips$obs ~ testingPred.mips$Cancer), col = alpha("brown4", alpha = 0.5), asp = FALSE)
par(new=TRUE)
plot(y, col = alpha("orange", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 1, asp=FALSE)
legend("bottomright", legend = c("29 transcript - AUC:0.90", "PSA - AUC:0.65",
                           "Retrained Mips: - AUC:0.80", "High-grade Mips - AUC:0.72"),
       col = c("red", "blue","brown4","orange"), lty=1, cex = .70)

dev.off()
###testing graph
#pdf(file = "/mnt/DATA6/kevin_recovery/urineRNA/20190915urineModelsTesting.pdf", onefile = TRUE, useDingbats = TRUE, pointsize = 18, width = 10, height = 10)
plot(rocGraph.psa, col = alpha("blue", alpha = 0.5), main = "ROCs: Testing n=36", ylab = "", xlab="", lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph.log, col = alpha("red", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 1, asp=FALSE)
par(new = TRUE)
plot(rocGraph.mips, col = alpha("brown4", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 1, asp=FALSE)
par(new = TRUE)
plot(z, col = alpha("orange", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 1, asp=FALSE)
legend("bottomright", legend = c("29 transcript - AUC:0.82", "PSA - AUC:0.69",
                           "Retrained Mips: - AUC:0.74", "High-grade Mips - AUC:0.69"),
       col = c("red", "blue","brown4","orange"), lty=1, cex = .70)
dev.off()



### below is for making the prob dataframes from above
transcriptPredTraning.predict <-  predict(object = fiveFoldModel.log,
                                          type=c("prob"), newdata = a.7[,2:ncol(a.7)])

extremeDfprobTrain <- data.frame(transcriptPredTraning.predict$Cancer)
rownames(extremeDfprobTrain) <- rownames(a.1)[-combined.subset]
psaPredTraning.predict <-  predict(object = fiveFoldModel.psa,
                                   type=c("prob"), newdata = a.5)
extremeDfprobTrain <- cbind(extremeDfprobTrain ,psaPredTraning.predict$Cancer)
rMipsPredTraning.predict <-  predict(object = fiveFoldModel.mips,
                                     type=c("prob"), newdata = mips.genes)
extremeDfprobTrain <- cbind(extremeDfprobTrain , rMipsPredTraning.predict$Cancer)
extremeDfprobTrain <- cbind(extremeDfprobTrain , fullHighMipHeldTrain.predict)
extremeDfprobTrain <- cbind(rownames(extremeDfprobTrain), extremeDfprobTrain)
colnames(extremeDfprobTrain) <- c("Samples","29_Transcript_Model", "PSA", "Retrained_Mips", "High-grade_Mips")
rownames(extremeDfprobTrain) <- NULL




extremeDfprobTesting <- data.frame(probROC.log$Cancer)
rownames(extremeDfprobTesting) <- rownames(a.1)[combined.subset]
extremeDfprobTesting <- cbind(extremeDfprobTesting, probROC.psa$Cancer)
extremeDfprobTesting <- cbind(extremeDfprobTesting, probROC.mips$Cancer)
extremeDfprobTesting <- cbind(extremeDfprobTesting, fullHighMipHeld.predict)
extremeDfprobTesting <- cbind(rownames(extremeDfprobTesting), extremeDfprobTesting)
colnames(extremeDfprobTesting) <- c("Samples","29_Transcript_Model", "PSA", "Retrained_Mips", "High-grade_Mips")
rownames(extremeDfprobTesting) <- NULL


write.table(extremeDfprobTrain,"/mnt/DATA6/kevin_recovery/urineRNA/20190628extremeTrainProbTable.txt", row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = "\t")
write.table(extremeDfprobTesting,"/mnt/DATA6/kevin_recovery/urineRNA/20190628extremeTestingProbTable.txt", row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = "\t")





