### 2021/03/21 Kevin Hu (kevhu@umich.edu)
### Code provided only includes code to use figures for paper
### Some of the code is abrdiged for clarity/brevity

### General setup: loading libraries and datasets 
library(VSURF)
library(parallel)
library(readxl)
library(ggplot2)
library(pROC)
library(reshape2)
library(gridExtra)
library(caret)
library(glmnet)
library(stringr)

allData <- read_xlsx("/mnt/DATA6/kevin_recovery/urineRNA/AC17.2ExtDes2_mastermatrix_allsamples_byKLK3.xlsx")
PSA <- read_xlsx("/path/to/file")
regModelCoeffs <- c(-4.466250739,0.357331671,0.577094327,0.194125163)
highRiskCoeffs <- c(-5.85877915,0.590375797,0.553164487,0.143707409)


### Removing GG2 tumors for extreme cohort design
### Creating dataframes (training/testing) downstream anlayses
allData <- allData[-which(allData$`Grade group` == 2),]
allData_transcripts <- allData[1:109,3:86]
rownames(allData_transcripts) <- as.vector(allData$Contig_ID[1:109])
ids <- allData$`Grade group bin`[1:109]
g1 <- which(ids == "Benign/GG1")
g2 <- which(ids == "GG2-5")
set.seed(2222)
g1_subset <- sample(g1, length(g1)/3)
set.seed(2222)
g2_subset <- sample(g2, length(g2)/3)
combined_subset <- c(g1_subset,g2_subset)

save(g1_subset, "/path/to/save/file")
save(g2_subset, "/path/to/save/file")

heldOut <- allData_transcripts[combined_subset,]
trainingSet <- allData_transcripts[-combined_subset,]

ids2 <- factor(ids[-combined_subset])
ids3 <- str_replace_all(ids2, "GG2-5", "Cancer")
ids3 <- str_replace_all(ids3, "Benign/GG1", "NoCancer")

### Use random forest in order to perform variable (transcript) selection
### Variable importance calculated from standard deviation
### Setting a specific type of seed for all machine-learning f(x)
### Regular seed types don't give reproducible results

testRandom <- list(allData_transcripts, ids2)
names(testRandom) <- c("x","y")
set.seed(2222, kind = "L'Ecuyer-CMRG")
vsurf.parallel <- VSURF(testRandom$x, testRandom$y,
                        mtry = 100, parallel = TRUE, ncores = 24,
                        clusterType = "FORK")
summary(vsurf.parallel)
colnames(a.1)[vsurf.parallel$varselect.interp]
plot(vsurf.parallel)

genesThres <- colnames(allData_transcripts)[vsurf.parallel$varselect.thres]
allData_genesThres <- allData_transcripts[,genesThres]
allData_genesThres <- allData_genesThres[, order(colnames(allData_genesThres))]

save(geneThres, "/path/to/save/file")

allData_caret <- cbind(ids3, log2(allData_genesThres[-combined_subset,] + 1))
allData_caret_m <- melt(data = allData_caret,id.vars = "ids3")

### Below I create weights for uneveness of the classes and parameters 
### used in for elasticnet (glmnet)


train_control <- trainControl(method = "cv", number = 5, savePredictions = TRUE,
                                       classProbs = TRUE, summaryFunction = twoClassSummary)
model_weights <- ifelse(allData_caret_m$ids3 == "Cancer",
                        (1/table(allData_caret_m$ids3)[1]) * 0.5,
                        (1/table(allData_caret_m$ids3)[2]) * 0.5)

set.seed(2222, kind = "L'Ecuyer-CMRG")
fiveFoldModel <- train(ids3 ~., data = allData_caret_m, trControl=train_control,
                       method="glmnet", metric = "ROC",
                       family = "binomial", weights = model_weights)

save(fiveFoldModel, "/path/to/save/file")


### Below I perform the same fittings but for the different models
### shown in the paper
### Note: For PSA, elasticnet should reduce to a simple logistic regression
### since it's only PSA + dummy variable (ones)

PSA_training <- PSA[-combined_subset]
PSA_caret <- data.frame(ids3,PSA_training)
PSA_caret$ones <- rep(1,nrow(PSA_caret))


set.seed(2222, kind = "L'Ecuyer-CMRG")
fiveFoldModel_psa <- train(ids3 ~., data = PSA_caret, trControl=train_control,
                           method="glmnet",metric = "ROC",
                           family = "binomial", weights = model_weights)

save(fiveFoldModel_psa, "/path/to/save/file")

mipsGenes <- c("PSA","PCA3","TMPRSS2:ERG")
mips_caret <- cbind(ids3, log2(allData_genesThres[-combined_subset, mipsGenes] + 1))

set.seed(2222, kind = "L'Ecuyer-CMRG")
fiveFoldModel_mips <- train(ids3 ~., data = mips_caret, trControl=train_control,
                            method="glmnet",metric = "ROC",
                            family = "binomial", weights = model_weights)
save(fiveFoldModel_mips, "/path/to/save/file")


### Regular glm models used for recreating selectMdx and exome
### Didn't have coefficients for these two tests

selectMdx <- c("HOXC6.E1E2", "DLX1.E1E2")
selectMdxDf <- cbind(ids3, log2(allData_genesThres[-combined_subset, selectMdx] + 1))
colnames(selectMdxf) <- c("hox", "dlx")
logReg_mdx <- glm(ids ~ selectMdxDf$hox + selectMdxDf$dlx , family = binomial())

selectMdxDf_test <- cbind(ids[combined_subset],log2(allData_genesThres[combined_subset, selectMdx] + 1))
colnames(selectMdxDf_test) <- c("hox", "dlx")

exome <- c("PCA3.E2E3", "ERG.E1E2")
exomeDf <- cbind(ids3, log2(allData_genesThres[-combined_subset, exome] + 1))
colnames(exomeDf) <- c("pca3", "erg")
logReg_mdx <- glm(ids ~ exomeDf$pca3 + exomeDf$erg , family = binomial())

exomeDf_test <- cbind(ids[combined_subset], log2(allData_genesThres[combined_subset, exome] + 1))
colnames(exomeDf_test) <- c("pca3", "erg")
### For Mips model that we have weights for, I train a dummy glm object
### and replace the coefficients with the known one
### Note: Mips created from elasticnet was a retrained version

mipsDf <- cbind(ids3, log2(allData_genesThres[-combined_subset, mipsGenes] + 1))
mipsDf_heldout <- cbind(ids[combined_subset], log2(allData_genesThres[combined_subset, mipsGenes] + 1))
dummyModel <- andiAnno.model <- glm(data = fullDf, formula = ids3 ~ .,
                                    family = binomial("logit"))

mips_cancer <- dummyModel
mips_cancer$coefficients <- regModelCoeffs
mips_highgrade <- dummyModel
mips_highgrade$coefficients <- highRiskCoeffs

### This portion is creating area under the reciever operator curve (AUROC)
### and confidene intervals (CI) - different models almost identical workflows
### as first instance


### Fivefold model training and testing ROCs
training_fivefold_probs <- fiveFoldModel$pred
training_fivefold_probs <- training_fivefold_probs[which(training_fivefold_probs$lambda == fiveFoldModel$bestTune$lambda &
                                                           training_fivefold_probs$alpha == fiveFoldModel$bestTune$alpha),]
roc_fivefold_train <- plot.roc(pROC::roc(training_fivefold_probs$obs ~ training_fivefold_probs$high))
pred_fivefold_test <- predict(fiveFoldModel, type = c("prob"),
                              newdata = log2(heldOut[,genesThres] + 1))
roc_fivefold_test <- roc(pred_fivefold_test$Grade ~ pred_fivefold_test$high)

ci.auc(roc_fivefold_train, conf.level = 0.95, method = "delong")
ci.auc(roc_fivefold_test, conf.level = 0.95, method = "delong")

### PSA

training_psa_probs <- fiveFoldModel_psa$pred
training_psa_probs <- training_psa_probs[which(training_psa_probs$lambda == fiveFoldModel_psa$bestTune$lambda &
                                                           training_psa_probs$alpha == fiveFoldModel_psa$bestTune$alpha),]
roc_psa_train <- plot.roc(pROC::roc(training_psa_probs$obs ~ training_psa_probs$high))
pred_psa_test <- predict(fiveFoldModel_psa, type = c("prob"),
                              newdata = PSA[combined_subset])
roc_psa_test <- roc(pred_psa_test$Grade ~ pred_psa_test$high)

ci.auc(roc_psa_train, conf.level = 0.95, method = "delong")
ci.auc(roc_psa_test, conf.level = 0.95, method = "delong")

### Retrained Mips variables
mipsDf_retrain_heldout <- log2(allData_genesThres[-combined_subset, mipsGenes] + 1)

training_mips_probs <- fiveFoldModel_mips$pred
training_mips_probs <- training_mips_probs[which(training_mips_probs$lambda == fiveFoldModel_mips$bestTune$lambda &
                                                           training_mips_probs$alpha == fiveFoldModel_mips$bestTune$alpha),]
roc_mips_train <- plot.roc(pROC::roc(training_mips_probs$obs ~ training_mips_probs$high))
pred_mips_test <- predict(fiveFoldModel_mips, type = c("prob"),
                              newdata = mipsDf_retrain_heldout)
roc_mips_test <- roc(pred_mips_test$Grade ~ pred_mips_test$high)

ci.auc(roc_mips_train, conf.level = 0.95, method = "delong")
ci.auc(roc_mips_test, conf.level = 0.95, method = "delong")


### Regular and high-grade Mips

trainReg_predict <- predict.glm(object = mips_cancer,
                                type=c("response"), newdata = mipsDf)
trainReg_roc <- roc(mips_genes$ids3 ~ trainReg_predict)
testingReg_predict <- predict.glm(object = mips_cancer, type=c("response"),
                                  newdata = mipsDf_heldout)
testingReg_roc <- roc(mips_cancer$ids ~ testingReg_predict)

ci.auc(trainReg_roc, conf.level = 0.95, method = "delong")
ci.auc(testingReg_roc, conf.level = 0.95, method = "delong")

trainHigh_predict <- predict.glm(object = mips_highgrade,
                                type=c("response"), newdata = mipsDf)
trainHigh_roc <- roc(mips_genes$ids3 ~ trainHigh_predict)
testingHigh_predict <- predict.glm(object = mips_highgrade, type=c("response"),
                                  newdata = mipsDf_heldout)
testingHigh_roc <- roc(mips_highgrade$ids ~ testingHigh_predict)

ci.auc(trainHigh_roc, conf.level = 0.95, method = "delong")
ci.auc(testingHigh_roc, conf.level = 0.95, method = "delong")



### Last portion is making AUROC figures
###
###


#### Training data graph
pdf(file = "/path/to/file", onefile = TRUE, useDingbats = TRUE,
    pointsize = 18, width = 10, height = 10)
plot.roc(roc_psa_train, col = alpha("blue", alpha = 0.5), asp = FALSE,
         main = "ROCs: Training n=73", ylab = "", xlab="")
par(new = TRUE)
plot.roc(roc_fivefold_train, col = alpha("red", alpha = 0.5), asp = FALSE)
par(new = TRUE)
plot.roc(roc_mips_train, col = alpha("brown4", alpha = 0.5), asp = FALSE)
par(new=TRUE)
plot(y, col = alpha("orange", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 1, asp=FALSE)
legend("bottomright", legend = c("29 transcript - AUC:0.90", "PSA - AUC:0.65",
                                 "Retrained Mips: - AUC:0.80", "High-grade Mips - AUC:0.72"),
       col = c("red", "blue","brown4","orange"), lty=1, cex = .70)

dev.off()

### Testing data graph
pdf(file = "/path/to/file", onefile = TRUE, useDingbats = TRUE,
    pointsize = 18,width = 10, height = 10)
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
