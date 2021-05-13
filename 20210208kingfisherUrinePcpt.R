
library(VSURF)
library(readxl)
library(ggplot2)
library(pROC)
library(reshape2)
library(gridExtra)
library(caret)
library(stringr)
library(Rcpp)
library(parallel)

getUrineIds <- function(filepath, page){
  tmpXls <- readxl::read_xlsx(path = as.character(filepath), sheet = as.numeric(page))
  listOfNames <- str_remove(string = tmpXls$LABEL, pattern = "CCGC_")
  listOfNames <- str_remove(listOfNames, pattern = "_.")
  listOfNames <- listOfNames[-which(is.na(listOfNames))]
}

getUrineIds_PSA <- function(filepath, page){
  tmpXls <- readxl::read_xlsx(path = as.character(filepath), sheet = as.numeric(page))
  listOfNames <- str_remove(string = tmpXls$LABEL, pattern = "CCGC_")
  listOfNames <- str_remove(listOfNames, pattern = "_.")
  listOfNames <- listOfNames[-which(is.na(listOfNames))]
  listOfPSA <- tmpXls$`Pre-Bx PSA Value`
  listOfPSA <- listOfPSA[-which(is.na(listOfPSA))]
  tmpDf <- cbind(listOfNames, listOfPSA)
}





rmBad <- function(dat, badSample){
  res <- dat[-which(dat %in% str_remove(badSample, "UR.*"))]
  if (length(res) == 0) {
    res <- dat
  }
  return(res)
}


seqStats <- read.csv("/mnt/DATA6/kevin_recovery/urineRNA/20200421_All urine samps pheno.csv",
                     stringsAsFactors = FALSE)

urineData <- read.csv("/mnt/DATA6/kevin_recovery/urineRNA/20200421_All urine samps hdat.csv", 
                      stringsAsFactors = FALSE)

sampleNames <- str_remove(colnames(urineData[,2:ncol(urineData)]), "X")
geneNames <- urineData[,1]
urineData_processed <- t(urineData[,2:ncol(urineData)])
rownames(urineData_processed) <- sampleNames
colnames(urineData_processed) <- geneNames
urineData_processed <- data.frame(urineData_processed, stringsAsFactors = FALSE)


### line below is for testing all fusions ERG 
fusionCols <- c(grep("\\.ERG", colnames(urineData_processed)),
                grep("\\.ETV", colnames(urineData_processed)))

fusionVar <- apply(urineData_processed[,fusionCols], 1, sum)
urineData_processed <- urineData_processed[,-fusionCols]
urineData_processed$fusions <- fusionVar



### deciding not to use klk3 b/c of how it doesn't seem to correlate all that well
#housekeeping <- urineData_processed$KLK3.E2E3 

### need to normalize


#urineData_processed_norm <- apply(urineData_processed, 1, function(x) log2(x/(sum(x)/1e6) + 1))
#urineData_processed_norm <- t(urineData_processed_norm)

housekeeping <- urineData_processed$KLK3.E2E3
urineData_processed_norm <- NULL
for (i in 1:nrow(urineData_processed)) {
  tmpDf <- log2(urineData_processed[i,]/housekeeping[i] * 100000 + 1)
  urineData_processed_norm <- rbind(urineData_processed_norm, tmpDf) 
}

urineData_processed_norm <- urineData_processed_norm[,order(colnames(urineData_processed_norm))]


#colnames(urineData_processed_norm) <- geneNames

### filter bad samples. tot mapped reads > 300000 & %e2e > .50 & %OT > .6
### in this not filtering by OT % because it wasnt in sheet, but checked it all maually (were good)

badSamples <- seqStats$Sample[-which(seqStats$Total.reads > 300000 & seqStats$percent_e2e > .50)]
badSamples <- str_remove(badSamples, "X")
badSamples <- c(badSamples, "91669", "91170", "90956", "91423", "91068", "90630", "91362",
                "92297", "90658", "91212", "91149", "91573", "90679", "91291", "91322", "91588", "91659", "89684")

urineData_processed_norm <- urineData_processed_norm[-which(rownames(urineData_processed_norm) %in% badSamples),]





### getting variables for PCPT variables to then use for calculator 
### b0 + b1 * log(PSA) + b2 * fam. history + b3 * DRE + b4 * prior biopsy
### need to create functions to pull these out from samples. test that my own logistic equation works. then I can try w.e with the score

path <- "/mnt/DATA6/kevin_recovery/urineRNA/DOD_Aim 1 Cohort_Urine Sample List_Random Lists_Benign_3+3_3+4_4+3_4+4_9 and 10WITH identifiers_FINAL.xlsx"
allSampPcpt <- NULL
pcptDatPull <- for (i in 1:6) {
  tmpXls <- readxl::read_xlsx(path = as.character(path), sheet = i)
  tmpXls <- tmpXls[2:nrow(tmpXls),]
  listOfNames <- str_remove(string = tmpXls$LABEL, pattern = "CCGC_")
  listOfNames <- str_remove(listOfNames, pattern = "_.")
  tmpPsa <- tmpXls$`Pre-Bx PSA Value`
  tmpFamHis <- tmpXls$FAM_HISTORY
  tmpDre <- tmpXls$DRE
  tmpPrior <- tmpXls$`# Prior Biopsies`
  tmpAge <- tmpXls$BIRTH_DATE
  tmpAge2 <- as.numeric(unlist(lapply(str_split(tmpAge, "-"), '[[', 1)))
  tmpAge2 <- 2021 - tmpAge2
  grade <- i
  
  
  tmpDf <- cbind("SampNames" = listOfNames, "PSA" = tmpPsa, "FamHis" = tmpFamHis,
                 "DRE" = tmpDre, "Prior" = tmpPrior, "Age" = tmpAge2,
                 "gradeHolder" = grade)
  allSampPcpt <- rbind(allSampPcpt, tmpDf)
}

allVarsPcpt <- data.frame(allSampPcpt, stringsAsFactors = FALSE)
allSampPcpt <- allSampPcpt[,1:5]

pcptCoeff <- c(-1.80, 0.85, 0.27, 0.91, -0.45)
allSampPcpt2 <- data.frame(allSampPcpt, stringsAsFactors = FALSE)



allSampPcpt2$dummyVar <- "Cancer"
allSampPcpt2$dummyVar[1:100] <- "NoCancer"
allSampPcpt2$dummyVar <- factor(allSampPcpt2$dummyVar)
allSampPcpt2$dummyVar <- relevel(allSampPcpt2$dummyVar, "NoCancer")
allSampPcpt2 <- allSampPcpt2[-which(is.na(allSampPcpt2$DRE) | is.na(allSampPcpt2$FamHis) | is.na(allSampPcpt2$Prior)),]
allSampPcpt2$Prior[1] <- 1
allSampPcpt2_1 <- allSampPcpt2[,2:6]
allSampPcpt2_1[,1:4] <- lapply(allSampPcpt2_1[,1:4], function(x) as.numeric(x))
dummyModel <- glm(formula = dummyVar ~ PSA + FamHis + DRE + Prior, data = allSampPcpt2_1, family = binomial(link = "logit"))


pcptToUse <- allSampPcpt2[,2:6]
pcptToUse$Prior[1] <- 0
pcptToUse[,1:4] <- lapply(pcptToUse[,1:4], function(x) as.numeric(x))
pcptToUse$PSA <- log(pcptToUse$PSA)

pcptModel <- dummyModel
pcptModel$coefficients
pcptModel$coefficients <- pcptCoeff

pcptToUse$DRE[which(pcptToUse$DRE == 1)] <- 0
pcptToUse$DRE[which(pcptToUse$DRE == 2)] <- 1

pcpt_predict <-  predict(object = pcptModel, type=c("response"), newdata = pcptToUse[,1:4])
pcpt_predict

exp(sum(pcptToUse[2,1:4] * pcptCoeff[2:5],pcptCoeff[1])) / (1 + exp(sum(pcptToUse[2,1:4] * pcptCoeff[2:5], pcptCoeff[1])))

allSampPcpt2$Pcpt <- pcpt_predict

### treating the variales like they are .. i.e numeric for them all
### questions, unknown variables - skip sample? also pcpt variable is just the percentage chance? how to verify?


### try treating things as categorical variables - below shows, doesn't change much to make numeric off binary variables or dummy vars for
### famhis or remove prior bx which is all 0 
pcptToUse2 <- pcptToUse[,c(1:3,5)]
pcptToUse2 <- pcptToUse2[-which(pcptToUse2$DRE == 3),]
pcptToUse2$FamHis <- factor(pcptToUse2$FamHis)
pcptToUse2$DRE <- factor(pcptToUse2$DRE)

dummyModel2 <- glm(formula = dummyVar ~ PSA + FamHis + DRE, data = pcptToUse2, family = binomial(link = "logit"))
pcptModel2 <- dummyModel2
pcptModel2$coefficients <- pcptCoeff[1:4]

pcpt_predict2 <-  predict(object = pcptModel2, type=c("response"), newdata = pcptToUse2[,1:3])
pcpt_predict2



### above is all the same results - so good that basic logistic regression gets the probability
### 


trainingDf_benign_g1 <- read.table("/mnt/DATA6/kevin_recovery/urineRNA/20200512trainingDf_benign_g1.txt",
                                   sep = "\t", stringsAsFactors = FALSE, header = TRUE)
trainingDf_benign_g1$Grade <- factor(trainingDf_benign_g1$Grade, levels = c("low", "high"))

testingDf_benign_g1 <- read.table("/mnt/DATA6/kevin_recovery/urineRNA/20200512testingDf_benign_g1.txt",
                                  sep = "\t", stringsAsFactors = FALSE, header = TRUE)
testingDf_benign_g1$Grade <- factor(testingDf_benign_g1$Grade, levels = c("low","high"))



trainingDf_benign_g2 <- read.table("/mnt/DATA6/kevin_recovery/urineRNA/20200512trainingDf_benign_g2.txt",
                                   sep = "\t", stringsAsFactors = FALSE, header = TRUE)
trainingDf_benign_g2$Grade <- factor(trainingDf_benign_g2$Grade, levels = c("low", "high"))


testingDf_benign_g2 <- read.table("/mnt/DATA6/kevin_recovery/urineRNA/20200512testingDf_benign_g2.txt",
                                  sep = "\t", stringsAsFactors = FALSE, header = TRUE)
testingDf_benign_g2$Grade <- factor(testingDf_benign_g2$Grade, levels = c("low","high"))


reseqedSamps <- str_remove_all(rownames(urineData_processed_norm)[grep("\\.1", rownames(urineData_processed_norm))], "\\.1")
urineData_processed_norm <- urineData_processed_norm[-which(rownames(urineData_processed_norm) %in% reseqedSamps),]

### getting rid of samples with no pcpt score

pcptFinalPredicts <- allSampPcpt2
pcptFinalPredicts <- allSampPcpt2[-which(allSampPcpt2$DRE == 3),]
tmpVar <- str_remove(str_remove(rownames(urineData_processed_norm), "UR.*"), "\\.1")

newExpressionNorms <-  data.frame(urineData_processed_norm[match(pcptFinalPredicts$SampNames, tmpVar),],
                                  stringsAsFactors = FALSE)

newExpressionNorms$Pcpt <- pcptFinalPredicts$Pcpt
newExpressionNorms <- newExpressionNorms[-which(is.na(newExpressionNorms$AC009478.1.E3E4)),]


tmpVar2 <- str_remove(str_remove(rownames(newExpressionNorms), "UR.*"), "\\.1")
### vsurf portion



bengin_g1_rf_df <- data.frame(newExpressionNorms[match(as.character(trainingDf_benign_g1$Sample),tmpVar2),],
                              stringsAsFactors = FALSE)
bengin_g1_rf_df <- bengin_g1_rf_df[-which(is.na(bengin_g1_rf_df$AC009478.1.E3E4)),]
### sanity check
notInG1 <- str_remove_all(str_remove_all(rownames(bengin_g1_rf_df), "\\.1"), "UR.*")
trainingDf_benign_g1 <- trainingDf_benign_g1[-which((!trainingDf_benign_g1$Sample %in% notInG1)),]

ids_benign_g1 <- factor(trainingDf_benign_g1$Grade)

bengin_g2_rf_df <- data.frame(newExpressionNorms[match(as.character(trainingDf_benign_g2$Sample),tmpVar2),],
                              stringsAsFactors = FALSE)

bengin_g2_rf_df <- bengin_g2_rf_df[-which(is.na(bengin_g2_rf_df$AC009478.1.E3E4)),]
notInG2 <- str_remove_all(str_remove_all(rownames(bengin_g2_rf_df), "\\.1"), "UR.*")
trainingDf_benign_g2 <- trainingDf_benign_g2[-which((!trainingDf_benign_g2$Sample %in% notInG2)),]

ids_benign_g2 <- factor(trainingDf_benign_g2$Grade)





set.seed(182)
vsurf_parallel_benign_g1 <- VSURF(bengin_g1_rf_df, ids_benign_g1, parallel = TRUE, ncores = 24, clusterType = "FORK")
summary(vsurf_parallel_benign_g1)
plot(vsurf_parallel_benign_g1)

genesThres_benign_g1 <- colnames(bengin_g1_rf_df)[vsurf_parallel_benign_g1$varselect.thres]
genesInterp_benign_g1 <- colnames(bengin_g1_rf_df)[vsurf_parallel_benign_g1$varselect.interp]




set.seed(82123)
vsurf_parallel_benign_g2 <- VSURF(bengin_g2_rf_df, ids_benign_g2, parallel = TRUE, ncores = 24, clusterType = "FORK")
summary(vsurf_parallel_benign_g2)
plot(vsurf_parallel_benign_g2)

genesThres_benign_g2 <- colnames(bengin_g2_rf_df)[vsurf_parallel_benign_g2$varselect.thres]
genesInterp_benign_g2 <- colnames(bengin_g2_rf_df)[vsurf_parallel_benign_g2$varselect.interp]



### setting up the models

ids_b_g1 <- trainingDf_benign_g1$Grade
ids_b_g2 <- trainingDf_benign_g2$Grade


bengin_g1_rf_df_genes <- bengin_g1_rf_df[, genesThres_benign_g1]
bengin_g2_rf_df_genes <- bengin_g2_rf_df[, genesThres_benign_g2]

bengin_g1_rf_df_genes <- cbind(ids_b_g1, bengin_g1_rf_df_genes)
bengin_g2_rf_df_genes <- cbind(ids_b_g2, bengin_g2_rf_df_genes)

bengin_g1_rf_df_genes$ids_b_g1 <- relevel(bengin_g1_rf_df_genes$ids_b_g1, ref = "low")
bengin_g2_rf_df_genes$ids_b_g2 <- relevel(bengin_g2_rf_df_genes$ids_b_g2, ref = "low")

model_weights_g1 <- ifelse(trainingDf_benign_g1$Grade == "high",
                           (1/table(trainingDf_benign_g1$Grade)[1]) * 0.5,
                           (1/table(trainingDf_benign_g1$Grade)[2]) * 0.5)
model_weights_g2 <- ifelse(trainingDf_benign_g2$Grade == "high",
                           (1/table(trainingDf_benign_g2$Grade)[1]) * 0.5,
                           (1/table(trainingDf_benign_g2$Grade)[2]) * 0.5)

train_control_rules <- trainControl(method = "cv", number = 5, savePredictions = TRUE, classProbs = TRUE,
                                    summaryFunction = twoClassSummary)



RNGkind("L'Ecuyer-CMRG")
set.seed(54763)
fiveFoldModel_g1 <- caret::train(ids_b_g1~., data = bengin_g1_rf_df_genes, trControl = train_control_rules,
                                 method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g1)
#save(fiveFoldModel_g1, file = "/mnt/DATA6/kevin_recovery/urineRNA/20210214fivefold_b_g1_combined.Robj")

RNGkind("L'Ecuyer-CMRG")
set.seed(6654763)
fiveFoldModel_g2 <- caret::train(ids_b_g2 ~., data = bengin_g2_rf_df_genes, trControl = train_control_rules,
                                 method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g2)
#save(fiveFoldModel_g2, file = "/mnt/DATA6/kevin_recovery/urineRNA/20210214fivefold_b_g2_combined.Robj")

### can save models above for later

genesThres_benign_g1 <- fiveFoldModel_g1$coefnames
genesThres_benign_g2 <- fiveFoldModel_g2$coefnames

ids_b_g1 <- trainingDf_benign_g1$Grade
ids_b_g2 <- trainingDf_benign_g2$Grade


bengin_g1_rf_df_genes <- bengin_g1_rf_df[, genesThres_benign_g1]
bengin_g2_rf_df_genes <- bengin_g2_rf_df[, genesThres_benign_g2]

bengin_g1_rf_df_genes <- cbind(ids_b_g1, bengin_g1_rf_df_genes)
bengin_g2_rf_df_genes <- cbind(ids_b_g2, bengin_g2_rf_df_genes)


testing_b_g1_df <- data.frame(newExpressionNorms[match(as.character(testingDf_benign_g1$Sample),tmpVar2),],
                              stringsAsFactors = FALSE)
testing_b_g1_df <- testing_b_g1_df[-which(is.na(testing_b_g1_df$AC009478.1.E3E4)),]
testing_b_g1_df <- testing_b_g1_df[, genesThres_benign_g1]
notInG1_test <- str_remove_all(str_remove_all(rownames(testing_b_g1_df), "\\.1"), "UR.*")
testingDf_benign_g1 <- testingDf_benign_g1[-which((!testingDf_benign_g1$Sample %in% notInG1_test)),]



testing_b_g2_df <- data.frame(newExpressionNorms[match(as.character(testingDf_benign_g2$Sample),tmpVar2),],
                              stringsAsFactors = FALSE)
testing_b_g2_df <- testing_b_g2_df[-which(is.na(testing_b_g2_df$AC009478.1.E3E4)),]
testing_b_g2_df <- testing_b_g2_df[, genesThres_benign_g2]
notInG2_test <- str_remove_all(str_remove_all(rownames(testing_b_g2_df), "\\.1"), "UR.*")
testingDf_benign_g2 <- testingDf_benign_g2[-which((!testingDf_benign_g2$Sample %in% notInG2_test)),]



predRes_b_g1 <- predict(fiveFoldModel_g1, type = c("prob"),newdata = testing_b_g1_df)
rocGraph_b_g1 <- roc(testingDf_benign_g1$Grade ~ predRes_b_g1$high)
predRes_b_g2 <- predict(fiveFoldModel_g2, type = c("prob"), newdata = testing_b_g2_df)
rocGraph_b_g2 <- roc(testingDf_benign_g2$Grade ~ predRes_b_g2$high)
training_b_g1_probs <- fiveFoldModel_g1$pred
training_b_g1_probs <- training_b_g1_probs[which(training_b_g1_probs$lambda == fiveFoldModel_g1$bestTune$lambda &
                                                   training_b_g1_probs$alpha == fiveFoldModel_g1$bestTune$alpha),]
roc_train_b_g1 <- plot.roc(pROC::roc(training_b_g1_probs$obs ~ training_b_g1_probs$high))
training_b_g2_probs <- fiveFoldModel_g2$pred
training_b_g2_probs <- training_b_g2_probs[which(training_b_g2_probs$lambda == fiveFoldModel_g2$bestTune$lambda &
                                                   training_b_g2_probs$alpha == fiveFoldModel_g2$bestTune$alpha),]
roc_train_b_g2 <- plot.roc(pROC::roc(training_b_g2_probs$obs ~ training_b_g2_probs$high))




### make simple logistic regression for pcpt
###
###

bengin_g1_rf_df_pcpt <- data.frame("ids_b_g1" = bengin_g1_rf_df_genes$ids_b_g1,"ones" = 1, "Pcpt" = bengin_g1_rf_df_genes$Pcpt)
rownames(bengin_g1_rf_df_pcpt) <- rownames(bengin_g1_rf_df_genes)
RNGkind("L'Ecuyer-CMRG")
set.seed(547003)
fiveFoldModel_g1_pcpt <- caret::train(ids_b_g1~., data = bengin_g1_rf_df_pcpt, trControl = train_control_rules,
                                 method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g1)



bengin_g2_rf_df_pcpt <- data.frame("ids_b_g1" = bengin_g2_rf_df_genes$ids_b_g2,"ones" = 1, "Pcpt" = bengin_g2_rf_df_genes$Pcpt)
rownames(bengin_g2_rf_df_pcpt) <- rownames(bengin_g2_rf_df_genes)
RNGkind("L'Ecuyer-CMRG")
set.seed(547003)
fiveFoldModel_g2_pcpt <- caret::train(ids_b_g1~., data = bengin_g2_rf_df_pcpt, trControl = train_control_rules,
                                 method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g1)


testing_b_g1_df_pcpt <- data.frame("ones" = 1, "Pcpt" = testing_b_g1_df[,1])
testing_b_g2_df_pcpt <- data.frame("ones" = 1, "Pcpt" = testing_b_g2_df[,1])


predRes_b_g1_pcpt <- predict(fiveFoldModel_g1_pcpt, type = c("prob"),newdata = testing_b_g1_df_pcpt)
rocGraph_b_g1_pcpt <- roc(testingDf_benign_g1$Grade ~ predRes_b_g1_pcpt$high)
predRes_b_g2_pcpt <- predict(fiveFoldModel_g2_pcpt, type = c("prob"), newdata = testing_b_g2_df_pcpt)
rocGraph_b_g2_pcpt <- roc(testingDf_benign_g2$Grade ~ predRes_b_g2_pcpt$high)
training_b_g1_probs_pcpt <- fiveFoldModel_g1_pcpt$pred
training_b_g1_probs_pcpt <- training_b_g1_probs_pcpt[which(training_b_g1_probs_pcpt$lambda == fiveFoldModel_g1_pcpt$bestTune$lambda &
                                                   training_b_g1_probs_pcpt$alpha == fiveFoldModel_g1_pcpt$bestTune$alpha),]

roc_train_b_g1_pcpt <- plot.roc(pROC::roc(training_b_g1_probs_pcpt$obs ~ training_b_g1_probs_pcpt$high))
training_b_g2_probs_pcpt <- fiveFoldModel_g2_pcpt$pred
training_b_g2_probs_pcpt <- training_b_g2_probs_pcpt[which(training_b_g2_probs_pcpt$lambda == fiveFoldModel_g2_pcpt$bestTune$lambda &
                                                   training_b_g2_probs_pcpt$alpha == fiveFoldModel_g2_pcpt$bestTune$alpha),]
roc_train_b_g2_pcpt <- plot.roc(pROC::roc(training_b_g2_probs_pcpt$obs ~ training_b_g2_probs_pcpt$high))

### plots


dev.off()
###testing graph
#pdf(file = "/mnt/DATA5/tmp/kev/misc/20210214benignGG1_combined.pdf", onefile = TRUE, useDingbats = TRUE, pointsize = 18, width = 10, height = 10)
plot(roc_train_b_g1, col = alpha("red"), main = "ROCs: BGG1 vs GG2GG5", ylab = "", xlab="", lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g1, col = alpha("red"), lty = 2, asp = FALSE)
par(new = TRUE)
plot(roc_train_b_g1_pcpt, col = alpha("black"), lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g1_pcpt, col = alpha("black"),ylab = "Sensitivity", xlab="Specificity", lty = 2, asp = FALSE)
legend("bottomright", legend = c("PCPT Training - AUC:0.66", "PCPT Testing - AUC: 0.57",
                                 "Model Training - AUC:0.66", "Model Testing - AUC: 0.61"),
       col = c("black", "black", "red", "red"), lty=c(1,2,1,2), cex = .70)
dev.off()



dev.off()
###testing graph
#pdf(file = "/mnt/DATA5/tmp/kev/misc/20210214benignGG2_combined.pdf", onefile = TRUE, useDingbats = TRUE, pointsize = 18, width = 10, height = 10)
plot(roc_train_b_g2, col = alpha("red"), main = "ROCs: BGG2 vs GG3GG5", ylab = "", xlab="", lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g2, col = alpha("red"), lty = 2, asp = FALSE)
par(new = TRUE)
plot(roc_train_b_g2_pcpt, col = alpha("black"), lty = 2, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g2_pcpt, col = alpha("black"),ylab = "Sensitivity", xlab="Specificity", lty = 2, asp = FALSE)

legend("bottomright", legend = c("PCPT Training - AUC:0.68", "PCPT Testing - AUC: 0.60",
                                 "Model Training - AUC:0.74", "Testing - AUC: 0.51"),
       col = c("black", "black", "red", "red"), lty=c(1,2,1,2), cex = .70)
dev.off()

### get weights
genes_b_g1_weights <- as.matrix(coef(fiveFoldModel_g1$finalModel, fiveFoldModel_g1$bestTune$lambda))
genes_b_g2_weights <- as.matrix(coef(fiveFoldModel_g2$finalModel, fiveFoldModel_g2$bestTune$lambda))




#m_weights_g1 <- genes_b_g1_weights
#m_weights_g2 <- genes_b_g2_weights
m_weights_g1
m_weights_g2
#combined_weights_g1 <- genes_b_g1_weights
#combined_weights_g2 <- genes_b_g2_weights
combined_weights_g1
combined_weights_g2

### try summing up all the T2 erg isoforms + list of good/poor quality samples
### summary of age, psa .... DRE etc for each grade group testing and training - maybe boxplots (maybe just Mips genes)


allVarsPcpt2 <- allVarsPcpt[-which(is.na(allVarsPcpt$DRE) | is.na(allVarsPcpt$FamHis) | is.na(allVarsPcpt$Prior)),]
allVarsPcpt2 <- allVarsPcpt2[-which(allVarsPcpt2$DRE == "3"),]
allVarsPcpt2$split <- "none"



allVarsPcpt2$split[which(allVarsPcpt2$SampNames %in% trainingDf_benign_g1$Sample)] <- "training"
allVarsPcpt2$split[which(allVarsPcpt2$SampNames %in% testingDf_benign_g1$Sample)] <- "testing"


allVarsPcpt2 <- allVarsPcpt2[-which(allVarsPcpt2$split == "none"),]

tmpGeneExpression <- urineData_processed_norm[match(as.character(allVarsPcpt2$SampNames),tmpVar),
                                              c("TMPRSS2.ERG.T1E4.COSF125", "PCA3.E2E3")]

allVarsPcpt2 <- cbind(allVarsPcpt2, tmpGeneExpression)
allVarsPcpt2$gradeHolder <- paste0("GG_", as.numeric(allVarsPcpt2$gradeHolder) - 1)
allVarsPcpt2$FamHis <- str_replace_all(allVarsPcpt2$FamHis, "0", "No")
allVarsPcpt2$FamHis <- str_replace_all(allVarsPcpt2$FamHis, "1", "Yes")
allVarsPcpt2$DRE <- str_replace_all(allVarsPcpt2$DRE, "1", "Negative")
allVarsPcpt2$DRE <- str_replace_all(allVarsPcpt2$DRE, "2", "Positive")

allVarsPcpt2$PSA <- as.numeric(allVarsPcpt2$PSA)
allVarsPcpt2$Age <- as.numeric(allVarsPcpt2$Age)
allVarsPcpt2$TMPRSS2.ERG.T1E4.COSF125 <- as.numeric(allVarsPcpt2$TMPRSS2.ERG.T1E4.COSF125)
allVarsPcpt2$PCA3.E2E3 <- as.numeric(allVarsPcpt2$PCA3.E2E3)

colnames(allVarsPcpt2)[9] <- "T2ERG"
colnames(allVarsPcpt2)[10] <- "PCA3"





###
###
### for diff splits
allVarsPcpt2$gradeHolder[which(allVarsPcpt2$gradeHolder %in% c("GG_0", "GG_1"))] <- "low"
allVarsPcpt2$gradeHolder[which(allVarsPcpt2$gradeHolder %in% c("GG_2", "GG_3","GG_4", "GG_5"))] <- "high"

allVarsPcpt2$gradeHolder[which(allVarsPcpt2$gradeHolder %in% c("GG_0", "GG_1", "GG_2"))] <- "low"
allVarsPcpt2$gradeHolder[which(allVarsPcpt2$gradeHolder %in% c("GG_3","GG_4", "GG_5"))] <- "high"




for (i in unique(allVarsPcpt2$gradeHolder)) {
  tmpDf <- allVarsPcpt2[allVarsPcpt2$gradeHolder == i,]
  tmpDf_plot <- tmpDf[,c(2, 8, 9,10)]
  tmpDf_plot_m <- melt(tmpDf_plot, id.vars = "split")
  
  ntest <- table(tmpDf$split)[1]
  ntrain <- table(tmpDf$split)[2]
  
  q <- ggplot(tmpDf_plot_m, aes(x=variable, y=value)) + geom_boxplot(aes(fill = factor(split))) + ggtitle(paste0(i,"(n = ", ntest,";", ntrain, ")"))
  q <- q + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + ylim(c(0,25))
    
  assign(i, q)
  for (j in unique(allVarsPcpt2$split)) {
    tmpDf2 <- tmpDf[which(tmpDf$split == j),]
    print(paste0(i, "_FamHis_", j, ":", "No-", table(tmpDf2$FamHis)[1], ";Yes-", table(tmpDf2$FamHis)[2]))
    print(paste0(i, "_DRE_", j, ":", "Negative-", table(tmpDf2$DRE)[1], ";Positive-", table(tmpDf2$DRE)[2]))
    print(paste0(i, "_AGE_", j, ":", "Median(Range) ", median(as.numeric(tmpDf2$Age)), " ",
                 min(as.numeric(tmpDf2$Age)), " ", max(as.numeric(tmpDf2$Age))))
  }
}


#grid.arrange(GG_0, GG_1, GG_2, GG_3, GG_4, GG_5, ncol = 2)
grid.arrange(low, high, ncol = 2)


table(allVarsPcpt2$gradeHolder)
grouping_1 <- melt(allVarsPcpt2[,c(2, 7, 9,10)], id.vars = "gradeHolder")
group_1plot <- ggplot(grouping_1, aes(x=variable, y=value)) + geom_boxplot(aes(fill = factor(gradeHolder))) + ggtitle(paste0("BGG2vsGG3GG5","(n = ", 54,";", 86, ")"))
group_1plot <- group_1plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + ylim(c(0,25))
group_1plot


