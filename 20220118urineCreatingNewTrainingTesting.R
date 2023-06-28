#lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

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
  listOfNames <- listOfNames[-which(is.na(listOfNames))[1]]
  listOfPSA <- tmpXls$`Pre-Bx PSA Value`
  listOfPSA <- listOfPSA[-which(is.na(listOfPSA))[1]]
  tmpDf <- cbind(listOfNames, listOfPSA)
}


rmBad <- function(dat, badSample){
  res <- dat[-which(dat %in% str_remove(badSample, "UR.*"))]
  if (length(res) == 0) {
    res <- dat
  }
  return(res)
}

# seqStats <- read.csv("/mnt/DATA4/kevhu/urineRNA/20200421_All urine samps pheno.csv",
#                      stringsAsFactors = FALSE)
# 
# urineData <- read.csv("/mnt/DATA4/kevhu/urineRNA/20200421_All urine samps hdat.csv", 
#                       stringsAsFactors = FALSE)


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


### new model is created to have low-grade be benign to 3-4 and 4-3 + is high-grade
# path <- "/mnt/DATA4/kevhu/urineRNA/DOD_Aim 1 Cohort_Urine Sample List_Random Lists_Benign_3+3_3+4_4+3_4+4_9 and 10WITH identifiers_FINAL.xlsx"
# path <- "/mnt/DATA6/kevin_recovery/urineRNA/DOD_Aim 1 Cohort_Urine Sample List_Random Lists_Benign_3+3_3+4_4+3_4+4_9 and 10WITH identifiers_FINAL.xlsx"
path <- "/mnt/DATA5/tmp/kev/misc/DOD_Aim 1 Cohort_Urine Sample List_Random Lists_Benign_3+3_3+4_4+3_4+4_9 and 10WITH identifiers_FINAL_with updates RSM_12.10.2021.xlsx"


urine_benign <- getUrineIds(path, 1)
urine_3_3 <- getUrineIds(path, 2)
urine_3_4 <- getUrineIds(path, 3)
urine_4_3 <- getUrineIds(path, 4)
urine_4_4 <- getUrineIds(path, 5)
urine_9_10 <- getUrineIds(path, 6)


urine_benign <- rmBad(urine_benign, badSamples)
urine_3_3 <- rmBad(urine_3_3, badSamples)
urine_3_4 <- rmBad(urine_3_4, badSamples)
urine_4_3 <- rmBad(urine_4_3, badSamples)
urine_4_4 <- rmBad(urine_4_4, badSamples)
urine_9_10 <- rmBad(urine_9_10, badSamples)
urine_9_10 <- urine_9_10[-which(duplicated(urine_9_10))]


newGgTable <- NULL
for (i in 1:6) {
  tmpXls <- readxl::read_xlsx(path = as.character(path), sheet = as.numeric(i))
  listOfNames <- str_remove(string = tmpXls$LABEL, pattern = "CCGC_")
  listOfNames <- str_remove(listOfNames, pattern = "_.")
  listOfNames <- listOfNames[-which(is.na(listOfNames))]
  
  newGg <- tmpXls$GG_New
  newGg  <- newGg[-which(is.na(newGg))[1]]
  
  tmpTable <- data.frame("sample" = listOfNames, "GG" = newGg)
  newGgTable <- rbind(newGgTable, tmpTable)
}


allUrineSamples <- c(urine_benign, urine_3_3, urine_3_4,
                         urine_4_3, urine_4_4, urine_9_10)
# check newGgTable$sample[match(allUrineSamples, newGgTable$sample)]
allUrineGrades <-  newGgTable$GG[match(allUrineSamples, newGgTable$sample)]

allUrineSamplesDf <- data.frame("samples" = allUrineSamples,
                                "grades" = allUrineGrades)


### combined all of them together and then assign their designated Grade groups
## split after back into their repspective groups



### i think it'd be better to randomly pick from this by grade group ... then combine into high and low-grade - so we know it's balanced in
### in that regard too
### below don't run because it's made already


# lowGrade_benign_g1 <- c(urine_benign, urine_3_3)
# lowGrade_benign_g2 <- c(urine_benign, urine_3_3, urine_3_4)
# highGrade_benign_g1 <- c(urine_3_4, urine_4_3, urine_4_4, urine_9_10)
# highGrade_benign_g2 <- c(urine_4_3, urine_4_4, urine_9_10)


urine_benign <- allUrineSamplesDf$samples[which(allUrineSamplesDf$grades == 0)]
urine_3_3 <-  allUrineSamplesDf$samples[which(allUrineSamplesDf$grades == 1)]
urine_3_4 <-  allUrineSamplesDf$samples[which(allUrineSamplesDf$grades == 2)]
urine_4_3 <-  allUrineSamplesDf$samples[which(allUrineSamplesDf$grades == 3)]
urine_4_4 <-  allUrineSamplesDf$samples[which(allUrineSamplesDf$grades == 4)]
urine_9_10 <- allUrineSamplesDf$samples[which(allUrineSamplesDf$grades == 5)]



lowGrade_benign_g1 <- c(urine_benign, urine_3_3)
lowGrade_benign_g2 <- c(urine_benign, urine_3_3, urine_3_4)
lowGrade_benign_g3 <- c(urine_benign)
highGrade_benign_g1 <- c(urine_3_4, urine_4_3, urine_4_4, urine_9_10)
highGrade_benign_g2 <- c(urine_4_3, urine_4_4, urine_9_10)
highGrade_benign_g3 <- c(urine_3_3, urine_3_4, urine_4_3, urine_4_4, urine_9_10)


set.seed(1003)
urine_benign_test <- sample(urine_benign, length(urine_benign)/3)
set.seed(1003)
urine_3_3_test <- sample(urine_3_3, length(urine_3_3)/3)
set.seed(1003)
urine_3_4_test <- sample(urine_3_4, length(urine_3_4)/3)
set.seed(1003)
urine_4_3_test <- sample(urine_4_3, length(urine_4_3)/3)
set.seed(1003)
urine_4_4_test <- sample(urine_4_4, length(urine_4_4)/3)
set.seed(1003)
urine_9_10_test <- sample(urine_9_10, length(urine_9_10)/3)


lowGrade_test_benign_g1 <- c(urine_benign_test, urine_3_3_test)
lowGrade_test_benign_g2 <- c(urine_benign_test, urine_3_3_test, urine_3_4_test)
lowGrade_test_benign_g3 <- c(urine_benign_test)
lowGrade_training_benign_g1 <- lowGrade_benign_g1[-which(lowGrade_benign_g1 %in% lowGrade_test_benign_g1)]
lowGrade_training_benign_g2 <- lowGrade_benign_g2[-which(lowGrade_benign_g2 %in% lowGrade_test_benign_g2)]
lowGrade_training_benign_g3 <- lowGrade_benign_g3[-which(lowGrade_benign_g3 %in% lowGrade_test_benign_g3)]


highGrade_test_benign_g1 <- c(urine_3_4_test, urine_4_3_test, urine_4_4_test, urine_9_10_test)
highGrade_test_benign_g2 <- c(urine_4_3_test, urine_4_4_test, urine_9_10_test)
highGrade_test_benign_g3 <- c(urine_3_3_test, urine_3_4_test,urine_4_3_test, urine_4_4_test, urine_9_10_test)
highGrade_training_benign_g1 <- highGrade_benign_g1[-which(highGrade_benign_g1 %in% highGrade_test_benign_g1)]
highGrade_training_benign_g2 <- highGrade_benign_g2[-which(highGrade_benign_g2 %in% highGrade_test_benign_g2)]
highGrade_training_benign_g3 <- highGrade_benign_g3[-which(highGrade_benign_g3 %in% highGrade_test_benign_g3)]


trainingDf_benign_g1 <- rbind(data.frame("Sample" =  lowGrade_training_benign_g1, "Grade" = "low",stringsAsFactors = FALSE),
                              data.frame("Sample" = highGrade_training_benign_g1, "Grade" = "high", stringsAsFactors = FALSE))

testingDf_benign_g1 <- rbind(data.frame("Sample" =  lowGrade_test_benign_g1, "Grade" = "low",stringsAsFactors = FALSE),
                             data.frame("Sample" = highGrade_test_benign_g1, "Grade" = "high", stringsAsFactors = FALSE))

trainingDf_benign_g2 <- rbind(data.frame("Sample" =  lowGrade_training_benign_g2, "Grade" = "low",stringsAsFactors = FALSE),
                              data.frame("Sample" = highGrade_training_benign_g2, "Grade" = "high", stringsAsFactors = FALSE))

testingDf_benign_g2 <- rbind(data.frame("Sample" =  lowGrade_test_benign_g2, "Grade" = "low",stringsAsFactors = FALSE),
                             data.frame("Sample" = highGrade_test_benign_g2, "Grade" = "high", stringsAsFactors = FALSE))

trainingDf_benign_g3 <- rbind(data.frame("Sample" =  lowGrade_training_benign_g3, "Grade" = "low",stringsAsFactors = FALSE),
                              data.frame("Sample" = highGrade_training_benign_g3, "Grade" = "high", stringsAsFactors = FALSE))

testingDf_benign_g3 <- rbind(data.frame("Sample" =  lowGrade_test_benign_g3, "Grade" = "low",stringsAsFactors = FALSE),
                             data.frame("Sample" = highGrade_test_benign_g3, "Grade" = "high", stringsAsFactors = FALSE))



# write.table(trainingDf_benign_g1, "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118trainingDf_benign_g1.txt", sep = "\t",
#             col.names = TRUE, quote = FALSE, row.names = TRUE)
# 
# write.table(trainingDf_benign_g2, "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118trainingDf_benign_g2.txt", sep = "\t",
#             col.names = TRUE, quote = FALSE, row.names = TRUE)
#
# write.table(trainingDf_benign_g3, "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118trainingDf_benign_g3.txt", sep = "\t",
#             col.names = TRUE, quote = FALSE, row.names = TRUE)
# 
# 
# write.table(testingDf_benign_g1, "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118testingDf_benign_g1.txt", sep = "\t",
#             col.names = TRUE, quote = FALSE, row.names = TRUE)
# 
# write.table(testingDf_benign_g2, "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118testingDf_benign_g2.txt", sep = "\t",
#             col.names = TRUE, quote = FALSE, row.names = TRUE)
#
# write.table(testingDf_benign_g3, "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118testingDf_benign_g3.txt", sep = "\t",
#             col.names = TRUE, quote = FALSE, row.names = TRUE)
# 
### below are the new models


trainingDf_benign_g1 <- read.table("/mnt/DATA5/tmp/kev/2022urineSimpa/20220118trainingDf_benign_g1.txt",
                                   sep = "\t", stringsAsFactors = FALSE, header = TRUE)
trainingDf_benign_g1$Grade <- factor(trainingDf_benign_g1$Grade, levels = c("low", "high"))

testingDf_benign_g1 <- read.table("/mnt/DATA5/tmp/kev/2022urineSimpa/20220118testingDf_benign_g1.txt",
                                  sep = "\t", stringsAsFactors = FALSE, header = TRUE)
testingDf_benign_g1$Grade <- factor(testingDf_benign_g1$Grade, levels = c("low","high"))



trainingDf_benign_g2 <- read.table("/mnt/DATA5/tmp/kev/2022urineSimpa/20220118trainingDf_benign_g2.txt",
                                   sep = "\t", stringsAsFactors = FALSE, header = TRUE)
trainingDf_benign_g2$Grade <- factor(trainingDf_benign_g2$Grade, levels = c("low", "high"))


testingDf_benign_g2 <- read.table("/mnt/DATA5/tmp/kev/2022urineSimpa/20220118testingDf_benign_g2.txt",
                                  sep = "\t", stringsAsFactors = FALSE, header = TRUE)
testingDf_benign_g2$Grade <- factor(testingDf_benign_g2$Grade, levels = c("low","high"))


trainingDf_benign_g3 <- read.table("/mnt/DATA5/tmp/kev/2022urineSimpa/20220118trainingDf_benign_g3.txt",
                                   sep = "\t", stringsAsFactors = FALSE, header = TRUE)
trainingDf_benign_g3$Grade <- factor(trainingDf_benign_g3$Grade, levels = c("low", "high"))


testingDf_benign_g3 <- read.table("/mnt/DATA5/tmp/kev/2022urineSimpa/20220118testingDf_benign_g3.txt",
                                  sep = "\t", stringsAsFactors = FALSE, header = TRUE)
testingDf_benign_g3$Grade <- factor(testingDf_benign_g3$Grade, levels = c("low","high"))




reseqedSamps <- str_remove_all(rownames(urineData_processed_norm)[grep("\\.1", rownames(urineData_processed_norm))], "\\.1")
urineData_processed_norm <- urineData_processed_norm[-which(rownames(urineData_processed_norm) %in% reseqedSamps),]
tmpVar <- str_remove(str_remove(rownames(urineData_processed_norm), "UR.*"), "\\.1")
### vsurf portion

bengin_g1_rf_df <- data.frame(urineData_processed_norm[match(as.character(trainingDf_benign_g1$Sample),tmpVar),],
                              stringsAsFactors = FALSE)
### sanity check
str_remove_all(str_remove_all(rownames(bengin_g1_rf_df), "\\.1"), "UR.*")

ids_benign_g1 <- factor(trainingDf_benign_g1$Grade)

bengin_g2_rf_df <- data.frame(urineData_processed_norm[match(as.character(trainingDf_benign_g2$Sample),tmpVar),],
                              stringsAsFactors = FALSE)
str_remove_all(str_remove_all(rownames(bengin_g2_rf_df), "\\.1"), "UR.*")

ids_benign_g2 <- factor(trainingDf_benign_g2$Grade)


bengin_g3_rf_df <- data.frame(urineData_processed_norm[match(as.character(trainingDf_benign_g3$Sample),tmpVar),],
                              stringsAsFactors = FALSE)
str_remove_all(str_remove_all(rownames(bengin_g3_rf_df), "\\.1"), "UR.*")

ids_benign_g3 <- factor(trainingDf_benign_g3$Grade)



### quickly removing first pass samples. this is all assuming the ".1" samples ares the ones sequenced first


set.seed(182, kind = "L'Ecuyer-CMRG")
vsurf_parallel_benign_g1 <- VSURF(bengin_g1_rf_df, ids_benign_g1, parallel = TRUE, ncores = 24, clusterType = "FORK")
summary(vsurf_parallel_benign_g1)
plot(vsurf_parallel_benign_g1)

genesThres_benign_g1 <- colnames(bengin_g1_rf_df)[vsurf_parallel_benign_g1$varselect.thres]
genesInterp_benign_g1 <- colnames(bengin_g1_rf_df)[vsurf_parallel_benign_g1$varselect.interp]




set.seed(82123, kind = "L'Ecuyer-CMRG")
vsurf_parallel_benign_g2 <- VSURF(bengin_g2_rf_df, ids_benign_g2, parallel = TRUE, ncores = 24, clusterType = "FORK")
summary(vsurf_parallel_benign_g2)
plot(vsurf_parallel_benign_g2)

genesThres_benign_g2 <- colnames(bengin_g2_rf_df)[vsurf_parallel_benign_g2$varselect.thres]
genesInterp_benign_g2 <- colnames(bengin_g2_rf_df)[vsurf_parallel_benign_g2$varselect.interp]



set.seed(99993, kind = "L'Ecuyer-CMRG")
vsurf_parallel_benign_g3 <- VSURF(bengin_g3_rf_df, ids_benign_g3, parallel = TRUE, ncores = 24, clusterType = "FORK")
summary(vsurf_parallel_benign_g3)
plot(vsurf_parallel_benign_g3)

genesThres_benign_g3 <- colnames(bengin_g3_rf_df)[vsurf_parallel_benign_g3$varselect.thres]
genesInterp_benign_g3 <- colnames(bengin_g3_rf_df)[vsurf_parallel_benign_g3$varselect.interp]





### setting up the models

ids_b_g1 <- trainingDf_benign_g1$Grade
ids_b_g2 <- trainingDf_benign_g2$Grade
ids_b_g3 <- trainingDf_benign_g3$Grade

### 2022 testing to see if i use the inter RF genes

bengin_g1_rf_df_genes <- bengin_g1_rf_df[, genesThres_benign_g1]
bengin_g2_rf_df_genes <- bengin_g2_rf_df[, genesThres_benign_g2]
bengin_g3_rf_df_genes <- bengin_g3_rf_df[, genesThres_benign_g3]


bengin_g1_rf_df_genes <- cbind(ids_b_g1, bengin_g1_rf_df_genes)
bengin_g2_rf_df_genes <- cbind(ids_b_g2, bengin_g2_rf_df_genes)
bengin_g3_rf_df_genes <- cbind(ids_b_g3, bengin_g3_rf_df_genes)

bengin_g1_rf_df_genes$ids_b_g1 <- relevel(factor(bengin_g1_rf_df_genes$ids_b_g1), ref = "low")
bengin_g2_rf_df_genes$ids_b_g2 <- relevel(factor(bengin_g2_rf_df_genes$ids_b_g2), ref = "low")
bengin_g3_rf_df_genes$ids_b_g3 <- relevel(factor(bengin_g3_rf_df_genes$ids_b_g3), ref = "low")


model_weights_g1 <- ifelse(trainingDf_benign_g1$Grade == "high",
                           (1/table(trainingDf_benign_g1$Grade)[1]) * 0.5,
                           (1/table(trainingDf_benign_g1$Grade)[2]) * 0.5)
model_weights_g2 <- ifelse(trainingDf_benign_g2$Grade == "high",
                           (1/table(trainingDf_benign_g2$Grade)[1]) * 0.5,
                           (1/table(trainingDf_benign_g2$Grade)[2]) * 0.5)
model_weights_g3 <- ifelse(trainingDf_benign_g3$Grade == "high",
                           (1/table(trainingDf_benign_g3$Grade)[1]) * 0.5,
                           (1/table(trainingDf_benign_g3$Grade)[2]) * 0.5)


train_control_rules <- trainControl(method = "cv", number = 5, savePredictions = TRUE, classProbs = TRUE,
                                    summaryFunction = twoClassSummary)




set.seed(54763, kind = "L'Ecuyer-CMRG")
fiveFoldModel_g1 <- caret::train(ids_b_g1~., data = bengin_g1_rf_df_genes, trControl = train_control_rules,
                                 method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g1)
# save(fiveFoldModel_g1, file = "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118oldNormfiveFoldModel_b_g1.Robj")

set.seed(6654763, kind = "L'Ecuyer-CMRG")
fiveFoldModel_g2 <- caret::train(ids_b_g2 ~., data = bengin_g2_rf_df_genes, trControl = train_control_rules,
                                 method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g2)
# save(fiveFoldModel_g2, file = "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118oldNormFiveFoldModel_b_g2.Robj")


set.seed(66878, kind = "L'Ecuyer-CMRG")
fiveFoldModel_g3 <- caret::train(ids_b_g3 ~., data = bengin_g3_rf_df_genes, trControl = train_control_rules,
                                 method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g3)
# save(fiveFoldModel_g3, file = "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118oldNormFiveFoldModel_b_g3.Robj")


genes_b_g1_weights <- as.matrix(coef(fiveFoldModel_g1$finalModel, fiveFoldModel_g1$bestTune$lambda))
genes_b_g2_weights <- as.matrix(coef(fiveFoldModel_g2$finalModel, fiveFoldModel_g2$bestTune$lambda))
genes_b_g3_weights <- as.matrix(coef(fiveFoldModel_g3$finalModel, fiveFoldModel_g3$bestTune$lambda))



### looking at the testing data after
###
###


load("/mnt/DATA5/tmp/kev/2022urineSimpa/20220118oldNormfiveFoldModel_b_g1.Robj")
load("/mnt/DATA5/tmp/kev/2022urineSimpa/20220118oldNormFiveFoldModel_b_g2.Robj")
load("/mnt/DATA5/tmp/kev/2022urineSimpa/20220118oldNormFiveFoldModel_b_g3.Robj")


genesThres_benign_g1 <- fiveFoldModel_g1$coefnames
genesThres_benign_g2 <- fiveFoldModel_g2$coefnames
genesThres_benign_g3 <- fiveFoldModel_g3$coefnames



ids_b_g1 <- trainingDf_benign_g1$Grade
ids_b_g2 <- trainingDf_benign_g2$Grade
ids_b_g3 <- trainingDf_benign_g3$Grade

bengin_g1_rf_df_genes <- bengin_g1_rf_df[, genesThres_benign_g1]
bengin_g2_rf_df_genes <- bengin_g2_rf_df[, genesThres_benign_g2]
bengin_g3_rf_df_genes <- bengin_g3_rf_df[, genesThres_benign_g3]


bengin_g1_rf_df_genes <- cbind(ids_b_g1, bengin_g1_rf_df_genes)
bengin_g2_rf_df_genes <- cbind(ids_b_g2, bengin_g2_rf_df_genes)
bengin_g3_rf_df_genes <- cbind(ids_b_g3, bengin_g3_rf_df_genes)


testing_b_g1_df <- data.frame(urineData_processed_norm[match(as.character(testingDf_benign_g1$Sample),tmpVar),],
                              stringsAsFactors = FALSE)
testing_b_g1_df <- testing_b_g1_df[, genesThres_benign_g1]

### sanity check - works
#str_remove_all(str_remove_all(rownames(testing_b_g1_df), "\\.1"), "UR.*")

testing_b_g2_df <- data.frame(urineData_processed_norm[match(as.character(testingDf_benign_g2$Sample),tmpVar),],
                              stringsAsFactors = FALSE)
testing_b_g2_df <- testing_b_g2_df[, genesThres_benign_g2]

testing_b_g3_df <- data.frame(urineData_processed_norm[match(as.character(testingDf_benign_g3$Sample),tmpVar),],
                              stringsAsFactors = FALSE)
testing_b_g3_df <- testing_b_g3_df[, genesThres_benign_g3]



predRes_b_g1 <- predict(fiveFoldModel_g1, type = c("prob"),newdata = testing_b_g1_df)
rocGraph_b_g1 <- roc(testingDf_benign_g1$Grade ~ predRes_b_g1$high)
predRes_b_g2 <- predict(fiveFoldModel_g2, type = c("prob"), newdata = testing_b_g2_df)
rocGraph_b_g2 <- roc(testingDf_benign_g2$Grade ~ predRes_b_g2$high)
predRes_b_g3 <- predict(fiveFoldModel_g3, type = c("prob"), newdata = testing_b_g3_df)
rocGraph_b_g3 <- roc(testingDf_benign_g3$Grade ~ predRes_b_g3$high)


training_b_g1_probs <- fiveFoldModel_g1$pred
training_b_g1_probs <- training_b_g1_probs[which(training_b_g1_probs$lambda == fiveFoldModel_g1$bestTune$lambda &
                                                   training_b_g1_probs$alpha == fiveFoldModel_g1$bestTune$alpha),]
roc_train_b_g1 <- plot.roc(pROC::roc(training_b_g1_probs$obs ~ training_b_g1_probs$high))

training_b_g2_probs <- fiveFoldModel_g2$pred
training_b_g2_probs <- training_b_g2_probs[which(training_b_g2_probs$lambda == fiveFoldModel_g2$bestTune$lambda &
                                                   training_b_g2_probs$alpha == fiveFoldModel_g2$bestTune$alpha),]
roc_train_b_g2 <- plot.roc(pROC::roc(training_b_g2_probs$obs ~ training_b_g2_probs$high))

training_b_g3_probs <- fiveFoldModel_g3$pred
training_b_g3_probs <- training_b_g3_probs[which(training_b_g3_probs$lambda == fiveFoldModel_g3$bestTune$lambda &
                                                   training_b_g3_probs$alpha == fiveFoldModel_g3$bestTune$alpha),]
roc_train_b_g3 <- plot.roc(pROC::roc(training_b_g3_probs$obs ~ training_b_g3_probs$high))


### creating other models

path <- "/mnt/DATA5/tmp/kev/misc/DOD_Aim 1 Cohort_Urine Sample List_Random Lists_Benign_3+3_3+4_4+3_4+4_9 and 10WITH identifiers_FINAL_with updates RSM_12.10.2021.xlsx"

urine_benign <- getUrineIds_PSA(path, 1)
urine_3_3 <- getUrineIds_PSA(path, 2)
urine_3_4 <- getUrineIds_PSA(path, 3)
urine_4_3 <- getUrineIds_PSA(path, 4)
urine_4_4 <- getUrineIds_PSA(path, 5)
urine_9_10 <- getUrineIds_PSA(path, 6)



newGgTable <- NULL
for (i in 1:6) {
  tmpXls <- readxl::read_xlsx(path = as.character(path), sheet = as.numeric(i))
  listOfNames <- str_remove(string = tmpXls$LABEL, pattern = "CCGC_")
  listOfNames <- str_remove(listOfNames, pattern = "_.")
  listOfNames <- listOfNames[-which(is.na(listOfNames))]
  
  newGg <- tmpXls$GG_New
  newGg  <- newGg[-which(is.na(newGg))[1]]
  
  tmpTable <- data.frame("sample" = listOfNames, "GG" = newGg)
  newGgTable <- rbind(newGgTable, tmpTable)
}


allUrineSamples <- c(urine_benign[,1], urine_3_3[,1], urine_3_4[,1],
                     urine_4_3[,1], urine_4_4[,1], urine_9_10[,1])
allUrinePSA <- c(urine_benign[,2], urine_3_3[,2], urine_3_4[,2],
                 urine_4_3[,2], urine_4_4[,2], urine_9_10[,2])
# check newGgTable$sample[match(allUrineSamples, newGgTable$sample)]
allUrineGrades <-  newGgTable$GG[match(allUrineSamples, newGgTable$sample)]

allUrineSamplesDf <- data.frame("samples" = allUrineSamples,
                                "PSA" = allUrinePSA,
                                "grades" = allUrineGrades)


urine_benign <- allUrineSamplesDf[which(allUrineSamplesDf$grades == 0), 1:2]
urine_3_3 <-  allUrineSamplesDf[which(allUrineSamplesDf$grades == 1), 1:2]
urine_3_4 <-  allUrineSamplesDf[which(allUrineSamplesDf$grades == 2), 1:2]
urine_4_3 <-  allUrineSamplesDf[which(allUrineSamplesDf$grades == 3), 1:2]
urine_4_4 <-  allUrineSamplesDf[which(allUrineSamplesDf$grades == 4), 1:2]
urine_9_10 <- allUrineSamplesDf[which(allUrineSamplesDf$grades == 5), 1:2]



g1PsaMat_low <- data.frame(rbind(urine_benign,urine_3_3), stringsAsFactors = FALSE)
g1PsaMat_low$grade <- "low"
g1PsaMat_high <- data.frame(rbind(urine_3_4, urine_4_3, urine_4_4, urine_9_10), stringsAsFactors = FALSE)
g1PsaMat_high$grade <- "high"


g2PsaMat_low <- data.frame(rbind(urine_benign,urine_3_3, urine_3_4), stringsAsFactors = FALSE)
g2PsaMat_low$grade <- "low"
g2PsaMat_high <- data.frame(rbind(urine_4_3, urine_4_4, urine_9_10), stringsAsFactors = FALSE)
g2PsaMat_high$grade <- "high"

g3PsaMat_low <- data.frame(rbind(urine_benign), stringsAsFactors = FALSE)
g3PsaMat_low$grade <- "low"
g3PsaMat_high <- data.frame(rbind(urine_3_3, urine_3_4, urine_4_3, urine_4_4, urine_9_10), stringsAsFactors = FALSE)
g3PsaMat_high$grade <- "high"


g1PsaMat_all <- rbind(g1PsaMat_low, g1PsaMat_high)
g2PsaMat_all <- rbind(g2PsaMat_low, g2PsaMat_high)
g3PsaMat_all <- rbind(g3PsaMat_low, g3PsaMat_high)


colnames(g1PsaMat_all)[1:2] <- c("listOfNames", "listOfPSA")
colnames(g2PsaMat_all)[1:2] <- c("listOfNames", "listOfPSA")
colnames(g3PsaMat_all)[1:2] <- c("listOfNames", "listOfPSA")


g1PsaMat_all$listOfPSA <- as.numeric(g1PsaMat_all$listOfPSA)
g2PsaMat_all$listOfPSA <- as.numeric(g2PsaMat_all$listOfPSA)
g3PsaMat_all$listOfPSA <- as.numeric(g3PsaMat_all$listOfPSA)


g1Psa_training <- g1PsaMat_all[match(trainingDf_benign_g1$Sample, g1PsaMat_all$listOfNames),]
rownames(g1Psa_training) <- g1Psa_training$listOfNames
g1Psa_training <- g1Psa_training[,2:3]
colnames(g1Psa_training) <- c("PSA", "grade")
g1Psa_training$ones <- 1


g2Psa_training <- g2PsaMat_all[match(trainingDf_benign_g2$Sample, g2PsaMat_all$listOfNames),]
rownames(g2Psa_training) <- g2Psa_training$listOfNames
g2Psa_training <- g2Psa_training[,2:3]
colnames(g2Psa_training) <- c("PSA", "grade")
g2Psa_training$ones <- 1

g3Psa_training <- g3PsaMat_all[match(trainingDf_benign_g3$Sample, g3PsaMat_all$listOfNames),]
rownames(g3Psa_training) <- g3Psa_training$listOfNames
g3Psa_training <- g3Psa_training[,2:3]
colnames(g3Psa_training) <- c("PSA", "grade")
g3Psa_training$ones <- 1


g1Psa_testing <- g1PsaMat_all[match(testingDf_benign_g1$Sample, g1PsaMat_all$listOfNames),]
rownames(g1Psa_testing) <- g1Psa_testing$listOfNames
g1Psa_testing <- g1Psa_testing[,2:3]
colnames(g1Psa_testing) <- c("PSA", "grade")
g1Psa_testing$ones <- 1

g2Psa_testing <- g2PsaMat_all[match(testingDf_benign_g2$Sample, g2PsaMat_all$listOfNames),]
rownames(g2Psa_testing) <- g2Psa_testing$listOfNames
g2Psa_testing <- g2Psa_testing[,2:3]
colnames(g2Psa_testing) <- c("PSA", "grade")
g2Psa_testing$ones <- 1


g3Psa_testing <- g3PsaMat_all[match(testingDf_benign_g3$Sample, g3PsaMat_all$listOfNames),]
rownames(g3Psa_testing) <- g3Psa_testing$listOfNames
g3Psa_testing <- g3Psa_testing[,2:3]
colnames(g3Psa_testing) <- c("PSA", "grade")
g3Psa_testing$ones <- 1


model_weights_g1 <- ifelse(trainingDf_benign_g1$Grade == "high",
                           (1/table(trainingDf_benign_g1$Grade)[1]) * 0.5,
                           (1/table(trainingDf_benign_g1$Grade)[2]) * 0.5)


model_weights_g2 <- ifelse(trainingDf_benign_g2$Grade == "high",
                           (1/table(trainingDf_benign_g1$Grade)[1]) * 0.5,
                           (1/table(trainingDf_benign_g1$Grade)[2]) * 0.5)

model_weights_g3 <- ifelse(trainingDf_benign_g3$Grade == "high",
                           (1/table(trainingDf_benign_g1$Grade)[1]) * 0.5,
                           (1/table(trainingDf_benign_g1$Grade)[2]) * 0.5)



train_control_rules <- trainControl(method = "cv", number = 5, savePredictions = TRUE, classProbs = TRUE,
                                    summaryFunction = twoClassSummary)




set.seed(340280, kind = "L'Ecuyer-CMRG")
fiveFoldModel_g1_psa <- caret::train(grade~., data = g1Psa_training, trControl = train_control_rules,
                                     method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g1)


set.seed(45335, kind = "L'Ecuyer-CMRG")
fiveFoldModel_g2_psa <- caret::train(grade~., data = g2Psa_training, trControl = train_control_rules,
                                     method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g2)


set.seed(454, kind = "L'Ecuyer-CMRG")
fiveFoldModel_g3_psa <- caret::train(grade~., data = g3Psa_training, trControl = train_control_rules,
                                     method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g3)






mipsGenes <- c("TMPRSS2.ERG.T1E4.COSF125", "PCA3.E2E3")
mips_g1 <- bengin_g1_rf_df[,mipsGenes]
mips_g1 <- cbind(g1Psa_training[,1:2], mips_g1)

mips_g2 <- bengin_g2_rf_df[, mipsGenes]
mips_g2 <- cbind(g2Psa_training[,1:2], mips_g2)

mips_g3 <- bengin_g3_rf_df[, mipsGenes]
mips_g3 <- cbind(g3Psa_training[,1:2], mips_g3)


mips_g1_test <- data.frame(urineData_processed_norm[match(as.character(testingDf_benign_g1$Sample),tmpVar), mipsGenes],
                           stringsAsFactors = FALSE)
mips_g1_test <- cbind(g1Psa_testing[,1:2], mips_g1_test)

mips_g2_test <- data.frame(urineData_processed_norm[match(as.character(testingDf_benign_g2$Sample),tmpVar), mipsGenes],
                           stringsAsFactors = FALSE)
mips_g2_test <- cbind(g2Psa_testing[,1:2], mips_g2_test)

mips_g3_test <- data.frame(urineData_processed_norm[match(as.character(testingDf_benign_g3$Sample),tmpVar), mipsGenes],
                           stringsAsFactors = FALSE)
mips_g3_test <- cbind(g3Psa_testing[,1:2], mips_g3_test)




set.seed(345, kind = "L'Ecuyer-CMRG")
fiveFoldModel_g1_mips <- caret::train(grade~., data = mips_g1, trControl = train_control_rules,
                                      method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g1)


set.seed(254353, kind = "L'Ecuyer-CMRG")
fiveFoldModel_g2_mips <- caret::train(grade~., data = mips_g2, trControl = train_control_rules,
                                      method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g2)

set.seed(2353, kind = "L'Ecuyer-CMRG")
fiveFoldModel_g3_mips <- caret::train(grade~., data = mips_g3, trControl = train_control_rules,
                                      method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g3)






### graphs for PSA

predRes_b_g1_psa <- predict(fiveFoldModel_g1_psa, type = c("prob"),newdata = g1Psa_testing)
rocGraph_b_g1_psa <- roc(testingDf_benign_g1$Grade ~ predRes_b_g1_psa$high)

predRes_b_g2_psa <- predict(fiveFoldModel_g2_psa, type = c("prob"), newdata = g2Psa_testing)
rocGraph_b_g2_psa <- roc(testingDf_benign_g2$Grade ~ predRes_b_g2_psa$high)

training_b_g1_probs_psa <- fiveFoldModel_g1_psa$pred
training_b_g1_probs_psa <- training_b_g1_probs_psa[which(training_b_g1_probs_psa$lambda == fiveFoldModel_g1_psa$bestTune$lambda &
                                                           training_b_g1_probs_psa$alpha == fiveFoldModel_g1_psa$bestTune$alpha),]
roc_train_b_g1_psa <- plot.roc(pROC::roc(training_b_g1_probs_psa$obs ~ training_b_g1_probs$high))


training_b_g2_probs_psa <- fiveFoldModel_g2_psa$pred
training_b_g2_probs_psa <- training_b_g2_probs_psa[which(training_b_g2_probs_psa$lambda == fiveFoldModel_g2_psa$bestTune$lambda &
                                                           training_b_g2_probs_psa$alpha == fiveFoldModel_g2_psa$bestTune$alpha),]
roc_train_b_g2_psa <- plot.roc(pROC::roc(training_b_g2_probs_psa$obs ~ training_b_g2_probs$high))




### graphs for MIPs

predRes_b_g1_mips <- predict(fiveFoldModel_g1_mips, type = c("prob"),newdata = mips_g1_test)
rocGraph_b_g1_mips <- roc(testingDf_benign_g1$Grade ~ predRes_b_g1_mips$high)

predRes_b_g2_mips <- predict(fiveFoldModel_g2_mips, type = c("prob"), newdata = mips_g2_test)
rocGraph_b_g2_mips <- roc(testingDf_benign_g2$Grade ~ predRes_b_g2_mips$high)

predRes_b_g3_mips <- predict(fiveFoldModel_g3_mips, type = c("prob"), newdata = mips_g3_test)
rocGraph_b_g3_mips <- roc(testingDf_benign_g3$Grade ~ predRes_b_g3_mips$high)



training_b_g1_probs_mips <- fiveFoldModel_g1_mips$pred
training_b_g1_probs_mips <- training_b_g1_probs_mips[which(training_b_g1_probs_mips$lambda == fiveFoldModel_g1_mips$bestTune$lambda &
                                                             training_b_g1_probs_mips$alpha == fiveFoldModel_g1_mips$bestTune$alpha),]
roc_train_b_g1_mips <- plot.roc(pROC::roc(training_b_g1_probs_mips$obs ~ training_b_g1_probs$high))


training_b_g2_probs_mips <- fiveFoldModel_g2_mips$pred
training_b_g2_probs_mips <- training_b_g2_probs_mips[which(training_b_g2_probs_mips$lambda == fiveFoldModel_g2_mips$bestTune$lambda &
                                                             training_b_g2_probs_mips$alpha == fiveFoldModel_g2_mips$bestTune$alpha),]
roc_train_b_g2_mips <- plot.roc(pROC::roc(training_b_g2_probs_mips$obs ~ training_b_g2_probs$high))


training_b_g3_probs_mips <- fiveFoldModel_g3_mips$pred
training_b_g3_probs_mips <- training_b_g3_probs_mips[which(training_b_g3_probs_mips$lambda == fiveFoldModel_g3_mips$bestTune$lambda &
                                                             training_b_g3_probs_mips$alpha == fiveFoldModel_g3_mips$bestTune$alpha),]
roc_train_b_g3_mips <- plot.roc(pROC::roc(training_b_g3_probs_mips$obs ~ training_b_g3_probs$high))



ci.auc(roc_train_b_g1_psa, conf.level = 0.95, method = "delong")
ci.auc(rocGraph_b_g1_psa, conf.level = 0.95, method = "delong")
ci.auc(roc_train_b_g1_mips, conf.level = 0.95, method = "delong")
ci.auc(rocGraph_b_g1_mips, conf.level = 0.95, method = "delong")
ci.auc(roc_train_b_g1, conf.level = 0.95, method = "delong")
ci.auc(rocGraph_b_g1, conf.level = 0.95, method = "delong")



dev.off()
###testing graph - 
#pdf(file = "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118benignGG1_psa_mips.pdf", onefile = TRUE, useDingbats = TRUE, pointsize = 18, width = 10, height = 10)
plot(roc_train_b_g1_psa, col = alpha("red"), main = "ROCs: BGG1 vs GG2GG5", ylab = "", xlab="", lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g1_psa, col = alpha("red"),ylab = "Sensitivity", xlab="Specificity", lty = 2, asp = FALSE)
par(new = TRUE)
plot(roc_train_b_g1_mips, col = alpha("blue"),lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g1_mips, col = alpha("blue"), lty = 2, asp = FALSE)
par(new = TRUE)
plot(roc_train_b_g1, col = alpha("black"),  lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g1, col = alpha("black"),lty = 2, asp = FALSE)

legend("bottomright", legend = c("Training PSA - AUC:0.616 (0.490-0.730)", "Testing PSA - AUC: 0.763 (0.612-0.911)",
                                 "Training Mips - AUC:0.622 (0.507-0.737)", "Testing Mips - AUC: 0.634 (0.466-0.803)",
                                 "Training New Model - AUC: 0.616 (0.500-732)", "Testing New Model - AUC: 0.550 (0.374-0.725)"),
       col = c("red","red", "blue", "blue", "black", "black", "green"), lty=c(1,2,1,2,1,2), cex = .70)

dev.off()





ci.auc(roc_train_b_g2_psa, conf.level = 0.95, method = "delong")
ci.auc(rocGraph_b_g2_psa, conf.level = 0.95, method = "delong")
ci.auc(roc_train_b_g2_mips, conf.level = 0.95, method = "delong")
ci.auc(rocGraph_b_g2_mips, conf.level = 0.95, method = "delong")
ci.auc(roc_train_b_g2, conf.level = 0.95, method = "delong")
ci.auc(rocGraph_b_g2, conf.level = 0.95, method = "delong")





dev.off()
###testing graph - 
#pdf(file = "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118benignGG2_psa_mips.pdf", onefile = TRUE, useDingbats = TRUE, pointsize = 18, width = 10, height = 10)
plot(roc_train_b_g2_psa, col = alpha("red"), main = "ROCs: BGG2 vs GG3GG5", ylab = "", xlab="", lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g2_psa, col = alpha("red"),ylab = "Sensitivity", xlab="Specificity", lty = 2, asp = FALSE)
par(new = TRUE)
plot(roc_train_b_g2_mips, col = alpha("blue"),lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g2_mips, col = alpha("blue"), lty = 2, asp = FALSE)
par(new = TRUE)
plot(roc_train_b_g2, col = alpha("black"), lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g2, col = alpha("black"), lty = 2, asp = FALSE)

legend("bottomright", legend = c("Training PSA - AUC:0.683 (0.555-0.811)", "Testing PSA - AUC: 0.776 (0.626-0.926)",
                                 "Training Mips - AUC:0.643 (0.506-0.708)", "Testing Mips - AUC: 0.700 (0.532-0.867)",
                                 "Training New Model - AUC: 0.683 (0.555-811)", "Testing New Model - AUC: 0.466 (0.267-0.642)"),
       col = c("red","red", "blue", "blue", "black", "black", "green"), lty=c(1,2,1,2,1,2), cex = .70)

dev.off()


### PSA is such a strong signal in this dataset ... so I just added it to the genes
### do it with RF and just altogether without
###



bengin_g1_rf_df_psa <- bengin_g1_rf_df
bengin_g2_rf_df_psa <- bengin_g2_rf_df

### these filter are for minimum expression made retrospectively after boxplots below
### without  the gene expression filter

bengin_g1_rf_df_psa[bengin_g1_rf_df_psa < 5] <- 0
bengin_g1_rf_df_psa[bengin_g1_rf_df_psa < 5] <- 0

bengin_g1_rf_df_psa$PSA <- as.numeric(allUrineSamplesDf$PSA[match(str_remove(rownames(bengin_g1_rf_df_psa), "UR.*"), allUrineSamplesDf$samples)])
bengin_g2_rf_df_psa$PSA <- as.numeric(allUrineSamplesDf$PSA[match(str_remove(rownames(bengin_g2_rf_df_psa), "UR.*"), allUrineSamplesDf$samples)])

bengin_g1_rf_df_psa$PSA[which(bengin_g1_rf_df_psa$PSA < 10)] <- 0
bengin_g1_rf_df_psa$PSA[which(bengin_g1_rf_df_psa$PSA > 0)] <- 1
bengin_g1_rf_df_psa$PSA <- as.character(bengin_g1_rf_df_psa$PSA)

bengin_g2_rf_df_psa$PSA[which(bengin_g2_rf_df_psa$PSA < 10)] <- 0
bengin_g2_rf_df_psa$PSA[which(bengin_g2_rf_df_psa$PSA > 0)] <- 1
bengin_g2_rf_df_psa$PSA <- as.character(bengin_g2_rf_df_psa$PSA)



bengin_g1_df_genes_all <- bengin_g1_rf_df_psa
tmp <- 2^bengin_g1_df_genes_all[,grep("\\.ERG|\\.ETV|\\.ETS|\\.BRAF|\\.RAF1", colnames(bengin_g1_rf_df_psa))]
fusCombined <- log2(apply(tmp, 1, sum))
bengin_g1_df_genes_all <- bengin_g1_df_genes_all[,-grep("\\.ERG|\\.ETV|\\.ETS|\\.BRAF|\\.RAF1", colnames(bengin_g1_rf_df_psa))]
bengin_g1_df_genes_all$fusion <- fusCombined

bengin_g2_df_genes_all <- bengin_g2_rf_df_psa
tmp <- 2^bengin_g2_df_genes_all[,grep("\\.ERG|\\.ETV|\\.ETS|\\.BRAF|\\.RAF1", colnames(bengin_g2_rf_df_psa))]
fusCombined <- log2(apply(tmp, 1, sum))
bengin_g2_df_genes_all <- bengin_g2_df_genes_all[,-grep("\\.ERG|\\.ETV|\\.ETS|\\.BRAF|\\.RAF1", colnames(bengin_g2_rf_df_psa))]
bengin_g2_df_genes_all$fusion <- fusCombined

### new edition of creating cutoffs for genes and PSA. PSA should be 10


training_g1_low_df <- bengin_g1_df_genes_all[which(ids_benign_g1 == "low"),]
training_g1_high_df <- bengin_g1_df_genes_all[which(ids_benign_g1 == "high"),]
training_g1_low_df_m <- melt(training_g1_low_df)
training_g1_high_df_m <- melt(training_g1_high_df)
a <- ggplot(training_g1_low_df_m, aes(x = variable, y = value)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 16), limits = c(0,60)) + ggtitle("low-grade benign-gg1")
b <- ggplot(training_g1_high_df_m, aes(x = variable, y = value)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 16), limits = c(0,60)) + ggtitle("high-grade gg2-gg5")

training_g2_low_df <- bengin_g2_df_genes_all[which(ids_benign_g2 == "low"),]
training_g2_high_df <- bengin_g2_df_genes_all[which(ids_benign_g2 == "high"),]
training_g2_low_df_m <- melt(training_g2_low_df)
training_g2_high_df_m <- melt(training_g2_high_df)
c <- ggplot(training_g2_low_df_m, aes(x = variable, y = value)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 16), limits = c(0,60)) + ggtitle("low-grade benign-gg2")
d <- ggplot(training_g2_high_df_m, aes(x = variable, y = value)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 16), limits = c(0,60)) + ggtitle("high-grade gg3-gg5")

### using log2 cutoff of 5 for minimum expression 

gridExtra::grid.arrange(a,b,c,d, ncol = 2)
summary(training_g1_low_df$fusion)
summary(training_g1_high_df$fusion)
summary(training_g2_low_df$fusion)
summary(training_g2_high_df$fusion)


set.seed(808, kind = "L'Ecuyer-CMRG")
vsurf_parallel_benign_g1_all <- VSURF(bengin_g1_df_genes_all, ids_benign_g1, parallel = TRUE, ncores = 24, clusterType = "FORK")
summary(vsurf_parallel_benign_g1_all)
plot(vsurf_parallel_benign_g1_all)

genesThres_benign_g1_all <- colnames(bengin_g1_df_genes_all)[vsurf_parallel_benign_g1_all$varselect.thres]
genesInterp_benign_g1_all <- colnames(bengin_g1_df_genes_all)[vsurf_parallel_benign_g1_all$varselect.interp]


set.seed(512, kind = "L'Ecuyer-CMRG")
vsurf_parallel_benign_g2_all <- VSURF(bengin_g2_df_genes_all, ids_benign_g2, parallel = TRUE, ncores = 24, clusterType = "FORK")
summary(vsurf_parallel_benign_g2_all)
plot(vsurf_parallel_benign_g2_all)

genesThres_benign_g2_all <- colnames(bengin_g2_df_genes_all)[vsurf_parallel_benign_g2_all$varselect.thres]
genesInterp_benign_g2_all <- colnames(bengin_g2_df_genes_all)[vsurf_parallel_benign_g2_all$varselect.interp]


bengin_g1_df_genes_all <- bengin_g1_df_genes_all[,genesThres_benign_g1_all]
bengin_g1_df_genes_all <- cbind(ids_b_g1, bengin_g1_df_genes_all)
bengin_g1_df_genes_all$ids_b_g1 <- relevel(bengin_g1_df_genes_all$ids_b_g1, ref = "low")

bengin_g2_df_genes_all <- bengin_g2_df_genes_all[,genesThres_benign_g2_all]
bengin_g2_df_genes_all <- cbind(ids_b_g2, bengin_g2_df_genes_all)
bengin_g2_df_genes_all$ids_b_g2 <- relevel(bengin_g2_df_genes_all$ids_b_g2, ref = "low")


set.seed(89873, kind = "L'Ecuyer-CMRG")
fiveFoldModel_g1_all <- caret::train(ids_b_g1~., data = bengin_g1_df_genes_all, trControl = train_control_rules,
                                 method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g1)
# save(fiveFoldModel_g1_all, file = "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118oldNormfiveFoldModel_b_g1_all.Robj")

set.seed(993, kind = "L'Ecuyer-CMRG")
fiveFoldModel_g2_all <- caret::train(ids_b_g2 ~., data = bengin_g2_df_genes_all, trControl = train_control_rules,
                                 method = "glmnet", metric = "ROC", family = "binomial", weights = model_weights_g2)
# save(fiveFoldModel_g2_all, file = "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118oldNormFiveFoldModel_b_g2_all.Robj")


genes_b_g1_weights_all <- as.matrix(coef(fiveFoldModel_g1_all$finalModel, fiveFoldModel_g1_all$bestTune$lambda))
genes_b_g2_weights_all <- as.matrix(coef(fiveFoldModel_g2_all$finalModel, fiveFoldModel_g2_all$bestTune$lambda))



load("/mnt/DATA5/tmp/kev/2022urineSimpa/20220118oldNormfiveFoldModel_b_g1_all.Robj")
load("/mnt/DATA5/tmp/kev/2022urineSimpa/20220118oldNormFiveFoldModel_b_g2_all.Robj")

genesThres_benign_g1_all <- fiveFoldModel_g1_all$coefnames
genesThres_benign_g2_all <- fiveFoldModel_g2_all$coefnames
genesThres_benign_g2_all[12] <- "PSA"

ids_b_g1 <- trainingDf_benign_g1$Grade
ids_b_g2 <- trainingDf_benign_g2$Grade


bengin_g1_df_genes_all <- bengin_g1_rf_df_psa
bengin_g2_df_genes_all <- bengin_g2_rf_df_psa 

bengin_g1_df_genes_all <- cbind(ids_b_g1, bengin_g1_df_genes_all)
bengin_g2_df_genes_all <- cbind(ids_b_g2, bengin_g2_df_genes_all)


testing_b_g1_df_all <- data.frame(urineData_processed_norm[match(as.character(testingDf_benign_g1$Sample),tmpVar),],
                              stringsAsFactors = FALSE)
testing_b_g1_df_all[testing_b_g1_df_all < 5] <- 0
testing_b_g1_df_all$PSA <- as.numeric(allUrineSamplesDf$PSA[match(str_remove(rownames(testing_b_g1_df), "UR.*"), allUrineSamplesDf$samples)])
testing_b_g1_df_all <- testing_b_g1_df_all[, genesThres_benign_g1_all]
# tmp <- testing_b_g1_df_all[,grep("\\.ERG|\\.ETV|\\.ETS|\\.BRAF|\\.RAF1", colnames(testing_b_g1_df_all))]
# fusCombined <- apply(tmp, 1, sum)
# testing_b_g1_df_all <- testing_b_g1_df_all[,-grep("\\.ERG|\\.ETV|\\.ETS|\\.BRAF|\\.RAF1", colnames(testing_b_g1_df_all))]
# testing_b_g1_df_all$fusion <- fusCombined

testing_b_g2_df_all <- data.frame(urineData_processed_norm[match(as.character(testingDf_benign_g2$Sample),tmpVar),],
                              stringsAsFactors = FALSE)
testing_b_g2_df_all[testing_b_g2_df_all < 5] <- 0
testing_b_g2_df_all$PSA <- as.numeric(allUrineSamplesDf$PSA[match(str_remove(rownames(testing_b_g2_df_all), "UR.*"), allUrineSamplesDf$samples)])
tmp <- 2^testing_b_g2_df_all[,grep("\\.ERG|\\.ETV|\\.ETS|\\.BRAF|\\.RAF1", colnames(testing_b_g2_df_all))]
fusCombined <- log2(apply(tmp, 1, sum))
testing_b_g2_df_all <- testing_b_g2_df_all[,-grep("\\.ERG|\\.ETV|\\.ETS|\\.BRAF|\\.RAF1", colnames(testing_b_g2_df_all))]
testing_b_g2_df_all$fusion <- fusCombined
testing_b_g2_df_all <- testing_b_g2_df_all[, genesThres_benign_g2_all]
testing_b_g2_df_all$PSA[which(testing_b_g2_df_all$PSA < 10)] <- 0
testing_b_g2_df_all$PSA[which(testing_b_g2_df_all$PSA > 0)] <- 1
testing_b_g2_df_all$PSA <- as.character(testing_b_g2_df_all$PSA)





predRes_b_g1_all <- predict(fiveFoldModel_g1_all, type = c("prob"),newdata = testing_b_g1_df_all)
rocGraph_b_g1_all <- roc(testingDf_benign_g1$Grade ~ predRes_b_g1_all$high)
predRes_b_g2_all <- predict(fiveFoldModel_g2_all, type = c("prob"), newdata = testing_b_g2_df_all)
rocGraph_b_g2_all <- roc(testingDf_benign_g2$Grade ~ predRes_b_g2_all$high)

training_b_g1_probs_all <- fiveFoldModel_g1_all$pred
training_b_g1_probs_all <- training_b_g1_probs_all[which(training_b_g1_probs_all$lambda == fiveFoldModel_g1_all$bestTune$lambda &
                                                   training_b_g1_probs_all$alpha == fiveFoldModel_g1_all$bestTune$alpha),]
roc_train_b_g1_all <- plot.roc(pROC::roc(training_b_g1_probs_all$obs ~ training_b_g1_probs_all$high))
training_b_g2_probs_all <- fiveFoldModel_g2_all$pred
training_b_g2_probs_all <- training_b_g2_probs_all[which(training_b_g2_probs_all$lambda == fiveFoldModel_g2_all$bestTune$lambda &
                                                   training_b_g2_probs_all$alpha == fiveFoldModel_g2_all$bestTune$alpha),]
roc_train_b_g2_all <- plot.roc(pROC::roc(training_b_g2_probs_all$obs ~ training_b_g2_probs_all$high))


ci.auc(roc_train_b_g1_all, conf.level = 0.95, method = "delong")
ci.auc(rocGraph_b_g1_all, conf.level = 0.95, method = "delong")


dev.off()
###testing graph - 
pdf(file = "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118benignGG1_all_psa.pdf", onefile = TRUE, useDingbats = TRUE, pointsize = 18, width = 10, height = 10)
plot(roc_train_b_g1_all, col = alpha("black"), main = "ROCs: BGG1 vs GG2GG5 all panel genes + PSA", ylab = "", xlab="", lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g1_all, col = alpha("black"), ylab = "Sensitivity", xlab="Specificity", lty = 2, asp = FALSE)

legend("bottomright", legend = c("Training New Model - AUC: 0.652 (0.538-765)",
                                 "Testing New Model - AUC: 0.585 (0.411-0.759)"),
       col = c("black", "black"), lty=c(1,2), cex = .70)

dev.off()


ci.auc(roc_train_b_g2_all, conf.level = 0.95, method = "delong")
ci.auc(rocGraph_b_g2_all, conf.level = 0.95, method = "delong")

dev.off()
###testing graph - 
pdf(file = "/mnt/DATA5/tmp/kev/2022urineSimpa/20220118benignGG2_all_psa.pdf", onefile = TRUE, useDingbats = TRUE, pointsize = 18, width = 10, height = 10)
plot(roc_train_b_g2_all, col = alpha("black"),  main = "ROCs: BGG2 vs GG3GG5", ylab = "", xlab="", lty = 1, asp = FALSE)
par(new = TRUE)
plot(rocGraph_b_g2_all, col = alpha("black"), ylab = "Sensitivity", xlab="Specificity", lty = 2, asp = FALSE)

legend("bottomright", legend = c("Training New Model - AUC: 0.700 (0.584-0.814)",
                                 "Testing New Model - AUC: 0.620 (0.423-0.817)"),
       col = c("black", "black"), lty=c(1,2), cex = .70)

dev.off()

