lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)


library(VSURF)
library(parallel)
library(readxl)
library(ggplot2)
library(pROC)
library(reshape2)
library(gridExtra)
library(caret)
library(stringr)
library(Rcpp)

#load("/mnt/DATA4/kevhu/urineRNA/20180808urineCombinedModelGG2.Robj")
load("/mnt/DATA4/kevhu/urineRNA/20180808urineMipsModelGG2.Robj")
load("/mnt/DATA4/kevhu/urineRNA/20180808urineLogGeneModelGG2.Robj")
load("/mnt/DATA4/kevhu/urineRNA/20180808urinePsaModelGG2.Robj")


a <- read_xlsx("/mnt/DATA4/kevhu/urineRNA/AC17.2ExtDes2_mastermatrix_allsamples_byKLK3.xlsx")
b <- read_xlsx("/mnt/DATA4/kevhu/urineRNA/20189_07_26_AC17.2ExtDes2_mastermatrix.xlsx")

a.1 <- a[1:120,3:86]
rownames(a.1) <- as.vector(a$Contig_ID[1:120])
ids <- a$`Grade group bin`[1:120]

b <- b[match(rownames(a.1), b$Contig_ID),]
PSA <- b$LAB_RESNUM

set.seed(4321)
g1 <- which(ids == "Benign/GG1")
g2 <- which(ids == "GG2-5")

#set.seed(123)
g1.subset <- sample(g1, length(g1)/3)
#set.seed(123)
g2.subset <- sample(g2, length(g2)/3)
combined.subset <- c(g1.subset,g2.subset)

held.out <- a.1[combined.subset,]
a.2 <- a.1[-combined.subset,]

ids.2 <- factor(ids[-combined.subset])
ids.3 <- str_replace_all(ids.2, "GG2-5", "Cancer")
ids.3 <- str_replace_all(ids.3, "Benign/GG1", "NoCancer")


load("/mnt/DATA4/kevhu/urineRNA/20180807randoForestGenes.Robj")
#load("/mnt/DATA4/kevhu/urineRNA/20180808genesWGG226.Robj")


a.2.genes30 <- a.2[,genes30]
a.2.genes30 <- a.2.genes30[, order(colnames(a.2.genes30))]

a.2.log <- log2(a.1[-combined.subset,]+1)
a.2.log.30 <- cbind(ids.3,a.2.log[,genes30])

###logged data
### note the initial ordering does seem to matter, but may need ordering for case weights

a.7 <- a.2.log.30
a.7 <- a.7[order(a.7$ids.3),]
ids.4 <- ids.3
ids.3 <- ids.3[order(ids.3)]


train_control.combined <- trainControl(method = "cv", number = 5, savePredictions = TRUE,
                                       classProbs = TRUE, summaryFunction = twoClassSummary)
model_weights <- ifelse(a.7$ids.3 == "Cancer",
                        (1/table(a.7$ids.3)[1]) * 0.5,
                        (1/table(a.7$ids.3)[2]) * 0.5)

set.seed(4444)
fiveFoldModel.log <- train(ids.3 ~., data = a.7, trControl=train_control.combined, method="glmnet",
                           metric = "ROC", family = "binomial", weights = model_weights)
fiveFoldModel.log

####alpha = 0.55 and lambda = 0.04765286.

### makign auc for testing
testingPred.log <- fiveFoldModel.log$pred
testingPred.log <- testingPred.log[which(testingPred.log$lambda > 0.004),]
testingPred.log <- testingPred.log[which(testingPred.log$alpha == .55),]
plot.roc(pROC::roc(testingPred.log$obs ~ testingPred.log$Cancer))
### testing set for log

ids.held.out <- ids[combined.subset]
ids.held.out <- str_replace_all(ids.held.out, "GG2-5", "Cancer")
ids.held.out <- str_replace_all(ids.held.out, "Benign/GG1", "NoCancer")
ids.held.out <- as.factor(ids.held.out)
held.out.log <- log2(held.out[,genes30] + 1)

held.out.log2 <- cbind(ids.held.out, held.out.log)
colnames(held.out.log2)[1] <- "ids.3"
### don't need to order these like training when the weights were needed
#held.out.log <- held.out.log[order(ids.held.out),]
#ids.held.out <- ids.held.out[order(ids.held.out)]


prob.log <- predict(fiveFoldModel.log, newdata = held.out.log2)
confusionMatrix(data = factor(prob.log), ids.held.out, positive = "Cancer")
probROC.log <- predict(fiveFoldModel.log, type = c("prob"),newdata = held.out.log2)
rocGraph.log <- roc(ids.held.out ~ probROC.log$`Cancer`)
plot.roc(rocGraph.log, main = "PSA with model ROC")
rocGraph.log$auc

#save(fiveFoldModel.log, file = "/mnt/DATA4/kevhu/urineRNA/20180808urineLogGeneModelGG2.Robj")

#### PSA 
####
####

PSA.training <- PSA[-combined.subset]
PSA.training <- PSA.training[order(ids.4)]
a.5 <- data.frame(ids.3,PSA.training)
a.5$ones <- rep(1,nrow(a.5))


set.seed(4444)
fiveFoldModel.psa <- train(ids.3 ~., data = a.5, trControl=train_control.combined, method="glmnet",
                           metric = "ROC", family = "binomial", weights = model_weights)
fiveFoldModel.psa

### makign auc for testing
testingPred.psa <- fiveFoldModel.psa$pred
testingPred.psa <- testingPred.psa[which(testingPred.psa$lambda > 0.025),]
testingPred.psa <- testingPred.psa[which(testingPred.psa$alpha == 0.1),]
plot.roc(pROC::roc(testingPred.psa$obs ~ testingPred.psa$Cancer))
###

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

#save(fiveFoldModel.psa, file = "/mnt/DATA4/kevhu/urineRNA/20180808urinePsaModelGG2.Robj")



### combined
###
###

a.6 <- cbind(a.7[,1],PSA.training, a.7[,2:ncol(a.7)])
colnames(a.6)[1:2] <- c("ids.3","PSA")

set.seed(4444)
fiveFoldModel.combined <- train(ids.3 ~., data = a.6, trControl=train_control.combined, method="glmnet",
                                metric = "ROC", family = "binomial", weights = model_weights)
fiveFoldModel.combined

### making auc for testing
testingPred.combined <- fiveFoldModel.combined$pred
testingPred.combined <- testingPred.combined[which(testingPred.combined$lambda > 0.047),]
testingPred.combined <- testingPred.combined[which(testingPred.combined$alpha == 0.55),]
plot.roc(pROC::roc(testingPred.combined$obs ~ testingPred.combined$Cancer))
###

held.out.combined <- cbind(ids.held.out,held.out.30.psa[,2],held.out.log)
colnames(held.out.combined) <- colnames(a.6)
prob.combined <- predict(fiveFoldModel.combined, newdata = held.out.combined)
confusionMatrix(data = factor(prob.combined), ids.held.out, positive = "Cancer")
probROC.combined <- predict(fiveFoldModel.combined, type = c("prob"),newdata = held.out.combined)
rocGraph.combined <- roc(ids.held.out ~ probROC.combined$`Cancer`)
plot(rocGraph.combined)
rocGraph.combined$auc


#save(fiveFoldModel.combined, file = "/mnt/DATA4/kevhu/urineRNA/20180808urineCombinedModelGG2.Robj")

### Surrogate Mips
###
### ASK SCOTT IF THERE IS A WAY TO FIND THE FREQUENCY FOR THE DIFFERENT FUSIONS. Can then do a weighted mean for TMPRSS2 VAR

tmprss2.fusion.training <- log2(rowSums(a.2[,(which(grepl("TMPRSS2-ERG", colnames(a.2))))]) + 1)
tmprss2.fusion.training <- tmprss2.fusion.training[order(ids.4)]
PCA3 <- as.vector(log2(a.2[,which(grepl("PCA3", colnames(a.2)))] + 1))
PCA3 <- PCA3$PCA3.E2E3[order(ids.4)]
mips.genes <- data.frame(cbind(PSA.training, PCA3, tmprss2.fusion.training),stringsAsFactors = FALSE)
colnames(mips.genes) <- c("PSA","PCA3","TMPRSS2:ERG")
mips.genes <- cbind(ids.3, mips.genes) 

set.seed(4444)
fiveFoldModel.mips <- train(ids.3 ~., data = mips.genes, trControl=train_control.combined, method="glmnet",
                            metric = "ROC", family = "binomial", weights = model_weights)
fiveFoldModel.mips

### making auc for testing
testingPred.mips <- fiveFoldModel.mips$pred
testingPred.mips <- testingPred.mips[which(testingPred.mips$lambda > 0.003),]
testingPred.mips <- testingPred.mips[which(testingPred.mips$alpha == 0.1),]
plot.roc(pROC::roc(testingPred.mips$obs ~ testingPred.mips$Cancer))
###


held.out.tmrpss <- log2(rowSums(held.out[,(which(grepl("TMPRSS2-ERG", colnames(held.out))))])+1)
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


#save(fiveFoldModel.mips, file = "/mnt/DATA4/kevhu/urineRNA/20180808urineMipsModelGG2.Robj")


### creating overall graph with training and testing 
###
###
#plot(rocGraph.combined, main = "ROCs: Training n=73; Testing n=36", col = alpha("green",alpha = 0.5), xlim = c(1,0), ylim = c(0,1), lty = 2, asp = FALSE)
#par(new = TRUE)
#plot.roc(pROC::roc(testingPred.combined$obs ~ testingPred.combined$Cancer), col = alpha("green", alpha = 0.5), asp = FALSE)
#par(new = TRUE)

pdf(file = "/mnt/DATA4/kevhu/urineRNA/20180816urineVsGG2to5.pdf", onefile = TRUE, useDingbats = TRUE, pointsize = 15)
plot(rocGraph.psa, col = alpha("red", alpha = 0.5), main = NULL, ylab = "", xlab="", lty = 2, asp = FALSE)
par(new = TRUE)
plot.roc(pROC::roc(testingPred.psa$obs ~ testingPred.psa$Cancer), col = alpha("red", alpha = 0.5), asp = FALSE)
par(new = TRUE)
plot(rocGraph.log, col = alpha("blue", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 2, asp=FALSE)
par(new = TRUE)
plot.roc(pROC::roc(testingPred.log$obs ~ testingPred.log$Cancer), col = alpha("blue", alpha = 0.5), asp = FALSE)
par(new = TRUE)
plot(rocGraph.mips, col = alpha("orange", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 2, asp=FALSE)
par(new = TRUE)
plot.roc(pROC::roc(testingPred.mips$obs ~ testingPred.mips$Cancer), col = alpha("orange", alpha = 0.5), asp = FALSE)
legend("right", legend = c("29 genes - AUC:0.861 ; AUC:0.732", "PSA - AUC:0.640 ; AUC:0.351",
                           "Mips: - AUC:0.697 ; AUC:0.664",
                           "Training","Testing"),
       col = c("blue", "red","orange","black","black"), lty=c(rep(1,4),2), cex = .5)

dev.off()

a <- plot(rocGraph.log, col = alpha("blue", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 2, asp=FALSE)
b <- plot.roc(pROC::roc(testingPred.log$obs ~ testingPred.log$Cancer), col = alpha("blue", alpha = 0.5), asp = FALSE)
c <- plot(rocGraph.mips, col = alpha("orange", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 2, asp=FALSE)
d <- plot.roc(pROC::roc(testingPred.mips$obs ~ testingPred.mips$Cancer), col = alpha("orange", alpha = 0.5), asp = FALSE)
e <- plot(rocGraph.psa, col = alpha("red", alpha = 0.5), main = NULL, ylab = "", xlab="", lty = 2, asp = FALSE)
f <- plot.roc(pROC::roc(testingPred.psa$obs ~ testingPred.psa$Cancer), col = alpha("red", alpha = 0.5), asp = FALSE)
#save(a,file = "/mnt/DATA4/kevhu/urineRNA/tmp/a.Robj")
#save(b,file = "/mnt/DATA4/kevhu/urineRNA/tmp/b.Robj")
#save(c,file = "/mnt/DATA4/kevhu/urineRNA/tmp/c.Robj")
#save(d,file = "/mnt/DATA4/kevhu/urineRNA/tmp/d.Robj")
#save(e,file = "/mnt/DATA4/kevhu/urineRNA/tmp/e.Robj")
#save(f,file = "/mnt/DATA4/kevhu/urineRNA/tmp/f.Robj")

load("/mnt/DATA4/kevhu/urineRNA/tmp/a.Robj")
load("/mnt/DATA4/kevhu/urineRNA/tmp/b.Robj")
load("/mnt/DATA4/kevhu/urineRNA/tmp/c.Robj")
load("/mnt/DATA4/kevhu/urineRNA/tmp/d.Robj")
load("/mnt/DATA4/kevhu/urineRNA/tmp/e.Robj")
load("/mnt/DATA4/kevhu/urineRNA/tmp/f.Robj")


pdf(file = "/mnt/DATA4/kevhu/urineRNA/20180816urineVsGG2to5.pdf", onefile = TRUE, useDingbats = TRUE, pointsize = 18,width = 9)
plot(e, col = alpha("red", alpha = 0.5), main = "ROCs: Training n=80; Testing n=40", ylab = "", xlab="", lty = 2, asp = FALSE)
par(new = TRUE)
plot.roc(f, col = alpha("red", alpha = 0.5), asp = FALSE)
par(new = TRUE)
plot(a, col = alpha("blue", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 2, asp=FALSE)
par(new = TRUE)
plot.roc(b, col = alpha("blue", alpha = 0.5), asp = FALSE)
par(new = TRUE)
plot(c, col = alpha("orange", alpha = 0.5),main = NULL, ylab = "", xlab="", lty = 2, asp=FALSE)
par(new = TRUE)
plot.roc(d, col = alpha("orange", alpha = 0.5), asp = FALSE)

legend("right", legend = c("29 genes - AUC:0.861 ; AUC:0.732", "PSA - AUC:0.687 ; AUC:0.616",
                           "Mips: - AUC:0.697 ; AUC:0.664",
                           "Training","Testing"),
       col = c("blue", "red","orange","black","black"), lty=c(rep(1,4),2), cex = .5)

dev.off()

