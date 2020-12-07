lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

library(VSURF)
library(parallel)
library(readxl)
library(ggplot2)
library(pROC)
library(reshape2)
library(gridExtra)
library(stringr)
library(caret)



a <- read_xlsx("/mnt/DATA4/kevhu/urineRNA/20189_07_26_AC17.2ExtDes2_mastermatrix.xlsx")

a <- a[-which(is.na(a$`SAT FIRST ANALYSIS`)),]
a <- a[which(a$e2e > 300000 & a$`%e2e` > .50 & a$`%OT` > .6),]
PSA <- a$LAB_RESNUM
ids <- a$`SAT FIRST ANALYSIS`
a <- a[,7:92]
a <- a[,-which(colnames(a) == "SUMT2ERG")]
a <- a[,-which(colnames(a) == "AVGT2ERG")]

set.seed(8889)
g1 <- which(ids == 0)
g2 <- which(ids == 1)

g1.subset <- sample(g1, length(g1)/3)
g2.subset <- sample(g2, length(g2)/3)
combined.subset <- c(g1.subset,g2.subset)
held.out <- a[combined.subset,]
a.2 <- a[-combined.subset,]

ids.2 <- factor(ids[-combined.subset])


ids.3 <- str_replace_all(ids.2, "1", "Cancer")
ids.3 <- str_replace_all(ids.3, "0", "NoCancer")
ids.3 <- factor(ids.3)

testRandom <- list(a.2, ids.2)
names(testRandom) <- c("x","y")
set.seed(8889)
vsurf.parallel <- VSURF(a.2, ids.3, parallel = TRUE, ncores = 24, clusterType = "FORK")


summary(vsurf.parallel)
plot(vsurf.parallel)

genesThres <- colnames(a.2)[vsurf.parallel$varselect.thres]
genesInterp <- colnames(a.2)[vsurf.parallel$varselect.interp]


a.2.log <- log2(a.2 +1)
ids.4 <- relevel(ids.3, "NoCancer")
a.2.log.thres <- cbind(ids.4,a.2.log[,genesThres])
a.2.log.thres.melt <- melt(data = a.2.log.thres,id.vars = "ids.4")
colnames(a.2.log.thres.melt)[1] <- "Classes"


### need to switch to non-cancer first, also need to get original 30 genes found by random forest


ggplot(data = a.2.log.thres.melt, aes(x = variable, y = value, colour = Classes)) + geom_boxplot() +
  facet_wrap(~variable, scales = "free") + theme_bw() + ylab(label = "log2 normalized expresssion") + 
  ggtitle(label = "Urine panel: Selected genes") + xlab(label = "") +
  theme(axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5),panel.grid = element_blank(), legend.title = element_text(hjust = 0.5)) +
  theme(legend.key.size = unit(0.5, "inches"), title = element_text(face = "bold"), strip.text.x = element_text(size = 7.5))



a.3 <- cbind(ids.3,a.2[,genesThres])
train_control <- trainControl(method = "cv", number = 5,
                              classProbs = TRUE, summaryFunction = twoClassSummary)

set.seed(8889)
fiveFoldModel <- train(ids.3 ~., data = a.3, trControl=train_control, method="glmnet",
                       metric = "ROC", family = "binomial", maxit = 10000000)
fiveFoldModel


ids.held.out <- ids[combined.subset]
ids.held.out <- str_replace_all(ids.held.out, "1", "Cancer")
ids.held.out <- str_replace_all(ids.held.out, "0", "NoCancer")
ids.held.out <- as.factor(ids.held.out)
held.out.thres <- held.out[,genesThres]

prob <- predict(fiveFoldModel, newdata = held.out.thres)
confusionMatrix(data = factor(prob), ids.held.out, positive = "Cancer")

probROC <- predict(fiveFoldModel, type = c("prob"),newdata = held.out.thres)

rocGraph <- roc(ids.held.out ~ probROC$`Cancer`)
plot(rocGraph)
rocGraph$auc

### the model itself isn't too great so let's try a log2 transform first, then maybe adding PSA


testRandom.log <- list(a.2.log, ids.2)
names(testRandom) <- c("x","y")
set.seed(8889)
vsurf.parallel.log<- VSURF(a.2.log, ids.3, parallel = TRUE, ncores = 24, clusterType = "FORK")

genesThres.log <- colnames(a.2.log)[vsurf.parallel.log$varselect.thres]
genesInterp.log <- colnames(a.2.log)[vsurf.parallel.log$varselect.interp]


a.4 <- cbind(ids.3,a.2.log[,genesThres.log])

set.seed(8889)
fiveFoldModel.log <- train(ids.3 ~., data = a.4, trControl=train_control, method="glmnet",
                       metric = "ROC", family = "binomial", maxit = 10000000)
fiveFoldModel.log


held.out.thres.log <- held.out[,genesThres.log]
held.out.thres.log <- log2(held.out.thres.log + 1)
prob.log <- predict(fiveFoldModel.log, newdata = held.out.thres.log)
confusionMatrix(data = factor(prob.log), ids.held.out, positive = "Cancer")

probROC.log <- predict(fiveFoldModel.log, type = c("prob"),newdata = held.out.thres.log)

rocGraph.log <- roc(ids.held.out ~ probROC.log$`Cancer`)
plot(rocGraph.log)
rocGraph.log$auc

### model with log transformed data seeems a lot cleaner 
### let's try adding PSA too

PSA

PSA.trianing <- PSA[-combined.subset]
one.training <- rep(1, length(PSA.trianing))
a.5 <- data.frame(PSA.trianing, one.training)

set.seed(8889)
fiveFoldModel.psa <- train(ids.3 ~., data = a.5, trControl=train_control, method="glmnet",
                           metric = "ROC", family = "binomial", maxit = 10000000)
fiveFoldModel.psa

PSA.held.out <- data.frame(PSA[combined.subset])
PSA.held.out$ones <- rep(1, nrow(PSA.held.out))
colnames(PSA.held.out) <- colnames(a.5)

prob.psa <- predict(fiveFoldModel.psa, newdata = PSA.held.out)
confusionMatrix(data = factor(prob.psa), ids.held.out, positive = "Cancer")

probROC.psa <- predict(fiveFoldModel.psa, type = c("prob"),newdata = PSA.held.out)
rocGraph.psa <- roc(ids.held.out ~ probROC.psa$`Cancer`)
plot(rocGraph.psa)
rocGraph.psa$auc

### combining PSA with logged model

a.6 <- cbind(PSA.trianing,a.2.log[,genesThres.log])
set.seed(8889)
fiveFoldModel.combined <- train(ids.3 ~., data = a.6, trControl=train_control, method="glmnet",
                           metric = "ROC", family = "binomial", maxit = 10000000)
fiveFoldModel.combined

combined.held.out <- cbind(PSA.held.out[,1], held.out.thres.log)
colnames(combined.held.out) <- colnames(a.6)
prob.combined <- predict(fiveFoldModel.combined, newdata = combined.held.out)
confusionMatrix(data = factor(prob.combined), ids.held.out, positive = "Cancer")

probROC.combined <- predict(fiveFoldModel.combined, type = c("prob"),newdata = combined.held.out)
rocGraph.combined <- roc(ids.held.out ~ probROC.combined$Cancer)
plot(rocGraph.combined)
rocGraph.combined$auc



tableOfStats <- data.frame("Model + Data" = c("Training 30","Testing 30 Genes","Training PSA","Testing PSA",
                                              "Training Logged Data","Testing Logged Data",
                                              "Training Combined","Testing Combined"),
                           "AUC" = c("0.77","0.62","0.69", "0.67","0.84","0.76","0.86","0.78"),
                           "Sensitivity" = c("61%","47%","22","33%","62%","60%","65%","67%"),
                           "Specificity" = c("75%","79%","98%","95%","82%","86%","82%","82%"))



#tableOfStats.grob <- tableGrob(tableOfStats)
#png("/mnt/DATA4/kevhu/urineRNA/tableOfStats.png",width = 800,height = 400)
#grid.arrange(tableOfStats.grob)
#dev.off()
