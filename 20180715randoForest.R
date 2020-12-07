lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)


# install.packages("VSURF")
library(VSURF)
library(parallel)
library(readxl)
library(ggplot2)
library(pROC)
library(reshape2)
library(gridExtra)
library(caret)
library(stringr)
library(glmnet)
library(Rcpp)

a <- read_xlsx("/mnt/DATA4/kevhu/urineRNA/AC17.2ExtDes2_mastermatrix_allsamples_byKLK3.xlsx")
b <- read_xlsx("/mnt/DATA4/kevhu/urineRNA/20189_07_26_AC17.2ExtDes2_mastermatrix.xlsx")


a <- a[-which(a$`Grade group` == 2),]
a.1 <- a[1:109,3:86]
rownames(a.1) <- as.vector(a$Contig_ID[1:109])
ids <- a$`Grade group bin`[1:109]


###need to reorder b in same as a so PSA matches
b <- b[match(rownames(a.1), b$Contig_ID),]
PSA <- b$LAB_RESNUM


table(ids)
###need to subsample equally between the two groups
set.seed(2222)
g1 <- which(ids == "Benign/GG1")
g2 <- which(ids == "GG2-5")

g1.subset <- sample(g1, length(g1)/3)
g2.subset <- sample(g2, length(g2)/3)
combined.subset <- c(g1.subset,g2.subset)

held.out <- a.1[combined.subset,]
a.2 <- a.1[-combined.subset,]
ids.2 <- factor(ids[-combined.subset])
ids.3 <- str_replace_all(ids.2, "GG2-5", "Cancer")
ids.3 <- str_replace_all(ids.3, "Benign/GG1", "NoCancer")

ids.full <- factor(ids)
ids.full <- str_replace_all(ids.full, "GG2-5", "Cancer")
ids.full <- str_replace_all(ids.full, "Benign/GG1", "NoCancer")

testRandom <- list(a.2, ids.2)
names(testRandom) <- c("x","y")
set.seed(2222, kind = "L'Ecuyer-CMRG")
vsurf.parallel <- VSURF(testRandom$x, testRandom$y, parallel = TRUE, ncores = 24, clusterType = "FORK")

summary(vsurf.parallel)
plot(vsurf.parallel)

colnames(a.2)[vsurf.parallel$varselect.interp]
colnames(a.2)[vsurf.parallel$varselect.thres]

genes30 <- colnames(a.2)[vsurf.parallel$varselect.thres]
genes3 <- colnames(a.2)[vsurf.parallel$varselect.interp]


load("/mnt/DATA4/kevhu/urineRNA/20180807randoForestGenes.Robj")

a.2.genes30 <- a.2[,genes30]
a.2.genes30 <- a.2.genes30[, order(colnames(a.2.genes30))]


featurePlot(x = log2(a.2.genes30 + 1),
            y = ids.3,
            plot = "box",
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")),
            par.strip.text=list(cex=.6),
            par.strip.naames=list(FALSE,FALSE),
            pch= "|",
            main = "Boxplots of reduced gene set")

###prior was no log2 expression value - this portion uses the log2 transformed data


a.2.log <- log2(a.1[-combined.subset,]+1)


testRandom.log <- list(a.2.log, ids.2)
names(testRandom.log) <- c("x","y")
set.seed(2222, kind = "L'Ecuyer-CMRG")
vsurf.parallel.log <- VSURF(testRandom$x, testRandom$y, parallel = TRUE, ncores = 24, clusterType = "FORK")

summary(vsurf.parallel)

genes3.log <- colnames(a.2.log)[vsurf.parallel.log$varselect.interp]
genes30.log <- colnames(a.2.log)[vsurf.parallel.log$varselect.thres]

featurePlot(x = a.2.log[,genes30.log],
            y = ids.3,
            plot = "box",
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")),
            par.strip.text=list(cex=.6),
            xlab = "",
            pch= "|")


### better overall boxplot than featurePlot - mainly because I'm not too familiar with lattice graphing parameters which is used by featurePlot

#a.2.log <- log2(a.1[-combined.subset,]+1)
a.2.log <- log2(a.1 +1)
a.2.log.30 <- cbind(ids.full,a.2.log[,genes30])
a.2.log.30.melt <- melt(data = a.2.log.30,id.vars = "ids.full")
colnames(a.2.log.30.melt)[1] <- "Classes"

### need to switch to non-cancer first, also need to get original 30 genes found by random forest

pdf(file = "/mnt/DATA4/kevhu/urineRNA/20180831urineBoxplots.pdf", onefile = TRUE, useDingbats = TRUE, width = 10)
ggplot(data = a.2.log.30.melt, aes(x = variable, y = value, colour = Classes)) + geom_boxplot() +
  facet_wrap(~variable, scales = "free") + theme_bw() + ylab(label = "log2 normalized expresssion") + 
  ggtitle(label = "Urine panel: 29 selected genes") + xlab(label = "") +
  theme(axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5),panel.grid = element_blank(), legend.title = element_text(hjust = 0.5)) +
  theme(legend.key.size = unit(0.5, "inches"), title = element_text(face = "bold"), strip.text.x = element_text(size = 3))
dev.off()




### below is just building a custom model - not the best, just testing if it works

ids.3 <- str_replace_all(ids.2, "GG2-5", "Cancer")
ids.3 <- str_replace_all(ids.3, "Benign/GG1", "NoCancer")
#ids.2 <- as.character(ids.2)
#ids.3 <- str_replace_all(ids.2, "GG2-5", "1")
#ids.3 <- str_replace_all(ids.3, "Benign/GG1", "0")
#ids.3 <- as.numeric(ids.3)

out1 <- glm(ids.3 ~ a.2.genes8$AC009478.1.E3E4 + a.2.genes8$ERG.E1E2 + a.2.genes8$HOXC6.E1E2 + a.2.genes8$HPN.E10E11 + 
              a.2.genes8$`NKX3-1.E1E2` + a.2.genes8$PDLIM5.E6E7 + a.2.genes8$`RP11-314O13.1.E3E4` + a.2.genes8$SCHLAP1.E1E2 + 
              a.2.genes8$TDRD1.E16E17 + a.2.genes8$`TMPRSS2-ERG.EF194202` + a.2.genes8$`TMPRSS2-ERG.T1E4.COSF125` + 
              a.2.genes8$`TMPRSS2-ERG.T1EIIIc_4` + a.2.genes8$`TMPRSS2-ERG.T2E5.COSF129`, family = binomial("logit"))

#summary(out1)
#prob <- predict(out1, type = c("response"), data = ids.3)
#rocGraph <- roc(ids.3 ~ prob)
#plot(rocGraph)


###use caret for CV
###
###
###




ids.3 <- str_replace_all(ids.2, "GG2-5", "Cancer")
ids.3 <- str_replace_all(ids.3, "Benign/GG1", "NoCancer")

ids.3 <- factor(ids.3)
a.3 <- cbind(ids.3,a.2.genes30)
train_control <- trainControl(method = "cv", number = 5,
                              classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)

set.seed(2222)
fiveFoldModel <- train(ids.3 ~., data = a.3, trControl=train_control, method="glmnet",
                       metric = "ROC", family = "binomial")
fiveFoldModel

### alpha was tuned to 1 and lambda to 0.002680061

ids.held.out <- ids[combined.subset]
ids.held.out <- str_replace_all(ids.held.out, "GG2-5", "Cancer")
ids.held.out <- str_replace_all(ids.held.out, "Benign/GG1", "NoCancer")
ids.held.out <- as.factor(ids.held.out)
held.out.30 <- held.out[,genes30]

prob <- predict(fiveFoldModel, newdata = held.out.30)
confusionMatrix(data = factor(prob), ids.held.out, positive = "Cancer")

probROC <- predict(fiveFoldModel, type = c("prob"),newdata = held.out.30)

rocGraph <- pROC::roc(ids.held.out ~ probROC$`Cancer`)
plot(rocGraph)

rocGraph$auc
### next just using PSA
###
###
###

### adding 

PSA.heldout <- PSA[-combined.subset]
a.5 <- data.frame(ids.3,PSA.heldout)
a.5$ones <- rep(1,nrow(a.5))
train_control.psa <- trainControl(method = "cv", number = 5,
                              classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
set.seed(2222)
fiveFoldModel.psa <- train(ids.3 ~ ., data = a.5, trControl=train_control.psa, method="glmnet",
                       metric = "ROC", family = "binomial")
fiveFoldModel.psa
### alpha was tuned to .1 and lambda to 0.03105096 

held.out.30.psa <- data.frame(PSA[combined.subset])
held.out.30.psa$ones <- rep(1, nrow(held.out.30.psa))
colnames(held.out.30.psa)[1] <- "PSA.heldout"

prob.psa <- predict(fiveFoldModel.psa, newdata = held.out.30.psa)
confusionMatrix(data = factor(prob.psa), ids.held.out, positive = "Cancer")

probROC.psa <- predict(fiveFoldModel.psa, type = c("prob"),newdata = held.out.30.psa)
rocGraph.psa <- roc(ids.held.out ~ probROC.psa$`Cancer`)
plot(rocGraph.psa)
rocGraph.psa$auc



### after let's try and run the model with important genes + psa
###
###

a.6 <- cbind(ids.3,a.2.genes30,PSA.heldout)

train_control.combined <- trainControl(method = "cv", number = 5,
                              classProbs = TRUE, summaryFunction = twoClassSummary)
set.seed(2222)
fiveFoldModel.combined <- train(ids.3 ~., data = a.6, trControl=train_control.combined, method="glmnet",
                       metric = "ROC", family = "binomial")
fiveFoldModel.combined

held.out.combined <- held.out[,genes30]
held.out.combined <- cbind(held.out.combined, PSA[combined.subset])
colnames(held.out.combined)[ncol(held.out.combined)] <- "PSA.heldout"

prob.combinined <- predict(fiveFoldModel.combined, newdata = held.out.combined)
confusionMatrix(data = factor(prob.combinined), ids.held.out, positive = "Cancer")

probROC.combined <- predict(fiveFoldModel.combined, type = c("prob"),newdata = held.out.combined)
rocGraph.combined <- roc(ids.held.out ~ probROC.combined$`Cancer`)
plot.roc(rocGraph.combined, main = "PSA with model ROC")
rocGraph.combined$auc



### doing it with logged data ...


a.7 <- a.2.log.30

train_control.combined <- trainControl(method = "cv", number = 5,
                                       classProbs = TRUE, summaryFunction = twoClassSummary)
set.seed(2222)
fiveFoldModel.log <- train(ids.3 ~., data = a.7, trControl=train_control.combined, method="glmnet",
                                metric = "ROC", family = "binomial")
fiveFoldModel.log

held.out.log <- log2(held.out[,genes30] + 1)
prob.log <- predict(fiveFoldModel.log, newdata = held.out.log)
confusionMatrix(data = factor(prob.log), ids.held.out, positive = "Cancer")

probROC.log <- predict(fiveFoldModel.log, type = c("prob"),newdata = held.out.log)
rocGraph.log <- roc(ids.held.out ~ probROC.log$`Cancer`)
plot.roc(rocGraph.log, main = "PSA with model ROC")
rocGraph.log$auc


### Next is creating the ROC graphs for the 3 types of analysis, just gene, just PSA, PSA + gene

plot(rocGraph.combined, main = "ROCs:Testing Data n = 36", col = "green", xlim = c(1,0), ylim = c(0,1), asp = FALSE)
par(new = TRUE)
plot(rocGraph.psa, col = "red", main = NULL, ylab = "", xlab="", asp = FALSE)
par(new = TRUE)
plot(rocGraph, col = "blue",main = NULL, ylab = "", xlab="", asp=FALSE)
legend("right", legend = c("30 Genes     AUC:0.862", "PSA     AUC:0.625","30 Genes + PSA     AUC:0.79"),
       col = c("blue", "red","green"), lty = 1, cex = 0.65)



### making dataframe of summary stats for the different models

tableOfStats <- data_frame("Model + Data" = c("Training 30 Genes","Testing 30 Genes","Training PSA","Testing PSA","Training Combined","Testing Combined"),
                           "AUC" = c("0.84","0.86","0.68", "0.63","0.82","0.79"),
                           "Sensitivity" = c("76%","85.7%","18%","7%","68%","71%"),
                           "Specificity" = c("75%","68%","100%","90%","75%","73%"))

tableOfStats.grob <- tableGrob(tableOfStats)
png("/mnt/DATA4/kevhu/urineRNA/tableOfStats.png",width = 800,height = 400)
grid.arrange(tableOfStats.grob)
dev.off()
