# devtools::install_github("https://github.com/rvalieris/signeR")
library(signeR)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)

### looking at example of mut and opp file so i make sure the format is correct

example_mut <- read.table(system.file("extdata","21_breast_cancers.mutations.txt",
                              package="signeR"), header=TRUE, check.names=FALSE)
example_opp <- read.table(system.file("extdata","21_breast_cancers.opportunity.txt",
                              package="signeR"))
example_Pmat <- as.matrix(read.table(system.file("extdata","Cosmic_signatures_BRC.txt",
                                                package="signeR"), sep="\t", check.names=FALSE))


### creating opportunity martix 

mut <- read.table("/mnt/DATA5/tmp/kev/sigProfileExtractor/20221004concordMutMat_AllMut.txt",
                  sep = "\t", stringsAsFactors = FALSE, header = TRUE)
mut2 <- t(mut[, 2:ncol(mut)])
colnames(mut2) <- colnames(example_mut)

Pmatrix <- as.matrix(read.table("/mnt/DATA5/tmp/kev/tmpDbs/cosmic/COSMIC_v3.3_SBS_GRCh38.txt", sep="\t", check.names=FALSE, header = TRUE))
Pmatrix2 <- Pmatrix[,2:ncol(Pmatrix)]
Pmatrix2 <- apply(Pmatrix2, 2, as.numeric)
rownames(Pmatrix2) <- rownames(example_Pmat)

cn_bed <- read.table("/mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed", sep = "\t",
                     stringsAsFactors = FALSE, header = FALSE)
cn_bed2 <- cn_bed[, 1:3]
write.table(cn_bed2, "/mnt/DATA5/tmp/kev/misc/20221205_IAD202670_reduced.bed", sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)

target_regions <- import(con="/mnt/DATA5/tmp/kev/misc/20221205_IAD202670_reduced.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Mmusculus.UCSC.mm10,
                                target_regions, nsamples=10)


### could try a few versions of this, i.e fixedP true and false, no opportunity matrix, and then having nsig fixed to 1: total is 9 different runs
### getting rid of artifact signatures
artSigs <- paste0("SBS", c(27, 43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 95))
Pmatrix3 <- Pmatrix2[, -which(colnames(Pmatrix2) %in% artSigs)]

signatures.Pstart <- signeR(M=mut2, P=Pmatrix3, fixedP=TRUE, main_eval=100, EM_eval=50, EMit_lim=20)
BICboxplot(signatures.Pstart)

ExposureBarplot(signatures.Pstart$SignExposures)

exposures <- Median_exp(signatures.Pstart$SignExposures)
signatures <- Median_sign(signatures.Pstart$SignExposures)
apply(exposures, 1, sum)[order(apply(exposures, 1, sum))]

# highestSigs <- c("SBS32", "SBS22", "SBS37", "SBS54", "SBS12", "SBS9")
highestSigs <- c("SBS32", "SBS22", "SBS37", "SBS25", "SBS26", "SBS85", "SBS34", "SBS16")
exposuresHigh <- exposures[which(rownames(exposures) %in% highestSigs),]
exposuresHighDf <- data.frame(cbind(rownames(exposuresHigh), exposuresHigh))
colnames(exposuresHighDf)[1] <- "Signature"

exposuresHighDf_melt <- reshape2::melt(exposuresHighDf, id.vars = "Signature")
exposuresHighDf_melt$value <- as.numeric(exposuresHighDf_melt$value)
exposuresHighDf_melt$variable <- factor(exposuresHighDf_melt$variable,
                                        levels = c("PP22.002_Dx82", "PP22.008_Dx88", "PP22.005_Dx85", "PP22.006_Dx86", "PP22.007_Dx87",
                                                   "PP22.003_Dx83", "PP22.004_Dx84", "PP22.009_Dx89", "PP22.001_Dx81", "PP22.010_Dx90"))
# my_col <- c("yellow", "pink", "darkorange", "darkblue", "darkred", "darkgreen")
my_col <- c("yellow", "darkgreen", "pink", "darkorange", "darkblue", "darkred", "cyan", "lightgreen")
ggplot(exposuresHighDf_melt, aes(fill=Signature, y=value, x=variable)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values = my_col) + ggtitle("signeR activity plot") + theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("exposures") + xlab("samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


SignPlot(signatures.Pstart$SignExposures[,1])

### can recreate my own exposure plots by subsetting the 5 largest exposures then make histogram
### or i can try to refit to the 5 largest exposures, subset the Pmatrix ...to have 5 largest exposures

Pmatrix2_reduced <- Pmatrix2[, which(colnames(Pmatrix2) %in% highestSigs)]
signatures.Pstart2 <- signeR(M=mut2, P=Pmatrix2_reduced, fixedP=TRUE, main_eval=100, EM_eval=50, EMit_lim=20)
ExposureBarplot(signatures.Pstart2$SignExposures)

### all using same KL start for NMF initialization

signatures.Pstart_KL <- signeR(M=mut2, P=Pmatrix2, fixedP=TRUE, main_eval=100, EM_eval=50, EMit_lim=20, start = "KL")
exposures_KL <- Median_exp(signatures.Pstart_KL$SignExposures)

signatures.Pstart_KL_hyper <- signeR(M=mut2, P=Pmatrix2, fixedP=TRUE, estimate_hyper = TRUE,
                                     main_eval=100, EM_eval=50, EMit_lim=20, start = "KL")
exposures_KL_hyper <- Median_exp(signatures.Pstart_KL_hyper$SignExposures)

highestSigs_KL_hyper <- c("SBS32", "SBS22", "SBS37", "SBS54", "SBS85", "SBS30")

exposuresHigh_KL_hyper <- exposures_KL_hyper[which(rownames(exposures_KL_hyper) %in% highestSigs_KL_hyper),]
exposuresHighDf_KL_hyper <- data.frame(cbind(rownames(exposuresHigh_KL_hyper), exposuresHigh_KL_hyper))
colnames(exposuresHighDf_KL_hyper)[1] <- "Signature"

exposuresHighDf_KL_hyper_melt <- reshape2::melt(exposuresHighDf_KL_hyper, id.vars = "Signature")
exposuresHighDf_KL_hyper_melt$value <- as.numeric(exposuresHighDf_KL_hyper_melt$value)
ggplot(exposuresHighDf_KL_hyper_melt, aes(fill=Signature, y=value, x=variable)) + 
  geom_bar(position="stack", stat="identity")



mexposuresig <- readRDS("/mnt/DATA5/tmp/kev/lrivaExporsureStudy/mexposuresig.rds")
mexposure <- readRDS("/mnt/DATA5/tmp/kev/lrivaExporsureStudy/mexposure.rds")
mSBS <- readRDS("/mnt/DATA5/tmp/kev/lrivaExporsureStudy/mSBSs.rds")
mSBS_input <- mSBS
rownames(mSBS_input) <- rownames(Pmatrix2)

signatures.Pstart_mSBS <- signeR(M=mut2, P=mSBS_input, fixedP=TRUE, main_eval=100, EM_eval=50, EMit_lim=20)
ExposureBarplot(signatures.Pstart_mSBS$SignExposures)



signatures.Pstart_KL_hyper_nlim3 <- signeR(M=mut2, P=Pmatrix2, fixedP=TRUE, estimate_hyper = TRUE,
                                     main_eval=100, EM_eval=50, EMit_lim=20, start = "KL", nsig = 3)
ExposureBarplot(signatures.Pstart_KL_hyper_nlim3$SignExposures)
exposures_KL_hyper_nsig3 <- Median_exp(signatures.Pstart_KL_hyper_nlim3$SignExposures)

sigsFromAlexandrov <- c("SBS5", "SBS22", "SBS32")
Pmatrix2_reduced3 <- Pmatrix2[, which(colnames(Pmatrix2) %in% sigsFromAlexandrov)]
signatures.Pstart_alexandrov <- signeR(M=mut2, P=Pmatrix2_reduced3, fixedP=TRUE, main_eval=100, EM_eval=50, EMit_lim=20)
ExposureBarplot(signatures.Pstart_alexandrov$SignExposures)

exposures_top3 <- Median_exp(signatures.Pstart_alexandrov$SignExposures)
exposuresHighDf_top3 <- data.frame(cbind(rownames(exposures_top3), exposures_top3))
colnames(exposuresHighDf_top3)[1] <- "Signature"

exposuresHighDf_top3_melt <- reshape2::melt(exposuresHighDf_top3, id.vars = "Signature")
exposuresHighDf_top3_melt$value <- as.numeric(exposuresHighDf_top3_melt$value)
exposuresHighDf_top3_melt$variable <- factor(exposuresHighDf_top3_melt$variable,
                                        levels = c("PP22.002_Dx82", "PP22.008_Dx88", "PP22.005_Dx85", "PP22.006_Dx86", "PP22.007_Dx87",
                                                   "PP22.003_Dx83", "PP22.004_Dx84", "PP22.009_Dx89", "PP22.001_Dx81", "PP22.010_Dx90"))

my_col2 <- c( "plum2", "darkorange", "seagreen3")
ggplot(exposuresHighDf_top3_melt, aes(fill=Signature, y=value, x=variable)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values = my_col2) + ggtitle("signeR activity plot: fitted from sigprofilerextractor results") + theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("exposures") + xlab("samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 



### testing recreating the mutation group plots

allMutDecomposedMatrix <- read.table("/mnt/DATA5/tmp/kev/sigProfileExtractor/allMut/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/Decomposed_Mutation_Probabilities.txt",
                                     sep = "\t", header = TRUE, stringsAsFactors = FALSE)
allMutMatrix <- read.table("/mnt/DATA5/tmp/kev/sigProfileExtractor/20221004concordMutMat_AllMut.txt",
                            sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# allMutDecomposedMatrix_signer <- read.table("/mnt/DATA5/tmp/kev/sigProfilerAssignment/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities.txt",
#                                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)

allMutDecomposedMatrix_signer <- read.table("/mnt/DATA5/tmp/kev/sigProfilerAssignment/noArt/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities.txt",
                                            sep = "\t", header = TRUE, stringsAsFactors = FALSE)

### no real probability singal in sbs8 and sbs37 so getting rid of them
allMutDecomposedMatrix_signer <- allMutDecomposedMatrix_signer[, c("Sample.Names","MutationType", "SBS12", "SBS16","SBS22", "SBS32")]



i <- unique(allMutDecomposedMatrix$Sample.Names)[2]
finalActvityMatrix <- NULL
mutBySbsCat <- NULL
mutBySbsCat_signer <- NULL
for (i in unique(allMutDecomposedMatrix$Sample.Names)) {
  tmpMatrix <- round(allMutDecomposedMatrix[which(allMutDecomposedMatrix$Sample.Names == i), 3:6], digits = 2)
  rownames(tmpMatrix) <- allMutMatrix$Mutation.Types
  tmpRes <- apply(tmpMatrix, 2, function(x) x * allMutMatrix[, which(colnames(allMutMatrix) == i)])

  ### checking if total is same as graphs from output
  tmpRes2 <- apply(tmpRes, 2, sum)
  finalActvityMatrix <- rbind(finalActvityMatrix, c(i, tmpRes2))
  tmpMatrix2 <- cbind("sbs" = allMutMatrix$Mutation.Types, tmpMatrix)
  # mutBySbsCat <- rbind(mutBySbsCat, data.frame("sample" = rep(i, nrow(melt(tmpRes))), melt(tmpRes)))
  mutBySbsCat <- rbind(mutBySbsCat, data.frame(melt(tmpRes)))
  
  ### do the same thing but for signer samples
  # tmpMatrix_signer <- round(allMutDecomposedMatrix_signer[which(allMutDecomposedMatrix_signer$Sample.Names == i), 3:8], digits = 2)
  tmpMatrix_signer <- round(allMutDecomposedMatrix_signer[which(allMutDecomposedMatrix_signer$Sample.Names == i), 3:6], digits = 2)
  rownames(tmpMatrix_signer) <- allMutMatrix$Mutation.Types
  tmpRes_signer <- apply(tmpMatrix_signer, 2, function(x) x * allMutMatrix[, which(colnames(allMutMatrix) == i)])
  mutBySbsCat_signer <- rbind(mutBySbsCat_signer, data.frame(melt(tmpRes_signer)))
}

### create graphing function  to do ordred stacked histograms

# tmpGraphVector <- melt(tmpMatrix2)
# tmpGraphVector$sbs <- factor(tmpGraphVector$sbs, levels = allMutMatrix$Mutation.Types)
# ggplot(data = tmpGraphVector, aes(x = sbs, y = value, fill = variable)) + geom_bar(stat = "identity") +
#   theme(panel.background = element_blank(), axis.text.x=element_text(angle = 90, vjust = 0.5))

mycolVector1 <- my_col2 <- c("lightgreen", "seagreen3", "plum2", "darkorange")

tmpGraphVector <- melt(tmpRes)
tmpGraphVector$Var1 <- factor(tmpGraphVector$Var1, levels = allMutMatrix$Mutation.Types)
a <- ggplot(data = tmpGraphVector, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = my_col2) +
  theme(panel.background = element_blank(), axis.text.x=element_text(angle = 90, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("mutations") + xlab("SBS") + ggtitle("PP22.002_Dx82 mutataions (n = 183)")


tmpGraphVector_all <- melt(mutBySbsCat)
tmpGraphVector_all$Var1 <- factor(tmpGraphVector_all$Var1, levels = allMutMatrix$Mutation.Types)
a_whole <- ggplot(data = tmpGraphVector_all, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = my_col2) +
  theme(panel.background = element_blank(), axis.text.x=element_text(angle = 90, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("mutations") + xlab("SBS") + ggtitle("All mutations exposure signalProfilerExtractor (n = 183)")



# my_col3 <- c("darkred", "darkgreen", "plum2", "darkorange", "darkblue", "purple3")
my_col3 <- c("yellow", "darkgreen", "plum2", "darkorange")
tmpGraphVector_signer <- melt(tmpRes_signer)
tmpGraphVector_signer$Var1 <- factor(tmpGraphVector_signer$Var1, levels = allMutMatrix$Mutation.Types) 
b <- ggplot(data = tmpGraphVector_signer, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = my_col3) +
  theme(panel.background = element_blank(), axis.text.x=element_text(angle = 90, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("mutations") + xlab("SBS") + ggtitle("PP22.002_Dx82 mutataions (n = 183)")


tmpGraphVector_signer_all <- melt(mutBySbsCat_signer)
tmpGraphVector_signer_all$Var1 <- factor(tmpGraphVector_signer_all$Var1, levels = allMutMatrix$Mutation.Types) 
b_whole <- ggplot(data = tmpGraphVector_signer_all, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = my_col3) +
  theme(panel.background = element_blank(), axis.text.x=element_text(angle = 90, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("mutations") + xlab("SBS") + ggtitle("All mutataions signeR (n = 183)")



# gridExtra::grid.arrange(a, b, nrow = 2)
gridExtra::grid.arrange(a_whole, b_whole, nrow = 2)


### probability matrix for highest signer signatures - put this and sample matrix to get decomposed probabilities
# signerPhat <- signatures.Pstart$Phat
# colnames(signerPhat) <- colnames(Pmatrix2)
# signerPhat_reduced <- data.frame("MutationType" = allMutMatrix$Mutation.Types, signerPhat[, which(colnames(signerPhat) %in% highestSigs)])
# write.table(signerPhat_reduced, "/mnt/DATA5/tmp/kev/sigProfilerAssignment/20221220signerHighestExposureSigs.txt", sep = "\t",
#             quote = FALSE, row.names = FALSE, col.names = TRUE)

### making matrix for probability fitting
signatures.Pstart$Phat
signerPhat <- signatures.Pstart$Phat
colnames(signerPhat) <- colnames(Pmatrix3)
signerPhat_reduced <- data.frame("MutationType" = allMutMatrix$Mutation.Types, signerPhat[, which(colnames(signerPhat) %in% highestSigs)])
write.table(signerPhat_reduced, "/mnt/DATA5/tmp/kev/sigProfilerAssignment/20230102signerHighestExposureSigsNoArt.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)



