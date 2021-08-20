andiXl <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/New_fishplotFigure.xlsx")

andiBase <- andiXl[2:8,3:21]
andiBase <- data.frame(lapply(andiBase, as.numeric))
colnames(andiBase) <- NULL

andiTreat <- andiXl[2:8, 22:40]
andiTreat <- data.frame(lapply(andiTreat, as.numeric))
colnames(andiTreat) <- NULL

andiBase2 <- NULL
for (i in 1:nrow(andiBase)) {
  tmpRow <- as.numeric(andiBase[i,])
  tmpRow[which(is.na(tmpRow))] <- median(tmpRow[-which(is.na(tmpRow))])
  andiBase2  <- rbind(andiBase2 , tmpRow)
}

andiBase2_vec <- signif(apply(andiBase2, 1, function(x) mean(x) * 100), digits = 4)

andiTreat2 <- NULL
for (i in 1:nrow(andiTreat)) {
  tmpRow <- as.numeric(andiTreat[i,])
  tmpRow[which(is.na(tmpRow))] <- median(tmpRow[-which(is.na(tmpRow))])
  andiTreat2  <- rbind(andiTreat2 , tmpRow)
}

andiTreat2_vec <- signif(apply(andiTreat2, 1, function(x) mean(x) * 100), digits = 4)

fishFract <- t((as.matrix(cbind(andiBase2_vec, andiTreat2_vec))))
rownames(fishFract) <- NULL
colnames(fishFract) <- NULL
parents <- c(0,1)
timepoints <- c(0,75,150)

fishFract <- matrix(c(89, 0, 
                      89, 45,
                      89, 89),
                    ncol = length(timepoints))


fish = createFishObject(fishFract,parents,
                        timepoints=timepoints,
                        clone.annots=c("CDH1 p.L220fsdel", "ESR1 p.Y537C"),
                        clone.labels=as.character(1:2))
fish = layoutClones(fish)
fish = setCol(fish,col=c("#1B9E77","#D95F02"))
sample.times = c(0,75, 150)


fishPlot(fish,shape="spline",vlines=sample.times, vlab=sample.times,
         cex.title=0.5, bg.type="solid")



