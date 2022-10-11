### expected BAF is from ASCAT paper (PNAS)

### generalized expected b-allele frequency is: b_i = n_Bi / ( n_Ai + n_Bi)
### where b_i is b-allele frequency at i genomic location and n is allele specific copy numbers of of allele A and B at i location
### e.g diploid het assuming 100 tumor content is b_i = 1/ (1 + 1) = 0.5

### adding in tumor content
### (1 - p + p*n_bi)/ 2 - 2p + p(n_Ai + n_Bi); p is cell fraction
### p is tumor content while assuming 2 copies of gene in normal cell i.e contamination
### e.g in diploid with 90% tumor content (1 - 0.90 + 0.90(1))/ (2 - 2(0.90) + p(1 + 1))

expectedBafCalc <- function(p, na, nb){
  # p is tumor content, na and nb are number of alleles a and b
  # if (na  == 0) {
  #   stop("A allele can't have 0 copies")
  # }
  res <- (1 - p + p * nb) / (2 - 2 * p + p * (na + nb))
}


expectedLogRCalc <- function(p, na, nb, pl){
  # if (na  == 0) {
  #   stop("A allele can't have 0 copies")
  # }
  res <- log2( (2* (1-p) + p* (na + nb)) / pl )
}



### testing, works!
expectedBafCalc(1, 2, 1)
tumorContents <- seq(0.1, 1, 0.1)


### diploid

diploidTable <- data.frame("ploidy" = rep(2, 6 * length(tumorContents)),
                           "tc" = rep(tumorContents, 6),
                           "na" = c(rep(1,10), rep(1,10), rep(2,10), rep(2,10), rep(3,10), rep(0,10)),
                           "nb" = c(rep(1,10), rep(0,10), rep(1,10), rep(0,10), rep(1,10), rep(0,10)),
                           "change" = c(rep("no change", 10), rep("one loss", 10),
                                        rep("one gain", 10), rep("copy neutral", 10),
                                        rep("two gain", 10), rep("hom loss", 10)))

dipExpectRes <- NULL
for (i in 1:nrow(diploidTable)) {
  tmpRes <- round(expectedBafCalc(diploidTable$tc[i], diploidTable$na[i], diploidTable$nb[i]), digits = 2)
  tmpRes2 <- 1 - tmpRes
  tmpRes3 <- paste0(tmpRes, ",", tmpRes2)
  dipExpectRes <- c(dipExpectRes, tmpRes3)
}

diploidTable$baf <- dipExpectRes

### triploid

triploidTable <- data.frame("ploidy" = rep(3, 10 * length(tumorContents)),
                            "tc" = rep(tumorContents, 10),
                            "na" = c(rep(2,10), rep(2,10), rep(3,10), rep(1,10), rep(2,10), rep(3,10), rep(1,10), rep(0,10), rep(4,10), rep(3,10)),
                            "nb" = c(rep(1,10), rep(0,10), rep(1,10), rep(1,10), rep(2,10), rep(0,10), rep(0,10), rep(0,10), rep(1,10), rep(2,10)),
                            "change" = c(rep("no change", 10), rep("one loss", 10),
                                         rep("one gain", 10), rep("one loss", 10),
                                         rep("one gain", 10), rep("copy neutral", 10),
                                         rep("two loss", 10), rep("hom loss", 10),
                                         rep("two gain", 10), rep("two gain", 10)))

tripExpectRes <- NULL
for (i in 1:nrow(triploidTable)) {
  tmpRes <- round(expectedBafCalc(triploidTable$tc[i], triploidTable$na[i], triploidTable$nb[i]), digits = 2)
  tmpRes2 <- 1 - tmpRes
  tmpRes3 <- paste0(tmpRes, ",", tmpRes2)
  tripExpectRes <- c(tripExpectRes, tmpRes3)
}
triploidTable$baf <- tripExpectRes

### tetraploid

tetraploidTable <- data.frame("ploidy" = rep(4, 8 * length(tumorContents)),
                            "tc" = rep(tumorContents, 8),
                            "na" = c(rep(2,10), rep(2,10), rep(2,10), rep(3,10), rep(3,10), rep(4,10), rep(0,10), rep(4,10)),
                            "nb" = c(rep(2,10), rep(1,10), rep(0,10), rep(1,10), rep(2,10), rep(0,10), rep(0,10), rep(2,10)),
                            "change" = c(rep("no change", 10), rep("one loss", 10),
                                         rep("two loss", 10), rep("copy neutral", 10),
                                         rep("one gain", 10), rep("copy neutral", 10),
                                         rep("hom loss", 10), rep("two gain", 10)))

tetraExpectRes <- NULL
for (i in 1:nrow(tetraploidTable)) {
  tmpRes <- round(expectedBafCalc(tetraploidTable$tc[i], tetraploidTable$na[i], tetraploidTable$nb[i]), digits = 2)
  tmpRes2 <- 1 - tmpRes
  tmpRes3 <- paste0(tmpRes, ",", tmpRes2)
  tetraExpectRes <- c(tetraExpectRes, tmpRes3)
}

tetraploidTable$baf <- tetraExpectRes

### pentaploid

pentaploidTable <- data.frame("ploidy" = rep(5, 4 * length(tumorContents)),
                              "tc" = rep(tumorContents, 4),
                              "na" = c(rep(3,10), rep(3,10), rep(3,10), rep(4,10)),
                              "nb" = c(rep(2,10), rep(1,10), rep(0,10), rep(2,10)),
                              "change" = c(rep("no change", 10), rep("one loss",10),
                                           rep("two loss", 10), rep("one gain", 10)))

pentaExpectRes <- NULL
for (i in 1:nrow(pentaploidTable)) {
  tmpRes <- round(expectedBafCalc(pentaploidTable$tc[i], pentaploidTable$na[i], pentaploidTable$nb[i]), digits = 2)
  tmpRes2 <- 1 - tmpRes
  tmpRes3 <- paste0(tmpRes, ",", tmpRes2)
  pentaExpectRes <- c(pentaExpectRes, tmpRes3)
}

pentaploidTable$baf <- pentaExpectRes



### final table

allPloidyDf <- rbind(diploidTable, triploidTable, tetraploidTable, pentaploidTable)

log2Res <- NULL
for (i in 1:nrow(allPloidyDf)) {
  tmpRes <- round(expectedLogRCalc(allPloidyDf$tc[i], allPloidyDf$na[i], allPloidyDf$nb[i], allPloidyDf$ploidy[i]), digits = 2)
  log2Res <- c(log2Res, tmpRes)
}

allPloidyDf$log2R <- log2Res

write.table(allPloidyDf, "/mnt/DATA5/tmp/kev/misc/20220328ploidyChart.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

### based on aneuploidy detection script ... I can use it to help predict ploidy
### convert log2 ratios into integer based values  for 2n, 3n, 4n, 5n 

### I NEED TO REDO THE BAF PLOTS (1) get rid of TC plot (2) REDO VAR FILT (3n and 5n ratio is biased against for log2R = 0)
### (3) ADD MORE LINES TO DISTINGUISH



### ask scott if they use some type of scaling factor when doing cn-plot. cause sometimes the estimates seems a little off i.e underestimating losses
### kind of funny how in some cases lower tc i.e at 0.5 can be easier b/c you filter out subclonal noise
