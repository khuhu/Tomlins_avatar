library(RTCGAToolbox)
library(rtracklayer)

getFirehoseDatasets()
getFirehoseRunningDates(last = 3)

getFirehoseData(dataset="COAD", runDate="20160128",
                destdir = "/mnt/DATA5/tmp/kev/misc/",
                CNASNP = TRUE, CNASeq = TRUE, CNACGH = TRUE)

getFirehoseData(dataset="COAD", runDate="20150821",
                destdir = "/mnt/DATA5/tmp/kev/misc/",
                Mutation = TRUE)

mm10tohg19_chain <- import.chain("/mnt/DATA5/tmp/kev/tmpDbs/ucscChainFiles/mm10ToHg19.over.chain")
varResults <- read.table("/mnt/DATA6/mouseData/reportAnno/Auto_user_AUS5-120-MG_EFD4_BBN_334_304_anno.txt",
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)

fdpFilt <- which(varResults$FDP > 20)
faoFilt <- which(varResults$FAO > 5)
freqFilt <- which(varResults$AF > 0.05)
hrunFilt <- which(varResults$HRUN < 4)
strandRatio <- intersect(which(varResults$FSAF/varResults$FSAR > 0.2),
                         which(varResults$FSAF/varResults$FSAR < 5))
qualFilt <- which(varResults$QUAL > 60)
goodSamps <- Reduce(intersect, list(fdpFilt, faoFilt, freqFilt, strandRatio, hrunFilt, qualFilt))
varResults<- varResults[goodSamps,]
mgpVars <- c("hom", "het")
varResults <- varResults[-which(varResults$mm10_mpgpv6_Indels %in% mgpVars),]


hg19out <- NULL
for(i in 1:nrow(varResults)){
  tmpRange <- GRanges(seqnames = varResults$Chr[i],
                      IRanges(start = varResults$Start[i], end = varResults$End[i]))
  tmpOut <- data.frame(liftOver(tmpRange, mm10tohg19_chain), stringsAsFactors = FALSE)
  if(nrow(tmpOut) == 0){
    next()
    tmpOut <- data.frame("hg19chr" = NA,"hg19start" = NA,"hg19end" = NA)
  } else {
    tmpOut <- tmpOut[,c("seqnames", "start", "end")]
    colnames(tmpOut) <- c("hg19chr", "hg19start", "hg19end")
  }
  
  tmpOrig <- varResults[i,c(5:6, 1:4)]
  tmpOut <- cbind(tmpOut, tmpOrig)
  hg19out <- rbind(hg19out, tmpOut)
}

#hg19out_filt <- hg19out[-which(is.na(hg19out$hg19start)),]
write.table(hg19out, "/mnt/DATA5/tmp/kev/misc/20210608liftOver_hg19.avinput",
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


tmpConv <- read.table("/mnt/DATA5/tmp/kev/misc/20210608liftOver.hg19_multianno.txt",
                      sep = "\t", stringsAsFactors = FALSE, header = FALSE,
                      skip = 1)

colnames(tmpConv) <- c("hg19_chr", "hg19_start", "hg19_end", "Ref", "Alt",
                       "Func.refGene", "Gene.refGene", "GeneDetail.refGene",
                       "ExonicFunc.refGene", "AAChange.refGene", 
                       "SampleName", "mm10_chr", "mm10_start", "mm10_end")

### use hg19 annovar database in dir /home/kevhu/programs


write.table(tmpConv, "/mnt/DATA5/tmp/kev/misc/20210608liftOverColNames.hg19_multianno.txt",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
