
setwd("/mnt/DATA4/kevhu/tmp2/test/")
trimmedFiles <- system('find . -name "*anno_UM.full.trim.txt"',intern = TRUE)
trimmedFiles <- sub("\\./","/mnt/DATA4/kevhu/tmp2/test/",trimmedFiles)




listOfNames <- c("649-CTC_R3","RaoCS1_LftEye","RaoCS1_RtEye","646-CTC_R6","650-CTC_R16")
listOfBarcodes <- c("IonXpress73", "IonXpress_070","IonXpress_071","IonXpress_072","IonXpress_074")
OCPNames <- c("RaoCS1_LftEye","RaoCS1_RtEye")
pangGUnames <- c("649-CTC_R3","646-CTC_R6","650-CTC_R16")


newMat <- NULL 
finalMat <- NULL
for(i in seq_along(trimmedFiles)){
  a <-  read.table(trimmedFiles[i],sep = "\t", stringsAsFactors = FALSE)
  numRep <- nrow(a)
  numRep <- numRep - 1 
  binds <- c(listOfNames[i], listOfBarcodes[i])
  binds2 <- NULL
  labels <- c("Sample","Barcode")
  for(j in 1:numRep){
    binds2 <- rbind(binds2, binds)
  }
  if(i == 1){
    binds2 <- rbind(labels, binds2)
    newMat <- cbind(binds2,a)
    finalMat <- rbind(finalMat, newMat)
  }
  if(i > 1){
    binds2 <- rbind(labels, binds2)
    newMat <- cbind(binds2,a)
    newMat <- newMat[-1, ]
    finalMat <- rbind(finalMat, newMat)
  }
}

colnames(finalMat) <- finalMat[1,]
finalMat <- finalMat[-1,]
colnames(finalMat)[1:2] <- c("Sample","Barcode")

ocpHayesDat <- finalMat[which(finalMat$Sample %in% OCPNames),]
panGUHayesDat <- finalMat[which(finalMat$Sample %in% pangGUnames),]


write.table(ocpHayesDat ,"./20180806HayesLabVariantsOCP.txt", col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(panGUHayesDat ,"./20180806HayesLabVariantsPanGU.txt", col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)


#dummy <- read.table("/mnt/DATA4/kevhu/tmp2/649-CTC_R3_IonXpress_073.anno_UM.full.trim.txt",sep = "\t", stringsAsFactors = FALSE)


### Part two is for getting the bed staitsitics; main thing is extracting the relevant information from Dan's table
###
###

library(hashmap)
hayesTabOCP <- read.table("/mnt/DATA4/kevhu/tmp2/test/20180806HayesLabVariantsOCP.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(hayesTabOCP)[1:2] <-c("SAMPLE","BARCODE")
hayesTabPanGU <- read.table("/mnt/DATA4/kevhu/tmp2/test/20180806HayesLabVariantsPanGU.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(hayesTabPanGU )[1:2] <-c("SAMPLE","BARCODE")

wholeAnno <- read.table("/mnt/DATA4/kevhu/subsetVarTab.txt", header = TRUE, stringsAsFactors = FALSE)
wholeAnno <- wholeAnno[-which(wholeAnno$BED == "OCP_20130724_designed_noTrack"),]

wholeAnno.ocp <- wholeAnno[which(wholeAnno$BED == "OCP_20130724_designed"),]
wholeAnno.ocp <- rbind(hayesTabOCP[,which(colnames(hayesTabOCP) %in% colnames(wholeAnno.ocp))],
                       wholeAnno.ocp[,which(colnames(wholeAnno.ocp) %in% colnames(hayesTabOCP))])

wholeAnno.gu <- wholeAnno[which(wholeAnno$BED == "WG_IAD127899_20170720_designed"),]
wholeAnno.gu <- rbind(hayesTabPanGU[,which(colnames(hayesTabPanGU) %in% colnames(wholeAnno.gu))],
                      wholeAnno.gu[,which(colnames(wholeAnno.gu) %in% colnames(hayesTabPanGU))])

unique( wholeAnno$BED_TOT_NS)
bed_tot_ns.ocp <- 665 + 2
bed_tot_ns.gu <- 503 + 3

#below should be calculating the tot obs and pos pct

#function for hash counter
hash_counter <- function(x,y){
  for(i in seq_along(x)){
    if(y$has_key(x[i]) == TRUE){
      counter <- y$find(x[i])
      y$insert(x[i],c(counter + 1))
    }
  }
}


lookup <- paste0(wholeAnno.ocp$Chr,wholeAnno.ocp$Start, wholeAnno.ocp$End, wholeAnno.ocp$Ref)
key <- unique(lookup)
keyVals <- rep(0, length(key))
hash <- hashmap(keys = key, values = keyVals)

hash_counter(lookup, hash)

bed_tot_obs.ocp <- hash$data()
bed_tot_pos.ocp <- hash$data()/bed_tot_ns.ocp

#now for tot_var.pct for ocp
lookup2 <- paste0(wholeAnno.ocp$Chr,wholeAnno.ocp$Start, wholeAnno.ocp$End, wholeAnno.ocp$Ref, wholeAnno.ocp$Alt)
key2 <- unique(lookup2)
keyVals2 <- rep(0, length(key2))
hash2 <- hashmap(keys = key2, values = keyVals2)


hash_counter(lookup2, hash2)

bed_tot_var_pct <- hash2$data()

#creating table for ocp
match_add <- function(x,y){
  dummyVector <- NULL
  for(i in seq_along(x)){
    a <- y[which(names(y) == x[i])]
    print(a)
    dummyVector <- c(dummyVector, a)
  }
  return(dummyVector)
}

hayesOCPkeys.1 <- paste0(hayesTabOCP$Chr,hayesTabOCP$Start, hayesTabOCP$End, hayesTabOCP$Ref)
hayesOCPkeys.2 <- paste0(hayesTabOCP$Chr,hayesTabOCP$Start, hayesTabOCP$End, hayesTabOCP$Ref, hayesTabOCP$Alt)


hayesTabOCP$bedTotObs <- match_add(hayesOCPkeys.1, bed_tot_obs.ocp)
hayesTabOCP$bedTotNs <- rep(bed_tot_ns.ocp, length(hayesOCPkeys.1))
hayesTabOCP$bedTotPosPct <- round(hayesTabOCP$bedTotObs/hayesTabOCP$bedTotNs,digits = 3)
hayesTabOCP$bedTotVarPct <-  round(match_add(hayesOCPkeys.2, bed_tot_var_pct)/hayesTabOCP$bedTotNs, digits = 3)

write.table(hayesTabOCP ,"./20180806HayesLabVariantsOCP2.txt", col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)

### repeating same steps for PanGU
###
###

lookup <- paste0(wholeAnno.gu$Chr,wholeAnno.gu$Start, wholeAnno.gu$End, wholeAnno.gu$Ref)
key <- unique(lookup)
keyVals <- rep(0, length(key))
hash <- hashmap(keys = key, values = keyVals)

hash_counter(lookup, hash)

bed_tot_obs.ocp <- hash$data()
bed_tot_pos.ocp <- hash$data()/bed_tot_ns.ocp

#now for tot_var.pct for ocp
lookup2 <- paste0(wholeAnno.gu$Chr,wholeAnno.gu$Start, wholeAnno.gu$End, wholeAnno.gu$Ref, wholeAnno.gu$Alt)
key2 <- unique(lookup2)
keyVals2 <- rep(0, length(key2))
hash2 <- hashmap(keys = key2, values = keyVals2)
hash_counter(lookup2, hash2)
bed_tot_var_pct <- hash2$data()


hayesGukeys.1 <- paste0(hayesTabPanGU$Chr,hayesTabPanGU$Start, hayesTabPanGU$End, hayesTabPanGU$Ref)
hayesGukeys.2 <- paste0(hayesTabPanGU$Chr,hayesTabPanGU$Start, hayesTabPanGU$End, hayesTabPanGU$Ref, hayesTabPanGU$Alt)

hayesTabPanGU$bedTotObs <- match_add(hayesGukeys.1, bed_tot_obs.ocp)
hayesTabPanGU$bedTotNs <- rep(bed_tot_ns.gu, length(hayesGukeys.1))
hayesTabPanGU$bedTotPosPct <- round(hayesTabPanGU$bedTotObs/hayesTabPanGU$bedTotNs,digits = 3)
hayesTabPanGU$bedTotVarPct <-  round(match_add(hayesGukeys.2, bed_tot_var_pct)/hayesTabPanGU$bedTotNs, digits = 3)

write.table(hayesTabPanGU ,"./20180806HayesLabVariantsPanGU2.txt", col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)



