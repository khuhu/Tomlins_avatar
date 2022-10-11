library("QDNAseq")
library("matrixStats")
library("ggplot2")
library("R.cache")
#library(Biobase)
#install.packages("phenoData")
source("http://bioconductor.org/biocLite.R")
#biocLite("CGHregions")
#biocLite("BiocUpgrade")
#biocLite("QDNAseq")
library("CGHregions")


# capture error output
emsg <- file("/mnt/DATA/analysis/thruPlex/manuscript/smoothed/errors.all_alpha.txt",open="wt")
#emsg <- file("/mnt/DATA/analysis/thruPlex/umuc/in_silico/dilution/dil_ds/errors.all_alpha.txt",open="wt")
sink(emsg,type="message")
#sink()
#bvec <- c(1,5,10,15,50,100,500,1000)
#bvec <- c(1000,500,100,50,15,10,5,1)
bvec <- c(50,15)
#bvec <- c(1000,500,100)
#bvec <- c(15)
for(i in 1:length(bvec)) {
  #i = 1
  tryCatch ({
    ow<-options(warn=2)
    bsize <- bvec[i]

    ### NORMAL SAMPLES
    # Process normal samples
    path <- "/mnt/DATA/analysis/thruPlex/20141219/bams/normals/male"

    #autos <- c(22)
    #autos <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
    bins <- getBinAnnotations(binSize=bsize)

    # count reads in each bin for each included samples
    readCounts2 <- binReadCounts(bins,path=path,cache=TRUE)
    #highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)
    
    # Apply filters to initial log2 bin counts
    readCountsFiltered <- applyFilters(readCounts2,chromosomes=c("Y"))
    
    # filter without blacklist; keep chrX
    #readCountsFilteredX <- applyFilters(readCounts, residual=FALSE,chromosomes=c("Y"))
    #readCountsFilteredX <- applyFilters(readCounts, residual=FALSE,chromosomes=c("Y"))
    
    # Manually add in chrX data to final data frame
    #fData(readCountsFiltered)$use[fData(readCountsFiltered)$chromosome=="X"] <- 
    #  fData(readCountsFilteredX)$use[fData(readCountsFilteredX)$chromosome=="X"]

    # plot filtered read counts
    #isobarPlot(readCountsFiltered)
    
    # estimate necessary GC content and mappability correction
    readCountsFiltered <- estimateCorrection(readCountsFiltered)
    #noisePlot(readCountsFiltered)
    
    # apply GC Content and mappability correction
    copyNumbers <- correctBins(readCountsFiltered)
    
    # apply median normalization
    copyNumbersNormalized <- normalizeBins(copyNumbers)
    
    # smooth outlier bins
    copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
    normCompCNsmooth <- poolRuns(copyNumbersSmooth,rep("male_compNorm",ncol(copyNumbersSmooth)))
    
    #plot(copyNumbersSmooth)
    
    ########################################

    #tums <- c("vcap_100","vcap_50","vcap_10","vcap_5","vcap_1")
    #atums <- c("v0_n100","v1_n99","v5_n95","v10_n90","v15_n85","v20_n80","v25_n75","v30_n70","v50_n50","v75_n25","v100_n0")
    #tums <- c("v100_n0")
    #tums <- c("vcap_1","vcap_5")
    
    # all patient samples
    #atums <- c("TP1176.IonXpress_013","TP1175.IonXpress_012","TP1193.IonXpress_004","TP1117.IonXpress_005",
    #           "TP1171.IonXpress_012","TP1120.IonXpress_006","TP1173.IonXpress_003","TP1123.IonXpress_007",
    #           "TP1147.IonXpress_009","TP1325.IonXpress_006","TP1330.IonXpress_007","TP1336.IonXpress_008",
    #           "TP1342.IonXpress_012","TP1346.IonXpress_013","TP1349.IonXpress_014","TP1353.IonXpress_015",
    #           "TP1354.IonXpress_016","TP1302.IonXpress_003","TP1318.IonXpress_004","TP1317.IonXpress_005",
    #           "TP1337.IonXpress_009","TP1338.IonXpress_010","TP1340.IonXpress_011")
    #tums <- c("TP1241.IonXpress_001","TP1242.IonXpress_002","TP1244.IonXpress_003","TP1270.IonXpress_005",
    #          "TP1271.IonXpress_006","TP1275.IonXpress_007","TP1283.IonXpress_008","TP1284.IonXpress_009",
    #          "TP1285.IonXpress_010","TP1291.IonXpress_011","TP1292.IonXpress_012","TP1293.IonXpress_013",
    #          "TP1299.IonXpress_014","TP1303.IonXpress_015","TP1319.IonXpress_016")
    
    #tums <- c("TP1319.IonXpress_016")
    #tums <- c("umuc5_100.IonXpress_012","umuc5_50.IonXpress_013","umuc5_10.IonXpress_014","umuc5_5.IonXpress_015")
    #atums <- c("umuc5.IonXpress_012")
    # subset of patient samples (all high-tum content & 2 'no tumor')
    #if(FALSE) {
    #atums <- c("TP1171.IonXpress_012","TP1173.IonXpress_003","TP1330.IonXpress_007","TP1349.IonXpress_014",
    #           "TP1353.IonXpress_015","TP1354.IonXpress_016","TP1337.IonXpress_009","TP1338.IonXpress_010",
    #           "TP1340.IonXpress_011")
    #atums <- c("vcap100.IonXpress_004")
    #atums <- c("TP1242.IonXpress_002","TP1241.IonXpress_001","TP1244.IonXpress_003","TP1270.IonXpress_005",
    #           "TP1271.IonXpress_006","TP1275.IonXpress_007","TP1283.IonXpress_008","TP1284.IonXpress_009",
    #           "TP1285.IonXpress_010","TP1291.IonXpress_011","TP1292.IonXpress_012","TP1293.IonXpress_013",
    #           "TP1299.IonXpress_014","TP1303.IonXpress_015")
    #atums <- c("TP1346.IonXpress_013")

    if (FALSE) {
      #btums <- c("0_001x","0_005x","0_025x","0_01x","0_05x","0_1x")
      btums <- c("0_005x","0_01x","0_1x")
      #btums <- c("0_01x","0_05x")
      tums <- c(NA)
      for (a in 1:length(atums)) {
        for (b in 1:length(btums)) {
          #newtum <- paste(atums[a],"umuc5_inSil",btums[b],sep=".")
          newtum <- paste(atums[a],btums[b],sep=".")
          tums <- c(tums,newtum)
        }
      }
      tums <- tums[!is.na(tums)]
    }
    #a <- seq(0,100,1)
    #b <- seq(100,0,1)
    #if(FALSE) {
    if(FALSE) {
      btums <- c("0_1x","0_05x","0_005x")
      tums <- c("v0_n100.umuc5_dil.0_005x","v0_n100.umuc5_dil.0_01x")
      for(b in 1:length(btums)) {
        ds <- btums[b]
        for(s in 1:100) {
          n <- paste("v",s,"_n",100-s,".umuc5_dil.",ds,sep="")
          tums <- c(tums,n)
        }
      }
    }
    #}
    #tums <- c("TP1337.IonXpress_009.0_01x","TP1337.IonXpress_009.0_05x","TP1337.IonXpress_009.0_1x","TP1337.IonXpress_009.0_25x","TP1337.IonXpress_009.0_5x")
    #tums <- c("TP1337.test")
    #}
    #tums <- c("vcap100.IonXpress_004.0_005x","vcap100.IonXpress_004.0_025x","vcap100.IonXpress_004.0_1x")
    #tums <- c("PDL1002_1.IonXpress_013","PDL1005_1.IonXpress_010","PDL1005_4.IonXpress_011",
    #          "PDL1006_1.IonXpress_014","PDL1007_1.IonXpress_015","RDART006_1.IonXpress_007",
    #          "RDART007_1.IonXpress_008","RDART008_1.IonXpress_009","ULMC118_1.IonXpress_001",
    #          "ULMC185_1.IonXpress_002","ULMC185_2.IonXpress_005","ULMC185_3.IonXpress_006",
    #          "ULMC242_1.IonXpress_003","ULMC244_1.IonXpress_004","ULMC245_1.IonXpress_012")
    
    # Yoda: SATProton-152  
    #tums <- c("TP1034.IonXpress_014","TP1069.IonXpress_015","TP1070.IonXpress_016",
    #          "TP1083.IonXpress_001","TP1103.IonXpress_002","TP1105.IonXpress_003",
    #          "TP1143.IonXpress_004","TP1201.IonXpress_005","TP1205.IonXpress_006",
    #          "TP1216.IonXpress_007","TP1221.IonXpress_008","TP1295.IonXpress_009",
    #          "TP1358.IonXpress_010","TP1367.IonXpress_011","TP1384.IonXpress_012",
    #          "TP1387.IonXpress_013")
    
    # Zeus: 1Proton-34
    #tums <- c("TP1052.IonXpress_009","TP1084.IonXpress_001","TP1085.IonXpress_008",
    #          "TP1096.IonXpress_002","TP1108.IonXpress_007","TP1113.IonXpress_003",
    #          "TP1153.IonXpress_004","TP1210.IonXpress_005","TP1225.IonXpress_006",
    #          "TP1281.IonXpress_013","TP1282.IonXpress_014","TP1320.IonXpress_016",
    #          "TP1371.IonXpress_010","TP1377.IonXpress_015","TP1393.IonXpress_011",
    #          "TP1400.IonXpress_012")
    
    # quick check (2017/05/19, DHH)
    tums <- c("TP1281.IonXpress_013","TP1282.IonXpress_014")
    
    
    #Yoda: SATProton-155
    #tums <- c("TP1031.IonXpress_001","TP1042.IonXpress_002","TP1139.IonXpress_003",
    #          "TP1226.IonXpress_005","TP1228.IonXpress_004","TP1246.IonXpress_006","TP1267.IonXpress_007",
    #          "TP1300.IonXpress_008","TP1318.IonXpress_009","TP1357.IonXpress_010","TP1359.IonXpress_011",
    #          "TP1361.IonXpress_012","TP1388.IonXpress_013","TP1401.IonXpress_014","TP1405.IonXpress_015",
    #          "TP1408.IonXpress_016")
    
    #Yoda: SATProton-156
    #tums <- c("TP1019.IonXpress_001","TP1081.IonXpress_003","TP1098.IonXpress_004",
    #          "TP1151.IonXpress_005","TP1152.IonXpress_002","TP1154.IonXpress_006","TP1182.IonXpress_007",
    #          "TP1183.IonXpress_008","TP1184.IonXpress_009","TP1364.IonXpress_010","TP1378.IonXpress_011",
    #          "TP1398.IonXpress_012","TP1398-pr.IonXpress_013","TP1413-pax.IonXpress_015",
    #          "TP1413-pe24hr.IonXpress_016","TP1413-pe.IonXpress_014")
    
    for(j in 1:length(tums)) {
      #j = 1
      tryCatch ({
        
        # ds_dil bams
        #j = 1
        bambase <- strsplit(tums[j],"\\.")
        dsval <- bambase[[1]][3]
          
        ## TUMOR Samples
        # count reads in each bin for each included samples
        
        #filepath <- paste("/mnt/DATA/analysis/thruPlex/in_silico/mixbam/sort/",tums[j],".mix.sort.bam",sep="")
        filepath <- paste("/mnt/DATA/analysis/thruPlex/sample_bams/Tomlins_1Proton_34/",tums[j],"_rawlib.bam",sep="")
        #filepath <- paste("/mnt/DATA/analysis/thruPlex/umuc/in_silico/dilution/mixam/downsample/",tums[j],".bam",sep="")
				#filepath <- paste("/mnt/DATA/analysis/thruPlex/umuc/in_silico/dilution/dil_ds/",dsval,"/",tums[j],".bam",sep="")
				#filepath <- paste("/mnt/DATA/analysis/thruPlex/vcap/bams/vcap100_ds/",tums[j],".bam",sep="")
        #filepath <- paste("/mnt/DATA/analysis/thruPlex/umuc/in_silico/dilution/mixbam/",tums[j],".mix.sort.bam",sep="")
        
        #tumorpath <- paste("/mnt/DATA/analysis/thruPlex/20150211/bams/",tums[j],sep="")
        #tumorpath <- paste("/mnt/DATA/analysis/thruPlex/in_silico/mixbam/sort",tums[j],"mix.sort",sep="")
        #tumorpath <- paste("/mnt/DATA/analysis/thruPlex/sample_bams/TP1337_ds/",tums[j],sep="")
        #tumorpath <- paste("/mnt/DATA/analysis/thruPlex/sample_bams/TP1337_ds/",tums[j],sep="")
        #readCounts <- binReadCounts(bins,path=tumorpath,cache=TRUE)
        readCounts <- binReadCounts(bins,bamfiles=filepath,cache=TRUE)
        #highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)
        
        # Apply filters to initial log2 bin counts
        #readCountsFiltered <- applyFilters(readCounts)
				readCountsFiltered <- applyFilters(readCounts,chromosomes=c("Y"))
        
        # filter without blacklist; keep chrX
        #readCountsFilteredX <- applyFilters(readCounts, residual=FALSE,chromosomes=c("Y"))
        #readCountsFilteredX <- applyFilters(readCounts, residual=FALSE,chromosomes=c("Y"))
        
        # Manually add in chrX data to final data frame
        #fData(readCountsFiltered)$use[fData(readCountsFiltered)$chromosome=="X"] <- 
        #  fData(readCountsFilteredX)$use[fData(readCountsFilteredX)$chromosome=="X"]
        
        # plot filtered read counts
        #isobarPlot(readCountsFiltered)
        
        # estimate necessary GC content and mappability correction
        readCountsFiltered <- estimateCorrection(readCountsFiltered)
        #noisePlot(readCountsFiltered)
        
        # apply GC Content and mappability correction
        copyNumbers <- correctBins(readCountsFiltered)
        
        # apply median normalization
        copyNumbersNormalized <- normalizeBins(copyNumbers)
        
        # smooth outlier bins
        copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
        #plot(copyNumbersSmooth)
        
        # COMBINE WITH NORMALS
        cnCombined <- combine(copyNumbersSmooth,normCompCNsmooth)
        
        # tumorVNormal
        tumorVNormal <- compareToReference(cnCombined,c(2,FALSE))
        
        # try different alpha values
        #al <- c(0.01,0.001,0.00001,0.0000001,0.0000000001)
        #altxt <- c("1e2","1e3","1e5","1e7","1e10")
        al <- c(0.01,0.00001,0.0000000001)
        altxt <- c("1e2","1e5","1e10")
				#al <- c(0.01,0.0000000001)
				#altxt <- c("1e2","1e10")
        for (k in 1:length(al)) {
          #k = 1
          tryCatch({
            alphset <- al[k]
            
            # Segmentation and plot
            copyNumbersSegmented <- segmentBins(tumorVNormal,alpha=alphset)
            #copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
            copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
            b <- bvec[i]
            t <- tums[j]
            alph <- altxt[k]
            #pdfname1 <- paste("/mnt/DATA/analysis/thruPlex/umuc/in_silico/dilution/mixbam/downsample/",t,".binsize_",b,".alpha_",alph,".segment.png",sep="")
            #pdfname1b <- paste("/mnt/DATA/analysis/thruPlex/umuc/in_silico/dilution/mixbam/downsample/",t,".binsize_",b,".alpha_",alph,".segment.pdf",sep="")
            #pdfname2 <- paste("/mnt/DATA/analysis/thruPlex/umuc/in_silico/dilution/mixbam/downsample/",t,".binsize_",b,".alpha_",alph,".call.png",sep="")
            #pdfname1 <- paste("/mnt/DATA/analysis/thruPlex/manuscript/ds_analysis/umuc_dil_ds_plots/bin",b,"kb/",dsval,"/",t,".binsize_",b,".alpha_",alph,".segment.png",sep="")
            #pdfname1b <- paste("/mnt/DATA/analysis/thruPlex/manuscript/ds_analysis/umuc_dil_ds_plots/bin",b,"kb/",dsval,"/",t,".binsize_",b,".alpha_",alph,".segment.pdf",sep="")
            #pdfname2 <- paste("/mnt/DATA/analysis/thruPlex/manuscript/ds_analysis/umuc_dil_ds_plots/bin",b,"kb/",dsval,"/",t,".binsize_",b,".alpha_",alph,".call.png",sep="")
            pdfname1 <- paste("/mnt/DATA/analysis/thruPlex/manuscript/smoothed/",t,".binsize_",b,".alpha_",alph,".segment.png",sep="")
            pdfname1b <- paste("/mnt/DATA/analysis/thruPlex/manuscript/smoothed/",t,".binsize_",b,".alpha_",alph,".segment.pdf",sep="")
            pdfname1 <- paste("/mnt/DATA/analysis/thruPlex/manuscript/smoothed/",t,".binsize_",b,".alpha_",alph,".call.png",sep="")
            #pdf(width=11,height=5,pdfname1)
#            par(mar=c(0.25,0.25,0.25,0.25))
					  png(filename=pdfname1,height=400,width=1200)
            plot(copyNumbersSegmented)
            dev.off()
            pdf(pdfname1b,height=5,width=10)
            plot(copyNumbersSegmented)
            dev.off()
            
            # make segmented data file
            #cgh <- makeCgh(copyNumbersSegmented)
            #saveRDS(cgh,"/mnt/DATA/analysis/thruPlex/20141219/cghSeg.rds")
            
            # Call aberrations & plot - cannot get to work reliably (always fails)!!! (2015/09/28)
            #sink(file="/mnt/DATA/analysis/thruPlex/20141219/callBins_log.txt",type="output")
            if(FALSE) {
              copyNumbersCalled <- callBins(copyNumbersSegmented)
              #pdf(width=11,height=5,pdfname2)
              png(filename=pdfname1,height=400,width=1200)
              plot(copyNumbersCalled)
              dev.off()
            }
            
            # Export to file
            segstxt <- paste("/mnt/DATA/analysis/thruPlex/manuscript/smoothed/",t,".binsize_",b,".alpha_",alph,".segs.txt",sep="")
            cnstxt <- paste("/mnt/DATA/analysis/thruPlex/manuscript/smoothed/",t,".binsize_",b,".alpha_",alph,".cns.txt",sep="")
            
            #calltxt <- paste("/mnt/DATA/analysis/thruPlex/manuscript/smoothed/",t,".binsize_",b,".alpha_",alph,".calls.txt",sep="")
            #segstxt <- paste("/mnt/DATA/analysis/thruPlex/manuscript/ds_analysis/umuc_dil_ds_plots/bin",b,"kb/",dsval,"/",t,".binsize_",b,".alpha_",alph,".segs.txt",sep="")
            #cnstxt <- paste("/mnt/DATA/analysis/thruPlex/manuscript/ds_analysis/umuc_dil_ds_plots/bin",b,"kb/",dsval,"/",t,".binsize_",b,".alpha_",alph,".cns.txt",sep="")
            #segstxt <- paste("/mnt/DATA/analysis/thruPlex/umuc/in_silico/dilution/cn_calls/",t,".binsize_",b,".alpha_",alph,".segs.txt",sep="")
            #cnstxt <- paste("/mnt/DATA/analysis/thruPlex/umuc/in_silico/dilution/cn_calls/",t,".binsize_",b,".alpha_",alph,".cns.txt",sep="")
            #exportBins(copyNumbersCalled,calltxt,type="calls",format="tsv")
            exportBins(copyNumbersSegmented,segstxt,type="segments",format="igv")
            exportBins(copyNumbersSegmented,cnstxt,type="copynumber",format="igv")
          },error=function(e){options(ow)
                              #sink(type="message")
                              cat("ERROR1 :",conditionMessage(e),"\n")})
        }
      },error=function(e){options(ow)
                          #sink(type="message")
                          cat("ERROR2 :",conditionMessage(e),"\n")
                          })
    }
  
  },error=function(e){options(ow)
                      #sink(type="message")
                      cat("ERROR3 :",conditionMessage(e),"\n")})
}

# reset sink
sink()


# annotation of copy number calls
#cgh <- makeCgh(copyNumbersCalled)
#regions <- CGHregions(cgh)

# annotation info 
#r <- as.data.frame(regions(regions))
#c <- as.data.frame(chromosomes(regions),ncol=1)
#s <- as.data.frame(bpstart(regions),ncol=1)
#e <- as.data.frame(bpend(regions),ncol=1)

#anno <- cbind(c,s,e,r)



