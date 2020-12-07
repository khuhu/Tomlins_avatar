###Generating the necessary input files for Dan's script
###Mainly used to test out the CN calling on Cho data
###Note to self: I need to make my own copy of both scripts so that I can edit the parts of the scripts for mouse and then redirect them .........may prove to be a slight problem b/c of the paths . we'll see


###Dan's script needs the following arguments: path to the script is /home/hovelson/scripts/runCNA.wAmps.forAuto.pl
# Usage: ./runCNA.pl --idx <coverage_Index> --gc <gc_bed_file> --out <output_dir> --normals <comma-separated-list_of_normals||or||file_wList_of_NormIDs> --flag <optional: set to 1 if gc bed has ampID in col4>\n"

###From above, I know I jsut need to generate the idx file and gc bedfile

###The idx file is not too bad the as it needs 3 fields that are tab-delimited
#(1)the sample name, (2)barcode and (3)path to the amplicon cov file - use full not and relative
###Now the second portion, the GC file ... I think ..... 



library(jsonlite)

### cho sample path: /mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-166-Cho_mouse_1_353_410/plugin_out/coverageAnalysis_out.721
### /mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-130-Hayes_CTC_AC05_Run1_OCPv3_295_328/plugin_out/coverageAnalysis_out.621
###2017/09/27: testing whether or not the two UM1001 libraries are consistent through different sequencing runs
###report1: Auto_user_1Proton-12-OCP_20151027_DABD25812_51_025
###report2: Auto_user_1Proton-13-OCP_20151027_DABD25844_52_027 


###making set for testing outlier algorithm for OCPv3 samples
# /mnt/DATA3/Yoda_data_dump/Auto_user_SAT-94-ARNegPR_OCPV3_170_233/plugin_out/coverageAnalysis_out.410
# /mnt/DATA3/Yoda_data_dump/Auto_user_SAT-98-Young_PR_new_177_242/plugin_out/coverageAnalysis_out.433



###new set of test samples - OCP.20150630.designed.bed
#/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-61-DNA_Aaronssamples_PR_AD_HN_LU_ML_132_186/plugin_out/coverageAnalysis_out.323/ #check
#/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-42-20170324_Intraductal_DNA_RNA_88_100/plugin_out/coverageAnalysis_out.191 #check
#/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-151-AS_UT_PR_OCP_Rerun_335_379/plugin_out/coverageAnalysis_out.672 #check
#/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-145-UT_Cohort_OCP_328_367/plugin_out/coverageAnalysis_out.650/local_beds #check
#/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-146-PR_AS_OCPDNA_PRRNATISSUE_329_369/plugin_out/coverageAnalysis_out.652 #check

###only for mouse files:
#mouseNames <- c("mouse","15788","15790","14827","14829","15104","15737","12260","12118","12119","B6","3867","15622")


#directory <- "/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-96-SHlotanPgu-6_171_265/plugin_out/coverageAnalysis_out.442"
#directory <- "/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-86-UT_DNA_PanGU_161_241/plugin_out/coverageAnalysis_out.401"
directory <- "/mnt/DATA3/Zeus_data_dump/Auto_user_1Proton-117-KI_DNA_panGU_193_323/plugin_out/coverageAnalysis_out.517"

setwd(directory)

barcodes <- system('find . -type d -name "IonXpress*"', intern = TRUE)
#system('find . -type d -name "IonXpress*"')
for(i in seq_along(barcodes)){
  barcodes[i] <- gsub("./","",barcodes[i])
}


jsonFile <- fromJSON("./results.json")



###making one large index file 
finalTable <- NULL
finalTable <- c("SampleName","Barocde","path")
barcodes <- sort(barcodes)

for(i in seq_along(barcodes)){
  #dummyVar <- file.path("jsonFile$barcodes$",barcodes[i],"$`Sample Name`", sep = "")
  #sampleName <- dummyVar
  if(length(jsonFile$barcodes[[i]]$`Uniformity of amplicon coverage`) == 0){
    next()
  }
  sampleName <- jsonFile$barcodes[[i]]$`Sample Name`
  print(sampleName)
  path <- paste(directory, "/",barcodes[i], "/",sep = "")
  print(path)
  setwd(path)
  ampliCovFile <- system('find . -name "*.amplicon.cov.xls" | grep IonXpress',intern = TRUE)
  ampliCovFile <- sub("./","",ampliCovFile)
  ampliCovFile <- paste(path, ampliCovFile, sep = "")
  #print(ampliCovFile)
  combined <- c(sampleName, barcodes[i], ampliCovFile)
  finalTable <- rbind(finalTable, combined)
}

rownames(finalTable) <- NULL
colnames(finalTable) <- NULL
#setwd("/mnt/DATA4/kevhu/testData")
finalTable <- finalTable[-1,]
finalTable <- data.frame(finalTable, stringsAsFactors = FALSE)
#finalTable <- finalTable[-c(4,9),]
#normalIDlist <- finalTable[,1]

write.table(x = finalTable, file = "/home/kevhu/data/201712201AaronKiPanGU.idx.txt",
            sep = '\t',row.names = FALSE, quote = FALSE, col.names = FALSE)


#write.table(x = normalIDlist, file = "/home/kevhu/data/normals/PanGU/ampliconCovIdx.PanGU_UT.n23.IDlist.txt",
#            sep = '\t',row.names = FALSE, quote = FALSE, col.names = FALSE)



