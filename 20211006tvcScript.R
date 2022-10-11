### 20211005: script should be used to invoke the variantAnalysis.py script
### variant_caller_pipeline.py -b /home/mouseData/bedFiles/IAD202670_167_Designed.gc.bed -s /home/mouseData/mgp.v6.combinedMouseFilt.vcf.gz 
### -z testSample -i /mnt/DATA3/eros_tmp/Auto_user_AUS5-142-MG_cho_20210701_357_353/IonXpress_033_rawlib.bam -r /home/reference/mm10/mm10_amp.fa

library(stringr)

dir <- "/mnt/DATA3/eros_tmp/Auto_user_AUS5-120-MG_EFD4_BBN_334_304/"
setwd(dir)
bamFiles <- system('find . -name "*.bam*" | grep -v "plugin" | grep -v "barcodes"| grep -v "bai"', intern = TRUE);
summaryFile <- read.table(system('find . -name "*.bc_summary*"', intern = TRUE), sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)

setwd("/mnt/DATA6/mouseData/tvcOut/Auto_user_AUS5-120-MG_EFD4_BBN_334_304")
for (i in seq_along(summaryFile$Barcode.ID)) {
  tmpBamLocation <- str_replace(bamFiles[grep(summaryFile$Barcode.ID[i], bamFiles)], "\\./", dir)
  cmd <- paste("variant_caller_pipeline.py -b /mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed",
               "-s /mnt/DATA6/mouseData/mgp.v6.combinedMouseFilt.vcf.gz -z",
               summaryFile$Sample.Name[i], "-i", tmpBamLocation,
               "-r /home/reference/mm10/mm10_amp.fa")
  system(eval(cmd))
}

i