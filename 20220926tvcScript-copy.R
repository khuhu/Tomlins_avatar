### 20211005: script should be used to invoke the variantAnalysis.py script
### variant_caller_pipeline.py -b /home/mouseData/bedFiles/IAD202670_167_Designed.gc.bed -s /home/mouseData/mgp.v6.combinedMouseFilt.vcf.gz 
### -z testSample -i /mnt/DATA3/eros_tmp/Auto_user_AUS5-142-MG_cho_20210701_357_353/IonXpress_033_rawlib.bam -r /home/reference/mm10/mm10_amp.fa

library(stringr)
library(optparse)

option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL,
              help="location of bam files", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output location", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dir <- read.table(opt$directory, stringsAsFactors = FALSE)$V1
setwd(dir)
bamFiles <- system('find . -name "*.bam*" | grep -v "plugin" | grep -v "barcodes"| grep -v "bai"', intern = TRUE);
summaryFile <- read.table(system('find . -name "*.bc_summary*"', intern = TRUE), sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)

outdir <-read.table(opt$outdir, stringsAsFactors = FALSE)$V1
setwd(outdir)
for (i in seq_along(summaryFile$Barcode.ID)) {
  tmpBamLocation <- str_replace(bamFiles[grep(summaryFile$Barcode.ID[i], bamFiles)], "\\./", dir)
  sampleName <- str_remove_all(summaryFile$Sample.Name[i], " ")
  if (file.exists(paste0(summaryFile$Sample.Name[i],".vcf"))) {
    next()
  } else{
    cmd <- paste("variant_caller_pipeline.py -b /mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed",
                 "-s /mnt/DATA6/mouseData/mgp.v6.combinedMouseFilt.vcf.gz -z",
                 sampleName, "-i", tmpBamLocation,
                 "-r /home/reference/mm10/mm10_amp.fa -N 10")
    system(eval(cmd))
  }
}


setwd(outdir)
for (i in seq_along(summaryFile$Barcode.ID)) {
  tmpBamLocation <- str_replace(bamFiles[grep(summaryFile$Barcode.ID[i], bamFiles)], "\\./", dir)
  if (file.exists(paste0(summaryFile$Sample.Name[i],"_geno.vcf"))) {
    next()
  } else{
    cmd <- paste("variant_caller_pipeline.py -b /mnt/DATA6/mouseData/bedFiles/IAD202670_167_Designed.gc.bed",
                 "-s /mnt/DATA6/mouseData/20211122_IAD202670_geno.vcf.gz -z",
                 paste0(summaryFile$Sample.Name[i],"_geno") , "-i", tmpBamLocation,
                 "-r /home/reference/mm10/mm10_amp.fa -N 10")
    system(eval(cmd))
  }
}



### only using this as an endpoint so i can process the individual vcf files by themselves
if(!file.exists("summary.txt")){
  write.table(summaryFile, "./summary.txt", sep = "\t", col.names = TRUE, row.names = FALSE,
              quote = FALSE)
}


if(!file.exists("check.txt")){
  system(sprintf("echo 'done' > check.txt"))
}



