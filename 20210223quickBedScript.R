#!/usr/bin/env Rscript
.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))
library(optparse)
library(stringr)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="path to bed from cov out in the form of /mnt/DATA3/eros_tmp/Auto_user_AUS5-92-P53_HD_297_240/plugin_out/coverageAnalysis_out.431/local_beds/IAH203804_170_Designed.gc.bed", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Need input file (input file).n", call.=FALSE)
}

if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Need output directory (input file).n", call.=FALSE)
}

tmpBed <- read.table(opt$input, sep = "\t", header = FALSE, skip = 1)
#tmpBed <- read.table("/mnt/DATA3/eros_tmp/Auto_user_AUS5-92-P53_HD_297_240/plugin_out/coverageAnalysis_out.431/local_beds/IAH203804_170_Designed.gc.bed",
#                     sep = "\t", header = FALSE, skip = 1)


tmpBed$V7 <- tmpBed$V3 - tmpBed$V2
tmpBed$V8 <- tmpBed$V6/tmpBed$V7
geneColumn <- str_remove(tmpBed$V5, pattern = "\\;.*")
geneColumn <- str_remove(geneColumn, pattern = "^.*\\=")
geneColumn2 <- geneColumn

geneColumn <- str_remove(geneColumn, pattern = "exon.*")
geneColumn <- str_remove(geneColumn, pattern = "Exon.*")
geneColumn <- str_remove(geneColumn, pattern = "^.*\\,")


tmpBed$V9 <- geneColumn
newGcBed <- tmpBed[,c(1:4,6:9)]

fname <- str_remove(opt$input, ".*local_beds\\/")


write.table(newGcBed, file = paste0(opt$outdir, "/",fname),
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)



