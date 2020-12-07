### This script will take a single input fastq file named "crispr.fastq" and sort it based on barcode 
### information which must be supplied in a csv file named "Barcode table.csv". Both files should be placed in 
### the working directory. "Barcode table.csv" is a 3 column file which requires the headers "Clone", "FWD" and "REV"
### for: the name of the clone, the sequence of the forward barcode and the sequence of the reverse barcode, respectively
### For example, one row may look like this: A3,CTTGACACCGC,CCTGGTTGTC.

# load ShortRead package. If not installed then run: 
# source("http://bioconductor.org/biocLite.R")
# biocLite("ShortRead")

library(ShortRead)

# read in fastq file
reads <- readFastq("crispr.fastq")

# create dataframe of plate IDs
Barcode_Table <- read.csv(file="Barcode table.csv")
Barcode_Table$FWD <- as.character(Barcode_Table$FWD)
Barcode_Table$REV <- as.character(Barcode_Table$REV)

for (i in seq_len(nrow(Barcode_Table)))
{
  Barcode1 <- Barcode_Table$FWD[i]
  Barcode2 <- Barcode_Table$REV[i]
  name <- Barcode_Table$Clone[i]
  
  # create filters
  RowFilter <- srFilter(function(x){
    substr(sread(x),1,nchar(Barcode1))==Barcode1
  },name="Row Filter")
  
  ColumnFilter <- srFilter(function(x){
    substr(sread(x),1,nchar(Barcode2))==Barcode2
  },name="Column Filter")
  
  # subset row and column reads
  RowReads <- reads[RowFilter(reads)]
  ColumnReads <- reads[ColumnFilter(reads)]
  
  # create reverse complement objects of row and column reads
  RowReads_Reverse <- ShortReadQ(reverseComplement(sread(RowReads)), FastqQuality(reverse(quality(quality(RowReads)))), id(RowReads))
  ColumnReads_Reverse <- ShortReadQ(reverseComplement(sread(ColumnReads)), FastqQuality(reverse(quality(quality(ColumnReads)))), id(ColumnReads))
  
  # Search Row Reads Rev for Column Reads and vice versa
  Clone_1 <- RowReads_Reverse[ColumnFilter(RowReads_Reverse)]
  Clone_2 <- ColumnReads_Reverse[RowFilter(ColumnReads_Reverse)]
  
  # Append two directions
  Clone <- append(Clone_1,Clone_2)
  
  # trim barcodes; change the number to the length of the longest barcode
  Clone <- narrow(Clone, 13)
  Clone <- ShortReadQ(reverseComplement(sread(Clone)), FastqQuality(reverse(quality(quality(Clone)))), id(Clone))
  Clone <- narrow(Clone, 13)
  
  # write fastq
  writeFastq(Clone, file=(sprintf("%s.fastq", name)))
}


