library(GenomicRanges)
library(vcfR)
tmpBed <- read.table("/home/kevhu/data/bedFiles/IAD202670_167_Designed.gc.noChr.bed", sep = "\t",
                     header = FALSE, stringsAsFactors = FALSE)

black6 <- read.vcfR("/mnt/DATA6/kevin_recovery/ComprehensiveMouse/vcfs/6J_100_500.vcf")


snpBed <- tmpBed[grep("SNP", tmpBed$V8),]
snpBed2 <- snpBed
snpBed2_GR <- GRanges(seqnames = snpBed2$V1,
                      IRanges(start = snpBed2$V2, end = snpBed2$V3))

snpPos_GR <- GRanges(seqnames = black6@fix[,1],
                     IRanges(start = as.numeric(black6@fix[,2]), end = as.numeric(black6@fix[,2])))


keepIdx <- queryHits(findOverlaps(snpPos_GR, snpBed2_GR))

black6_2 <- black6
black6_2@gt <- black6_2@gt[keepIdx, ]
black6_2@fix <- black6_2@fix[keepIdx, ]

vcfR::write.vcf(black6_2, file = "/mnt/DATA5/tmp/kev/misc/20211122_IAD202670_geno.vcf.gz")

### renamed to 20211120_IAD202670_geno.vcf.gz for brevity
