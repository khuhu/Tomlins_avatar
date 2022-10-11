#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Define script & dataset locations where appropriate
my $pardir = "/home/hovelson/";
my $annodir = "/srv/annovar";
my $convert = "$annodir/convert2annovar.pl";
my $anno = "$annodir/annotate_variation.pl";
my $annotab = "$annodir/table_annovar.pl";
my $multi_allelic = "$pardir/scripts/fixMultiAllelicIndels.pl";
my $parseAnnovar = "$pardir/scripts/parseAnnovarOutput.pl";
my $datadir = "$annodir/data/";
my $cfilec = "$pardir/data/ref/cosmic/CosmicCodingMuts_v68.vcf.gz";
my $cfilenc = "$pardir/data/ref/cosmic/CosmicNonCodingMuts_v68.vcf.gz";
my $dbfile = "$pardir/data/ref/dbSNP/dbSNP.all.v135.vcf.gz";
my $ref = "/mnt/DATA/reference/hg19/hg19.fasta";
my ($in,$opre);
my $params = GetOptions("in=s" => \$in,
												"opre=s" => \$opre);

unless ($in) {
	die "Usage: ./runAnnovar.pl --in input.vcf --opre output.prefix\n";
}
my $pre;
if (length($opre) > 0) {
	$pre = $opre;
}
else {
	if ($in =~ /\S+gz$/) {
		$pre = substr($in,0,-7);
	}
	elsif ($in =~ /\S+vcf$/) {
		$pre = substr($in,0,-4);
	}
}

# Fix multi-allelic indel entries - remove unnecessary fix on 2016/04/12, with updated version of annovar
#my $fixMA = "perl $multi_allelic $in > $pre.fixMA.vcf; cp $in $pre.origMA.vcf.gz; mv $pre.fixMA.vcf $pre.vcf";
#print "$fixMA\n";`$fixMA`;

# bgzip and tabix files for bcftools processing; split to one variant per row
my $split;
if ($in =~ /\S+gz$/) {
	$split = "gunzip -c $in > $pre.orig.vcf; bgzip -f $pre.orig.vcf; tabix $pre.orig.vcf.gz; bcftools norm -m-both -o $pre.split.vcf $pre.orig.vcf.gz 2> $pre.log";
} elsif ($in =~ /\S+vcf$/) {
	$split = "cp $pre.orig.vcf; bgzip -f $pre.orig.vcf; tabix $pre.orig.vcf.gz; bcftools norm -m-both -o $pre.split.vcf $pre.orig.vcf.gz 2> $pre.log";
}
print "$split\n"; `$split`;

# left-align
my $leftalign = "bcftools norm -f $ref -o $pre.split.lalign.vcf $pre.split.vcf; cp $pre.split.lalign.vcf $pre.vcf";
print "$leftalign\n"; `$leftalign`;

=cut
my $fixMAvcf = "cp $in $pre.origMA.vcf.gz; cp $in $pre.vcf";
my $fixMAgz = "cp $in $pre.origMA.vcf.gz; gunzip -c $in > $pre.vcf";
if ($in =~ /\S+gz$/) {
	print "$fixMAgz\n";`$fixMAgz`;
} else {
	print "$fixMAvcf\n";`$fixMAvcf`;
}
=cut

# Convert to annovar format
my $avinput = "$pre.avinput";
my $c1 = "perl $convert -format vcf4 -outfile $avinput -include -withzyg $pre.vcf 2> $pre.log";
print "$c1\n";`$c1`;

# Annotate with refGene, Cosmic, ESP/1000G
my $c2 = "perl $annotab $avinput $datadir/ -build hg19 -out $pre -protocol refGene,esp6500si_all,1000g2012apr_all,snp135,ljb_all,cosmic68,clinvar_20150330,genomicSuperDups,phastConsElements46way -operation g,f,f,f,f,f,f,r,r -nastring NA -remove -otherinfo 2> $pre.log";
print "$c2\n";`$c2`;

# Parse Annovar output
my $annoTbl = "$pre.hg19_multianno.txt";
my $c3 = "perl $parseAnnovar $pre.vcf $annoTbl 2> $pre.log";
print "$c3\n";`$c3`;

