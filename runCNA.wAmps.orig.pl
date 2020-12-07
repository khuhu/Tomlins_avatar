#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $gcProg = "/home/hovelson/scripts/GCcalc/GCcalculator.py";
#my $cnProg = "/home/hovelson/scripts/cna-analysis-17Dec13.R";
my $cnProg = "/home/hovelson/scripts/cna-analysis-17Dec13.outputAmps.R";
#my $cnProg = "/home/hovelson/scripts/cna-analysis-20140401.R";
#my $cnProg = "/home/hovelson/scripts/cna-analysis-20140425.expandWidth_forCCP.R";
#my $cnProg = "/home/hovelson/scripts/cna-analysis-17Dec13.expandYaxis.R";
my ($index,$conf,$gcbed,$odir,$normtxt);
my $flag = 0;
my $params = GetOptions("idx=s" => \$index,
												"gc=s" => \$gcbed,
												"out=s" => \$odir,
												"flag=i" => \$flag,
												"normals=s" => \$normtxt);

#unless ($index && $gcbed && $odir && $normtxt) {
unless ($index && $gcbed && $odir) {
	die "Usage: ./runCNA.pl --idx <coverage_Index> --gc <gc_bed_file> --out <output_dir> --normals <comma-separated-list_of_normals||or||file_wList_of_NormIDs> --flag <optional: set to 1 if gc bed has ampID in col4>\n";
}
#$gcbed = $ARGV[0];
#$index = $ARGV[1];
#$odir = "./";
#$flag = 0;

my $c1 = "mkdir -p $odir";print "$c1\n";`$c1`;
my $amp_pre = "AMP";
my $ampGCinput = "$odir/amplicon.GCinput.txt";
my $ampCOVinput = "$odir/amplicon.combinedCoverage.input.txt";
my $ampIdx = "$odir/amplicon.index.txt";
my $sampleInfo = "$odir/sampleInfo.input.txt";
my $log = "$odir/copyNumberAnalysis.log";
#my @tumors = ("NEW_PR_pool","BL_195","BL_176");
my @norms;
if ($normtxt =~ /,/) {
	@norms = split(",",$normtxt);
}
else {
	open(NF,"$normtxt") || die "cannot open normal ID list file!\n";
	while(<NF>) {
		chomp;
		push(@norms,$_);
	}
}

# set tumor hash
=cut
my %tumh;
foreach my $t (@tumors) {
	$tumh{$t} = 1;
}	
=cut

# set tumor hash
my %normh;
foreach my $n (@norms) {
	$normh{$n} = 1;
}

# amplicon hash
my %amph;

# Create amplicon GC file for CNA input
open(BIN,"$gcbed") || die "cannot open input bed file!\n";
open(OUT,">$ampGCinput") || die "cannot open amplicon GC for writing!\n";
open(AMPIDX,">$ampIdx") || die "cannot open amplicon index for writing!\n";
my %amphash;
my $ampcount = 0;
print OUT "AmpliconId\tAmpliconIndex\tChromNum\tStartPos\tEndPos\tGene\tNumGC\tLength\tGC\n";
print AMPIDX "Attributes\tAmpliconId\tAmpliconIndex\tChromNum\tStartPos\tEndPos\tGene\tNumGC\tLength\tGC\n";
while(<BIN>) {
	chomp;
	my $l = $_;
	my @f = split("\t",$l);
	my ($chrraw,$spos,$epos,$amp,$gcbase,$totbase,$gcpct,$generaw);
	unless ($#f == 7) {
		my $msg = "\nError: input bed file does not contain required number of columns at line $ampcount:\n$_\n";
		$msg .= "\n\tDoes your bed file have track info in the header?\n\tHave you generated GC content info for this bed file?";
		$msg .= "\n\nTo obtain GC content info for target amplicons, run:\n\tpython $gcProg your_bed_file\n\n";
		die "$msg";
	}
	elsif ($f[5] =~ /[+-]+/) {
		my $msg = "\n\tError: input bed file does not contain total reads for amplicon:\n\t$l\n\n";
		$msg .= "\n\nTo obtain GC content info for target amplicons, run:\n\tpython $gcProg your_bed_file\n\n";
		die "$msg";
	}
	else {
		$ampcount++;
		($chrraw,$spos,$epos,$amp,$gcbase,$totbase,$gcpct,$generaw) = @f;
		my ($gene,$ampID,$ampName) = ("NA","NA","NA");
		# parse generaw field
		if ($generaw =~ /GENE_ID=(.+)_[a-z]+$/) {
			$gene = $1;
		}
		elsif ($generaw =~ /GENE_ID=(.+)$/) {
			$gene = $1;
		}
		else {
			$gene = $generaw;
		}
		# check to see if set indicating col4 is amplicon ID
		$ampID = "$amp_pre\_$ampcount";
		if ($flag == 1) {
			$ampID = $amp;
		}
		$amphash{$ampID} = 1;
		# strip 'chr' string from chr
		my $chr = $chrraw;
		if ($chrraw =~ /chr(.+)/) {
			$chr = $1;
		}
		if ($chr eq "X") {
			$chr = 23;
			$chrraw = "chr23";
		}
		# Ion Torrent amp coverage output uses 1-based amplicon coverage output!
		my $spos1 = $spos+1;
		$amph{"$chrraw:$spos1:$epos:$amp"} = $ampID;
		print OUT "$ampID\t$ampcount\t$chr\t$spos\t$epos\t$gene\t$gcbase\t$totbase\t".sprintf("%0.4f",$gcpct)."\n";
		print AMPIDX "$amp\t$ampID\t$ampcount\t$chr\t$spos\t$epos\t$gene\t$gcbase\t$totbase\t".sprintf("%0.4f",$gcpct)."\n";
	}
}
print "amplicon ($ampIdx) & GC ($ampGCinput) file created for CNA input!\n";

# get total number of samples to parse
my $t = "wc -l $index";
my $string = substr(`$t`,0,-1);
my $ns = 0;
if ($string =~ /([0-9]+)\s+.+/) {
	$ns = $1;
}

# set key/array-of-length-ns pair in hash of arrays for each amplicon
my %ampRdCts;
foreach my $k (keys %amph) {
	@{$ampRdCts{$k}}[$ns-1] = 0;
}

my @tumors;
# parse amplicon coverage for each sample
open(SAMPINFO,">$sampleInfo") || die "cannot open sample info file for writing!\n";
print SAMPINFO "DNA\tSample\tSampleClass\n";
my @samps;
my $sct = 0;
open(IDX,"$index") || die "cannot open amplicon index file!\n";
while(<IDX>) {
	chomp;
	my ($samp,$bc,$ampfile) = split("\t",$_);
	push(@samps,$samp);
	# determine whether sample is tumor or normal for sample info file output
	if (exists $normh{$samp}) {
		print SAMPINFO "$samp\t$samp\tNormal\n";
	}
	else {
		print SAMPINFO "$samp\t$samp\tTumor\n";
		push(@tumors,$samp);
	}	
	print "Processing amplicon coverage for $samp: $ampfile...";
	open(TMP,"$ampfile") || die "cannot open amplicon file ($samp:$ampfile)!\n";
	while(<TMP>) {
		chomp;
		if ($_ =~ /^contig_id/) {}
		else {
			my ($contig_id,$contig_srt,$contig_end,$region_id,$attributes,$gc,$overlaps,$fwd_e2e,$rev_e2e,$total_reads,$fwd_reads,$rev_reads) = split("\t",$_);
			if ($contig_id eq "chrX") {
				$contig_id = "chr23";
			}			
			my $tk = "$contig_id:$contig_srt:$contig_end:$region_id";		
			my $ampID = "NA";
			# Note: Ion Torrent amp coverage info uses 1-based start positions!
			# Start-positions here are 1-based
			if (exists $amph{$tk}) {
				$ampID = $amph{$tk};
				# Read counts (fwd + rev) corresponding to this amplicon for this sample
				@{$ampRdCts{$tk}}[$sct] = $total_reads;
			}
		}
	}
	close(TMP);
	$sct++;
	print "done!\n";
}

# Print out combined coverage file	
open(OUTAMP,">$ampCOVinput") || die "cannot open combined coverage output file!\n";
print OUTAMP "AmpliconId,".join(",",@samps)."\n";
foreach my $k2 (keys %ampRdCts) {
	my @dat = @{$ampRdCts{$k2}};
	if (scalar @dat < 17) {
		#print "$k2\n";
	}
	my $ampID = $amph{$k2};
	print OUTAMP "$ampID,".join(",",@dat)."\n";
}
print "combined amplicon coverage file ($ampCOVinput) created!\n";


# Run CN analysis
#	Usage: [Rscript] AMPLICON_INFO SAMPLE_INFO READ_COUNTS [SNP_FILE] [OPTIONS]	
#my $c = "Rscript $cnProg $ampGCinput $sampleInfo $ampCOVinput 2> $log";
my $c = "Rscript $cnProg $ampGCinput $sampleInfo $ampCOVinput --min-amplicons-per-gene=3 2> $log";
#my $c = "Rscript $cnProg $ampGCinput $sampleInfo $ampCOVinput --min-amplicons-per-gene=10 > $log";
print "\n\nRunning '$c' ...\n";
`$c`;

# Move output files to correct output directory (temporary - need to fix R script to output correctly
my $output1 = "mv out.*pdf $odir";
my $output2 = "mv calls.*txt $odir";
my $output3 = "mv ampLevel*txt $odir";
print "$output1\n$output2\n$output3\n";
`$output1`;
`$output2`;
`$output3`;

# combine CN calls into single doc
my $combine = "(head -1 $odir/calls.$tumors[0].txt; cat $odir/calls.*txt | grep -v ^Gene) > $odir/combined.allSamples.CNcalls.txt";
`$combine`;

