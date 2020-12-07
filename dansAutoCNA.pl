#!/usr/bin/perl
use strict;
use warnings;
use File::stat;

# Check last mod time of alldone.OK
my $alldone = "/mnt/DATA2/share/\@CN_CALLS/alldone.OK";
my $aref = stat($alldone);
my @a = @$aref;
my $ltime = localtime($a[9]);
my ($lsec,$lmin,$lhour,$lmday,$lmon,$lyear,$lwday,$lyday,$lisdst) = localtime($a[9]);
$lyear += 1900;
$lmon += 1;
my $last_mod_time = "$lyear$lmon$lmday";

# Current date/time
my $ctime = localtime(time);
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon += 1;
my $curtime = "$hour:$min:$sec";
my $datestamp = "$year$mon$mday";

# Assign paths/dirs
my $cnadir = "/mnt/DATA2/share/\@CN_CALLS";
my $logdir = "/mnt/DATA2/share/\@CN_CALLS/logs";
my $annodir = "/mnt/DATA/annotated_vcfs/anno";
my $beddir = "/mnt/DATA/analysis/data/copy_number/normals";
my $gcdir = "/home/hovelson/lists/bed";
my $datdir = "/mnt/DATA/croncopy/";
my $toskip = "/home/hovelson/lists/RunsToSkipForAnno.20150812.txt";
my $annoscript = "/home/hovelson/scripts/runCNA.wAmps.forAuto.pl";
my $CCPwg = "/home/hovelson/scripts/CCP.geneLevel.forAuto.R";
my $annolog = "$logdir/cnaLog.$datestamp.log";
my $pipelinelog = "$logdir/cnaMakeLog.$datestamp.log";
my $masterAmpCovIndex = "$cnadir/master.ampCovIndex.txt";
my $oldMasterAmpCovIndex = "$logdir/master.ampCovIndex.$datestamp.txt";

if (-e $masterAmpCovIndex) {
	my $c = "cp $masterAmpCovIndex $oldMasterAmpCovIndex";
	`$c`;
}

# set bed hashes for all DNA panels
my %normFhash = (
		"OCP_20150630_designed" => "$beddir/OCP1c/ampliconCovIdx.OCP1c_normals.n28.txt",
		"OCP_20130724_designed" => "$beddir/OCP1/normals.n2.20140811.txt",
		"OCP2_20131028_designed" => "$beddir/OCP2/normals.n4.20140317.index",
		"OCP3_20140522_designed" => "$beddir/OCP3/n4.OCP3_normals.20140915.txt",
		"CCP" => "$beddir/CCP/n11.CCPnormals.20160226.txt",
);

my %normLhash = (
		"OCP_20150630_designed" => "$beddir/OCP1c/ampliconCovIdx.OCP1c_normals.IDlist.txt",
		"OCP_20130724_designed" => "$beddir/OCP1/normals.n2.20140811.IDlist.txt",
		"OCP2_20131028_designed" => "$beddir/OCP2/normals.n4.20140317.IDlist.txt",
		"OCP3_20140522_designed" => "$beddir/OCP3/n4.OCP3_normals.20140915.IDlist.txt",
		"CCP" => "$beddir/CCP/n11.CCPnormals.IDlist.txt",
);

my %bedhash = (
    "OCP_20150630_designed" => "$gcdir/OCP.20150630.designed.noTrack.GCgenes.bed",
    "OCP_20130724_designed" => "$gcdir/OCP.20130724.designed.noTrack.GC.geneName.bed",
    "OCP2_20131028_designed" => "$gcdir/OCP2.20131028.designed.noTrack.GC.bed",
    "OCP3_20140522_designed" => "$gcdir/OCP3.20140506.designed.noTrack.GC.minusNF1probe.20140619.bed",
    "CCP" => "$gcdir/CCP.noTrack.GC.bed",
);

# open log file
open(LOG,">>$annolog") || die "cannot open log file!\n";
open(STDERR,">>$annolog") || die "cannot open log file!\n";
print LOG "\n------------------------------\n";
print LOG "(last alldone.OK timestamp: $ltime)\n";
print LOG "(current time: $ctime)\n";

# open master amp cov index file
open(MASTER,">$masterAmpCovIndex") || die "cannot open master amp cov index file!\n";
print MASTER "SAMPLE\tBARCODE\tAMP_COV_FILE\tBED\tMACHINE\tREPORT\tDIRECTORY\tSAT_REP_DIR\tRUN\n";
# create Makefile for copynumber analysis
my %cnhash;
open(IN,"$annodir/annoIndex.master.txt") || die "cannot open master anno Index!\n";
while(<IN>) {
	chomp;
	if ($_ =~ /^SAMPLE/) {}
	else {
		my $cmd = "";
		my ($sample_,$bc,$run,$report,$dir,$bed,$machine) = split("\t",$_);
	
		# fix CCP bed file name
		if ($bed =~ /CCP_20131001_designed/) {
			$bed = "CCP";
		}

		### (1) First transform id's to 8-digit suffixes
		$sample_ = uc($sample_);
		my $sample = $sample_;
		if ($sample_ =~ /^([A-Z]{2})[-_]{1}([1-9]+.+)/) {
			my $prelet = $1;
			my $dig = $2;
			my $ndig = length($dig);
			my $nzer = 5 - $ndig;
			my $zertxt = "0" x $nzer;
			$sample = "$prelet-$zertxt"."$dig";
		}
		elsif ($sample_ =~ /^([A-Z]{2})([1-9]+.+)/) {
			my $prelet = $1;
			my $dig = $2;
			my $ndig = length($dig);
			my $nzer = 5 - $ndig;
			my $zertxt = "0" x $nzer;
			$sample = "$prelet-$zertxt"."$dig";
		}
		elsif ($sample_ =~ /^([A-Z]{2})([1-9]+)/) {
			my $prelet = $1;
			my $dig = $2;
			my $ndig = length($dig);
			my $nzer = 5 - $ndig;
			my $zertxt = "0" x $nzer;
			$sample = "$prelet-$zertxt"."$dig";
		}
		elsif ($sample_ =~ /^([0-9]+)/) {
			$sample = "SAT_CNA-$1";
		}
		elsif ($sample_ =~ /.+DNA/) {
			$sample = $sample_;
		}
		else {
			#print "$sample__ || $sample_ || $sample\n";
		}
		# for R output
		my $sample2 = $sample;
		$sample2 =~ s/-/_/g;

		
		### (2) Next, check if already analyzed -- if not, cue for analysis
		my $repdir = "$cnadir/$machine/$report";
		my $odir = "$cnadir/$machine/$report/$sample";
		my $okay = "$odir/$sample.$bed.cna.OK";
		my $idx = "$odir/$sample.$bed.idx.txt";
		my $nidx = "$odir/$sample.$bed.idx.wNorms.txt";
		my $machdir = "/mnt/DATA/croncopy/$machine/$report/coverage";
		my $samplog = "$odir/$sample.$bed.log";
		#open(STDERR,">>$odir/$sample.$bed.$datestamp.log") || die "cannot open sample-dir log file!\n";
		
		# print out to master amp cov index file
		my $oampcovcmd = "ls $machdir/$bc*amplicon*";
    my $oampcov = `$oampcovcmd`;
    $oampcov = substr($oampcov,0,-1);
		if (-e $oampcov) {
			print MASTER "$sample2\t$bc\t$oampcov\t$bed\t$machine\t$report\t$dir\t$repdir\t$run\n";
		}

		unless (-e $okay) {
			if (-e $normFhash{$bed}) {
				unless (-d $repdir) {
					my $mdcmd = "mkdir -p $repdir; chmod 777 $repdir";
					`$mdcmd`;
				}
				unless (-d $odir) {
					my $mdcmd2 = "mkdir -p $odir; chmod 777 $odir";#print "$mdcmd\n";
					`$mdcmd2`;
				}
				my $ampcovcmd = "ls $machdir/$bc*amplicon*";
				my $ampcov = `$ampcovcmd`;
				$ampcov = substr($ampcov,0,-1); 
				unless (length($ampcov) <= 0) {
					if (-e $ampcov) {
						$cmd .= "$okay:";
						$cmd .= "\n\tawk \'BEGIN \{print \"$sample2\\t$bc\\t$ampcov\"\}\' > $idx";
						$cmd .= "\n\tcat $idx $normFhash{$bed} > $nidx";
						$cmd .= "\n\tperl $annoscript --idx $nidx --gc $bedhash{$bed} --out $odir --normals $normLhash{$bed} 2> $samplog";
						if ($bed == "CCP") {
							my $callf = "$odir/calls.$sample2.txt";
							if (-e $callf) {
								$cmd .= "\n\tRscript $CCPwg $sample calls.$sample2.txt $odir 2> $samplog";
							}
						}
						$cmd .= "\n\ttouch $okay";
						$cnhash{$okay} = $cmd;
#					print "$cmd\n";
					}
				}
			}	
		}
	}
}
open(MAKE,">$cnadir/Makefile") || die "cannot open cna makefile!\n";
print MAKE "all:";
foreach my $k (sort keys %cnhash) {
	print MAKE " $k";
}
print MAKE "\n\ttouch $alldone\n\n";
foreach my $k (sort keys %cnhash) {
	print MAKE "$cnhash{$k}\n\n";
}

# Run batch CNA
print LOG "($curtime) Running CNA batch calling...\n";
my $mfile = "/mnt/DATA2/share/\@CN_CALLS/Makefile";
if (-e $mfile) {
  my $cmd = "make -f $mfile all -j 20 2> $annolog";
  sleep(5);
  print LOG "($curtime) $cmd\n";
  `$cmd`;
}
		

							

