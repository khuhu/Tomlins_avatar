#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
#my $tvcdir = "/srv/tvc4.2.3/tvc-4.2.3";
#my $tvcdir = "/srv/tvc4.2.3/tvc-4.2.3-Ubuntu_12.04_x86_64-binary";
#my $tvc = "/srv/tvc4.2.3/tvc-4.2.3-Ubuntu_12.04_x86_64-binary/bin/variant_caller_pipeline.py";
#my $ref = "/srv/tvc4.2.3/ref/hg19/hg19.fasta";
#my $bdir = "/srv/tvc4.2.3/ref/bed/";
#my $pbdir = "/srv/tvc4.2.3/ref/ptrim_bed/";
my $tvcdir = "/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary";
my $tvc = "/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary/bin/variant_caller_pipeline.py";
my $ref = "/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary/ref/hg19/hg19.fasta";
my $bdir = "/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary/ref/bed/";
my $pbdir = "/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary/ref/ptrim_bed/";
my $proton = "/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary/ref/json/somatic_low_stringency_proton.json";

#my $proton = "/srv/tvc4.2.3/ref/json/somatic_low_stringency_proton.json";
#my $proton = "/srv/tvc4.2.3/ref/json/somatic_high_stringency_proton.json";
#my $pgm = "/srv/tvc4.2.3/ref/json/somatic_low_stringency_pgm.json";
#my $pgm = "/srv/tvc4.2.3/ref/json/somatic_high_stringency_pgm.json";
my $pgm = "/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary/ref/json/somatic_low_stringency_pgm.json";
#my $hvcf1 = "/mnt/DATA/analysis/Hardiman/illumina_vars.inIAD34847_Hotspot.preprocess.vcf";
#my $hvcf2 = "/mnt/DATA/analysis/Hardiman/illumina_vars.inIAD41516_Hotspot.preprocess.vcf";
#my $adrenal_hs = "/mnt/DATA/analysis/Adrenal/hs_matrix/adrenal.hotspot_panel_variants.processed.vcf";
my $hvcf1 = "/mnt/DATA/analysis/PCa_multifocal/allVcalls.update.CCP.preprocess.vcf";
my $hvcf2 = "/mnt/DATA/analysis/PCa_multifocal/allVcalls.update.OCP1.preprocess.vcf";
my $hvcf3 = "/mnt/DATA/analysis/PCa_multifocal/allVcalls.update.OCP3.preprocess.vcf";
my $adrenal_hs = "/mnt/DATA/analysis/Adrenal/hs_matrix/v2/adrenal.v2.hotspots.processed.vcf";
my $adrenal_hs_master = "/mnt/DATA/analysis/Adrenal/apcc_apa/bed/adrenal.FDP150_BEDpt1.uniqCalls.sorted.inV1panel.pp.vcf";
my $ut_hs = "/mnt/DATA/analysis/UT/UT.added_vars.preprocess.vcf";
my $cfdna_mo_hs = "/mnt/DATA/analysis/thruPlex/manuscript/OCP/MiOncoSamples_TVC/MiOnco_calls.inOCPreg.preprocess.vcf";
#my $proton = "/srv/tvc4.2.3/ref/json/somatic_high_stringency_proton.json";
#my $pgm = "/srv/tvc4.2.3/ref/json/somatic_high_stringency_pgm.json";
#my $tbed = "/home/hovelson/lists/Horizon_Diagnostics.wPos.20141031.OCPv2.pos_only.bed";
#my $tbed = "/home/hovelson/lists/Horizon_Diagnostics.wPos.20141014.v2.pos_only.bed";
#my $hvcf = "/mnt/DATA/project_data/lifetech_seq_run/DNA_MO_100--VT3ZO/R_2014_06_04_10_55_44_user_C01-1067-Ex715_OCP_MO_100--91780--dna_mo_100--IonXpress_060--55--TVC.vcf";
#my $hvcf = "/home/hovelson/lists/vcfs/OCP_hs2.inHorizonDxregions.vcf";
#my $hvcf = "/home/hovelson/lists/vcfs/OCP_hs2.plusUMforcedCalls.inHorizonDxregions.vcf";
#my $hvcf = "/home/hovelson/lists/vcfs/AcrometrixOncologyHotspot.wheader.vcf";
#my $hvcf = "/home/hovelson/lists/vcfs/MO_Cohort.allVCalls.wheader.pre_processed.vcf";
my $hvcf = "";

my ($ix,$odir,$pfile,$machine);
my $params = GetOptions("idx=s" => \$ix,
                        "odir=s" => \$odir,
												"machine=s" => \$machine); #,
#												"pfile=s" => \$pfile);
unless ($ix && $odir && $machine) {
	die "Usage: ./runTVC.pl --idx <index.txt> --odir <output_dir> --machine <pgm,proton>\n";
}

my %ih;
open(IX,"$ix") || die "cannot open index file $ix!\n";
while(<IX>) {
	chomp;
	my ($samp,$bc,$bam,$panel) = split(" ",$_);
	$ih{"$bc;$samp;$panel"} = "$samp;$bam;$panel";
}


my @o;
my $otxt;
open(MAKE,">Makefile") || die "cannot open makefile!\n";
foreach my $i (sort keys %ih) {
	if ($ih{$i} =~ /^None/) {}
	else {
		my ($bc,$samp,$panel) = split(";",$i);
		my ($samp_,$bam,$panel_) = split(";",$ih{$i});
		my $ibam = $bam;
		my $pre = substr($ibam,0,-4);
		my $ovcf = "$odir/$samp.$pre.vcf";
		my $odir2 = "$odir/$samp.$bc.$panel";
		my $sampID = "$samp.$bc.$panel";
		my $ok = "$odir2/$sampID.OK";
		push(@o,$ok);
		my $bedroot = "";
		#if ($panel =~ /^OCPv1/) {$bedroot = "OCP.20130724.designed";}
		#elsif ($panel =~ /^OCPv2/) {$bedroot = "OCP2.20131028.designed";}
		#elsif ($panel =~ /^OCPv3/) {$bedroot = "OCP3.20140506.designed";}
		#elsif ($panel =~ /^OCP3/) {$bedroot = "OCP3.20140506.designed";}
		#elsif ($panel =~ /^CCP/) {$bedroot = "CCP";}

		if ($panel =~ /^OCPv1/) {$bedroot = "OCP.20130724.designed";}
    elsif ($panel =~ /^OCP_20130724/) {$bedroot = "OCP.20130724.designed";
      # for multifocal PCa - 2016/10/27 DHH
      $hvcf = $hvcf2;
    }
    elsif ($panel =~ /^OCPv2/) {$bedroot = "OCP2.20131028.designed";}
    elsif ($panel =~ /^OCPv3/) {$bedroot = "OCP3.20140506.designed";}
    elsif ($panel =~ /^OCP3/) {$bedroot = "OCP3.20140506.designed";
      # for multifocal PCa - 2016/10/27 DHH
      $hvcf = $hvcf3;
    }
    elsif ($panel =~ /^OCP1c/) {
			$bedroot = "OCP.20150630.designed";
			$hvcf = "$cfdna_mo_hs";
		}
    elsif ($panel =~ /^CCP/) {$bedroot = "CCP";
      # for multifocal PCa - 2016/10/27 DHH
      #$hvcf = $hvcf1;
			$hvcf = "$ut_hs";
    }
 

		elsif ($panel =~ /^IAD_41516/) {
			$bedroot = "IAD41516_47_Submitted-1";
			$hvcf = $hvcf2;
		}
		elsif ($panel =~ /^IAD_34847/) {
			$bedroot = "IAD34847_Submitted-1";
			$hvcf = $hvcf1;
		}
		elsif ($panel =~ /^IAD46903/) {
      $bedroot = "IAD46903_31_Designed";
			#$hvcf = $adrenal_hs;
      $hvcf = $adrenal_hs_master;
    }
    elsif ($panel =~ /^IAD79597_173/) {
      $bedroot = "IAD79597_173_Designed";
      $hvcf = $adrenal_hs_master;
    }

		elsif ($panel =~ /^CHP/) {
			$bedroot = "CHP2_designed_20120806";
		}
		elsif ($panel =~ /^TargetSeq-Exome-50Mb-hg19_revA/) {
			$bedroot = "TargetSeq-Exome-50Mb-hg19_revA";
		}

		my $c1 = "$ok:\n\t";
		my $hlen = length $hvcf;
		if ($hlen > 0) {
			if ($machine =~ /proton/) {	
		# command for running tvc
	#			$c1 = "mkdir -p $odir/$samp; chmod 777 $odir/$samp; $tvc -b $bdir/$bedroot.bed -i $ibam -r $ref -o $odir/$samp -p $proton -B $tvcdir --primer-trim-bed $pbdir/$bedroot.bed / 2>&1 | tee $odir/$samp/$samp.variantCalling.log";
				$c1 .= "mkdir -p $odir2; chmod 777 $odir2; $tvc -i $ibam -r $ref -o $odir2 -p $proton -B $tvcdir -b $bdir/$bedroot.bed -s $hvcf -N 8 / 2>&1 | tee $odir2/$sampID.variantCalling.log";
				$c1 .= "\n\ttouch $ok\n\n";
			}
			else {
	#			$c1 = "mkdir -p $odir/$samp; chmod 777 $odir/$samp; $tvc -b $bdir/$bedroot.bed -i $ibam -r $ref -o $odir/$samp -p $pgm -B $tvcdir --primer-trim-bed $pbdir/$bedroot.bed / 2>&1 | tee $odir/$samp/$samp.variantCalling.log";
				$c1 .= "mkdir -p $odir2; chmod 777 $odir2; $tvc -i $ibam -r $ref -o $odir2 -p $pgm -B $tvcdir -b $bdir/$bedroot.bed -s $hvcf -N 8 / 2>&1 | tee $odir2/$sampID.variantCalling.log";
				$c1 .= "\n\ttouch $ok\n\n";
			}
			$otxt .= $c1;
		}
		else {
			if ($machine =~ /proton/) {
				$c1 .= "mkdir -p $odir2; chmod 777 $odir2; $tvc -b $bdir/$bedroot.bed -i $ibam -r $ref -o $odir2 -p $proton -B $tvcdir --primer-trim-bed $pbdir/$bedroot.bed 2>&1 | tee $odir2/$sampID.variantCalling.log";
				$c1 .= "\n\ttouch $ok\n\n";
			}
			else {
				$c1 .= "mkdir -p $odir2; chmod 777 $odir2; $tvc -b $bdir/$bedroot.bed -i $ibam -r $ref -o $odir2 -p $pgm -B $tvcdir --primer-trim-bed $pbdir/$bedroot.bed 2>&1 | tee $odir2/$sampID.variantCalling.log";
				$c1 .= "\n\ttouch $ok\n\n";
			}
			$otxt .= $c1;
		}
	}
}

print MAKE "all: ".join(" ",@o)."\n\ttouch all.OK\n\n";
print MAKE $otxt;

=pod
  -h, --help            show this help message and exit
  -b BEDFILE, --region-bed=BEDFILE
                        Limit variant calling to regions in this BED file
                        (optional)
  -s HOTSPOT_VCF, --hotspot-vcf=HOTSPOT_VCF
                        VCF file specifying exact hotspot positions (optional)
  -i BAMFILE, --input-bam=BAMFILE
                        BAM file containing aligned reads (required)
  -r REFERENCE, --reference-fasta=REFERENCE
                        FASTA file containing reference genome (required)
  -o OUTDIR, --output-dir=OUTDIR
                        Output directory (default: current)
  -p PARAMFILE, --parameters-file=PARAMFILE
                        JSON file containing variant calling parameters
                        (recommended)
  -B TVCROOTDIR, --tvc-root-dir=TVCROOTDIR
                        Directory path to location of variant caller programs.
                        Defaults to the directory this script is located
  -n NUMTHREADS, --num-threads=NUMTHREADS
                        Set TVC number of threads (default: 12)
  --primer-trim-bed=PTRIM_BED
                        Perform primer trimming using provided BED file.
                        (optional)
  --postprocessed-bam=POSTPROCESSED_BAM
                        Perform primer trimming, storing the results in
                        provided BAM file name (optional)

