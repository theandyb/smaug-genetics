#!/usr/local/bin/perl

##############################################################################
# This script annotates tab-delimited vcf summary files with sequence motif
# and mutation category information, and counts the number of mutable motifs
# in the reference genome
##############################################################################

##############################################################################
#Initialize inputs, options, and defaults and define errors
##############################################################################
use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Path qw(make_path);
use Benchmark;
use FindBin;
# use lib "$FindBin::Bin";
use FaSlice;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

print "Script will run with the following parameters:\n";
for (sort keys %{$config}) {
    say "$_: $config->{$_}";
}

# my $adj = $config->{adj};
my $adj=1;
my $mac = $config->{mac};
my $binw = $config->{binw};
my $data = $config->{data};
# my $bin_scheme = $config->{bin_scheme};
my $bin_scheme = "all";
my $parentdir = $config->{parentdir};
my $count_motifs = $config->{count_motifs};
my $expand_summ = $config->{expand_summ};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef getMotif);

my $chr=$ARGV[0];

##############################################################################
#Process inputs
##############################################################################
my $subseq=2;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

my $bw=$binw/1000;

##############################################################################
# Read in files and initialize outputs
##############################################################################
my $in_path = "/net/bipolar/jedidiah/testpipe/summaries";
my $out_path = "$parentdir/output/${subseq}bp_${bw}k_${mac}_${data}2";
make_path("$out_path");

# print "Getting reference for chr$chr...\n";
# my $f_fasta;
# if($data eq "mask"){
#   $f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
# } else {
#   $f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
# }
#
# my $seq=getRef($f_fasta, $chr);
#
# my $seqlength=length($seq);
# print "chr$chr seqlength: $seqlength\n";
# print "Done\n";

##############################################################################
# Counts possible mutable sites per bin for 6 main categories
# and for local sequence analysis if selected
##############################################################################
if($count_motifs eq "TRUE"){
	my $start_time=new Benchmark;
	print "Counting motifs...\n";
	# print "seqlength: $length\n";

	# my $bin_flag;
	my $bin_out;
	if($bin_scheme eq "fixed"){
		$bin_out = "$out_path/chr$chr.motif_counts_fixed.txt";
	} elsif($bin_scheme eq "band") {
		$bin_out = "$out_path/chr$chr.motif_counts_band.txt";
	} else {
		$bin_out = "$out_path/chr$chr.motif_counts_all.txt";
	}

	open(my $outFH, '>', $bin_out) or die "can't write to $bin_out: $!\n";

  my $fname = "$parentdir/reference_data/human_g1k_v37/chr$chr.fasta.gz";
  if ( -e "$fname$chr.fasta.gz" ) { $fname = "$fname$chr.fasta.gz"; }
  my $fa = FaSlice->new(file=>$fname, size=>5_000_000);

  # print $binFH "CHR\tBIN\tMOTIF\tCOUNT\n";
  my $startpos;
  my $endpos;
	my @motifs;

	if($bin_scheme eq "fixed"){
    my $fixedfile = "$parentdir/reference_data/genome.${bw}kb.sorted.bed";
    open my $fixedFH, '<', $fixedfile or die "$fixedfile: $!";
    readWindows($fixedFH, $outFH, $fname);
    # my $bandno = 1;
    # while(<$fixedFH>){
    #   chomp;
    #   my @line=split(/\t/, $_);
    #   my $chrind=$line[0];
    #
    #   if($chrind eq "chr$chr"){
    #     $startpos = $line[1]+1;
    #     $endpos = $line[2];
    #
    #     my $binseq = $fa->get_slice($chr, $startpos, $endpos);
    #
    #     @motifs = ($binseq =~ /(?=([ACGT]{$subseq}))/g);
    #
    #     writeCounts($_, $bandno, \@motifs, $binFH);
    #     $bandno = $bandno+1;
    #   }
    # }
  } elsif($bin_scheme eq "band") {
      my $bandfile = "$parentdir/reference_data/cytoBand.txt";
      open my $bandFH, '<', $bandfile or die "$bandfile: $!";
      readWindows($bandFH, $outFH, $fname);
      # my $bandno=1;
      # while(<$bandFH>){
      #   chomp;
      #   my @line=split(/\t/, $_);
      #   my $chrind=$line[0];
      #   if($chrind eq "chr$chr"){
      #     $startpos = $line[1]+1;
      #     $endpos = $line[2];
      #
      #     my $binseq = $fa->get_slice($chr, $startpos, $endpos);
      #     my $length=length($binseq);
      #     print "$bandno length: $length\n";
      #     @motifs = ($binseq =~ /(?=([ACGT]{$subseq}))/g);
      #     # print join(", ", @motifs);
      #     writeCounts($_, $bandno, \@motifs, $binFH);
      #     # print BIN "$_\t$countstr\n";
      #     $bandno = $bandno+1;
      #   }
      # }
    }	else {
      my $genome = "$parentdir/reference_data/hg19.genome";
      open my $gFH, '<', $genome or die "can't open $genome: $!";
      readWindows($gFH, $outFH, $fname);
      # my $length;
      # while(<$gFH>){
      #   chomp;
      #   my @line=split(/\t/, $_);
      #   if ($line[0] eq "chr$chr") { $length = $line[1]; last;}
      # }
      #
      #
      # $startpos=1;
      # $endpos=$length;
      # my $binseq = $fa->get_slice($chr, $startpos, $endpos);
  		# @motifs = ($binseq =~ /(?=([ACGT]{$subseq}))/g);
  		# my $bin = 0;
      # writeCounts($chr, $bin, \@motifs, $outFH);
      # print BIN "$chr\t$countstr\n";
	}

	my $end_time=new Benchmark;
	my $difference = timediff($end_time, $start_time);
	print "Done. ";
	print "Runtime: ", timestr($difference), "\n";
}

sub readWindows {
  my $windowFH = shift;
  my $outFH = shift;
  my $fname = shift;
  my $fa = FaSlice->new(file=>$fname, size=>5_000_000);

  my $startpos;
  my $endpos;
  my @motifs;
  my $bandno = 1;
  while(<$windowFH>){
    chomp;
    my @line=split(/\t/, $_);
    my $chrind=$line[0];

    if($chrind eq "chr$chr"){
      $startpos = $line[1]+1;
      $endpos = $line[2];

      my $binseq = $fa->get_slice($chr, $startpos, $endpos);

      @motifs = ($binseq =~ /(?=([ACGT]{$subseq}))/g);

      writeCounts($_, $bandno, \@motifs, $outFH);
      $bandno = $bandno+1;
    }
  }
}

##############################################################################
# Read motif counts from hash table, sum counts symmetric motifs and write out
# counting strategy modified from https://www.biostars.org/p/5143/
##############################################################################
sub writeCounts {
  my $first = $_[0];
	my $bin = $_[1];
	my @motifs = @{$_[2]};
  my $outFH = $_[3];
	my %tri_count=();
	$tri_count{$_}++ for @motifs;

	foreach my $motif (sort keys %tri_count) {
		my $altmotif = $motif;
		$altmotif =~ tr/ACGT/TGCA/;
		$altmotif = reverse $altmotif;

		my $seqp = "$motif\($altmotif\)";

		my $sum;
		if(exists($tri_count{$motif}) && exists($tri_count{$altmotif})){
			$sum=$tri_count{$motif}+$tri_count{$altmotif};
		} elsif(exists($tri_count{$motif}) && !exists($tri_count{$altmotif})) {
			$sum=$tri_count{$motif};
		} elsif(!exists($tri_count{$motif}) && exists($tri_count{$altmotif})) {
			$sum=$tri_count{$altmotif};
		}

    print $outFH "$first\t$bin\t$seqp\t$sum\n";
	}
}
