#!/usr/local/bin/perl

##############################################################################
# Used to obtain full data from logistic regression model
# loops through reference genome and outputs 1 line per base, as long as
# valid covariate data exists
##############################################################################

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils 'pairwise';
use Cwd;
use Benchmark;
use Tie::File;

my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";
my $categ="AT_CG";

# initialize singleton file
my $f_mlist = "$parentdir/output/logmod_data/${categ}_mlist.txt";
open my $mlist, '<', $f_mlist or die "can't open $f_mlist: $!";

our %hash=();
my @fh;
while(<$mlist>){
  chomp;
  #my @line=split(/\t/, $_);
  #my $motif=$line[3];
  #my $key=join("\t", @line[0 .. 1]);
  #my $pcs=join("\t", @line[2 .. $#line]);

  local *FILE;
  open(FILE, ">${categ}_tmp_{$_}.txt") || die;
  push(@fh, *FILE);
  $hash{$_}=*FILE;
}

# initialize singleton file
my $f_positions = "$parentdir/output/logmod_data/${categ}_full.txt";
open my $positions, '<', $f_positions or die "can't open $f_positions: $!";



while(<$positions>){
  chomp;
  my @line=split(/\t/, $_);
  my $motif=$line[3];

  my $file=$hash{$motif};
  print $file "$_\n";

}
