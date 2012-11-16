#!/usr/bin/env perl
# noise_reduction.pl
# Mike Covington
# created: 2012-09-21
#
# Description:
#
use strict;
use warnings;
use Getopt::Long;

use genotyping_commander;

my $usage = <<USAGE_END;

USAGE:
genotype_parents+nr.pl
  --id1        Parent 1 ID
  --id2        Parent 2 ID
  --bam1       Parent 1 alignment file (.bam)
  --bam2       Parent 2 alignment file (.bam)
  --fasta      Reference file (.fasta/.fa)
  --out_dir    Output directory [current]
  --threads    Number of threads [1]
  --verbose
  --help

USAGE_END

my (
    $par1,    $par2,    $par1_bam, $par2_bam, $fasta_file,
    $out_dir, $threads, $verbose,  $help
);

my $options = GetOptions(
    "id1=s"     => \$par1,
    "id2=s"     => \$par2,
    "bam1=s"    => \$par1_bam,
    "bam2=s"    => \$par2_bam,
    "fasta=s"   => \$fasta_file,
    "out_dir=s" => \$out_dir,
    "threads=i" => \$threads,
    "verbose"   => \$verbose,
    "help"      => \$help,
);

die $usage unless $options;
die $usage if $help;
die $usage
  unless defined $par1
  && defined $par2
  && defined $par1_bam
  && defined $par2_bam
  && defined $fasta_file;

my $geno = genotyping_commander->new(
    par1    => $par1,
    par2    => $par2,
    fasta   => $fasta_file,
    out_dir => $out_dir,
    threads => $threads,
    verbose => $verbose,
);

my %parents = (
    "par1" => {
        "id"  => $par1,
        "bam" => $par1_bam,
    },
    "par2" => {
        "id"  => $par2,
        "bam" => $par2_bam,
    },
);

for ( sort keys %parents ) {
    $geno->id( $parents{$_}->{id} );
    $geno->bam( $parents{$_}->{bam} );
    $geno->before_noise_reduction(1);
    $geno->extract_mpileup;
    $geno->genotype;
}

$geno->noise_reduction;

exit;
