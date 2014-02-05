#!/usr/bin/env perl
# genotype_parents+nr.pl
# Mike Covington
# created: 2012-09-21
#
# Description:
#
use strict;
use warnings;
use Getopt::Long;

use FindBin qw($Bin);
use lib "$Bin/../../lib";
use SNPtools::Genotype;

my $usage = <<USAGE_END;

USAGE:
$0
  --id1         Parent 1 ID
  --id2         Parent 2 ID
  --bam1        Parent 1 alignment file (.bam)
  --bam2        Parent 2 alignment file (.bam)
  --fasta       Reference file (.fasta/.fa)
  --seq_list    OPTIONAL: Comma-delimted list of sequence IDs to analyze
                (By default, this list is generated from the bam file header.)
  --out_dir     Output directory [current]
  --threads     Number of threads [1]
  --nr_ratio    Noise Reduction Ratio [0.7]
  --no_nr       Disable Noise Reduction
  --verbose
  --help

USAGE_END

my (
    $par1,    $par2,    $par1_bam, $par2_bam, $fasta_file, $seq_list,
    $out_dir, $threads, $nr_ratio, $no_nr,    $verbose,    $help
);

my $options = GetOptions(
    "id1=s"      => \$par1,
    "id2=s"      => \$par2,
    "bam1=s"     => \$par1_bam,
    "bam2=s"     => \$par2_bam,
    "fasta=s"    => \$fasta_file,
    "seq_list=s" => \$seq_list,
    "out_dir=s"  => \$out_dir,
    "threads=i"  => \$threads,
    "nr_ratio=f" => \$nr_ratio,
    "no_nr"      => \$no_nr,
    "verbose"    => \$verbose,
    "help"       => \$help,
);

die $usage unless $options;
die $usage if $help;
die $usage
  unless defined $par1
  && defined $par2
  && defined $par1_bam
  && defined $par2_bam
  && defined $fasta_file;

my $geno = SNPtools::Genotype->new(
    par1     => $par1,
    par2     => $par2,
    fasta    => $fasta_file,
    seq_list => $seq_list,
    out_dir  => $out_dir,
    threads  => $threads,
    nr_ratio => $nr_ratio,
    verbose  => $verbose,
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

$geno->noise_reduction unless $no_nr;

exit;
