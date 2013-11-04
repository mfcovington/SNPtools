#!/usr/bin/env perl
# genotype_all.pl
# Mike Covington
# created: 2011-12-05
#
# Description:
#
use strict;
use warnings;
use Getopt::Long;

use SNPtools::Genotype;

my $usage = <<USAGE_END;

USAGE:
$0
  --fasta      Reference file (.fasta/.fa)
  --out_dir    Output directory [current]
  --threads    Number of threads [1]
  --verbose
  --help

USAGE_END

my ( $id, $par1, $par2, $bam_file, $fasta_file, $out_dir, $threads, $verbose, $help );
my $options = GetOptions(
    "fasta=s"   => \$fasta_file,
    "out_dir=s" => \$out_dir,
    "threads=i" => \$threads,
    "verbose"   => \$verbose,
    "help"      => \$help,
);

die $usage unless $options;
die $usage if $help;
die $usage unless defined $fasta_file;

my $geno = SNPtools::Genotype->new(
    fasta   => $fasta_file,
    out_dir => $out_dir,
    threads => $threads,
    verbose => $verbose,
);

my $map_dir = $out_dir . "/mapped/";
opendir my $map_dh, $map_dir;
my @ids = readdir $map_dh;

for (@ids) {
    my $bam_file = $map_dir . $_ . "/merged/" . $_ . "REPLACE_WITH_STANDARDIZED_NAME.bam";
    my $geno = SNPtools::Genotype->new(
        id      => $_,
        bam     => $bam_file,
        fasta   => $fasta_file,
        out_dir => $out_dir,
        threads => $threads,
        verbose => $verbose,
    );
    $geno->extract_mpileup;
    $geno->genotype;
}

exit;
