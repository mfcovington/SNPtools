#!/usr/bin/env perl
# coverage_calculator.v2.pl
# Mike Covington
# created: 2011-12-05
#
# Description:
#
use strict;
use warnings;
use autodie;
use Getopt::Long;
use Parallel::ForkManager;

use snp_commander;

my $usage = <<USAGE_END;

USAGE:
snp_finder.pl
  --id           Sample identifier
  --bam          Sample alignment file (.bam)
  --fasta        Reference file (.fasta/.fa)
  --out_dir      Output directory
  --cov_min      Minimum coverage [4]
  --snp_min      Minimum fraction of reads matching reference (for snps) [0.33]
  --indel_min    Minimum fraction of reads matching reference (for indels) [0.66]
  --threads      Number of threads [1]
  --verbose
  --help

USAGE_END

my $threads = 1;
my $out_dir = "./";
my ( $id, $bam_file, $fasta_file, $cov_min, $snp_min, $indel_min, $verbose,
    $help );
my $options = GetOptions(
    "id=s"        => \$id,
    "bam=s"       => \$bam_file,
    "fasta=s"     => \$fasta_file,
    "out_dir=s"   => \$out_dir,
    "cov_min=i"   => \$cov_min,
    "snp_min=f"   => \$snp_min,
    "indel_min=f" => \$indel_min,
    "threads=i"   => \$threads,
    "verbose"     => \$verbose,
    "help"        => \$help,
);

die $usage unless $options;
die $usage if $help;
die $usage
  unless defined $id
  && defined $bam_file
  && defined $fasta_file;

my $snps = snp_commander->new(
    bam     => $bam_file,
    fasta   => $fasta_file,
    verbose => $verbose,
);

$snps->cov_min($cov_min)     if defined $cov_min;
$snps->snp_min($snp_min)     if defined $snp_min;
$snps->indel_min($indel_min) if defined $indel_min;

my @chromosomes = $snps->get_seq_names;

# put this inside module??
my $pm = new Parallel::ForkManager($threads);
foreach my $chr (@chromosomes) {
    $pm->start and next;

    my $cov_out = "$out_dir/$id.$chr.snps";
    $snps->chromosome($chr);
    $snps->out_file($cov_out);
    $snps->identify_snps;

    $pm->finish;
}
$pm->wait_all_children;
exit;
