#!/usr/bin/env perl
# snp_finder.pl
# Mike Covington
# created: 2011-12-05
#
# Description:
#
use strict;
use warnings;
use autodie;
use Getopt::Long;

use SNPtools::SNPfinder;
use SNPtools::Coverage;

###TODO:
#move looping through chromosomes into module - FINISHED, UNTESTED
#add subsequent steps
#make so out_dir contains several directories (one for each piece of analysis) (allow this to be over-riden)
#port threading to coverage_commander

my $usage = <<USAGE_END;

USAGE:
snp_finder.pl
  --id           Sample identifier
  --bam          Sample alignment file (.bam)
  --fasta        Reference file (.fasta/.fa)
  --seq_list     OPTIONAL: Comma-delimted list of sequence IDs to analyze
                 (By default, this list is generated from the bam file header.)
  --out_dir      Output directory [current]
  --cov_min      Minimum coverage [4]
  --snp_min      Minimum fraction of reads matching reference (for snps) [0.33]
  --indel_min    Minimum fraction of reads matching reference (for indels) [0.66]
  --threads      Number of threads [1]
  --verbose
  --help

USAGE_END

my ( $id, $bam_file, $fasta_file, $seq_list, $out_dir, $cov_min, $snp_min, $indel_min, $threads, $verbose,
    $help );
my $options = GetOptions(
    "id=s"        => \$id,
    "bam=s"       => \$bam_file,
    "fasta=s"     => \$fasta_file,
    "seq_list=s"  => \$seq_list,
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

my $snps = SNPtools::SNPfinder->new(
    id        => $id,
    bam       => $bam_file,
    fasta     => $fasta_file,
    seq_list  => $seq_list,
    out_dir   => $out_dir,
    cov_min   => $cov_min,
    snp_min   => $snp_min,
    indel_min => $indel_min,
    threads   => $threads,
    verbose   => $verbose,
);

my $coverage = SNPtools::Coverage->new(
    id       => $id,
    bam      => $bam_file,
    seq_list => $seq_list,
    out_dir  => $out_dir,
    threads  => $threads,
    verbose  => $verbose,
);

$snps->identify_snps;
$coverage->get_coverage_db;
$snps->flanking_cov;

exit;
