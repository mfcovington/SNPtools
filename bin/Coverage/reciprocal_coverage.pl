#!/usr/bin/env perl
# reciprocal_coverage.pl
# Mike Covington
# created: 2013-11-08
#
# Description:
#
use strict;
use warnings;
use Getopt::Long;
use SNPtools::Coverage;

my $usage = <<USAGE_END;

USAGE:
$0
  --bam          Sample alignment file (.bam) [Unused, but required for now.]
  --par1         Parent 1 ID
  --par2         Parent 2 ID
  --par1_bam     Parent 1 alignment file (.bam)
  --par2_bam     Parent 2 alignment file (.bam)
  --seq_list     OPTIONAL: Comma-delimted list of sequence IDs to analyze
                 (By default, this list is generated from the bam file header.)
  --out_dir      Output directory [current]
  --threads      Number of threads [1]
  --verbose
  --help

USAGE_END

my (
    $bam_file, $par1,    $par2,    $par1_bam, $par2_bam,
    $seq_list, $out_dir, $threads, $verbose,  $help
);
my $options = GetOptions(
    "bam=s"      => \$bam_file, # Need temporarily
    "par1=s"     => \$par1,
    "par2=s"     => \$par2,
    "par1_bam=s" => \$par1_bam,
    "par2_bam=s" => \$par2_bam,
    "seq_list=s" => \$seq_list,
    "out_dir=s"  => \$out_dir,
    "threads=i"  => \$threads,
    "verbose"    => \$verbose,
    "help"       => \$help,
);

die $usage unless $options;
die $usage if $help;
die $usage
  unless
     defined $par1
  && defined $par2
  && defined $bam_file
  && defined $par1_bam
  && defined $par2_bam;

my $coverage = SNPtools::Coverage->new(
    par1     => $par1,
    par2     => $par2,
    par1_bam => $par1_bam,
    par2_bam => $par2_bam,
    bam      => $bam_file,
    seq_list => $seq_list,
    out_dir  => $out_dir,
    threads  => $threads,
    verbose  => $verbose,
);

$coverage->reciprocal_coverage;
exit;
