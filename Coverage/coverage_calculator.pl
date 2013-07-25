#!/usr/bin/env perl
# coverage_calculator.pl
# Mike Covington
# created: 2011-12-05
#
# Description:
#
use strict;
use warnings;
use autodie;
use Getopt::Long;
use coverage_commander;

my $usage = <<USAGE_END;

USAGE:
coverage_calculator.pl
  --id          Sample identifier
  --bam         Sample alignment file (.bam)
  --seq_list    OPTIONAL: Comma-delimted list of sequence IDs to analyze
                (By default, this list is generated from the bam file header.)
  --out_dir     Output directory [current]
  --threads     Number of threads [1]
  --verbose
  --help

USAGE_END

my ( $id, $bam_file, $seq_list, $out_dir, $threads, $verbose, $help );
my $options = GetOptions(
    "id=s"       => \$id,
    "bam=s"      => \$bam_file,
    "seq_list=s" => \$seq_list,
    "out_dir=s"  => \$out_dir,
    "threads=i"  => \$threads,
    "verbose"    => \$verbose,
    "help"       => \$help,
);

die $usage if $help;
die $usage unless defined $bam_file && defined $id;

my $coverage = coverage_commander->new(
    id       => $id,
    bam      => $bam_file,
    seq_list => $seq_list,
    out_dir  => $out_dir || "./",
    threads  => $threads,
    verbose  => $verbose,
);
$coverage->get_coverage_all;

exit;

