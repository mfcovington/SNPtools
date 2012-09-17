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
  --id         Sample identifier
  --bam        Sample alignment file (.bam)
  --out_dir    Output directory [current]
  --threads    Number of threads [1]
  --verbose
  --help

USAGE_END

my ( $bam_file, $id, $out_dir, $threads, $verbose, $help );
my $options = GetOptions(
    "bam=s"     => \$bam_file,
    "id=s"      => \$id,
    "out_dir=s" => \$out_dir,
    "threads=i" => \$threads,
    "verbose"   => \$verbose,
    "help"      => \$help,
);

die $usage if $help;
die $usage unless defined $bam_file && defined $id;

my $coverage = coverage_commander->new(
    bam     => $bam_file,
    id      => $id,
    out_dir => $out_dir || "./",
    threads => $threads,
    verbose => $verbose,
);
$coverage->get_coverage_all;

exit;

