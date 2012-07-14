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
use feature 'say';
use Getopt::Long;
use IO::File;
use coverage_commander;

my $usage = <<USAGE_END;

USAGE:
coverage_calculator.v2.pl
  --bam     </PATH/TO/file.bam>
  --id      <sample identifier for file.bam>
  --out_dir </PATH/TO/DESTINATION/DIRECTORY/>
  --verbose
  --help

USAGE_END

my ( $bam_file, $id, $verbose, $help );
my $out_dir = "./";
my $options = GetOptions(
    "bam=s"     => \$bam_file,
    "id=s"      => \$id,
    "out_dir=s" => \$out_dir,
    "verbose"   => \$verbose,
    "help"      => \$help,
);

die $usage if $help;
die $usage unless defined $bam_file && defined $id;

my $coverage = coverage_commander->new(
    bam     => $bam_file,
    verbose => $verbose,
);
my @chromosomes = $coverage->get_seq_names;

foreach my $chr (@chromosomes) {
    my $cov_out = "$out_dir/$id.$chr.coverage";
    $coverage->chromosome($chr);
    $coverage->out_file($cov_out);
    $coverage->get_coverage;
}
exit;

