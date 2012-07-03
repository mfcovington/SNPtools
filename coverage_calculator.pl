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
use Bio::DB::Sam;
use IO::File;
use coverage_commander;

my $usage = <<USAGE_END;

USAGE:
coverage_calculator.v2.pl
  --ref </PATH/TO/reference.fa>
  --bam </PATH/TO/file.bam>
  --id  <sample identifier for file.bam>
  --out </PATH/TO/DESTINATION/DIRECTORY/>
  --help

USAGE_END

my ( $bam_file, $id, $help );
my $out_dir = "./";
my $options = GetOptions(
    "bam=s" => \$bam_file,
    "id=s"  => \$id,
    "out=s" => \$out_dir,
    "help"  => \$help,
);

die $usage if $help;
die $usage unless defined $bam_file && defined $id;

my $sam = Bio::DB::Sam->new(
    -bam   => $bam_file,
    # -fasta => $ref_fa
);
my @chromosomes = $sam->seq_ids;

# open my $log_fh, ">", "$out_dir/$id.log";
# $log_fh->autoflush(1);

my $coverage = coverage_commander->new(
    bam     => "/Users/mfc/sandbox/genotyping/bwa_tophat_M82-Slyc.sorted.dupl_rm.bam",
    verbose => 1,
    # log     => "$out_dir/$id.log",
);

foreach my $chr (@chromosomes) {
    my $cov_out = "$out_dir/$id.$chr.coverage";
    $coverage->{chromosome} = $chr;
    $coverage->{out_file}   = $cov_out;
    $coverage->get_coverage();
}
# close $log_fh;
exit;

