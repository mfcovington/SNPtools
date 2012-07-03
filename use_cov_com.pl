#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2012-07-02
#
# Description: 
#
use Modern::Perl;
use Data::Printer;
use coverage_commander;


my $coverage = coverage_commander->new(
    bam        => "/Users/mfc/sandbox/genotyping/bwa_tophat_M82-Slyc.sorted.dupl_rm.bam",
    chromosome => "SL2.40ch00",
    pos_start  => 570800,
    pos_end    => 570900,
    include_intron_gaps => 0,
    verbose => 1,
  );
# p $coverage;  
say $coverage->samtools_cmd();
# system( $coverage->samtools_cmd() );

$coverage->get_coverage();