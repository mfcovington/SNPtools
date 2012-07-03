#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2012-07-02
#
# Description: 
#
use Modern::Perl;
# use Data::Printer;
use coverage_commander;


my $coverage = coverage_commander->new(
    bam        => "/Users/mfc/sandbox/genotyping/bwa_tophat_M82-Slyc.sorted.dupl_rm.bam",
    chromosome => "SL2.40ch01",
    # pos_start  => 570800,
    # pos_end    => 570900,
    out_file   => "chr01.test",
    verbose    => 1,
  );
# p $coverage;  
say $coverage->samtools_cmd_gaps;
say $coverage->samtools_cmd_nogaps;
# say $coverage->pos_start;
# $coverage->{pos_start} = 570850;
# say $coverage->pos_start;

$coverage->get_coverage();