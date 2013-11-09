#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2012-07-02
#
# Description:
#
use Modern::Perl;
use Data::Printer;
use SNPtools::Coverage;


my $coverage = SNPtools::Coverage->new(
    bam        => "/Users/mfc/sandbox/genotyping/bwa_tophat_M82-Slyc.sorted.dupl_rm.bam",
    chromosome => "SL2.40ch01",
    # pos_start  => 570800,
    # pos_end    => 570900,
    out_file   => "chr01.test",
    verbose    => 1,
    # nogap => 0,
  );
# p $coverage;
say $coverage->samtools_cmd_gaps;
say $coverage->samtools_cmd_nogaps;
# say $coverage->pos_start;
# $coverage->{pos_start} = 570850;
# say $coverage->pos_start;
say $coverage->samtools_cmd_header;

my @lengths = $coverage->get_seq_lengths;
my @names = $coverage->get_seq_names;
p @lengths;
p @names;
# $coverage->get_coverage();