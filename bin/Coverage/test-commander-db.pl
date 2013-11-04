#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2013-06-13
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use SNPtools::Coverage;
my $seq_list = 'A01,A02,A03,A04,A05,A06,A07,A08,A09,A10';
my $out_dir = '/Users/mfc/git.repos/SNPtools/sample-files/output';
my $coverage = SNPtools::Coverage->new(
    id       => 'R500',#'test',
    bam      => "sample-files/bam/R500.10kb.bam",#"sample-files/bam/R500.trunc.bam",#"test/test.bam",
    seq_list => $seq_list,
    out_dir  => $out_dir || "./",
    threads  => 2,#4,
    verbose  => 1,
    # _chromosome => 'SL2.40ch01',
);
# $coverage->get_coverage_all;

$coverage->get_coverage_db;

# $coverage->add_positions;
use Data::Printer;
my $cov_pos_ref = $coverage->cov_pos;
# print scalar keys $$cov_pos_ref{'SL2.40ch01'}, "\n";

# $coverage->populate_CoverageDB_by_chr;

exit;
