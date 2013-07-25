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
use coverage_commander;
# my $seq_list = 'SL2.40ch01';
my $coverage = coverage_commander->new(
    id       => 'test',
    bam      => "test/test.bam",
    # seq_list => $seq_list,
    # out_dir  => $out_dir || "./",
    threads  => 4,
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
