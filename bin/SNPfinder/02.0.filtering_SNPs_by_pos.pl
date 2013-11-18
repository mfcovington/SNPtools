#!/usr/bin/env perl
# 02.0.filtering_SNPs_by_pos.pl
# Mike Covington
# created: 2011-12-12
#
# Description:
#
use strict;
use warnings;
use autodie;
use List::Util qw( sum );
use File::Basename;

my $ratio_threshold    = 1.2;
my $coverage_threshold = 4;

my @files = @ARGV;
my ( $counter, $counter_passed );
foreach my $file (@files) {
    my $snp_file = $file;
    my ( $filename, $directories, $suffix ) = fileparse( $snp_file, ".csv" );
    my $out_file = "$directories$filename.FILTERED.csv";
    my %bases_n;

    open my $snp_fh, "<", $snp_file;
    open my $out_fh, ">", $out_file;

    my $header = <$snp_fh>;
    print $out_fh $header;

    while (<$snp_fh>){
        $counter++;

        my @elements    = split /,/, $_;
        my @acgt_counts = @elements[ 3 .. 6 ];
        my $del_count   = $elements[7];

        # skip SNPs/indels supported by fewer than $coverage_threshold
        next if ( sum( @acgt_counts ) < $coverage_threshold
                  &&   $del_count     < $coverage_threshold );

        my ( $nogap_lt, $nogap_cov, $nogap_rt, $gap_lt, $gap_cov, $gap_rt ) =
          @elements[ 9 .. 14 ];

        my $nogap_rt_ratio = $nogap_rt / $nogap_cov;
        my $nogap_lt_ratio = $nogap_lt / $nogap_cov;
        my $gap_rt_ratio   = $gap_rt / $gap_cov;
        my $gap_lt_ratio   = $gap_lt / $gap_cov;

        # next if $gap_rt_ratio == 0;
        # next if $gap_lt_ratio == 0;
        if ( $gap_rt_ratio == 0 || $gap_lt_ratio == 0 ) {
            chomp @elements;
            say $out_fh join ",", @elements, "BADCOV";
            next;
        }

        # print SNPs/indels to output if they pass the flanking coverage test
        unless (
            $nogap_cov == 0 || $gap_cov == 0    #avoid illegal division by zero
            # || (   $nogap_rt_ratio > $ratio_threshold
            #     && $gap_rt_ratio < $ratio_threshold ) #intron-exon junction
            # || (   $nogap_lt_ratio > $ratio_threshold
            #     && $gap_lt_ratio < $ratio_threshold ) #exon-intron junction
            || $nogap_rt_ratio / $gap_rt_ratio > $ratio_threshold #intron-exon junction
            || $nogap_lt_ratio / $gap_lt_ratio > $ratio_threshold #exon-intron junction
          )
        {
            print $out_fh join ",", @elements;
            $counter_passed++;
        }

    }

    close $snp_fh;
    close $out_fh;

    my $counter_diff = $counter - $counter_passed;
    print "Filtered from $counter down to $counter_passed ("
      . $counter_diff
      . " SNPs/indels removed ("
      . ( $counter - $counter_passed ) * 100 / $counter
      . " %)) \n---Done!\n";
}

exit;