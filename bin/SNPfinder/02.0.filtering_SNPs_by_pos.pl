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

my $ratio_1   = 1.2;
my $ratio_2   = 0.8;
my $alt_ratio = 2;
my $cov_min   = 4;

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
    chomp $header;
    say $out_fh join ",", $header, "lt_ratio", "rt_ratio";

    while (<$snp_fh>){
        $counter++;

        my @elements    = split /,/, $_;
        my @acgt_counts = @elements[ 3 .. 6 ];
        my $del_count   = $elements[7];

        # skip SNPs/indels supported by fewer than $cov_min
        next if ( sum( @acgt_counts ) < $cov_min
                  &&   $del_count     < $cov_min );

        my ( $alt, $nogap_lt, $nogap_cov, $nogap_rt, $gap_lt, $gap_cov, $gap_rt ) =
          @elements[ 8 .. 14 ];

        # Avoid illegal division by zero (is this necessary w/ coverage check?)
        next if $nogap_cov == 0 || $gap_cov == 0;

        my $nogap_lt_ratio = $nogap_lt / $nogap_cov;
        my $nogap_rt_ratio = $nogap_rt / $nogap_cov;
        my $gap_lt_ratio   = $gap_lt / $gap_cov;
        my $gap_rt_ratio   = $gap_rt / $gap_cov;

        my $lt_ratio;
        my $rt_ratio;
        my $filter = 0;

        if ( $gap_lt_ratio > 0 && $gap_rt_ratio > 0 ) {
            $lt_ratio = sprintf( '%.3f', $nogap_lt_ratio / $gap_lt_ratio );
            $rt_ratio = sprintf( '%.3f', $nogap_rt_ratio / $gap_rt_ratio );

            # 1: exon-intron junction
            # 2: intron-exon junction
            $filter++
              if ( $lt_ratio >= $ratio_1 && $rt_ratio <= $ratio_2 )
              || ( $rt_ratio >= $ratio_1 && $lt_ratio <= $ratio_2 );
        }
        else { # Alternate filter to avoid illegal division by 0
            $lt_ratio = 'DIVby0';
            $rt_ratio = 'DIVby0';

            # 1: exon-intron junction
            # 2: intron-exon junction
            $filter++
              if ( $nogap_lt_ratio > $alt_ratio && $gap_lt_ratio < $alt_ratio )
              || ( $nogap_rt_ratio > $alt_ratio && $gap_rt_ratio < $alt_ratio );
        }

        # Don't apply filter for deletions
        $filter = 0 if $alt eq 'del';

        # Skip if filtered
        next if $filter;

        push @elements, $lt_ratio, $rt_ratio;
        chomp @elements;
        say $out_fh join ",", @elements;
        $counter_passed++;

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