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

my $ratio_threshold    = 2;
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
        my @elements = split ",", $_;
        my ( $scaff_name, $snp_pos, $reference, $snp_base ) = @elements[ 0 .. 2, 8 ];

        #REMOVE THE DOT DECIMAL FROM THE INSERTIONS: 1001.1 => 1001
        my $short_snp_pos = $snp_pos;
        $short_snp_pos =~ s/\..*//;

        # skip SNPs/indels supported by fewer than $coverage_threshold
        next if ( sum( @elements[ 3 .. 6 ] ) < $coverage_threshold       #ACGT
                 &&    $elements[7]          < $coverage_threshold );    #deletion

        # print SNPs/indels to output unless they don't pass the flanking coverage test
        unless (   $elements[10] == 0 || $elements[13] == 0                  #avoid illegal division by zero
                || (   $elements[11] / $elements[10] > $ratio_threshold
                    && $elements[14] / $elements[13] < $ratio_threshold )    #intron-exon junction
                || (   $elements[9]  / $elements[10] > $ratio_threshold
                    && $elements[12] / $elements[13] < $ratio_threshold )    #exon-intron junction
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