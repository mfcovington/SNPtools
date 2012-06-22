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
# use Bio::DB::Sam;
# use Data::Dumper;
use List::Util qw( sum );
use File::Basename;

my $ratio_threshold = 2;
my $coverage_threshold = 4;
my ( $counter, $counter_passed );

my @files = @ARGV;
foreach my $file (@files) {
	my $SNP_file = $file;
    my ( $filename, $directories, $suffix ) = fileparse( $SNP_file, ".csv" );
	my $out_file = "$directories$filename.FILTERED.csv";
    my %bases_n;
	
    open SNPs, "<", $SNP_file;
    open OUT,  ">", $out_file;

    my $header = <SNPs>;
    print OUT $header;
	
	while (<SNPs>){
		$counter++;
        my @elements = split ",", $_;
		my ( $scaff_name, $snp_pos, $reference, $snp_base ) = @elements[ 0 .. 2, 8 ];

		#REMOVE THE DOT DECIMAL FROM THE INSERTIONS: 1001.1 => 1001
		my $short_snp_pos = $snp_pos;
		$short_snp_pos =~ s/\..*//;
	
        # skip SNPs/indels supported by fewer than $coverage_threshold
        if ( sum( @elements[ 3 .. 6 ] ) < $coverage_threshold    #ACGT
            &&    $elements[7]          < $coverage_threshold    #deletion
        ) {
            next;
        }

		# print SNPs/indels to output unless they don't pass the flanking coverage test
        unless (   $elements[26] == 0                                        #avoid illegal division by zero
                || (   $elements[27] / $elements[26] > $ratio_threshold
                    && $elements[30] / $elements[29] < $ratio_threshold )    #intron-exon junction
                || (   $elements[25] / $elements[26] > $ratio_threshold
                    && $elements[28] / $elements[29] < $ratio_threshold )    #exon-intron junction
               )
        {
            print OUT join ",", @elements;
            $counter_passed++;
        }

	}
close OUT;

my $counter_diff = $counter - $counter_passed;
print "Filtered from $counter down to $counter_passed ("
  . $counter_diff
  . " SNPs/indels removed ("
  . ( $counter - $counter_passed ) * 100 / $counter
  . " %)) \n---Done!;\n";

}
exit;

