#!/usr/bin/env perl
# flanking_coverage_calculator.pl
# Mike Covington
# created: 2011-12-11
#
# Description: 
#
use strict;
use warnings;
use autodie;
use feature 'say';
use File::Basename;

# new version of coverage data is tab-delimited and only has data for regions with coverage

my $usage = <<USAGE_END;

USAGE:
flanking_coverage_calculator.pl
  --snp_file   </PATH/TO/snp_file.csv>
  --cov_prefix </PATH/TO/cov_file_prefix_only -- w/o .cov_(no)gap>
  --help

USAGE_END

my ( $snp_in, $cov_prefix, $help );
my $out_dir = "./";
my $options = GetOptions(
    "snp_file=s"   => \$snp_in,
    "cov_prefix=s" => \$cov_prefix,
    "help"         => \$help,
);

die $usage if $help;
die $usage unless defined $snp_in && defined $cov_prefix;

my ( $cov_nogaps, $cov_gaps ) = ( $cov_prefix ) x 2 ;
$cov_nogaps .= ".cov_nogaps";
$cov_gaps   .= ".cov_gaps";
my ( $filename, $directories, $suffix ) = fileparse( $snp_in, ".csv" );
my $snp_out = $directories . $filename . ".nogap.gap.csv";

#open files
open $snp_in_fh,     "<", $snp_in;
open $snp_out_fh,    ">", $snp_out;
open $cov_nogaps_fh, "<", $cov_nogaps;
open $cov_gaps_fh,   "<", $cov_gaps;


#read in first 17 lines
my (@nogaps_chr, @nogaps_pos, @nogaps_cov, @gaps_chr, @gaps_pos, @gaps_cov, $nogaps_line, $gaps_line);
for my $count ( 0 .. 16 ) {
    $nogaps_line  = <$cov_nogaps_fh>;
    $gaps_line = <$cov_gaps_fh>;
    chomp( $nogaps_line, $gaps_line );
    ( $nogaps_chr[$count], $nogaps_pos[$count], $nogaps_cov[$count] )    = split( /\t/, $nogaps_line );
    ( $gaps_chr[$count], $gaps_pos[$count], $gaps_cov[$count] ) = split( /\t/, $gaps_line );
}


my $header = <$snp_in_fh>;
chomp $header;
say $snp_out_fh join( ",", $header, "nogap_pos-8", "nogap_pos", "nogap_pos+8", "gap_pos-8", "gap_pos", "gap_pos+8" );
while ( my $snp_line = <$snp_in_fh> ) {
    chomp $snp_line;
    my ( $snp_chr, $snp_pos_unsplit, $snp_remainder ) = split( /,/, $snp_line, 3 );
    my ( $snp_pos, $snp_pos_index ) = split( /\./, $snp_pos_unsplit );
    while ( $snp_pos > $nogaps_pos[8] ) {
        shift @nogaps_chr;
        shift @nogaps_pos;
        shift @nogaps_cov;
        shift @gaps_chr;
        shift @gaps_pos;
        shift @gaps_cov;
        $nogaps_line = <$cov_nogaps_fh>;
        $gaps_line   = <$cov_gaps_fh>;
        chomp( $nogaps_line, $gaps_line );
        my @nogaps_line_elements = split( /\t/, $nogaps_line );
        my @gaps_line_elements   = split( /\t/, $gaps_line );
        push @nogaps_chr, $nogaps_line_elements[0];
        push @nogaps_pos, $nogaps_line_elements[1];
        push @nogaps_cov, $nogaps_line_elements[2];
        push @gaps_chr,   $gaps_line_elements[0];
        push @gaps_pos,   $gaps_line_elements[1];
        push @gaps_cov,   $gaps_line_elements[2];
    }

    if ( $snp_chr eq $nogaps_chr[8] && $snp_pos == $nogaps_pos[8] ) {
        say $snp_out_fh join ( ",", $snp_chr, $snp_pos_unsplit, $snp_remainder, $nogaps_cov[0], $nogaps_cov[8], $nogaps_cov[16], $gaps_cov[0], $gaps_cov[8], $gaps_cov[16] );
    }
    else {
        die "Something strange going on here...";
    }
}

close $snp_in_fh;
close $cov_nogaps_fh;
close $cov_gaps_fh;
close $snp_out_fh;


