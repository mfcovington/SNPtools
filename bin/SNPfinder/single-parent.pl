#!/usr/bin/env perl
# Mike Covington
# created: 2014-04-16
#
# Description: Creates polyDB master SNP files from filter SNP files
#              for SNPs between a single parent and the reference
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use File::Basename;
use File::Path 'make_path';

# TODO: Incorporate this functionality into module

my ( $par1, @snps_file_list ) = @ARGV;

for my $snps_file (@snps_file_list) {

    my ($expected_chr)
        = $snps_file =~ /$par1\.(.+)\.snps\.nogap\.gap\.FILTERED\.csv/;

    die "Unable to extract $par1 chromosome ID for $snps_file."
        unless defined $expected_chr;

    my $snps_dir = ( fileparse($snps_file) )[1];
    my $out_dir  = "$snps_dir/../snp_master";
    make_path($out_dir);

    my $polydb_file = "$out_dir/polyDB.$expected_chr";
    die "Output file ($polydb_file) already exists" if -e $polydb_file;
    open my $polydb_fh, ">", $polydb_file;
    say $polydb_fh join "\t", 'chr', 'pos', 'ref_base', 'snp_base',
        'genotype', 'insert_position', 'SNP_CLASS';

    open my $snps_fh, "<", $snps_file;
    <$snps_fh>;
    for (<$snps_fh>) {
        my ( $chr, $pos_raw, $ref, $alt ) = ( split /,/, $_ )[ 0 .. 2, 8 ];
        die "Chromosome ID ($chr) does not match expected ID ($expected_chr)"
            if $chr ne $expected_chr;
        my ( $pos, $ins_pos ) = split /\./, $pos_raw;
        $ins_pos //= 'NA';
        say $polydb_fh join "\t", $chr, $pos, $ref, $alt, $par1, $ins_pos,
            'SNP';
    }
    close $snps_fh;
    close $polydb_fh;
}

exit;
