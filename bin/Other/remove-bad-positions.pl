#!/usr/bin/env perl
# Mike Covington
# created: 2014-05-28
#
# Description: Remove positions from SNP & genotyped files
#
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use File::Path 'make_path';

# TODO: Make output directories customizable at CLI
# TODO: Use Getopt::Long

die "Specify 'bad positions file' and 'source directory'"
    unless scalar @ARGV == 2;

my $bad_positions_file = $ARGV[0];
my $source_dir         = $ARGV[1];
my $output_dir         = "$source_dir/clean";

my $snp_dir  = "$source_dir/snp_master";
my $geno_dir = "$source_dir/genotyped";

my $snp_out_dir  = "$output_dir/snp_master";
my $geno_out_dir = "$output_dir/genotyped";
make_path $snp_out_dir, $geno_out_dir;

open my $bad_positions_fh, "<", $bad_positions_file;
my %bad_positions;
while (<$bad_positions_fh>) {
    chomp;
    my ( $chr, $pos ) = split;
    $bad_positions{$chr}{$pos} = 1;
}
close $bad_positions_fh;

opendir(my $snp_dh, $snp_dir);
my @snp_file_list
    = grep { /^polyDB\..+\.nr/ && -f "$snp_dir/$_" } readdir($snp_dh);
closedir $snp_dh;

opendir(my $geno_dh, $geno_dir);
my @geno_file_list
    = grep { /.+\.genotyped\.nr/ && -f "$geno_dir/$_" } readdir($geno_dh);
closedir $geno_dh;

for my $snp_file (@snp_file_list) {
    open my $snp_fh,     "<", "$snp_dir/$snp_file";
    open my $snp_out_fh, ">", "$snp_out_dir/$snp_file";
    my $header = <$snp_fh>;
    print $snp_out_fh $header;
    while (<$snp_fh>) {
        my ( $chr, $pos ) = split;
        next if exists $bad_positions{$chr}{$pos};
        print $snp_out_fh $_;
    }
    close $snp_fh;
    close $snp_out_fh;

}

for my $geno_file (@geno_file_list) {
    open my $geno_fh,     "<", "$geno_dir/$geno_file";
    open my $geno_out_fh, ">", "$geno_out_dir/$geno_file";
    while (<$geno_fh>) {
        my ( $chr, $pos ) = split;
        next if exists $bad_positions{$chr}{$pos};
        print $geno_out_fh $_;
    }
    close $geno_fh;
    close $geno_out_fh;
}

exit;
