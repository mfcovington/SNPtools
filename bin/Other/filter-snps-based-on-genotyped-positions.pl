#!/usr/bin/env perl
# Mike Covington
# created: 2014-05-02
#
# Description: Filters a set of polymorphism files based on positions observed
#              in a set of genotyped files
#
use strict;
use warnings;
use autodie;
use feature 'say';
use File::Path 'make_path';
use Getopt::Long;

my ( $snp_in_dir, $snp_out_dir );

my $options = GetOptions(
    "snp_in_dir=s"  => \$snp_in_dir,
    "snp_out_dir=s" => \$snp_out_dir,
);

my @filtered_genotype_file_list = @ARGV;

die "Must define SNP input directory.\n"  unless defined $snp_in_dir;
die "Must define SNP output directory.\n" unless defined $snp_out_dir;
die "Must provide list of filtered genotyped files.\n"
    unless scalar @ARGV > 0;

make_path $snp_out_dir;

my %snps_filtered;

for my $genotyped_file (@filtered_genotype_file_list) {
    open my $genotyped_fh, "<", $genotyped_file;
    while (<$genotyped_fh>) {
        my ( $chr, $pos ) = split;
        $snps_filtered{$chr}{$pos} = 1;
    }
    close $genotyped_fh;
}

my @snp_files = `ls $snp_in_dir/*nr`;
chomp @snp_files;
$_ =~ s|$snp_in_dir/|| for @snp_files;

for my $file (@snp_files) {
    my %snps_in;
    my $chr;

    open my $snp_in_fh, "<", "$snp_in_dir/$file";
    my $header = <$snp_in_fh>;
    while (<$snp_in_fh>) {
        ( $chr, my $pos ) = split;
        $snps_in{$pos} = $_;
    }
    close $snp_in_fh;

    open my $snp_out_fh, ">", "$snp_out_dir/$file";
    print $snp_out_fh $header;
    for my $pos ( sort { $a <=> $b } keys %snps_in ) {
        print $snp_out_fh $snps_in{$pos} if exists $snps_filtered{$chr}{$pos};
    }
    close $snp_out_fh;
}
