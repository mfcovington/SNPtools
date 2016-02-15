#!/usr/bin/env perl
# Mike Covington
# created: 2014-05-02
#
# Description: Filters a set of genotyped files based on positions observed
#              in a set of polymorphism files.
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use Getopt::Long;


my ( $id, $no_nr, $out_dir );

my $options = GetOptions(
    "id=s"      => \$id,
    "no_nr"     => \$no_nr,
    "out_dir=s" => \$out_dir,
);

die
    "Must define an output directory that contains 'snp_master/' and 'genotyped/'.\n"
    unless defined $out_dir;
die "Must define the ID of the sample to filter.\n" unless defined $id;

my @snp_file_list;
if ($no_nr) {
    @snp_file_list = `ls $out_dir/snp_master/polyDB.* | grep -v nr\$`;
}
else {
    @snp_file_list = `ls $out_dir/snp_master/polyDB.*.nr`;
}
chomp @snp_file_list;

my %snps_filtered;

for my $snp_file (@snp_file_list) {
    open my $snp_fh, "<", $snp_file;
    my $header = <$snp_fh>;
    while (<$snp_fh>) {
        my ( $chr, $pos ) = split;
        $snps_filtered{$chr}{$pos} = 1;
    }
    close $snp_fh;
}

for my $chr ( sort keys %snps_filtered ) {
    my $geno_in_file = "$out_dir/genotyped/$id.$chr.genotyped";
    $geno_in_file .= '.nr' unless $no_nr;

    my $geno_out_file = "$out_dir/genotyped/$id.filtered.$chr.genotyped";
    $geno_out_file .= '.nr' unless $no_nr;

    open my $geno_in_fh, "<", $geno_in_file;
    open my $geno_out_fh, ">", $geno_out_file;

    my $header = <$geno_in_fh>;
    print $geno_out_fh $header;
    while (<$geno_in_fh>) {
        my ( $chr, $pos ) = split;
        print $geno_out_fh $_ if exists $snps_filtered{$chr}{$pos};
    }

    close $geno_in_fh;
    close $geno_out_fh;
}
