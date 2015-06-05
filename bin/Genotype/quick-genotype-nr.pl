#!/usr/bin/env perl
# Mike Covington
# created: 2015-06-04
#
# Description: Quickly filter genotyped files for SNPs that passed noise-reduction
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use Getopt::Long;

my $snp_dir  = 'snp_master';
my $geno_dir = 'genotyped';
my ( $force, $help );

my $options = GetOptions (
    "snp_dir=s"  => \$snp_dir,
    "geno_dir=s" => \$geno_dir,
    "force"      => \$force,
    "help"       => \$help,
);

my @snp_file_list  = get_files( $snp_dir,  qr/^polyDB\..+\.nr$/ );
my @geno_file_list = get_files( $geno_dir, qr/^.+\.genotyped$/ );

my $positions = {};
for (@snp_file_list) {
    get_snp_positions( "$snp_dir/$_", $positions );
}

for (@geno_file_list) {
    my $geno_file = "$geno_dir/$_";

    if ( !-e "$geno_file.nr" || $force ) {
        filter_genotyped_file( $geno_file, $positions );
    }
}

exit;

sub get_files {
    my ( $dir, $file_regex ) = @_;

    opendir(my $dh, $dir);
    my @file_list
        = grep { /$file_regex/ && -f "$dir/$_" } readdir($dh);
    closedir $dh;

    return @file_list;
}

sub get_snp_positions {
    my ( $snp_file, $positions ) = @_;

    open my $snp_fh, "<", $snp_file;
    <$snp_fh>;
    while (<$snp_fh>) {
        my ( $seqid, $pos ) = split;
        $$positions{$seqid}{$pos} = 1;
    }
    close $snp_fh;
}

sub filter_genotyped_file {
    my ( $geno_file, $positions ) = @_;

    open my $geno_fh, "<", $geno_file;
    open my $geno_nr_fh, ">", "$geno_file.nr";
    while (<$geno_fh>) {
        my ( $seqid, $pos ) = split;
        print $geno_nr_fh $_ if exists $$positions{$seqid}{$pos};
    }
    close $geno_fh;
    close $geno_nr_fh;
}
