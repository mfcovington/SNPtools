#!/usr/bin/env perl
use strict;
use warnings;
use feature 'say';
use CoverageDB::Main;
use autodie;
my $schema = CoverageDB::Main->connect('dbi:SQLite:db/coverage.db');

my $flank_dist = 8;
my @chromosomes =
  qw( SL2.40ch01 SL2.40ch02 SL2.40ch03 SL2.40ch04 SL2.40ch05 SL2.40ch06 SL2.40ch07 SL2.40ch08 SL2.40ch09 SL2.40ch10 );

for my $chr (@chromosomes) {
    my %cov_pos;
    add_positions( $chr, \%cov_pos );
    populate_CoverageDB_by_chr( $chr, \%cov_pos );
}

sub add_positions {
    my ( $chr, $cov_pos_ref ) = @_;
    open my $snps_fh, "<", "../genotyping/snp_master/polyDB.$chr.nr";
    <$snps_fh>;
    while (<$snps_fh>) {
        my $snp_pos = [ split /\t/ ]->[1];
        $$cov_pos_ref{$chr}{$snp_pos}                 = 1;
        $$cov_pos_ref{$chr}{ $snp_pos - $flank_dist } = 1;
        $$cov_pos_ref{$chr}{ $snp_pos + $flank_dist } = 1;
    }
    close $snps_fh;
}

sub populate_CoverageDB_by_chr {
    my ( $chromosome, $cov_pos_ref ) = @_;
    my $bam_file = "test/test.bam";

    my $sam_gap_cmd = "samtools mpileup -r $chromosome $bam_file | cut -f1-2,4";
    my $sam_nogap_cmd = "samtools depth -r $chromosome $bam_file";

    open my $gap_fh,   "-|", $sam_gap_cmd;
    open my $nogap_fh, "-|", $sam_nogap_cmd;

    my $sample_id = "test";
    my $count     = 1;
    my @cov_data;

    while ( my $gap_line = <$gap_fh> ) {
        my $nogap_line = <$nogap_fh>;
        chomp( $gap_line, $nogap_line );
        my ( $chr, $pos, $gap_cov ) = split /\t/, $gap_line;
        my $nogap_cov = [ split /\t/, $nogap_line ]->[2];

        if ( exists $$cov_pos_ref{$chr}{$pos} ) {
            $count++;
            push @cov_data, [ $sample_id, $chr, $pos, $gap_cov, $nogap_cov ];
        }

        populate_and_reset( \$count, \@cov_data ) if $count % 100000 == 0;
    }
    close $gap_fh;
    close $nogap_fh;

    populate_and_reset( \$count, \@cov_data ) if scalar @cov_data;
}

sub populate_and_reset {
    my ( $count_ref, $cov_data_ref ) = @_;
    say $$count_ref++;
    $schema->populate(
        'Coverage',
        [
            [qw/sample_id chromosome position gap_cov nogap_cov/],
            @$cov_data_ref
        ]
    );
    @$cov_data_ref = ();
}
