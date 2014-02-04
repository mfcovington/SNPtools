#!/usr/bin/env perl
# snp_finder.pl
# Mike Covington
# created: 2014-02-04
#
# Description: Temporary script for use while extracting
#              common code to SNPtools.pm from other modules
#
use strict;
use warnings;
use autodie;

use SNPtools::SNPfinder;
use SNPtools::Coverage;
use SNPtools::Genotype;
use SNPtools::Plot;
use feature 'say';


my $id       = "test_id";
my $base_dir = "/Users/mfc/git.repos";
my $bam_file = "$base_dir/snp_identification/1.1.bam";
my $fa_file  = "$base_dir/sample-files/fa/S_lycopersicum_chromosomes.2.40.fa";

my $snps = SNPtools::SNPfinder->new(
    id    => $id,
    bam   => $bam_file,
    fasta => $fa_file,
);

my $coverage = SNPtools::Coverage->new(
    id  => $id,
    bam => $bam_file,
);

my $genotype = SNPtools::Genotype->new(
    id    => $id,
    bam   => $bam_file,
    fasta => $fa_file,
);

my $plot = SNPtools::Plot->new(
    id  => $id,
    bam => $bam_file,
);

report( "SNP", $snps);
report( "COV", $coverage);
report( "GENO", $genotype);
report( "PLOT", $plot);

sub report {
    my ( $class_id, $self ) = @_;

    say $class_id;
    say "ID: ", $self->id;
    say "OUT_DIR: ", $self->out_dir;
    say $_ for $self->get_seq_names;
    say "----------"
}

exit;
