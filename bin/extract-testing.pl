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
    say "bam: ",      $self->bam      // "MISSING";
    say "fasta: ",    $self->fasta    // "MISSING";
    say "id: ",       $self->id       // "MISSING";
    say "out_dir: ",  $self->out_dir  // "MISSING";
    say "out_file: ", $self->out_file // "MISSING";
    say "par1: ",     $self->par1     // "MISSING";
    say "par2: ",     $self->par2     // "MISSING";
    say "seq_list: ", $self->seq_list // "MISSING";
    say "threads: ",  $self->threads  // "MISSING";
    say "verbose: ",  $self->verbose  // "MISSING";

    say "_chromosome: ",    $self->_chromosome    // "MISSING";
    say "_genotyped_dir: ", $self->_genotyped_dir // "MISSING";
    say "_plot_dir: ",      $self->_plot_dir      // "MISSING";
    say "_mpileup_dir: ",   $self->_mpileup_dir   // "MISSING";
    say "_snp_dir: ",       $self->_snp_dir       // "MISSING";

    say $_ for $self->get_seq_names;
    say "----------"
}

exit;
