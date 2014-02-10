#!/usr/bin/env perl
# Mike Covington
# created: 2014-02-09
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use List::Util 'sum';

use Data::Printer;
$|++;
my $bam_file = "sample-files/bam/R500.10kb.bam";
my $fasta_ref = "sample-files/fa/B.rapa_genome_sequence_0830.fa";
my $chromosome = "A01";

my $min_cov = 4;
my $min_snp_ratio = 0.66;
my $min_ins_ratio = 0.5;

open my $mpileup_fh, "-|", "samtools mpileup -A -r $chromosome -f $fasta_ref $bam_file";
# open my $mpileup_fh, "-|", "samtools mpileup -A -r $chromosome:19001-19300 -f $fasta_ref $bam_file";
# open my $mpileup_fh, "-|", "samtools mpileup -A -r $chromosome:19262-19262 -f $fasta_ref $bam_file";
# open my $mpileup_fh, "-|", "samtools mpileup -A -r $chromosome:10197-10197 -f $fasta_ref $bam_file";

say "seq_id,pos,ref,a,c,g,t,del,consensus";

while (<$mpileup_fh>) {
    my ( $seqid, $pos, $ref, $depth, $read_bases, $read_quals ) = split;

    $read_bases =~ tr/acgt/ACGT/;

    my ( $inserts, $read_bases_no_ins ) = get_inserts($read_bases);
    my $ins_counts = get_insert_counts($inserts);

    clean_deletions(\$read_bases_no_ins);

    my ( $total_counts, $counts ) = count_bases($ref, $read_bases_no_ins);

    my $consensus = get_consensus_base($counts);

    output_snp( $seqid, $pos, $ref, $counts, $consensus, $total_counts,
        $min_cov, $min_snp_ratio );

    output_insert( $seqid, $pos, $ref, $ins_counts, $counts, $min_ins_ratio );
}

close $mpileup_fh;

exit;

sub get_inserts {    # Capture sequences of variable length inserts
                     # and remove them from $read_bases
    my $read_bases = shift;

    my %inserts;
    for my $ins_len ( $read_bases =~ m/\+(\d+)/g ) {
        $inserts{$1}++ if $read_bases =~ s/\+(?:$ins_len)([ACGT]{$ins_len})//;
    }
    return \%inserts, $read_bases;
}

sub get_insert_counts {    # Get insert counts by position and nucleotide
    my $inserts = shift;

    my %ins_counts;
    for my $insert ( keys $inserts ) {
        my $ins_pos = sprintf "%02d", 1;
        for my $nt ( split //, $insert ) {
            $ins_counts{$ins_pos}{$nt} += $$inserts{$insert};
            $ins_pos++;
        }
    }
    return \%ins_counts;
}

sub clean_deletions {    # Keep '*' deletion indicator, but remove '-1A'-type
    my $pileup_ref = shift;

    for my $del_len ( $$pileup_ref =~ m/\-(\d+)/g ) {
        $$pileup_ref =~ s/\-(?:$del_len)[ACGT]{$del_len}//;
    }
}

sub count_bases {
    my ( $ref, $read_bases_no_ins ) = @_;

    my %counts;

    for my $base (qw(A C G T del)) {
        $counts{$base} = 0;
    }

    $counts{A}++    for $read_bases_no_ins =~ m/A/ig;
    $counts{C}++    for $read_bases_no_ins =~ m/C/ig;
    $counts{G}++    for $read_bases_no_ins =~ m/G/ig;
    $counts{T}++    for $read_bases_no_ins =~ m/T/ig;
    $counts{$ref}++ for $read_bases_no_ins =~ m/[.,]/g;
    $counts{del}++  for $read_bases_no_ins =~ m/\*/ig;

    my $total_counts = sum values %counts;

    return $total_counts, \%counts;
}

sub get_consensus_base {
    my $counts = shift;

    my ($consensus)
        = scalar keys $counts == 1
        ? keys $counts
        : sort { $$counts{$b} <=> $$counts{$a} } keys $counts;

    return $consensus;
}

sub output_snp {
    my ( $seqid, $pos, $ref, $counts, $consensus, $total_counts, $min_cov,
        $min_snp_ratio )
        = @_;

    say join ",", $seqid, $pos, $ref, $$counts{A}, $$counts{C}, $$counts{G},
        $$counts{T}, $$counts{del}, $consensus
        if ( $ref ne $consensus
        && $total_counts >= $min_cov
        && $$counts{$consensus} >= $min_snp_ratio * $total_counts );
}

sub output_insert {
    my ( $seqid, $pos, $ref, $ins_counts, $counts, $min_ins_ratio ) = @_;

    for my $ins_pos ( sort { $a <=> $b } keys $ins_counts ) {
        my ($ins_base)
            = sort { $$ins_counts{$ins_pos}{$b} <=> $$ins_counts{$ins_pos}{$a} }
            keys $$ins_counts{$ins_pos};
        next
            unless $$ins_counts{$ins_pos}{$ins_base}
            >= $$counts{$ref} * $min_ins_ratio;
        say join ",", $seqid, "$pos.$ins_pos", "INS",
            $$ins_counts{$ins_pos}{A} // 0, $$ins_counts{$ins_pos}{C} // 0,
            $$ins_counts{$ins_pos}{G} // 0, $$ins_counts{$ins_pos}{T} // 0,
            0, $ins_base;
    }
}
