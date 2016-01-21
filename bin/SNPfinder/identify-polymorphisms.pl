#!/usr/bin/env perl
# Mike Covington
# created: 2014-02-09
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use List::Util 'sum';
use Getopt::Long;

my $bam_file   = "../../sample-files/bam/R500.10kb.bam";
my $fasta_ref  = "../../sample-files/fa/B.rapa_genome_sequence_0830.fa";
my $chromosome = "A01";

my $min_cov       = 4;
my $min_snp_ratio = 0.66;
my $min_ins_ratio = 0.33;

my $outputfile = "$chromosome.snps";    # Temporary file name

my $options = GetOptions(
    "chromosome=s"    => \$chromosome,
    "outputfile=s"    => \$outputfile,
    "min_cov=i"       => \$min_cov,
    "min_snp_ratio=f" => \$min_snp_ratio,
    "min_ins_ratio=f" => \$min_ins_ratio,
    "fasta_ref=s"     => \$fasta_ref,
    "bam_file=s"      => \$bam_file,
);

$outputfile = "$outputfile.csv";

open my $mpileup_fh, "-|", "samtools mpileup -A -r $chromosome -f $fasta_ref $bam_file";
open my $out_fh, ">", $outputfile;

say $out_fh "seq_id,pos,ref,a,c,g,t,del,consensus";

while (<$mpileup_fh>) {
    my ( $seqid, $pos, $ref, $depth, $read_bases, $read_quals ) = split;

    $ref =~ tr/acgtn/ACGTN/;
    next if $ref eq "N";

    clean_pileup(\$read_bases);

    my ( $inserts, $top_ins, $read_bases_no_ins ) = get_inserts($read_bases);

    clean_deletions(\$read_bases_no_ins);

    my ( $total_counts, $counts ) = count_bases($ref, $read_bases_no_ins);

    my $consensus = get_consensus_base($counts);

    output_snp( $seqid, $pos, $ref, $counts, $consensus, $total_counts,
        $min_cov, $min_snp_ratio, $out_fh );

    output_deletions( $seqid, $pos, $ref, $counts, $total_counts, $min_cov,
        $min_ins_ratio, $out_fh );

    output_insert( $seqid, $pos, $ref, $inserts, $top_ins, $counts, $min_cov,
        $min_ins_ratio, $out_fh );
}

close $mpileup_fh;
close $out_fh;

exit;

sub clean_pileup {
    my $read_bases = shift;

    $$read_bases =~ tr/acgt/ACGT/;
    $$read_bases =~ s/\^.|\$//g;    # Start of read + MapQ score | end of read
}

sub get_inserts {    # Capture sequences of variable length inserts
                     # and remove them from $read_bases
    my $read_bases = shift;

    my %inserts;
    for my $ins_len ( $read_bases =~ m/\+(\d+)/g ) {
        $inserts{$1}++ if $read_bases =~ s/\+(?:$ins_len)([ACGT]{$ins_len})//;
    }

    my ($top_ins) = sort { $inserts{$b} <=> $inserts{$a} } keys %inserts;

    return \%inserts, $top_ins, $read_bases;
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
        = scalar keys %$counts == 1
        ? keys %$counts
        : sort { $$counts{$b} <=> $$counts{$a} } keys %$counts;

    return $consensus;
}

sub output_snp {
    my ( $seqid, $pos, $ref, $counts, $consensus, $total_counts, $min_cov,
        $min_snp_ratio, $out_fh )
        = @_;

    return
        if $consensus eq "del"
        || $ref eq $consensus
        || $total_counts < $min_cov
        || $$counts{$consensus} < $min_snp_ratio * $total_counts;

    say $out_fh join ",", $seqid, $pos, $ref, $$counts{A}, $$counts{C},
        $$counts{G}, $$counts{T}, $$counts{del}, $consensus;
}

sub output_deletions {
    my ( $seqid, $pos, $ref, $counts, $total_counts,
        $min_cov, $min_ins_ratio, $out_fh )
        = @_;

    my $del_counts = $$counts{del};

    return
        if $del_counts == 0
        || $del_counts < $$counts{$ref} * $min_ins_ratio
        || $del_counts < $min_cov;

    say $out_fh join ",", $seqid, $pos, $ref, $$counts{A}, $$counts{C},
        $$counts{G}, $$counts{T}, $$counts{del}, "del";
}

sub output_insert {
    my ( $seqid, $pos, $ref, $inserts, $top_ins, $counts, $min_cov,
        $min_ins_ratio, $out_fh )
        = @_;

    return if scalar keys %$inserts == 0;

    my $top_count = $$inserts{$top_ins};
    my $total_ins_counts = sum values $inserts;

    return
        if $top_count < $$counts{$ref} * $min_ins_ratio
        || $top_count <= 0.5 * $total_ins_counts
        || $total_ins_counts < $min_cov;

    my @ins_bases = split //, $top_ins;

    my %ins_counts;
    my $ins_pos = sprintf "%02d", 1;
    for my $nt ( split //, $top_ins ) {
        $ins_counts{$ins_pos}{$nt} += $top_count;
        $ins_pos++;
    }

    for my $ins_pos ( sort { $a <=> $b } keys %ins_counts ) {

        my ($ins_base)
            = sort { $ins_counts{$ins_pos}{$b} <=> $ins_counts{$ins_pos}{$a} }
            keys %{ $ins_counts{$ins_pos} };

        say $out_fh join ",", $seqid, "$pos.$ins_pos", "INS",
            $ins_counts{$ins_pos}{A} // 0, $ins_counts{$ins_pos}{C} // 0,
            $ins_counts{$ins_pos}{G} // 0, $ins_counts{$ins_pos}{T} // 0,
            0, $ins_base;
    }
}
