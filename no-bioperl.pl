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

open my $mpileup_fh, "-|", "samtools mpileup -A -r $chromosome:11562-11572 -f $fasta_ref $bam_file";
# my $stop = 0;

while (<$mpileup_fh>) {
    my ( $seqid, $pos, $ref, $depth, $read_bases, $read_quals ) = split;
    my %counts;

    for my $base (qw(A C G T del)) {
        $counts{$base} = 0;
    }

say "------------------------------------";
say "depth: $depth";
$read_bases =~ tr/acgt/ACGT/;
say $read_bases;

my %inserts;
$inserts{$_}++ for $read_bases =~ m/\+\d+([ACGT]+)/ig;
# $inserts{$_}++ for $read_bases =~ m/[.,]\+\d+([ACGT]+)/ig;
p %inserts;
$read_bases =~ s/\+\d+([ACGT]+)//ig;
# $read_bases =~ s/[.,]\+\d+([ACGT]+)//ig;

# tr/acgt/ACGT/ for @inserts;

my ($insert) = sort { $inserts{$b} <=> $inserts{$a} } keys %inserts;
my $insert_count = defined $insert ? $inserts{$insert} : 0;
# my ( $insert ) = keys %inserts;
say "insert: ", defined $insert ? $insert : "NA";
say "insert count: ", defined $insert ? $inserts{$insert} : "NA";


say $read_bases;
    $counts{A}++ for $read_bases =~ m/A/ig;
    $counts{C}++ for $read_bases =~ m/C/ig;
    $counts{G}++ for $read_bases =~ m/G/ig;
    $counts{T}++ for $read_bases =~ m/T/ig;
    $counts{$ref}++ for $read_bases =~ m/[.,]/g;
    $counts{del}++ for $read_bases =~ m/\*/ig;
    # $counts{ins} += $inserts{$_} for keys %inserts;
    # $counts{ins} = sum values %inserts // 0;
    # $counts{ins} = ( sum values %inserts ) // 0;
p %counts;
my @vals = values %counts;
say "vals: @vals";
# my $consensus;
my ($consensus) = scalar keys %counts == 1 ? keys %counts : sort { $counts{$b} <=> $counts{$a} } keys %counts;
say "consensus: $consensus";
# exit;


    say join ",", $seqid, $pos, $ref, $counts{A}, $counts{C}, $counts{G}, $counts{T}, $counts{del}, $consensus;

    if ( $insert_count > $counts{$ref} * 0.5 ){
        my ( $insert ) = keys %inserts;
        my @ins_seq = split //, $insert;
        my $count = 0;
        my %ins_counts;
        for my $ins_base (@ins_seq) {
            ++$count;
            $ins_counts{$ins_base} = $insert_count;
            say join ",", $seqid, "$pos.$count", "INS", $ins_counts{A} // 0,$ins_counts{C} // 0, $ins_counts{G} // 0, $ins_counts{T} // 0, 0, $ins_base;
            delete $ins_counts{$ins_base};
        }
    }


    # say join ",", $seqid, $pos, $ref, $counts{A}, $counts{C}, $counts{G}, $counts{T}, $counts{del}, $counts{ins}, "consensus";

    # $stop++;
    # die if $stop == 10;
}

close $mpileup_fh;


