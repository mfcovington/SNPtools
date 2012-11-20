#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2012-01-17
#
# Description:
#
use strict;
use warnings;
use autodie;
use List::Util qw[max];
use Getopt::Long;
use File::Basename;

my ( $pileup_file, $snp_file, $par1_id, $par2_id );
my $out_dir = "./";

my $options = GetOptions(
    "pileup=s"  => \$pileup_file,
    "snp=s"     => \$snp_file,
    "par1_id=s" => \$par1_id,
    "par2_id=s" => \$par2_id,
    "out_dir=s" => \$out_dir        ### make output file
);

my $out_file = fileparse($pileup_file);
$out_file =~ s/mpileup/genotyped/;

open my $pileup_fh, "<", $pileup_file;
open my $snp_fh,    "<", $snp_file;
open my $out_fh,    ">", "$out_dir/$out_file";

my $header = <$snp_fh>; # 0 = chr, 1 = pos, 2 = ref_base, 3 = snp_base, 4 = genotype, 5 = insert_position, 6 = SNP_CLASS
my @snp = ( 0, 0 );    # sets = 0
my $snp_line;

while ( my $pileup_line = <$pileup_fh> ) {
    chomp $pileup_line;
    my @pileup = split( /\t/ , $pileup_line ); # 0 = chr, 1 = pos, 2 = ref_base, 3 = filtered_coverage, 4 = genotype_per_read, 5 = base_quality

    #### extract genotype info from pileup ####
    #remove "junk" characters from $pileup[4]:
    $pileup[4] =~ s/<|>|\^.{1}|\$//g; # < and > are CIGAR Ns (or reference skips); ^ and $ denote that the base is at the start or end of a read; the character after ^ indicates mapping quality for the read
    #remove pileups with deletions marked with a - sign.  They make no sense when I look at IGV, etc
    $pileup[4] = "" if $pileup[4] =~ m/[\.,]-\d+([ACGT]+)/gi;
    # deal with insertions
    my @insertions;
    while ( $pileup[4] =~ m/[\.,]\+\d+([ACGT]+)/gi ) { #collect insertions
        push ( @insertions, $1 );
    }
    $pileup[4] =~ s/[\.,]\+\d+[ACGT]+//gi; #remove insertions
    #convert , and . to ref_base
    $pileup[4] =~ s/\.|,/$pileup[2]/gi;
    #count
    my $A_count   = $pileup[4] =~ tr/Aa//;
    my $C_count   = $pileup[4] =~ tr/Cc//;
    my $T_count   = $pileup[4] =~ tr/Tt//;
    my $G_count   = $pileup[4] =~ tr/Gg//;
    my $del_count = $pileup[4] =~ tr/*//;
    my $insert_count = scalar @insertions;
    my $total_count = $A_count + $C_count + $T_count + $G_count + $del_count + $insert_count;
    my $par1_count = 0;
    my $par2_count = 0;

    #read in snps until get a position match to the current pileup
    until ( $pileup[1] == $snp[1] ) { #read snps until position matches with pileup
        $snp_line = <$snp_fh>;
        chomp $snp_line;
        @snp = split( /\t/, $snp_line ); # 0 = chr, 1 = pos, 2 = ref_base, 3 = snp_base, 4 = genotype, 5 = insert_position, 6 = SNP_CLASS
    }
    my ( $m82_base, $pen_base, @insert_bases, $insert_parent, $is_insert, $insert_seq, $del_parent, $is_del, $skip );

    if ( $snp[5] =~ /^[0-9]+$/ && $snp[5] > 1 ) { #solve problem created by removing "NOT" snps that were actually "DIFF_SNPS" with Iinserts of varying size  ###PART 1###
        $skip = 1;
    }

    unless ( eof($snp_fh) ) {
        while ($pileup[1] == $snp[1]) { #while positions match
            if ( $snp[6] eq "DIFF_SNP" ) {
                $m82_base = $snp[3] if $snp[4] eq $par1_id;
                $pen_base = $snp[3] if $snp[4] eq $par2_id;
            }
            elsif ( $snp[2] eq "INS" ) {
                $insert_parent = $snp[4];
                push( @insert_bases, $snp[3] );
                $is_insert = 1;
            }
            elsif ( $snp[3] eq "del" ) {
                $del_parent = $snp[4];
                $is_del     = 1;
                $m82_base   = $snp[2] if $del_parent eq $par2_id;
                $pen_base   = $snp[2] if $del_parent eq $par1_id;
            }
            else {
                if ( $snp[4] eq $par1_id ) {
                    $m82_base = $snp[3];
                    $pen_base = $snp[2];
                }
                else {    #$snp[4] eq $par2_id
                    $m82_base = $snp[2];
                    $pen_base = $snp[3];
                }
            }

            $snp_line = <$snp_fh>;
            chomp $snp_line;
            @snp = split( /\t/, $snp_line );
            last if eof($snp_fh);
        }

        $par1_count = 0;
        $par2_count = 0;
        if ($is_insert) {
            $insert_seq = join( "", @insert_bases );
            my $pileup_inserts = join( " ", @insertions );
            my $matched_insert_count = $pileup_inserts =~ s/\b$insert_seq\b/\b$insert_seq\b/gi; ##added word boundaries so AA wouldn't match AAAA twice, for example
            if ( $insert_parent eq $par1_id ) {
                $par1_count = $matched_insert_count;
                $par2_count = max( $A_count, $C_count, $T_count, $G_count );
            }
            else {    #$insert_parent eq $par2_id
                $par1_count = max( $A_count, $C_count, $T_count, $G_count );
                $par2_count = $matched_insert_count;
            }
        }
        elsif ($is_del) {
            if ( $del_parent eq $par1_id ) {
                $par1_count = $del_count;
                $par2_count = $pileup[4] =~ s/$pen_base/$pen_base/gi;
            }
            else {    #$del_parent eq $par2_id
                $par1_count = $pileup[4] =~ s/$m82_base/$m82_base/gi;
                $par2_count = $del_count;
            }
        }
        else {
            $par1_count = $pileup[4] =~ s/$m82_base/$m82_base/gi;
            $par2_count = $pileup[4] =~ s/$pen_base/$pen_base/gi;
        }
    }

    if ($skip) { #solve problem created by removing "NOT" snps that were actually "DIFF_SNPS" with Iinserts of varying size
        $total_count = 0;
        $par1_count   = 0;
        $par2_count   = 0;
    }

    $par1_count = 0 unless $par1_count;
    $par2_count = 0 unless $par2_count;

    say $out_fh join( "\t", @pileup[0,1], $par1_count, $par2_count, $total_count );
}

close $pileup_fh;
close $snp_fh;
close $out_fh;
exit;