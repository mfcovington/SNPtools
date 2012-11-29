#!/usr/bin/env perl
# genotyping_pileups.pl
# Mike Covington
# created: 2012-01-17
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use List::Util 'max';
use Scalar::Util 'looks_like_number';
use Getopt::Long;

###TODO: ¿¿¿¿incorporate into commander as sub????

my ( $mpileup_file, $snp_file, $par1_id, $par2_id, $out_file, $no_indels );

my $options = GetOptions(
    "mpileup=s"  => \$mpileup_file,
    "snp=s"      => \$snp_file,
    "par1_id=s"  => \$par1_id,
    "par2_id=s"  => \$par2_id,
    "out_file=s" => \$out_file,
    "no_indels"  => \$no_indels,
);

#build mpileup hash
open my $mpileup_fh, "<", $mpileup_file;
my %mpileups = map {
    chomp;
    my @delim = split /\t/;
    (
        $delim[1] => {
            'ref_base' => $delim[2],
            'mpileup'  => $delim[4],
        }
      )
} <$mpileup_fh>;

# build SNP/indel hash
open my $snp_fh, "<", $snp_file;
my $header = <$snp_fh>; # 0 = chr, 1 = pos, 2 = ref_base, 3 = snp_base, 4 = genotype, 5 = insert_position, 6 = SNP_CLASS
my %snps;
my $chromosome;
while (<$snp_fh>) {
    chomp;
    my @delim_snp = split /\t/;
    push @{ $snps{ $delim_snp[1] } },
      {
        'ref_base'   => $delim_snp[2],
        'snp_base'   => $delim_snp[3],
        'genotype'   => $delim_snp[4],
        'insert_pos' => $delim_snp[5],
        'snp_class'  => $delim_snp[6],
      };
    $chromosome = $delim_snp[0] unless defined $chromosome;
}

open my $out_fh, ">", $out_file;
for my $position ( sort { $a <=> $b } keys %mpileups ) {

    # skip positions for which snp record does not exist
    # - these should actually indicate an error, right?
    # - since that is the case, test to see how much overhead this check uses
    next unless exists $snps{$position};

    # skip indels if ignoring indels
    next
      if $no_indels
      && ( ${ $snps{$position} }[0]->{'ref_base'} eq 'INS'
        || ${ $snps{$position} }[0]->{'snp_base'} eq 'del' );

    # remove "junk" characters from mpileups:
    # < and > are CIGAR Ns (or reference skips)
    # ^ and $ denote that the base is at the start or end of a read
    # the character after ^ indicates mapping quality for the read
    $mpileups{$position}->{'mpileup'} =~ s/ < | > | \^.{1} | \$ //gx;

    # remove mpileups with deletions marked with a - sign.  They make no sense when I look at IGV, etc
    $mpileups{$position}->{'mpileup'} = "" if $mpileups{$position}->{'mpileup'} =~ m/ [\.,] - \d+ [ACGT]+ /gix;

    # collect insertions
    my @insertions;
    push ( @insertions, $1 ) while $mpileups{$position}->{'mpileup'} =~ m/ [\.,] \+ \d+ ([ACGT]+) /gix;

    # remove insertions
    $mpileups{$position}->{'mpileup'} =~ s/ [\.,] \+ \d+ [ACGT]+ //gix;

    # convert , and . to ref_base
    $mpileups{$position}->{'mpileup'} =~ s/ \. | , /$mpileups{$position}->{'ref_base'}/gix;

    #count
    my $A_count   = $mpileups{$position}->{'mpileup'} =~ tr/Aa//;
    my $C_count   = $mpileups{$position}->{'mpileup'} =~ tr/Cc//;
    my $T_count   = $mpileups{$position}->{'mpileup'} =~ tr/Tt//;
    my $G_count   = $mpileups{$position}->{'mpileup'} =~ tr/Gg//;
    my $del_count = $mpileups{$position}->{'mpileup'} =~ tr/*//;
    my $insert_count = scalar @insertions;
    my $total_count = $A_count + $C_count + $T_count + $G_count + $del_count + $insert_count;
    my $par1_count = 0;
    my $par2_count = 0;


    my (
        $par1_base,     $par2_base, @insert_bases,
        $insert_parent, $is_insert, $insert_seq,
        $del_parent,    $is_del,    $skip
    );

    # at least for now, skip polymorphisms where
    # one parent is a SNP and the other is an insertion
    # (very rare and not labeled as DIFF_SNPs)
    next
      if scalar @{ $snps{$position} } > 1
      && ${ $snps{$position} }[0]->{'insert_pos'} eq 'NA';

    # skip positions with problematic artifact created by removing "NOT" polymorphisms
    # that were actually "DIFF_SNPS" with inserts of varying size
    next
      if looks_like_number( ${ $snps{$position} }[0]->{'insert_pos'} )
      && ${ $snps{$position} }[0]->{'insert_pos'} > 1;

    # separate snp classes and determine genotype of parents from snp hash
    if ( ${ $snps{$position} }[0]->{'snp_class'} eq 'DIFF_SNP' ) {
        $par1_base = ${ $snps{$position} }[0]->{'snp_base'}
          if ${ $snps{$position} }[0]->{'genotype'} eq $par1_id;
        $par2_base = ${ $snps{$position} }[0]->{'snp_base'}
          if ${ $snps{$position} }[0]->{'genotype'} eq $par2_id;
    }
    elsif ( ${ $snps{$position} }[0]->{'ref_base'} eq 'INS' ) {
        $is_insert     = 1;
        $insert_parent = ${ $snps{$position} }[0]->{'genotype'};
        push( @insert_bases, ${ $snps{$position} }[$_]->{'snp_base'} )
          for ( 0 .. $#{ $snps{$position} } );
    }
    elsif ( ${ $snps{$position} }[0]->{'snp_base'} eq 'del' ) {
        $is_del     = 1;
        $del_parent = ${ $snps{$position} }[0]->{'genotype'};
        $par1_base  = ${ $snps{$position} }[0]->{'ref_base'}
          if $del_parent eq $par2_id;
        $par2_base = ${ $snps{$position} }[0]->{'ref_base'}
          if $del_parent eq $par1_id;
    }
    elsif ( ${ $snps{$position} }[0]->{'genotype'} eq $par1_id ) {
        $par1_base = ${ $snps{$position} }[0]->{'snp_base'};
        $par2_base = ${ $snps{$position} }[0]->{'ref_base'};
    }
    elsif ( ${ $snps{$position} }[0]->{'genotype'} eq $par2_id ) {
        $par1_base = ${ $snps{$position} }[0]->{'ref_base'};
        $par2_base = ${ $snps{$position} }[0]->{'snp_base'};
    }
    else { die "SOMETHING IS WRONG!!!"; }

    # count stuff up
    $par1_count = 0;
    $par2_count = 0;
    if ($is_insert) {
        $insert_seq = join "", @insert_bases;
        my $mpileup_inserts = join " ", @insertions;
        my $matched_insert_count =
          $mpileup_inserts =~ s/\b$insert_seq\b/\b$insert_seq\b/gi; ##added word boundaries so AA wouldn't match AAAA twice, for example
        if ( $insert_parent eq $par1_id ) {
            $par1_count = $matched_insert_count;
            $par2_count = max( $A_count, $C_count, $T_count, $G_count );
        }
        elsif ( $insert_parent eq $par2_id ) {
            $par1_count = max( $A_count, $C_count, $T_count, $G_count );
            $par2_count = $matched_insert_count;
        }
        else { die "SOMETHING IS WRONG!!!"; }
    }
    elsif ($is_del) {
        if ( $del_parent eq $par1_id ) {
            $par1_count = $del_count;
            $par2_count =
              $mpileups{$position}->{'mpileup'} =~ s/$par2_base/$par2_base/gi;
        }
        elsif ( $del_parent eq $par2_id ) {
            $par1_count =
              $mpileups{$position}->{'mpileup'} =~ s/$par1_base/$par1_base/gi;
            $par2_count = $del_count;
        }
        else { die "SOMETHING IS WRONG!!!"; }
    }
    else {
        $par1_count =
          $mpileups{$position}->{'mpileup'} =~ s/$par1_base/$par1_base/gi
          if defined $par1_base; # reminder: parent base undefined for inserts
        $par2_count =
          $mpileups{$position}->{'mpileup'} =~ s/$par2_base/$par2_base/gi
          if defined $par2_base; # reminder: parent base undefined for inserts
    }

    $par1_count = 0 unless $par1_count;
    $par2_count = 0 unless $par2_count;

    say $out_fh join( "\t", $chromosome, $position, $par1_count, $par2_count, $total_count );

}

close $mpileup_fh;
close $snp_fh;
close $out_fh;
exit;