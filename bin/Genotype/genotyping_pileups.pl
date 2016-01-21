#!/usr/bin/env perl
# genotyping_pileups.pl
# Mike Covington
# created: 2012-01-17
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use List::Util 'max';
use Scalar::Util 'looks_like_number';
use Getopt::Long;

###TODO: ¿¿¿¿incorporate into commander as sub????

my ($snp_file,  $par1_id,   $par2_id, $out_file,
    $no_indels, $fasta_ref, $bam_file
);

my $options = GetOptions(
    "snp=s"       => \$snp_file,
    "par1_id=s"   => \$par1_id,
    "par2_id=s"   => \$par2_id,
    "out_file=s"  => \$out_file,
    "no_indels"   => \$no_indels,
    "fasta_ref=s" => \$fasta_ref,
    "bam_file=s"  => \$bam_file,
);

# build mpileup hash
my $mpileups = get_mpileups( $snp_file, $fasta_ref, $bam_file );

# build SNP/indel hash
my ( $snps, $chromosome ) = get_snps( $snp_file );

genotype( $mpileups, $snps, $chromosome, $out_file );

exit;


sub get_mpileups {
    my ( $snp_file, $fasta_ref, $bam_file ) = @_;

    my %mpileups;
    open my $mpileup_fh, "-|",
        "samtools mpileup -A -l $snp_file -f $fasta_ref $bam_file";
    for (<$mpileup_fh>) {
        my @delim = split /\t/;
        $mpileups{ $delim[1] } = {
            'ref_base' => $delim[2],
            'mpileup'  => $delim[4]
        };
    }
    close $mpileup_fh;

    return \%mpileups;
}

sub get_snps {
    my $snp_file = shift;

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
    close $snp_fh;

    return \%snps, $chromosome;
}

sub genotype {
    my ( $mpileups, $snps, $chromosome, $out_file ) = @_;

    open my $out_fh, ">", $out_file;
    for my $position ( sort { $a <=> $b } keys %$mpileups ) {

        # skip positions for which snp record does not exist
        # - these should actually indicate an error, right?
        # - since that is the case, test to see how much overhead this check uses
        next unless exists $$snps{$position};

        my $ref_base = ${ $$snps{$position} }[0]->{'ref_base'};
        my $snp_base = ${ $$snps{$position} }[0]->{'snp_base'};

        # skip indels if ignoring indels
        next
            if $no_indels
            && ( $ref_base eq 'INS' || $snp_base eq 'del' );

        my $mpileup = $$mpileups{$position}->{'mpileup'};

        # remove "junk" characters from mpileups:
        # < and > are CIGAR Ns (or reference skips)
        # ^ and $ denote that the base is at the start or end of a read
        # the character after ^ indicates mapping quality for the read
        $mpileup =~ s/<|>|\^.{1}|\$//g;

        # remove mpileups with deletions marked with a - sign.  They make no sense when I look at IGV, etc
        $mpileup = "" if $mpileup =~ m/[\.,]-\d+[ACGT]+/gi;

        # collect insertions
        my @insertions;
        push @insertions, $1 while $mpileup =~ m/[\.,]\+\d+([ACGT]+)/gi;

        # remove insertions
        $mpileup =~ s/[\.,]\+\d+[ACGT]+//gi;

        # convert , and . to ref_base
        $mpileup =~ s/\.|,/$$mpileups{$position}->{'ref_base'}/gi;

        my $insert_pos = ${ $$snps{$position} }[0]->{'insert_pos'};

        # at least for now, skip polymorphisms where
        # one parent is a SNP/del and the other is an insertion
        # (very rare and not labeled as DIFF_SNPs)
        next
            if scalar @{ $$snps{$position} } > 1
            && $insert_pos eq 'NA';

        # skip positions with problematic artifact created by removing "NOT" polymorphisms
        # that were actually "DIFF_SNPS" with inserts of varying size
        next
            if looks_like_number($insert_pos)
            && $insert_pos > 1;

        # separate snp classes and determine genotype of parents from snp hash
        my ($par1_base,     $par2_base, @insert_bases,
            $insert_parent, $is_insert, $insert_seq,
            $del_parent,    $is_del,    $skip
        );

        my $genotype = ${ $$snps{$position} }[0]->{'genotype'};

        if ( ${ $$snps{$position} }[0]->{'snp_class'} eq 'DIFF_SNP' ) {
            $par1_base = $snp_base
                if $genotype eq $par1_id;
            $par2_base = $snp_base
                if $genotype eq $par2_id;
        }
        elsif ( $ref_base eq 'INS' ) {
            $is_insert     = 1;
            $insert_parent = $genotype;
            push( @insert_bases, ${ $$snps{$position} }[$_]->{'snp_base'} )
                for ( 0 .. $#{ $$snps{$position} } );
        }
        elsif ( $snp_base eq 'del' ) {
            $is_del     = 1;
            $del_parent = $genotype;
            $par1_base  = $ref_base
                if $del_parent eq $par2_id;
            $par2_base = $ref_base
                if $del_parent eq $par1_id;
        }
        elsif ( $genotype eq $par1_id ) {
            $par1_base = $snp_base;
            $par2_base = $ref_base;
        }
        elsif ( $genotype eq $par2_id ) {
            $par1_base = $ref_base;
            $par2_base = $snp_base;
        }
        else { die "SOMETHING IS WRONG!!!"; }

        # count stuff up
        my $A_count   = $mpileup =~ tr/Aa//;
        my $C_count   = $mpileup =~ tr/Cc//;
        my $T_count   = $mpileup =~ tr/Tt//;
        my $G_count   = $mpileup =~ tr/Gg//;
        my $del_count = $mpileup =~ tr/*//;
        my $insert_count = scalar @insertions;
        my $total_count
            = $A_count
            + $C_count
            + $T_count
            + $G_count
            + $del_count
            + $insert_count;
        my $par1_count = 0;
        my $par2_count = 0;

        if ($is_insert) {
            $insert_seq = join "", @insert_bases;
            my $mpileup_inserts = join " ", @insertions;
            my $matched_insert_count
                = $mpileup_inserts =~ s/\b$insert_seq\b/\b$insert_seq\b/gi
                ; # added word boundaries so AA wouldn't match AAAA twice, for example
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
                $par2_count = $mpileup =~ s/$par2_base/$par2_base/gi;
            }
            elsif ( $del_parent eq $par2_id ) {
                $par1_count = $mpileup =~ s/$par1_base/$par1_base/gi;
                $par2_count = $del_count;
            }
            else { die "SOMETHING IS WRONG!!!"; }
        }
        else {
            # reminder: parent base undefined for inserts
            $par1_count = $mpileup =~ s/$par1_base/$par1_base/gi
                if defined $par1_base;
            $par2_count = $mpileup =~ s/$par2_base/$par2_base/gi
                if defined $par2_base;
        }

        $par1_count = 0 unless $par1_count;
        $par2_count = 0 unless $par2_count;

        say $out_fh join( "\t",
            $chromosome, $position, $par1_count,
            $par2_count, $total_count );

    }

    close $out_fh;
}
