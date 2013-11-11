#!/usr/bin/env perl
#perl 01.1.SNP_calling_homos.pl -fasta_ref /Volumes/SolexaRAID/Solexa_runs_Data/00.Downloaded_references_and_others/S_lycopersicum_chromosomes.2.40.fa -chrm_start 0 -chrm_end 0 -n_reads 4 -ref_freq 0.66 -indel_freq 0.33 -o SNP_table.PEN.final.Picard_and_GATK.1.1 -bam_file ../2.40.Chromosomes_alignments/BWA_merged_files/PEN.final.Picard_and_GATK.bam
#adapted from Pepe's script by mfc
#2012-01-06: CIGAR subroutine altered per Pepe's suggestion
use strict;
use warnings;
use feature 'say';
use autodie;
use Getopt::Long;
use Bio::DB::Sam;
use List::Util 'sum';
use List::MoreUtils 'any';

my (
    $chromosome,
    $outputfile,
    $min_ratio,        # threshold for fraction of reads matching ref
    $min_ratio_ins,    # threshold for fraction of reads matching ref for indels
    $min_depth,        # threshold for number of reads
    $fasta_ref,
    $bam_file,
);

GetOptions(
    'chromosome:s' => \$chromosome,
    'o:s'          => \$outputfile,
    'n_reads:i'    => \$min_depth,
    'ref_freq:f'   => \$min_ratio,
    'indel_freq:f' => \$min_ratio_ins,
    'fasta_ref:s'  => \$fasta_ref,
    'bam_file:s'   => \$bam_file,
);

my $sam = Bio::DB::Sam->new(
    -bam   => "$bam_file",
    -fasta => "$fasta_ref"
);

say "Looking for SNPs on $chromosome ($bam_file)";

my $out_filename = "$outputfile.csv";
open my $snp_out_fh, ">", $out_filename;

my $snp_header = "seq_id,pos,ref,a,c,g,t,del,consensus";
say $snp_out_fh "$snp_header";

$sam->fast_pileup( $chromosome, \&snp_caller );

close $snp_out_fh;

exit;

###################################################
#SUBROUTINES

sub snp_caller {
    my ( $seqid, $pos, $p ) = @_;
    my %positions;
    my %distances;

    #GET BASE FOR THE REFERENCE SEQUENCE
    my $refbase = $sam->segment( $seqid, $pos, $pos )->dna;
    $refbase = uc $refbase;

    return if $refbase eq "N";
    return if scalar @$p < $min_depth;

    my %ins_hash = (
        'inserted_bases'       => [],
        'ins_point_matches'    => 0,
        'ins_point_mismatches' => 0,
        'preINSbases' => { 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0, 'del' => 0 }
    );

    my $coverage = 0;

    for my $pileup (@$p) {

        next if $pileup->is_refskip;

        my $b = $pileup->alignment;

        my $qpos = $pileup->qpos;
        my $qseq = $b->qseq;
        next unless $b->qscore->[$qpos] > 25;

        my $qbase = substr( $qseq, $qpos, 1 );
        next if $qbase eq 'N';

        my $read_length = $b->query->length;

        #GETTING THE POSITION IN THE READS WHERE THE SNP IS LOCATED
        my $pos_in_read = $pileup->pos;

        #SKIP THE READ IF THE POSITION IS FOUND INSIDE A GAP
        my $cigar = $b->cigar_str;
        if ( $cigar =~ /N/ ) {
            my @gap_pos = parse_cigar($cigar);
            next if any { $_ == $pos_in_read } @gap_pos;
        }

 #IF THE INDEL IS AN INSERTION ADD BASES, QUALITIES, AND MORE TO %ins_hash'
 #LATER WE SEE IF THESE ARE ENOUGH AND MAKE ONE LINE OF OUTPUT PER BASE INSERTED
        if ( $pileup->indel > 0 ) {
            my $ins_base = substr( $qseq, ($qpos) + 1, $pileup->indel );

            #ADD THE INFO FROM THE BASE THAT BEFORE THE INSERTION "PREINS"
            $ins_hash{'preINSbases'}->{$qbase}++;

            #ADD THE POSITION IN THE READ WHERE THE INSERTION IS LOCATED
            $positions{'ins'}->{$qbase}->{$pos_in_read}++;
            $distances{'ins'}->{$qbase}->{ $read_length - $pos_in_read }++;

         #ADD BASES AND QUALITIES TO THE INSERTION HASH BY POSITION IN INSERTION
            my @in_nucs = split '', $ins_base;
            for my $ins_pos ( 0 .. $#in_nucs ) {
                my $printf_pos;
                if ( ( $ins_pos + 1 ) < 10 ) {
                    $printf_pos = "$pos.0" . ( $ins_pos + 1 );
                }
                else {
                    $printf_pos = "$pos." . ( $ins_pos + 1 );
                }
                $ins_hash{'bases'}->{$printf_pos}->{ $in_nucs[$ins_pos] }++;
            }

            #THIS IS GOING to hold the BASES OF THE INSERTION
            push @{ $ins_hash{'inserted_bases'} }, $ins_base;
            if ( $refbase eq $qbase ) {
                $ins_hash{'ins_point_matches'}++;
            }
            else {
                $ins_hash{'ins_point_mismatches'}++;
            }

        }

#IF IT IS A DELETION, JUST REMEMBER THE SIZE OF THE INSERTION AND IN WHICH READ IT WAS
        elsif ( $pileup->is_del ) {
            $positions{'del'}->{$pos_in_read}++;
            $distances{'del'}->{ $read_length - $pos_in_read }++;
        }

        #NON DELETION READS, JUST ADD TO HASHES
        #THERE ARE NS IN THE READS, WE DO NOT COUNT THEM
        else {
            $positions{$qbase}->{$pos_in_read}++;
            $distances{$qbase}->{ $read_length - $pos_in_read }++;
        }

        $coverage++;

    }

    ###########################################################################
    #PARSING OUTPUT TO WRITE IN FILE
    ###########################################################################

    #GET NUMBER OF READS FOR DELETIONS, INSERTIONS AND REFERENCE
    my ( $insert_reads, $deletion_reads, $ref_reads, $all_reads ) = (0) x 4;
    $ref_reads = sum values %{ $positions{$refbase} }
      if exists $positions{$refbase};

    for (qw(A C G T)) {
        $all_reads += sum values %{ $positions{$_} } if exists $positions{$_};
    }

    if ( exists $positions{'ins'} ) {
        $insert_reads += sum values %{ $positions{'ins'}->{$_} }
          for keys %{ $positions{'ins'} }
    }

    $deletion_reads = sum values %{ $positions{'del'} }
      if exists $positions{'del'};

    #I CHECK FIRST IF IT IS AN INSERTION< THE AN DELETION AND THEN A SNP
    if (   $insert_reads > 1
        && ( $insert_reads >= round( $ref_reads * $min_ratio_ins ) )
        && ( $insert_reads >= round( $all_reads * ( $min_ratio / 2 ) ) ) )
    {
        #FIRST CHECK IF THERE IS A SNP IN THE BASE BEFORE THE INS
        #FOR THIS I"VE MADE A FUNCTION THAT CALCULATES THE BASES, POSITIONS ETC
        my ( $preins_line, $preins_pos, $preins_most_abundant_base_reads,
            $preins_most_abundant_base, $preins_total )
          = calculate_snp_in_pre_insertion( \%positions, $refbase );
        my ( $rpreins_line, $rpreins_pos ) =
          calculate_snp_in_pre_insertion( \%distances, $refbase );

        #THEN CHECK AS ALWAYS FOR SNPS...
        if (   $preins_most_abundant_base ne $refbase
            && $preins_most_abundant_base_reads >=
            ( $preins_total * $min_ratio )
            && $preins_most_abundant_base_reads >= $min_depth )
        {
            my $line =
                "$seqid,$pos,$refbase,$preins_line,$preins_most_abundant_base";
            say $snp_out_fh "$line";

            #say "SNP: $line";
        }

#NOW ADD INSERTION LINES TO THE OUTFILE IN  A LOOP (INFO IN %ins_hash)
#FIRST GET ATTRIBUTE LINE FOR THE PREINS BASE (CAN BE DONE FOR EACH BASE IN THE INSERTION)
        my %new_hash_with_positions_for_ins   = ();
        my %r_new_hash_with_positions_for_ins = ();
        for my $nt ( qw( A C G T del ) ) {
            for ( keys %{ $positions{'ins'}->{$nt} } ) {
                $new_hash_with_positions_for_ins{$_} =
                  $positions{'ins'}->{$nt}->{$_};
            }
            for ( keys %{ $distances{'ins'}->{$nt} } ) {
                $r_new_hash_with_positions_for_ins{$_} =
                  $distances{'ins'}->{$nt}->{$_};
            }
        }

        #PASS HASH TO SUB TO CALCULATE ATTRBUTES (POSITIONS OF READS)

        for my $new_pos ( sort { $a <=> $b } keys( %{ $ins_hash{'bases'} } ) ) {
            my $line = "$seqid,$new_pos,INS,";

            #FIRST GET ABUNDANCES
            my ( $nt_line, $most_abundant_ins, $most_abundant_ins_reads ) =
              ( "", 0, 0 );
            for my $nt ( qw( A C G T del ) ) {
                if ( exists( $ins_hash{'bases'}->{$new_pos}->{$nt} ) ) {
                    $nt_line .= "$ins_hash{'bases'}->{$new_pos}->{$nt},";
                    if ( $most_abundant_ins_reads <
                        $ins_hash{'bases'}->{$new_pos}->{$nt} )
                    {
                        $most_abundant_ins = $nt;
                        $most_abundant_ins_reads =
                          $ins_hash{'bases'}->{$new_pos}->{$nt};
                    }
                }
                else {
                    $nt_line .= "0,";
                }
            }
            $line .= $nt_line . "$most_abundant_ins";
            say $snp_out_fh "$line";
        }

    }
    elsif ( $deletion_reads > 1
        && ( $deletion_reads >= round( $ref_reads * $min_ratio_ins ) ) )
    {
        my ( $nt_line, $most_abundant_base, $most_abundant_base_reads, $total )
          = do_quantification_loop( \%positions );
        my $line = "$seqid,$pos,$refbase,$nt_line,del";
        say $snp_out_fh "$line";
    }
    else {

        my ( $nt_line, $most_abundant_base, $most_abundant_base_reads, $total )
          = do_quantification_loop( \%positions );

        if (   $most_abundant_base ne $refbase
            && $most_abundant_base_reads >= ( $coverage * $min_ratio )
            && $most_abundant_base_reads >= $min_depth )
        {
            my $line = "$seqid,$pos,$refbase,$nt_line,$most_abundant_base";
            say $snp_out_fh "$line";
        }

    }
    %positions = ();
    %distances = ();

};

sub round {
    my ($number) = shift;
    return int( $number + .5 );
}

sub calculate_snp_in_pre_insertion {
    my ( $positions_in_reads_hash, $refbase ) = @_;

    #DO THE LOOP TO CALCULATE ABUNDANCES OF NTS, ETC
    my ( $nt_line, $total, $most_abundant_base, $most_abundant_base_reads ) =
      ( "", 0, "del", 0 );
    for my $nt ( qw( A C G T del ) ) {
        my $number_of_reads = 0;
        $number_of_reads += sum( values %{ $positions_in_reads_hash->{$nt} } )
          if ( exists $positions_in_reads_hash->{$nt} );
        $number_of_reads +=
          sum( values %{ $positions_in_reads_hash->{'ins'}->{$nt} } )
          if ( exists $positions_in_reads_hash->{'ins'}->{$nt} );
        $total += $number_of_reads;
        $nt_line .= "$number_of_reads,";
        if ( $most_abundant_base_reads < $number_of_reads ) {
            $most_abundant_base       = $nt;
            $most_abundant_base_reads = $number_of_reads;
        }
    }
    chop $nt_line;

#MERGE HASH WITH POSITIONS SEEN FOR THE PREINSERTION BASE (THE ONES WITH AND WITH OUT INSERTIONS)...
    my %new_merged_hash;
    for ( keys %{ $positions_in_reads_hash->{$most_abundant_base} } ) {
        $new_merged_hash{$_} =
          $positions_in_reads_hash->{$most_abundant_base}->{$_};
    }
    for ( keys %{ $positions_in_reads_hash->{'ins'}->{$most_abundant_base} } ) {
        $new_merged_hash{$_} =
          $positions_in_reads_hash->{'ins'}->{$most_abundant_base}->{$_};
    }

    return ( $nt_line, \%new_merged_hash, $most_abundant_base_reads,
        $most_abundant_base, $total );
}

sub do_quantification_loop {
    my $positions_in_reads_hash = shift;
    my ( $nt_line, $total, $most_abundant_base, $most_abundant_base_reads ) =
      ( "", 0, "del", 0 );
    for my $nt ( qw( A C G T del ) ) {
        my $number_of_reads = 0;
        $number_of_reads = sum( values %{ $positions_in_reads_hash->{$nt} } )
          if ( exists $positions_in_reads_hash->{$nt} );
        $total += $number_of_reads;
        $nt_line .= "$number_of_reads,";
        if ( $most_abundant_base_reads < $number_of_reads ) {
            $most_abundant_base       = $nt;
            $most_abundant_base_reads = $number_of_reads;
        }
    }
    chop $nt_line;
    return ( $nt_line, $most_abundant_base, $most_abundant_base_reads, $total );
}

#NEW CIGAR SUBROUTINE FROM PEPE ON 2011-12-20:
sub parse_cigar {
    my @parts = split /\d*N/, shift;
    pop @parts;
    my @gap_start_pos;
    my $pos_sum = 1;
    for (@parts) {
        my @numbers = split /[A-Z]/;
        my @letters = split /\d*/;
        shift @letters;
        for my $i ( 0 .. $#numbers ) {
            $pos_sum += $numbers[$i] if ( $letters[$i] =~ /M|I/ );
            push @gap_start_pos, $pos_sum;
        }
    }
    return @gap_start_pos;
}

