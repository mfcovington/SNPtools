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

my (
    $chromosome,
    $outputfile,
    $threshold_fraction_of_reads_matching_ref,
    $threshold_fraction_of_reads_matching_ref_for_indels,
    $threshold_number_of_reads,
    $fasta_ref,
    $sps_bam_file,
);

GetOptions(
    'chromosome:s' => \$chromosome,
    'o:s'          => \$outputfile,
    'n_reads:i'    => \$threshold_number_of_reads,
    'ref_freq:f'   => \$threshold_fraction_of_reads_matching_ref,
    'indel_freq:f' => \$threshold_fraction_of_reads_matching_ref_for_indels,
    'fasta_ref:s'  => \$fasta_ref,
    'bam_file:s'   => \$sps_bam_file,
);

say "Look for SNPs on $chromosome ($sps_bam_file)";
##OPEN OUTPUT FILE
my $out_filename = "$outputfile.csv";
open my $snp_out_fh, ">", $out_filename;

#MAKE HEADER LINE AND PRINT IT TO THE OUTPUT FILE
my $snp_header = "seq_id,pos,ref,a,c,g,t,del";
say $snp_out_fh "$snp_header";

my $sam = Bio::DB::Sam->new(
    -bam   => "$sps_bam_file",
    -fasta => "$fasta_ref"
);

my $snp_caller = sub {
    my ( $seqid, $pos, $p ) = @_;
    my %positions_in_reads_seen;
    my %distances_in_reads_seen;

    #       say "$seqid,$pos,$p";
    #GET BASE FOR THE REFERENCE SEQUENCE
    my $refbase = $sam->segment( $seqid, $pos, $pos )->dna;
    $refbase = uc($refbase);

    #   if($pos == 4197587){
    #   say "$pos";

    if ( $refbase ne "N" ) {
        my %ins_hash = (
            'inserted_bases'       => [],
            'ins_point_matches'    => 0,
            'ins_point_mismatches' => 0,
            'preINSbases' =>
              { 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0, 'del' => 0 }
        );

        my $alignment_size = 0;

        for my $pileup (@$p) {
            my $b = $pileup->alignment;

            #NAME OF READ AND ORIENTATION
            my $name = $b->qname;
            my $flag = $b->flag;
            $name .= "_$flag";
            $name =~ s/trimmed_//;

#       say "NAME: $name";
#GET BASE AND QUALIY SCORE FOR THE BASE...(THERE IS A PROBLEM WITH READ THAT HAVE DELETIONS
#SO I NEED TO CHECK IF THE BASE WE ARE LOOKING AT HAS A POS GREATER THAN THE END OF THE READ)
            my ($qbase);
            my $read_length = $b->query->length;
            if ( $pileup->qpos < $read_length ) {
                $qbase = substr( $b->qseq, $pileup->qpos, 1 );
            }
            else {
                $qbase = substr( $b->qseq, ( $read_length - 1 ), 1 );
            }

            #GETTING THE POSITION IN THE READS WHERE THE SNP IS LOCATED
            my $pos_in_read = $pileup->pos;

            #SKIP THE READ IF THE POSION IS FOUND INSIDE THE SPLICING
            my $cigar                           = $b->cigar_str;
            my $skip_this_read_it_is_a_splicing = 0;
            if ( $cigar =~ /N/ ) {

#           print "$pos - $name;  ref: $refbase; base: $qbase; length: $read_length; ";
#           CIGAR_parsing($cigar) if ($cigar=~ /N/);
                my @splic_pos = CIGAR_parsing($cigar);
                foreach (@splic_pos) {
                    $skip_this_read_it_is_a_splicing = 1
                      if ( $pos_in_read == $_ );
                }

                #           say "";
            }
            next if ( $skip_this_read_it_is_a_splicing == 1 );

#       if ($pileup->indel){
#           say "indel";
#           say "indel at pos: $pos\tread name: $name\tref base: $refbase\tbase: $qbase\tquality: $qscore";
#           say "pos in read: ",$pileup->qpos,"\tquery length:",$b->query->length,"\tindel length: ",$pileup->indel,"";
#IF THE INDEL IS AN INSERTION ADD BASES, QUALITIES, AND MORE TO %ins_hash'
#LATER WE SEE IF THESE ARE ENOUGH AND MAKE ONE LINE OF OUTPUT PER BASE INSERTED
            if ( $pileup->indel > 0 ) {
                my $ins_base =
                  substr( $b->qseq, ( $pileup->qpos ) + 1, $pileup->indel );

#               say "indel at pos: $pos\tread name: $name\tref base: $refbase\tbase: $qbase";
#               say "read pos: ",$pileup->qpos,"\tquery length:",$b->query->length,"\tindel length: ",$pileup->indel,"\t$ins_base";
#ADD THE INFO FROM THE BASE THAT BEFORE THE INSERTION "PREINS"
                $ins_hash{'preINSbases'}->{$qbase}++;

                #ADD THE POSITION IN THE READ WHERE THE INSERTION IS LOCATED
                $positions_in_reads_seen{'ins'}->{$qbase}->{$pos_in_read}++;
                $distances_in_reads_seen{'ins'}->{$qbase}
                  ->{ $read_length - $pos_in_read }++;
                $alignment_size++;
#### GET COVERAGE

         #ADD BASES AND QUALITIES TO THE INSERTION HASH BY POSITION IN INSERTION
                my @in_nucs = split '', $ins_base;
                for ( my $ins_pos = 0 ; $ins_pos < @in_nucs ; $ins_pos++ ) {
                    my $printf_pos;
                    if ( ( $ins_pos + 1 ) < 10 ) {
                        $printf_pos = "$pos.0" . ( $ins_pos + 1 );
                    }
                    else {
                        $printf_pos = "$pos." . ( $ins_pos + 1 );
                    }
                    $ins_hash{'bases'}->{$printf_pos}->{ $in_nucs[$ins_pos] }++;
                }

                #               say "\nscore: $score_string";
                #THIS IS GOING to hold the BASES OF THE INSERTION
                push @{ $ins_hash{'inserted_bases'} }, $ins_base;
                if ( $refbase eq $qbase ) {
                    $ins_hash{'ins_point_matches'}++;
                }
                else {
                    $ins_hash{'ins_point_mismatches'}++;
                }

            }

#IF IT IS A DELETION, JUST REMEMBE THE SIZE OF THE INSERTION AND IN WHICH READ IT WAS
            if ( $pileup->is_del ) {
                $positions_in_reads_seen{'del'}->{$pos_in_read}++;
                $distances_in_reads_seen{'del'}
                  ->{ $read_length - $pos_in_read }++;
                $alignment_size++;
#### GET COVERAGE
                next;
            }

            #NON DELETION READS, JUST ADD TO HASHES
            #THERE ARE NS IN THE READS, WE DO NOT COUNT THEM
            if ( $qbase !~ /[nN]/ ) {
                $positions_in_reads_seen{$qbase}->{$pos_in_read}++;
                $distances_in_reads_seen{$qbase}
                  ->{ $read_length - $pos_in_read }++;
                $alignment_size++;
#### GET COVERAGE
                #               say "NON INDEL: $pos\t$name\t$refbase\t$qbase";
            }
        }

        ###########################################################################
        #PARSING OUTPUT TO WRITE IN FILE
        ###########################################################################

        #GET NUMBER OF READS FOR DELETIONS, INSERTIONS AND REFERENCE
        my ( $insert_reads, $deletion_reads, $ref_reads, $all_reads ) =
          ( 0, 0, 0, 0 );
        $ref_reads = sum( values %{ $positions_in_reads_seen{$refbase} } )
          if ( defined $positions_in_reads_seen{$refbase} );
        $all_reads += sum( values %{ $positions_in_reads_seen{"A"} } )
          if ( defined $positions_in_reads_seen{"A"} );
        $all_reads += sum( values %{ $positions_in_reads_seen{"C"} } )
          if ( defined $positions_in_reads_seen{"C"} );
        $all_reads += sum( values %{ $positions_in_reads_seen{"G"} } )
          if ( defined $positions_in_reads_seen{"G"} );
        $all_reads += sum( values %{ $positions_in_reads_seen{"T"} } )
          if ( defined $positions_in_reads_seen{"T"} );

        if ( defined $positions_in_reads_seen{'ins'} ) {
            for ( keys %{ $positions_in_reads_seen{'ins'} } ) {
                $insert_reads +=
                  sum( values %{ $positions_in_reads_seen{'ins'}->{$_} } );
            }
        }
        $deletion_reads = sum( values %{ $positions_in_reads_seen{'del'} } )
          if ( defined $positions_in_reads_seen{'del'} );

#I CHECK FIRST IF IT IS AN INSERTION< THE AN DELETION AND THEN A SNP
#   say "$refbase\t$ref_reads\tall reads: $all_reads\tins: $insert_reads\tdel: $deletion_reads";
        if (
            $insert_reads > 1
            && (
                $insert_reads >= round(
                    $ref_reads *
                      $threshold_fraction_of_reads_matching_ref_for_indels
                )
            )
            && (
                $insert_reads >= round(
                    $all_reads *
                      ( $threshold_fraction_of_reads_matching_ref / 2 )
                )
            )
          )
        {
         #FIRST CHECK IF THERE IS A SNP IN THE BASE BEFORE THE INS
         #FOR THIS I"VE MADE A FUNCTION THAT CALCULATES THE BASES, POSITIONS ETC
            my ( $preins_line, $preins_pos_in_reads_seen,
                $preins_most_abundant_base_reads,
                $preins_most_abundant_base, $preins_total )
              = calculate_snp_in_pre_insertion( \%positions_in_reads_seen,
                $refbase );
            my ( $rpreins_line, $rpreins_pos_in_reads_seen ) =
              calculate_snp_in_pre_insertion( \%distances_in_reads_seen,
                $refbase );

#       say "$pos: $refbase\ttotal: $preins_total\tref: $ref_reads\tins: $insert_reads > ".round($ref_reads*$threshold_fraction_of_reads_matching_ref_for_indels)."\t\n$preins_line\n$rpreins_line";

            #THEN CHECK AS ALWAYS FOR SNPS...
            if (   $preins_most_abundant_base ne $refbase
                && $preins_most_abundant_base_reads >=
                ( $preins_total * $threshold_fraction_of_reads_matching_ref )
                && $preins_most_abundant_base_reads >=
                $threshold_number_of_reads )
            {
                my $line =
"$seqid,$pos,$refbase,$preins_line,$preins_most_abundant_base,"
                  . calculate_attributes($preins_pos_in_reads_seen) . ","
                  . calculate_attributes($rpreins_pos_in_reads_seen);
                say $snp_out_fh "$line";

                #say "SNP: $line";
            }

#NOW ADD INSERTION LINES TO THE OUTFILE IN  A LOOP (INFO IN %ins_hash)
#FIRST GET ATTRIBUTE LINE FOR THE PREINS BASE (CAN BE DONE FOR EACH BASE IN THE INSERTION)
            my %new_hash_with_positions_for_ins   = ();
            my %r_new_hash_with_positions_for_ins = ();
            foreach my $nt ( ( 'A', 'C', 'G', 'T', 'del' ) ) {
                foreach ( keys %{ $positions_in_reads_seen{'ins'}->{$nt} } ) {
                    $new_hash_with_positions_for_ins{$_} =
                      $positions_in_reads_seen{'ins'}->{$nt}->{$_};
                }
                foreach ( keys %{ $distances_in_reads_seen{'ins'}->{$nt} } ) {
                    $r_new_hash_with_positions_for_ins{$_} =
                      $distances_in_reads_seen{'ins'}->{$nt}->{$_};
                }
            }

            #PASS HASH TO SUB TO CALCULATE ATTRBUTES (POSITIONS OF READS)
            my $preins_attr_line =
                calculate_attributes( \%new_hash_with_positions_for_ins ) . ","
              . calculate_attributes( \%r_new_hash_with_positions_for_ins );

            for
              my $new_pos ( sort { $a <=> $b } keys( %{ $ins_hash{'bases'} } ) )
            {
                my $line = "$seqid,$new_pos,INS,";

                #FIRST GET ABUNDANCES
                my ( $nt_line, $most_abundant_ins, $most_abundant_ins_reads ) =
                  ( "", 0, 0 );
                foreach my $nt ( ( 'A', 'C', 'G', 'T', 'del' ) ) {
                    if ( defined( $ins_hash{'bases'}->{$new_pos}->{$nt} ) ) {
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
                $line .= $nt_line . "$most_abundant_ins,$preins_attr_line";
                say $snp_out_fh "$line";

                #           say "INS: $line";
            }

        }
        elsif (
            $deletion_reads > 1
            && (
                $deletion_reads >= round(
                    $ref_reads *
                      $threshold_fraction_of_reads_matching_ref_for_indels
                )
            )
          )
        {
            my ( $nt_line, $most_abundant_base, $most_abundant_base_reads,
                $total )
              = do_quantification_loop( \%positions_in_reads_seen );
            my $line = "$seqid,$pos,$refbase,$nt_line,del,";
            $line .=
                calculate_attributes( \%{ $positions_in_reads_seen{'del'} } )
              . ","
              . calculate_attributes( \%{ $distances_in_reads_seen{'del'} } );
            say $snp_out_fh "$line";

            #       say "DEL: $line";
        }
        else {

            my ( $nt_line, $most_abundant_base, $most_abundant_base_reads,
                $total )
              = do_quantification_loop( \%positions_in_reads_seen );

#       say "$most_abundant_base ne $refbase && $most_abundant_base_reads > ($alignment_size*$threshold_fraction_of_reads_matching_ref) && $most_abundant_base_reads >= $threshold_number_of_reads";
            if (   $most_abundant_base ne $refbase
                && $most_abundant_base_reads >=
                ( $alignment_size * $threshold_fraction_of_reads_matching_ref )
                && $most_abundant_base_reads >= $threshold_number_of_reads )
            {
                my $line = "$seqid,$pos,$refbase,$nt_line,$most_abundant_base,";
                $line .= calculate_attributes(
                    \%{ $positions_in_reads_seen{$most_abundant_base} } )
                  . ","
                  . calculate_attributes(
                    \%{ $distances_in_reads_seen{$most_abundant_base} } );
                say $snp_out_fh "$line";

                #           say "SNP: $line";
            }
        }
        %positions_in_reads_seen = ();
        %distances_in_reads_seen = ();

    }

};

# my @targets = $sam->seq_ids;
# @targets = sort {$a cmp $b} @targets;
#say "@targets";
# my $chromosome_counter=0;
# foreach my $chromosome (@targets[$chromosome..$chrm_end]){
# $chromosome_counter++;
#CAPITALIZE FOR COLOMA!!!
#   substr($chromosome,0,2,"SL");
# say "$chromosome_counter. $chromosome\tlength: $chr_lengths{$chromosome}";
# `sleep 1`;
$sam->fast_pileup( $chromosome, $snp_caller );

#   $sam->fast_pileup("$chromosome:55900..56900",$snp_caller);
# }
close $snp_out_fh;

exit;

###################################################
#SUBROUTINES

sub round {
    my ($number) = shift;
    return int( $number + .5 );
}

sub calculate_attributes {
    my ($positions_from_reads_hash) = shift;
    my %quarters;
    foreach my $p ( keys %{$positions_from_reads_hash} ) {
        $quarters{ ( int( $p / 16.3 ) + 1 ) } +=
          $positions_from_reads_hash->{$p};
    }
    my ( $quart_count, $out_string ) = ( 0, "" );
    foreach my $s ( ( 1 .. 5 ) ) {
        if ( defined $quarters{$s} ) {
            $quart_count++;
            $out_string .= "$quarters{$s},";
        }
        else {
            $out_string .= "0,";
        }
    }
    $out_string .=
        sum( values %{$positions_from_reads_hash} )
      . ",$quart_count,"
      . scalar( keys %{$positions_from_reads_hash} );
    return $out_string;
}

sub calculate_snp_in_pre_insertion {
    my ( $positions_in_reads_hash, $refbase ) = @_;

    #DO THE LOOP TO CALCULATE ABUNDANCES OF NTS, ETC
    my ( $nt_line, $total, $most_abundant_base, $most_abundant_base_reads ) =
      ( "", 0, "del", 0 );
    foreach my $nt ( ( 'A', 'C', 'G', 'T', 'del' ) ) {
        my $number_of_reads = 0;
        $number_of_reads += sum( values %{ $positions_in_reads_hash->{$nt} } )
          if ( defined $positions_in_reads_hash->{$nt} );
        $number_of_reads +=
          sum( values %{ $positions_in_reads_hash->{'ins'}->{$nt} } )
          if ( defined $positions_in_reads_hash->{'ins'}->{$nt} );
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
    foreach my $nt ( ( 'A', 'C', 'G', 'T', 'del' ) ) {
        my $number_of_reads = 0;
        $number_of_reads = sum( values %{ $positions_in_reads_hash->{$nt} } )
          if ( defined $positions_in_reads_hash->{$nt} );
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

#OLD CIGAR SUBROUTINE:
# sub CIGAR_parsing {
#   my $cigar_tr =  shift;
#   my @parts = split /\d+N/,$cigar_tr;
#   my @bad_pos;
#   my $cum_sum = 1;
#   foreach(@parts){
#       $_ =~ s/M//g;
#       $cum_sum += $_;
#       push @bad_pos, $cum_sum;
#   }
#   pop @bad_pos;
# # print "@bad_pos";
#   return @bad_pos;
# }

#NEW CIGAR SUBROUTINE FROM PEPE ON 2011-12-20:
sub CIGAR_parsing {
    my $cigar_tr = shift;
    my @parts = split /\d*N/, $cigar_tr;
    pop @parts;
    my @bad_pos;
    my $cum_sum = 1;
    foreach my $sub_cigar (@parts) {
        my @numbers = split( /[A-Z]/, $sub_cigar );
        my @letters = split( /\d*/,   $sub_cigar );
        shift @letters;
        for ( my $i = 0 ; $i < @numbers ; $i++ ) {
            $cum_sum += $numbers[$i] if ( $letters[$i] =~ /(M|I)/ );
            push @bad_pos, $cum_sum;
        }
    }
    return @bad_pos;
}

