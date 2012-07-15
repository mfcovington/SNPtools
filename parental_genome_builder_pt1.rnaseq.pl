#!/usr/bin/perl
# parental_genome_builder_pt1.rnaseq.pl
# Mike Covington
# created: 2012-01-03
#
# Description: 
#
# Change Log:
#	2012-01-15: adapted for use with RNAseq data only (SNPs + indels)
#	2012-02-23: Added wildcard to and altered structure of SNP files to account for difference caused by changing upstream scripts
#	2012-02-23: Changed capitalization of coverage files to account for difference caused by changing upstream scripts
#	2012-02-23: Changed `mkdir $output_dir` to `mkdir -p $output_dir`
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;

# TODO:
    # get rid of references to rna
    # change m82/pen to par1/par2
    # transition from chr # to chr ID
    # make write_to_master_snp subroutine

my $usage = <<USAGE_END;

USAGE:                                            #this usage statement doesn't reflect reality, yet
parental_genome_builder_pt1.rnaseq.pl
  --chr     <chromosome number>                   #need to change this to chr id eventually!!
  --snp1    <parent #1 snp file>                  #currently, everything is pointing to directories
  --snp2    <parent #2 snp file>
  --cov1    <parent #1 coverage.cov_nogap file>
  --cov2    <parent #2 coverage.cov_nogap file>
  --name1   <parent #1 name>
  --name2   <parent #2 name>
  --out     <output directory>
  --min_cov <minimum coverage [4]>
  --help

USAGE_END

#defaults
my $m82_rna_snp_dir = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/M82_SNPS";
my $pen_rna_snp_dir = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/PENN_SNPS";
my $cov_dir         = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/coverage_out";
my $output_dir      = "./";
my $par1_name       = "PAR1"
my $par1_name       = "PAR2"
my $min_cov         = 4;
my ( $chr_num, $help );

GetOptions(
    "chr=i"     => \$chr_num,
    "msnp=s"    => \$m82_rna_snp_dir,
    "psnp=s"    => \$pen_rna_snp_dir,
    "cov=s"     => \$cov_dir,
    "name1=s"   => \$par1_name,
    "name2=s"   => \$par2_name,
    "out=s"     => \$output_dir,
    "min_cov=i" => \$min_cov,
    "help"      => \$help,
);

die $usage if $help;
die $usage unless defined $chr_num;

my $chr_offset = $chr_num + 1;
my $chr_format;
if ( $chr_num < 10 ) {
    $chr_format = "0$chr_num";
}
else {
    $chr_format = $chr_num;
}

# input files
my $m82_rna_snps = glob($m82_rna_snp_dir . "/01.2.SNP_table.bwa_tophat_*.sorted.dupl_rm.$chr_offset.$chr_offset.nogap.gap.FILTERED.csv");
my $pen_rna_snps = glob($pen_rna_snp_dir . "/01.2.SNP_table.bwa_tophat_*.sorted.dupl_rm.$chr_offset.$chr_offset.nogap.gap.FILTERED.csv");
my $m82_rna_cov_file = $cov_dir . "/M82.SL2.40ch$chr_format.coverage.col";
my $pen_rna_cov_file = $cov_dir . "/PEN.SL2.40ch$chr_format.coverage.col";
open $m82_rna_cov_fh, "<", $m82_rna_cov_file;
open $pen_rna_cov_fh, "<", $pen_rna_cov_file;

# snp/cov build hashes
my %par1_snps = snp_hash_builder( $m82_rna_snps );
my %par2_snps = snp_hash_builder( $pen_rna_snps );
my %par1_cov  = cov_hash_builder( $m82_rna_cov_file );
my %par2_cov  = cov_hash_builder( $pen_rna_cov_file );

# get snp positions with good coverage
my @all_snp_pos = sort { $a <=> $b } keys %{ { %par1_snps, %par2_snps } };
my @good_cov_pos;
for my $pos (@all_snp_pos) {
    my $pos_par1_cov = $par1_cov{$pos} || next;    # || next bypasses if undef
    my $pos_par2_cov = $par2_cov{$pos} || next;
    next unless $pos_par1_cov > 4 && $pos_par1_cov > 4;
    push @good_cov_pos, $pos;
}

# write master snp file
`mkdir -p $output_dir`;
my $master_snp_file = "$output_dir/master_snp_list_rnaseq.chr$chr_format";
open $master_fh, ">", $master_snp_file;
say $master_fh join "\t", "chr", "pos", "ref_base", "snp_base", "genotype", "insert_position";    # insert_position will be decimal portion of position for insertions

for my $pos (@good_cov_pos) {
    if ( defined $par1_snps{$pos} ) {
        my $insert_length = $par1_snps{$pos}{insert} || 0;
        if ( $insert_length ) {
            for my $insert_idx ( 1 .. $insert_length ){
                $insert_idx = "0$insert_idx" if $insert_idx < 10;
                say $master_fh join "\t", "insert chr here", $pos, @{ $par1_snps{$pos}{$insert_idx} }[0,1], $par1_name, $insert_idx;
            }
        }
        else{
            say $master_fh join "\t", "insert chr here", $pos, @{ $par1_snps{$pos}{snp_del} }[0,1], $par1_name, "NA";
        }
    }
    if ( defined $par2_snps{$pos} ) {
        my $insert_length = $par2_snps{$pos}{insert} || 0;
        if ( $insert_length ) {
            for my $insert_idx ( 1 .. $insert_length ){
                $insert_idx = "0$insert_idx" if $insert_idx < 10;
                say $master_fh join "\t", "insert chr here", $pos, @{ $par2_snps{$pos}{$insert_idx} }[0,1], $par2_name, $insert_idx;
            }
        }
        else{
            say $master_fh join "\t", "insert chr here", $pos, @{ $par2_snps{$pos}{snp_del} }[0,1], $par2_name, "NA";
        }
    }
}

close ($master_fh);
exit;

###############
# subroutines #
###############

sub snp_hash_builder {
    my $par_snps_file = shift;
    open my $par_snps_fh, "<", $par_snps_file;
    my $header = <$par_snps_fh>;    # do I care about this?
    my %par_snps;
    while ( <$par_snps_fh> ) {
        chomp;
        my @par_snp = split /,/;
        my @par_pos = split /\./, $par_snp[1];
        if ( defined $par_pos[1] ) {
            $par_snps{$par_pos[0]}{$par_pos[1]} = [ @par_snp[ 2, 8 ], $par_pos[1] ];
            $par_snps{$par_pos[0]}{insert}      = $par_pos[1];
        }
        else {
            $par_snps{$par_pos[0]}{snp_del} = [ @par_snp[ 2, 8 ] ];
        }
    }
    close $par_snps_fh;
    return %par_snps;
}

sub cov_hash_builder {
    my $par_cov_file = shift;
    open my $par_cov_fh, "<", $par_cov_file;
    my %par_cov = map { chomp; @{ [ split /\t/ ] }[ 1 .. 2 ] } <$par_cov_fh>;
    close $par_cov_fh;
    return %par_cov;
}

__DATA__

Example of snp_hash:

{
    2068768   {
        snp_del   [
            [0] "A",
            [1] "C"
        ]
    },
    2069223   {
        01   [
            [0] "INS",
            [1] "A",
            [2] 01
        ],
        insert   01
    },
    2078525   {
        01   [
            [0] "INS",
            [1] "A",
            [2] 01
        ],
        02   [
            [0] "INS",
            [1] "C",
            [2] 02
        ],
        03   [
            [0] "INS",
            [1] "T",
            [2] 03
        ],
        insert   03
    },
    2119625   {
        snp_del   [
            [0] "del",
            [1] "C"
        ]
    },
    2119631   {
        snp_del   [
            [0] "T",
            [1] "A"
        ]
    }
}
