#!/usr/bin/env perl
# parental_genome_builder_pt1.rnaseq.pl
# Mike Covington
# created: 2012-01-03
#
# Description: 
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use File::Path 'make_path';

# TODO:
    # make write_master_snp subroutine

my $usage = <<USAGE_END;

USAGE:
parental_genome_builder_pt1.rnaseq.pl
  --chr_id     chromosome ID
  --snp1       parent #1 snp file
  --snp2       parent #2 snp file
  --cov1       parent #1 coverage.cov_nogap file
  --cov2       parent #2 coverage.cov_nogap file
  --name1      parent #1 name
  --name2      parent #2 name
  --out        output directory
  --min_cov    minimum coverage [4]
  --help       (displays this usage statement)

USAGE_END

# defaults/options
my $output_dir = "./";
my $par1_name  = "PAR1";
my $par2_name  = "PAR2";
my $min_cov    = 4;
my ( $chr_id, $par1_snp_file, $par2_snp_file, $par1_cov_file, $par2_cov_file, $help );

GetOptions(
    "chr_id=s"  => \$chr_id,
    "snp1=s"    => \$par1_snp_file,
    "snp2=s"    => \$par2_snp_file,
    "cov1=s"    => \$par1_cov_file,
    "cov2=s"    => \$par2_cov_file,
    "name1=s"   => \$par1_name,
    "name2=s"   => \$par2_name,
    "out=s"     => \$output_dir,
    "min_cov=i" => \$min_cov,
    "help"      => \$help,
);

die $usage if $help;
die $usage
  unless defined $chr_id
      && defined $par1_snp_file
      && defined $par2_snp_file
      && defined $par1_cov_file
      && defined $par2_cov_file;

# snp/cov build hashes
my %par1_snps = snp_hash_builder($par1_snp_file);
my %par2_snps = snp_hash_builder($par2_snp_file);
my %par1_cov  = cov_hash_builder($par1_cov_file);
my %par2_cov  = cov_hash_builder($par2_cov_file);

# get snp positions with good coverage
my @all_snp_pos = sort { $a <=> $b } keys %{ { %par1_snps, %par2_snps } };
my @good_cov_pos;
for my $pos (@all_snp_pos) {
    my $pos_par1_cov = $par1_cov{$pos} || next;    # || next bypasses if undef
    my $pos_par2_cov = $par2_cov{$pos} || next;
    next if $pos_par1_cov < $min_cov || $pos_par1_cov < $min_cov;
    push @good_cov_pos, $pos;
}

# write master snp file
make_path($output_dir);
my $master_snp_file = "$output_dir/master_snp_list.$par1_name.vs.$par2_name.$chr_id";
open my $master_fh, ">", $master_snp_file;
say $master_fh join "\t", "chr", "pos", "ref_base", "snp_base", "genotype", "insert_position";

for my $pos (@good_cov_pos) {
    if ( defined $par1_snps{$pos} ) {
        my $insert_length = $par1_snps{$pos}{insert} || 0;
        if ($insert_length) {
            for my $insert_idx ( 1 .. $insert_length ) {
                $insert_idx = "0$insert_idx" if $insert_idx < 10;
                say $master_fh join "\t", $chr_id, $pos, @{ $par1_snps{$pos}{$insert_idx} }[ 0, 1 ], $par1_name, $insert_idx;
            }
        }
        else{
            say $master_fh join "\t", $chr_id, $pos, @{ $par1_snps{$pos}{snp_del} }[ 0, 1 ], $par1_name, "NA";
        }
    }
    if ( defined $par2_snps{$pos} ) {
        my $insert_length = $par2_snps{$pos}{insert} || 0;
        if ($insert_length) {
            for my $insert_idx ( 1 .. $insert_length ) {
                $insert_idx = "0$insert_idx" if $insert_idx < 10;
                say $master_fh join "\t", $chr_id, $pos, @{ $par2_snps{$pos}{$insert_idx} }[ 0, 1 ], $par2_name, $insert_idx;
            }
        }
        else{
            say $master_fh join "\t", $chr_id, $pos, @{ $par2_snps{$pos}{snp_del} }[ 0, 1 ], $par2_name, "NA";
        }
    }
}

close $master_fh;
exit;

###############
# subroutines #
###############

sub snp_hash_builder {
    my $par_snps_file = shift;
    open my $par_snps_fh, "<", $par_snps_file;
    my $header = <$par_snps_fh>;    # do I care about this?
    my %par_snps;
    while (<$par_snps_fh>) {
        chomp;
        my @par_snp = split /,/;
        my @par_pos = split /\./, $par_snp[1];
        if ( defined $par_pos[1] ) {
            $par_snps{$par_pos[0]}{$par_pos[1]} = [ @par_snp[ 2, 8 ], $par_pos[1] ];
            $par_snps{$par_pos[0]}{insert}      = $par_pos[1];
        }
        else {
            $par_snps{ $par_pos[0] }{snp_del} = [ @par_snp[ 2, 8 ] ];
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
