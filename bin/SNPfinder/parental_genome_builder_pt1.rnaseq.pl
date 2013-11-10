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

# TODO: Make write_master_snp subroutine
# TODO: Get coverage data from DB
# TODO: (which requires:) Write coverage data to DB for par1/par2 combined sites of interest (in Coverage.pm)

my $usage = <<USAGE_END;

USAGE:
parental_genome_builder_pt1.rnaseq.pl
  --chr        chromosome ID
  --snp1       parent #1 snp file
  --snp2       parent #2 snp file
  --cov1       parent #1 coverage.cov_nogap file
  --cov2       parent #2 coverage.cov_nogap file
  --par1       parent #1 name
  --par2       parent #2 name
  --out_dir    output directory
  --min_cov    minimum coverage [4]
  --help       (displays this usage statement)

USAGE_END

# defaults/options
my $out_dir = "./";
my $par1    = "PAR1";
my $par2    = "PAR2";
my $min_cov = 4;
my ( $chr, $par1_snp_file, $par2_snp_file, $par1_cov_file, $par2_cov_file, $help );

GetOptions(
    "chr=s"     => \$chr,
    "snp1=s"    => \$par1_snp_file,
    "snp2=s"    => \$par2_snp_file,
    "cov1=s"    => \$par1_cov_file,
    "cov2=s"    => \$par2_cov_file,
    "name1=s"   => \$par1,
    "name2=s"   => \$par2,
    "out_dir=s" => \$out_dir,
    "min_cov=i" => \$min_cov,
    "help"      => \$help,
);

die $usage if $help;
die $usage
  unless defined $chr
      && defined $par1_snp_file
      && defined $par2_snp_file
      && defined $par1_cov_file
      && defined $par2_cov_file;

# snp/cov build hashes
my %par1_cov  = cov_hash_builder($par1_cov_file);
my %par2_cov  = cov_hash_builder($par2_cov_file);
my %snps;
snp_hash_builder( $par1, $par1_snp_file, \%snps);
snp_hash_builder( $par2, $par2_snp_file, \%snps);

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
make_path($out_dir);
my $master_snp_file = "$out_dir/master_snp_list.$par1.vs.$par2.$chr";
open my $master_fh, ">", $master_snp_file;
say $master_fh join "\t", "chr", "pos", "ref_base", "snp_base", "genotype", "insert_position";

for my $pos (@good_cov_pos) {

    next unless defined $snps{$pos};

    for my $par_id ( sort keys $snps{$pos} ) {

        # TODO: Check whether the next line should be '|| 0' or '// 0'
        my $insert_length = $snps{$pos}{$par_id}{insert} || 0;
        if ($insert_length) {
            for my $insert_idx ( 1 .. $insert_length ) {
                $insert_idx = "0$insert_idx" if $insert_idx < 10;
                next unless defined $snps{$pos}{$par_id}{$insert_idx}[2];    # there are occasional gaps due to insufficient coverage
                say $master_fh join "\t", $chr, $pos, @{ $snps{$pos}{$par_id}{$insert_idx} }[ 0, 1 ], $par_id, $insert_idx;
            }
        }
        else{
            say $master_fh join "\t", $chr, $pos, @{ $snps{$pos}{$par_id}{snp_del} }[ 0, 1 ], $par_id, "NA";
        }
    }
}

close $master_fh;
exit;

###############
# subroutines #
###############

sub snp_hash_builder {
    my ( $par_id, $par_snps_file, $snps_ref ) = @_;
    open my $par_snps_fh, "<", $par_snps_file;
    my $header = <$par_snps_fh>;    # do I care about this?
    # my %par_snps;
    while (<$par_snps_fh>) {
        chomp;
        my @par_snp = split /,/;
        my @par_pos = split /\./, $par_snp[1];
        if ( defined $par_pos[1] ) {
            $$snps_ref{$par_pos[0]}{$par_id}{$par_pos[1]} = [ @par_snp[ 2, 8 ], $par_pos[1] ];
            $$snps_ref{$par_pos[0]}{$par_id}{insert}      = $par_pos[1];
        }
        else {
            $$snps_ref{ $par_pos[0] }{$par_id}{snp_del} = [ @par_snp[ 2, 8 ] ];
        }
    }
    close $par_snps_fh;
}

sub cov_hash_builder {
    my $par_cov_file = shift;
    open my $par_cov_fh, "<", $par_cov_file;
    my %par_cov = map { chomp; @{ [ split /\t/ ] }[ 1 .. 2 ] } <$par_cov_fh>;
    close $par_cov_fh;
    return %par_cov;
}

}
