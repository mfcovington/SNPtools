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

my $usage = <<USAGE_END;

USAGE:                                         #this usage statement doesn't reflect reality, yet
parental_genome_builder_pt1.rnaseq.pl
  --chr    <chromosome number>                 #need to change this to chr id eventually!!
  --snp1 <parent #1 snp file>                  #currently, everything is pointing to directories
  --snp2 <parent #2 snp file>
  --cov1 <parent #1 coverage.cov_nogap file>
  --cov2 <parent #2 coverage.cov_nogap file>
  --out  <output directory>
  --help

USAGE_END

#defaults
my $chr_num         = -1;
my $m82_rna_snp_dir = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/M82_SNPS";
my $pen_rna_snp_dir = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/PENN_SNPS";
my $cov_dir         = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/coverage_out";
my $output_dir      = "./";
my $help;

GetOptions(
    "chr=i"  => \$chr_num,
    "msnp=s" => \$m82_rna_snp_dir,
    "psnp=s" => \$pen_rna_snp_dir,
    "cov=s"  => \$cov_dir,
    "out=s"  => \$output_dir,
    "help"   => \$help,
);

die $usage if $help;
die $usage unless $chr_num >= 0;

my $chr_offset = $chr_num + 1;
my $chr_format;
if ( $chr_num < 10 ) {
    $chr_format = "0$chr_num";
}
else {
    $chr_format = $chr_num;
}



#OPEN SNP FILES
my $m82_rna_snps = glob($m82_rna_snp_dir . "/01.2.SNP_table.bwa_tophat_*.sorted.dupl_rm.$chr_offset.$chr_offset.nogap.gap.FILTERED.csv");
my $pen_rna_snps = glob($pen_rna_snp_dir . "/01.2.SNP_table.bwa_tophat_*.sorted.dupl_rm.$chr_offset.$chr_offset.nogap.gap.FILTERED.csv");
open $m82_rna_snps_fh, "<", $m82_rna_snps;
open $pen_rna_snps_fh, "<", $pen_rna_snps;

#OPEN COVERAGE FILES
my $m82_rna_cov_file = $cov_dir . "/M82.SL2.40ch$chr_format.coverage.col";
my $pen_rna_cov_file = $cov_dir . "/PEN.SL2.40ch$chr_format.coverage.col";
open $m82_rna_cov_fh, "<", $m82_rna_cov_file;
open $pen_rna_cov_fh, "<", $pen_rna_cov_file;

#OPEN MASTER SNP FILE
`mkdir -p $output_dir`;
my $master_snp_file = "$output_dir/master_snp_list_rnaseq.chr$chr_format";
open $master_fh, ">", $master_snp_file;
say $master_fh join "\t", "chr", "pos", "ref_base", "snp_base", "genotype", "insert_position";    # insert_position will be decimal portion of position for insertions

#discard headers
<$m82_rna_snps_fh>;
<$pen_rna_snps_fh>;

#read in first line of each SNP file
my (@m82_rna_snp, @pen_rna_snp, @m82_rna_pos, @pen_rna_pos, $line);

$line = <$m82_rna_snps_fh>;
@m82_rna_snp = split /,/, $line;    # 0=chr, 1=pos, 2=ref_base , 8=snp_base
@m82_rna_pos = split /\./, $m82_rna_snp[1];

$line = <$pen_rna_snps_fh>;
@pen_rna_snp = split /,/, $line;    # 0=chr, 1=pos, 2=ref_base , 8=snp_base
@pen_rna_pos = split /\./, $pen_rna_snp[1];


#read in cov files until at correct chromosome
my ( @m82_rna_cov, @pen_rna_cov );
while ( $line = <$m82_rna_cov_fh> ) {
    @m82_rna_cov = split /,/, $line;    # 0=chr, 1=pos, 2=coverage
    @pen_rna_cov = split /,/, <$pen_rna_cov_fh>;
    last if $m82_rna_cov[0] eq "SL2.40ch" . $chr_format;
}

until ( $m82_rna_cov[0] ne "SL2.40ch" . $chr_format ) {
    my $pos = $m82_rna_cov[1];    ###NEED TO CHANGE TO 

  IFM82: if ( $m82_rna_pos[0] == $pos ) {
        my $is_insert      = 0;
        my $insert_pos     = $m82_rna_pos[0];    ###TEST
        my $decimal        = "NA";
        my $sufficient_cov = 0;
        chomp $pen_rna_cov[2];

        if ( $m82_rna_pos[1] ) {
            $is_insert = 1;
            $decimal   = $m82_rna_pos[1];
        }

        if ( $pen_rna_cov[2] >=4 ) {
            say $master_fh join "\t", $m82_rna_snp[0], $m82_rna_pos[0], @m82_rna_snp[2,8], "M82", $decimal;
            $sufficient_cov = 1;
        }

        while ( $line = <$m82_rna_snps_fh> ) {
            @m82_rna_snp = split /,/, $line;    # 0=chr, 1=pos, 2=ref_base , 8=snp_base
            @m82_rna_pos = split /\./, $m82_rna_snp[1]; 

            goto IFM82 if ( $is_insert == 0 && $insert_pos == $m82_rna_pos[0] );    ###TEST

            if ( $is_insert == 1 && $insert_pos == $m82_rna_pos[0] ) {
                next unless $sufficient_cov == 1;
                $decimal = $m82_rna_pos[1];
                say $master_fh join "\t", $m82_rna_snp[0], $m82_rna_pos[0], @m82_rna_snp[2,8], "M82", $decimal;
            }
            else { last; }
        }
    }


  IFPEN: if ( $pen_rna_pos[0] == $pos ) {
        my $is_insert      = 0;
        my $insert_pos     = $pen_rna_pos[0];    ###TEST
        my $decimal        = "NA";
        my $sufficient_cov = 0;
        chomp $m82_rna_cov[2];

        if ( $pen_rna_pos[1] ) {
            $is_insert = 1;
            $decimal   = $pen_rna_pos[1];
        }

        if ( $m82_rna_cov[2] >=4 ) {
            say $master_fh join "\t", $pen_rna_snp[0], $pen_rna_pos[0], @pen_rna_snp[2,8], "PEN", $decimal;
            $sufficient_cov = 1;
        }

        while ( $line = <$pen_rna_snps_fh> ) {
            @pen_rna_snp = split /,/, $line;    # 0=chr, 1=pos, 2=ref_base , 8=snp_base
            @pen_rna_pos = split /\./, $pen_rna_snp[1];

            goto IFPEN if ( $is_insert == 0 && $insert_pos == $pen_rna_pos[0] );    ###TEST

            if ( $is_insert == 1 && $insert_pos == $pen_rna_pos[0] ) {
                next unless $sufficient_cov == 1;
                $decimal = $pen_rna_pos[1];
                say $master_fh join "\t", $pen_rna_snp[0], $pen_rna_pos[0], @pen_rna_snp[2,8], "PEN", $decimal;
            }
            else { last; }
        }
    }


    last if eof($m82_rna_cov_fh);
    @m82_rna_cov = split /,/, <$m82_rna_cov_fh>;    # 0=chr, 1=pos, 2=coverage
    @pen_rna_cov = split /,/, <$pen_rna_cov_fh>;
}
close ($m82_rna_snps_fh);
close ($pen_rna_snps_fh);
close ($m82_rna_cov_fh);
close ($pen_rna_cov_fh);
close ($master_fh);
exit;
## incorporate actual writing of parental genomes in until loop?  NO, I think having a separate script would be better.  Easier to filter out SNPs in common between m82 and pen.  also easier to filter out disagreements between RNAseq and RESCAN
