#!/usr/bin/perl
# parental_genome_builder_rnaseq_only_pt1_v2.solo.pl
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
#	2012-04-19: Changed to work with single genotype instead of comparing two.
#
use strict; use warnings;
use Getopt::Long;

my $chr_num = -1;
my $snp_dir = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/M82_SNPS";
my $output_dir = "./";
my $offset;
my $id = "bwa_tophat_*.sorted.dupl_rm";

GetOptions (
	"chr=i"		=>	\$chr_num,
	"snp=s"		=>	\$snp_dir,
	"out=s"		=>	\$output_dir,
	"offset"	=>	\$offset,
	"id=s"		=>	\$id,
);

die "USAGE: parental_genome_builder_rnaseq_only_pt1_v2.solo.pl --chr <> --snp <> --out <> --offset (flag, used if # of snp files is offset) --id\n" unless $chr_num >= 0;

my $chr_offset;
if ($offset) {
	$chr_offset = $chr_num + 1;
}else{
	$chr_offset = $chr_num;
}

my $chr_format;
if ($chr_num < 10) {
	$chr_format = "0$chr_num";
} else {
	$chr_format = $chr_num;
}

#OPEN SNP FILES
my $snps = glob($snp_dir . "/01.2.SNP_table.$id.$chr_offset.$chr_offset.nogap.gap.FILTERED.csv");
open (SNPS, $snps) or die "Can't open $snps";

#OPEN MASTER SNP FILE
`mkdir -p $output_dir`;
my $master_snp_file = ">$output_dir/master_snp_list_solo.chr$chr_format";
open (MASTER, $master_snp_file) or die "Can't open $master_snp_file";
print MASTER join("\t", "chr", "pos", "ref_base", "snp_base", "insert_position"), "\n"; # insert_position will be decimal portion of position for insertions

#discard headers
<SNPS>;

while (my $line = <SNPS>) {
	my @snp = split(/,/, $line); # 0=chr, 1=pos, 2=ref_base , 8=snp_base
	my @pos = split(/\./, $snp[1]); # 0=pos, 1=insert_decimal
	
	my $decimal = "NA";
	$decimal = $pos[1] if $pos[1];
	
	print MASTER join("\t", $snp[0], $pos[0], @snp[2,8], $decimal), "\n";	
}

close (SNPS);
close (MASTER);
exit;