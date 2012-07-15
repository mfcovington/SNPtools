#!/usr/bin/perl
# parental_genome_builder_pt1_v2.pl
# Mike Covington
# created: 2012-01-03
# 2012-01-15: adapted for use with RNAseq data only (SNPs + indels)
#
# Description: 
#
use strict; use warnings;
use Getopt::Long;

my $chr_num = -1;
my $m82_rna_snp_dir = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/M82_SNPS";
my $pen_rna_snp_dir = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/PENN_SNPS";
my $cov_dir = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/coverage_out";
my $output_dir = "./";

GetOptions (
	"chr=i" => \$chr_num,    # numeric
	"msnp=s"   => \$m82_rna_snp_dir,      # string
	"psnp=s"   => \$pen_rna_snp_dir,      # string
	"cov=s"   => \$cov_dir,      # string
	"out=s"	=>	\$output_dir
);

die "USAGE: parental_genome_builder_pt1.pl --chr <> --msnp <> --psnp <> --cov <> --out <>\n" unless $chr_num >= 0;

my $chr_offset = $chr_num + 1;
my $chr_format;
if ($chr_num < 10) {
	$chr_format = "0$chr_num";
} else {
	$chr_format = $chr_num;
}


#OPEN SNP FILES
my $m82_rna_snps = $m82_rna_snp_dir . "/01.2.SNP_table.bwa_tophat_m82.sorted.dupl_rm.realigned.sorted.$chr_offset.$chr_offset.nogap.gap.FILTERED.csv";
my $pen_rna_snps = $pen_rna_snp_dir . "/01.2.SNP_table.bwa_tophat_penn.sorted.dupl_rm.realigned.sorted.$chr_offset.$chr_offset.nogap.gap.FILTERED.csv";
open (M82_RNA_SNPS, $m82_rna_snps) or die "Can't open $m82_rna_snps";
open (PEN_RNA_SNPS, $pen_rna_snps) or die "Can't open $pen_rna_snps";

#OPEN COVERAGE FILES
my $m82_rna_cov_file = $cov_dir . "/m82.sl2.40ch$chr_format.coverage.col";
my $pen_rna_cov_file = $cov_dir . "/penn.sl2.40ch$chr_format.coverage.col";
open (M82_RNA_COV, $m82_rna_cov_file) or die "Can't open $m82_rna_cov_file";
open (PEN_RNA_COV, $pen_rna_cov_file) or die "Can't open $pen_rna_cov_file";

#OPEN MASTER SNP FILE
`mkdir $output_dir`;
my $master_snp_file = ">$output_dir/master_snp_list_rnaseq.chr$chr_format";
open (MASTER, $master_snp_file) or die "Can't open $master_snp_file";
print MASTER join("\t", "chr", "pos", "ref_base", "snp_base", "genotype", "insert_position"), "\n"; # insert_position will be decimal portion of position for insertions

#discard headers
<M82_RNA_SNPS>;
<PEN_RNA_SNPS>;

#read in first line of each SNP file
my (@m82_rna_snp, @pen_rna_snp, @m82_rna_pos, @pen_rna_pos, $line);

$line = <M82_RNA_SNPS>;
@m82_rna_snp = split(/,/, $line); # 0=chr, 1=pos, 2=ref_base , 8=snp_base
@m82_rna_pos = split(/\./, $m82_rna_snp[1]); 

$line = <PEN_RNA_SNPS>;
@pen_rna_snp = split(/,/, $line); # 0=chr, 1=pos, 2=ref_base , 8=snp_base
@pen_rna_pos = split(/\./, $pen_rna_snp[1]); 
		


#read in cov files until at correct chromosome
my (@m82_rna_cov, @pen_rna_cov, @m82_rescan_cov, @pen_rescan_cov);
while ($line = <M82_RNA_COV>) {
	@m82_rna_cov = split(/,/, $line); # 0=chr, 1=pos, 2=coverage
	@pen_rna_cov = split(/,/, <PEN_RNA_COV>);
	last if $m82_rna_cov[0] eq "SL2.40ch" . $chr_format;
}

until ($m82_rna_cov[0] ne "SL2.40ch" . $chr_format) { 
	my $pos = $m82_rna_cov[1];

	IFM82: if ($m82_rna_pos[0] == $pos) {
		my $is_insert = 0;
# 		my $insert_pos;
		my $insert_pos = $m82_rna_pos[0]; ###TEST
		my $decimal = "NA";
		my $sufficient_cov = 0;
		chomp $pen_rna_cov[2];
		
		if ($m82_rna_pos[1]) {
			$is_insert = 1;
# 			$insert_pos = $m82_rna_pos[0];
			$decimal = $m82_rna_pos[1];
		}

		if ($pen_rna_cov[2] >=4) {
			print MASTER join("\t", $m82_rna_snp[0], $m82_rna_pos[0], @m82_rna_snp[2,8], "M82", $decimal), "\n";
			$sufficient_cov = 1;
		}
		
		while ($line = <M82_RNA_SNPS>) {
			@m82_rna_snp = split(/,/, $line); # 0=chr, 1=pos, 2=ref_base , 8=snp_base
			@m82_rna_pos = split(/\./, $m82_rna_snp[1]); 
			
			goto IFM82 if ($is_insert == 0 && $insert_pos == $m82_rna_pos[0]); ###TEST
			
			if ($is_insert == 1 && $insert_pos == $m82_rna_pos[0]) {
				next unless $sufficient_cov == 1;
				$decimal = $m82_rna_pos[1];
				print MASTER join("\t", $m82_rna_snp[0], $m82_rna_pos[0], @m82_rna_snp[2,8], "M82", $decimal), "\n";
			}else{
				last;
			}
		}		
	}


	IFPEN: if ($pen_rna_pos[0] == $pos) {
		my $is_insert = 0;
# 		my $insert_pos;
		my $insert_pos = $pen_rna_pos[0]; ###TEST
		my $decimal = "NA";
		my $sufficient_cov = 0;
		chomp $m82_rna_cov[2];
		
		if ($pen_rna_pos[1]) {
			$is_insert = 1;
# 			$insert_pos = $pen_rna_pos[0];
			$decimal = $pen_rna_pos[1];
		}

		if ($m82_rna_cov[2] >=4) {
			print MASTER join("\t", $pen_rna_snp[0], $pen_rna_pos[0], @pen_rna_snp[2,8], "PEN", $decimal), "\n";
			$sufficient_cov = 1;
		}
		
		while ($line = <PEN_RNA_SNPS>) {
			@pen_rna_snp = split(/,/, $line); # 0=chr, 1=pos, 2=ref_base , 8=snp_base
			@pen_rna_pos = split(/\./, $pen_rna_snp[1]); 
		
			goto IFPEN if ($is_insert == 0 && $insert_pos == $pen_rna_pos[0]); ###TEST
			
			if ($is_insert == 1 && $insert_pos == $pen_rna_pos[0]) {
				next unless $sufficient_cov == 1;
				$decimal = $pen_rna_pos[1];
				print MASTER join("\t", $pen_rna_snp[0], $pen_rna_pos[0], @pen_rna_snp[2,8], "PEN", $decimal), "\n";
			}else{
				last;
			}
		}		
	}


	last if eof(M82_RNA_COV);
	@m82_rna_cov = split(/,/, <M82_RNA_COV>); # 0=chr, 1=pos, 2=coverage
	@pen_rna_cov = split(/,/, <PEN_RNA_COV>);
}
close (M82_RNA_SNPS);
close (PEN_RNA_SNPS);
close (M82_RNA_COV);
close (PEN_RNA_COV);
close (MASTER);
exit;
## incorporate actual writing of parental genomes in until loop?  NO, I think having a separate script would be better.  Easier to filter out SNPs in common between m82 and pen.  also easier to filter out disagreements between RNAseq and RESCAN
