#!/usr/bin/perl
# parental_genome_builder_pt1_v2.pl
# Mike Covington
# created: 2012-01-03
#
# Description: 
#
# Change Log:
#	2012-02-23: Added wildcard w/ glob to and altered structure of SNP files to account for difference caused by changing upstream scripts
#	2012-02-23: Changed capitalization of coverage files to account for difference caused by changing upstream scripts
#	2012-02-23: Changed `mkdir $output_dir` to `mkdir -p $output_dir`
#
use strict; use warnings;
use Getopt::Long;

my $chr_num = -1;
my $m82_rna_snp_dir = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/M82_SNPS";
my $pen_rna_snp_dir = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/PENN_SNPS";
my $m82_rescan_snps = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/M82_RESCAN/M82vsHz_repeat_frequency_filtered_uniqueSNP.txt";
my $pen_rescan_snps = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/PENN_RESCAN/PENvsHz_repeat_frequency_filtered_uniqueSNP.txt";
my $cov_dir = "/Volumes/Runner_3B/Mike_temp_SolRAID/Mike_SNPs/coverage_out";
my $output_dir = "./";

GetOptions (
	"chr=i" => \$chr_num,    # numeric
	"msnp=s"   => \$m82_rna_snp_dir,      # string
	"psnp=s"   => \$pen_rna_snp_dir,      # string
	"mrescan=s"   => \$m82_rescan_snps,      # string
	"prescan=s"   => \$pen_rescan_snps,      # string
	"cov=s"   => \$cov_dir,      # string
	"out=s"	=>	\$output_dir
);

die "USAGE: parental_genome_builder_pt1_v2.pl --chr <> --msnp <> --psnp <> --mrescan <> --prescan <> --cov <> --out <>\n" unless $chr_num >= 0;

my $chr_offset = $chr_num + 1;
my $chr_format;
if ($chr_num < 10) {
	$chr_format = "0$chr_num";
} else {
	$chr_format = $chr_num;
}


#OPEN SNP FILES
my $m82_rna_snps = glob($m82_rna_snp_dir . "/01.2.SNP_table.bwa_tophat_*.sorted.dupl_rm.$chr_offset.$chr_offset.nogap.gap.FILTERED.csv");
my $pen_rna_snps = glob($pen_rna_snp_dir . "/01.2.SNP_table.bwa_tophat_*.sorted.dupl_rm.$chr_offset.$chr_offset.nogap.gap.FILTERED.csv");
open (M82_RNA_SNPS, $m82_rna_snps) or die "Can't open $m82_rna_snps";
open (PEN_RNA_SNPS, $pen_rna_snps) or die "Can't open $pen_rna_snps";
open (M82_RE_SNPS, $m82_rescan_snps) or die "Can't open $m82_rescan_snps";
open (PEN_RE_SNPS, $pen_rescan_snps) or die "Can't open $pen_rescan_snps";

#OPEN COVERAGE FILES
my $m82_rna_cov_file = $cov_dir . "/M82.SL2.40ch$chr_format.coverage.col";
my $pen_rna_cov_file = $cov_dir . "/PEN.SL2.40ch$chr_format.coverage.col";
my $m82_rescan_cov_file = $cov_dir . "/m82_rescan_no_repeats.SL2.40ch$chr_format.coverage.whit";
my $pen_rescan_cov_file = $cov_dir . "/penn_rescan_no_repeats.SL2.40ch$chr_format.coverage.whit";
open (M82_RNA_COV, $m82_rna_cov_file) or die "Can't open $m82_rna_cov_file";
open (PEN_RNA_COV, $pen_rna_cov_file) or die "Can't open $pen_rna_cov_file";
open (M82_RE_COV, $m82_rescan_cov_file) or die "Can't open $m82_rescan_cov_file";
open (PEN_RE_COV, $pen_rescan_cov_file) or die "Can't open $pen_rescan_cov_file";

#OPEN MASTER SNP FILE
`mkdir -p $output_dir`;
my $master_snp_file = ">$output_dir/master_snp_list.chr$chr_format";
open (MASTER, $master_snp_file) or die "Can't open $master_snp_file";
print MASTER join("\t", "chr", "pos", "ref_base", "snp_base", "genotype", "RESCANvsRNAseq"), "\n";

#discard headers
<M82_RNA_SNPS>;
<PEN_RNA_SNPS>;

#read in first line of each SNP file; for RNAseq, need to skip if insertion/decimal in position (at least for now); for RESCAN SNPs, need to skip to correct chromosome and skip SNPs where ref_base eq "N"
my (@m82_rna_snp, @pen_rna_snp, @m82_rescan_snp, @pen_rescan_snp);
while (my $line = <M82_RNA_SNPS>) {
	@m82_rna_snp = split(/,/, $line); # 0=chr, 1=pos, 2=ref_base , 8=snp_base
	last unless (split(/\./, $m82_rna_snp[1]))[1]; #keep unless there is a decimal in position (and is, therefore, an insertion)
}
while (my $line = <PEN_RNA_SNPS>) {
	@pen_rna_snp = split(/,/, $line); # 0=chr, 1=pos, 2=ref_base , 8=snp_base
	last unless (split(/\./, $pen_rna_snp[1]))[1]; #keep unless there is a decimal in position (and is, therefore, an insertion)
}
while (my $line = <M82_RE_SNPS>) {
	@m82_rescan_snp = split(/\t/, $line); # 0=chr, 1=pos, 2=ref_base , 3=snp_base
	last if $m82_rescan_snp[0] eq "SL2.40ch" . $chr_format;
}
while ($m82_rescan_snp[2] eq "N") {
	@m82_rescan_snp = split(/\t/, <M82_RE_SNPS>);
}
while (my $line = <PEN_RE_SNPS>) {
	@pen_rescan_snp = split(/\t/, $line);
	last if $pen_rescan_snp[0] eq "SL2.40ch" . $chr_format;
}
while ($pen_rescan_snp[2] eq "N") {
	@pen_rescan_snp = split(/\t/, <PEN_RE_SNPS>);
}

#read in cov files until at correct chromosome
my (@m82_rna_cov, @pen_rna_cov, @m82_rescan_cov, @pen_rescan_cov);
while (my $line = <M82_RNA_COV>) {
	@m82_rna_cov = split(/,/, $line); # 0=chr, 1=pos, 2=coverage
	@pen_rna_cov = split(/,/, <PEN_RNA_COV>);
	@m82_rescan_cov = split(/,/, <M82_RE_COV>);
	@pen_rescan_cov = split(/,/, <PEN_RE_COV>);
	last if $m82_rna_cov[0] eq "SL2.40ch" . $chr_format;
}

#advance through coverage lines one at a time.  for each position, check against SNP positions.  if there is a match, confirm sufficient coverage in corresponding partner.  if sufficient coverage, confirm SNP and not indel.  If SNP, write to master SNP file (include chr, pos, ref_base, snp_base, genotype and RESCANvsRNAseq).  Advance SNP file that was used to next SNP position.  Make sure to allow for multiple SNPs at same position. (Can compare later)
until ($m82_rna_cov[0] ne "SL2.40ch" . $chr_format) { ##what about at EOF? Do I also need to check for defined/exists?
	my $pos = $m82_rna_cov[1];

	if ($m82_rna_snp[1] == $pos) {
		chomp $pen_rna_cov[2];
		if ($pen_rna_cov[2] >=4 && $m82_rna_snp[8] ne "del") {
			print MASTER join("\t", @m82_rna_snp[0,1,2,8], "M82", "RNAseq"), "\n";
		}
		while (my $line = <M82_RNA_SNPS>) {
			@m82_rna_snp = split(/,/, $line); # 0=chr, 1=pos, 2=ref_base , 8=snp_base
			last unless (split(/\./, $m82_rna_snp[1]))[1]; #keep unless there is a decimal in position (and is, therefore, an insertion)
		}
	}

	if ($pen_rna_snp[1] == $pos) {
		chomp $m82_rna_cov[2];
		if ($m82_rna_cov[2] >=4 && $pen_rna_snp[8] ne "del") {
			print MASTER join("\t", @pen_rna_snp[0,1,2,8], "PEN", "RNAseq"), "\n";
		}
		while (my $line = <PEN_RNA_SNPS>) {
			@pen_rna_snp = split(/,/, $line); # 0=chr, 1=pos, 2=ref_base , 8=snp_base
			last unless (split(/\./, $pen_rna_snp[1]))[1]; #keep unless there is a decimal in position (and is, therefore, an insertion)
		}
	}

	if ($m82_rescan_snp[1] == $pos) {
		last unless $m82_rescan_snp[0] eq "SL2.40ch" . $chr_format;
		chomp $pen_rescan_cov[2];
		if ($pen_rescan_cov[2] >=4) {
			chomp @m82_rescan_snp;
			print MASTER join("\t", @m82_rescan_snp[0..3], "M82", "RESCAN"), "\n";
		}
		while (my $line = <M82_RE_SNPS>) {
			@m82_rescan_snp = split(/\t/, $line); # 0=chr, 1=pos, 2=ref_base , 3=snp_base
			last unless $m82_rescan_snp[2] eq "N"; #keep unless ref_base eq "N"
		}
	}

	if ($pen_rescan_snp[1] == $pos) {
		last unless $pen_rescan_snp[0] eq "SL2.40ch" . $chr_format;
		chomp $m82_rescan_cov[2];
		if ($m82_rescan_cov[2] >=4) {
			chomp @pen_rescan_snp;
			print MASTER join("\t", @pen_rescan_snp[0..3], "PEN", "RESCAN"), "\n";
		}
		while (my $line = <PEN_RE_SNPS>) {
			@pen_rescan_snp = split(/\t/, $line); # 0=chr, 1=pos, 2=ref_base , 3=snp_base
			last unless $pen_rescan_snp[2] eq "N"; #keep unless ref_base eq "N"
		}
	}


	last if eof(M82_RNA_COV);
	@m82_rna_cov = split(/,/, <M82_RNA_COV>); # 0=chr, 1=pos, 2=coverage
	@pen_rna_cov = split(/,/, <PEN_RNA_COV>);
	@m82_rescan_cov = split(/,/, <M82_RE_COV>);
	@pen_rescan_cov = split(/,/, <PEN_RE_COV>);
}
close (M82_RNA_SNPS);
close (PEN_RNA_SNPS);
close (M82_RE_SNPS);
close (PEN_RE_SNPS);
close (M82_RNA_COV);
close (PEN_RNA_COV);
close (M82_RE_COV);
close (PEN_RE_COV);
close (MASTER);
exit;
## incorporate actual writing of parental genomes in until loop?  NO, I think having a separate script would be better.  Easier to filter out SNPs in common between m82 and pen.  also easier to filter out disagreements between RNAseq and RESCAN
