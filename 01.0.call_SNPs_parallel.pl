#!/usr/bin/perl
# 01.0.call_SNPs_parallel.pl
# Mike Covington
# created: 2012-01-07
#
# Description: Adapted from Pepe's script (a second time, since original adaptation is lost in backup purgatory).
# - Moved all parameters to be set to one location
# - Enabled customization of output directory
# - Changed to specify --fa --bam and --out at CLI
# - Calls 01.1.SNP_calling_homos.pl (formerly named: 01.1.SNP_calling_homo_based_on_pepes_cigar_update.pl)
#
# 2012-04-02: cleaned up, expanded getoptions and usage statement
#
# To-do:
#	Re-work the way the thresholds work.  Want it to be >= instead of >
#	Change naming of output files
#	Change from system calls?
#
use strict; use warnings;
use POSIX qw(ceil floor);
use File::Basename;
use Getopt::Long;

#DEFAULT PARAMETERS
my $chrm_start = 0;
my $number_of_chrms = 13;
my $number_of_threads = 13;
my $threshold_number_of_reads = 4; #threshold of coverage to call snps (only look at sites with this number or more reads)
my $threshold_fraction_of_reads_matching_ref = 0.66; ##0.66 is good because when there are 3 bases only, all three have to be snps to be called a snp(sites with this fraction of bases matching reference bases or less will be counted)
my $threshold_fraction_of_reads_matching_ref_for_indels = 0.33; ##0.33 is good for indels (sites with this fraction of bases matching reference bases or less will be counted)
my $ref_fasta = "/Volumes/SolexaRAID/Solexa_runs_Data/00.Downloaded_references_and_others/S_lycopersicum_chromosomes.2.40.fa"; #reference fasta file
my ($bam_file, $out_dir, $help);

my $options = GetOptions (
	"fa=s"			=>	\$ref_fasta,
	"bam=s"			=>	\$bam_file,
	"out=s"			=>	\$out_dir,
	"start=i"		=>	\$chrm_start,
	"total=i"		=>	\$number_of_chrms,
	"threads=i"		=>	\$number_of_threads,
	"depth_cut=i"	=>	\$threshold_number_of_reads,
	"snp_cut=f"		=>	\$threshold_fraction_of_reads_matching_ref,
	"indel_cut=f"	=>	\$threshold_fraction_of_reads_matching_ref_for_indels,
	"help"			=>	\$help,
);


my $usage = "
	USAGE: 01.0.call_SNPs_parallel.pl
		--fa		</path/to/reference.fa>
		--bam		</path/to/align.bam>
		--out		</path/to/output/directory/>
		--start		<chromosome to start analyzing from>
		--total		<total # of chromosomes to analyze>
		--threads	<number of cores to use at a time>
		--depth_cut	<coverage threshold for calling SNPs>
		--snp_cut	<>
		--indel_cut	<>
		--help		**print this help**

";

die $usage if (
	$help
	or not $bam_file
	or not $out_dir
);


#Parallelization
my $step = ceil(($number_of_chrms-$chrm_start)/$number_of_threads);

for (my $i=$chrm_start;$i<$number_of_chrms;$i += $step){
	if (($i+$step-1)>=$number_of_chrms){
		my $end = ($number_of_chrms-1);
		my $basename = basename($bam_file,(".bam"));
		my $out = $out_dir . "/01.2.SNP_table.$basename." . ($i+1) . "." . ($end+1);
		system("perl 01.1.SNP_calling_homos.pl -fasta_ref $ref_fasta -chrm_start $i -chrm_end $end -n_reads $threshold_number_of_reads -ref_freq $threshold_fraction_of_reads_matching_ref -indel_freq $threshold_fraction_of_reads_matching_ref_for_indels -o $out -bam_file $bam_file &");
	}else{
		my $end = ($i+$step-1);
		my $basename = basename($bam_file,(".bam"));
		my $out = $out_dir . "/01.2.SNP_table.$basename." . ($i+1) . "." . ($end+1);
		system("perl 01.1.SNP_calling_homos.pl -fasta_ref $ref_fasta -chrm_start $i -chrm_end $end -n_reads $threshold_number_of_reads -ref_freq $threshold_fraction_of_reads_matching_ref -indel_freq $threshold_fraction_of_reads_matching_ref_for_indels -o $out -bam_file $bam_file &");
	}
}

exit;


