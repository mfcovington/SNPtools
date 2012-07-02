#!/usr/bin/perl 

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::DB::Sam;
use List::Util qw(sum);

die "USAGE: FILE_NAME.pl </PATH/TO/file.bam> </PATH/TO/file.csv>\n" unless (@ARGV >= 1);
my ($bam_file, @files) = @ARGV;

my $sam = Bio::DB::Sam->new(
	-bam  => $bam_file,
	-fasta => "/Volumes/SolexaRAID/Solexa_runs_Data/00.Downloaded_references_and_others/S_lycopersicum_chromosomes.2.40.fa"
); #my $sam

open COVERAGE_OUT,	 ">coverage_files-col";
my @chromosomes = $sam->seq_ids;
foreach my $chr (@chromosomes) {
	$chr =~ tr/sl/SL/; # caps in chr names aren't consisent
	my $cov_pos = 1;
	my $chr_length = $sam->length($chr);
	print $chr_length, "\n";
# 	while ($cov_pos<$chr_length) {
# 		$cov_pos++;
# 		my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr,-start=>$cov_pos,-end=>$cov_pos+999999);
# 		my ($cov_val) = $coverage->coverage;
# #		my @coveragedata = $coverage->coverage;
# #		foreach my $cov_val (@coveragedata) {
# 		print COVERAGE_OUT "$chr,$cov_pos,$cov_val\n";
# 		last if $cov_pos > 999999;
# #		}
# 	}


# 	while ($cov_pos<$chr_length) {
# 		if ($cov_pos+99999<$chr_length) {
# 			$cov_pos++;
# 			my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr,-start=>$cov_pos,-end=>$cov_pos+99999);
# 			my ($cov_val) = $coverage->coverage;
# 			print COVERAGE_OUT "$chr,$cov_pos,$cov_val\n";
# 		} else {
# 			$cov_pos++;
# 			my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr,-start=>$cov_pos,-end=>$chr_length);
# 			my ($cov_val) = $coverage->coverage;
# 			print COVERAGE_OUT "$chr,$cov_pos,$cov_val\n";
# 		}
# 	}

	my $chunk_size = 999999;  #99999 = 2 min and 16G still free; 999 is much longer (only 60% after 15 min) and worse on memory; 999999 = 2min and 23G left; 9999999 = 2 min 18G free
	while ($cov_pos<=$chr_length) {
		if ($cov_pos+$chunk_size<$chr_length) {
			my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr,-start=>$cov_pos,-end=>$cov_pos+$chunk_size);
			my @coveragedata = $coverage->coverage;
			foreach my $cov_val (@coveragedata) {
				print COVERAGE_OUT "$chr,$cov_pos,$cov_val\n";
				$cov_pos++;
			}
		} else {
			my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr,-start=>$cov_pos,-end=>$chr_length);
			my @coveragedata = $coverage->coverage;
			foreach my $cov_val (@coveragedata) {
				print COVERAGE_OUT "$chr,$cov_pos,$cov_val\n";
				$cov_pos++;
			}
		}
	}


# foreach my $cov_val (@coveragedata) {
# 	print COVERAGE_OUT "$chr,$cov_pos,$cov_val\n";
# 	$cov_pos++;
# }


}
#print COVERAGE_OUT map { "$_ \n" } @coveragedata;
close COVERAGE_OUT;
exit;

