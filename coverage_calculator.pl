#!/usr/bin/perl
# coverage_calculator.pl
# Mike Covington
# created: 2011-12-05
#
# Description: 
#
use strict; use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::DB::Sam;
use List::Util qw(sum);

###TODO: INCORPORATE GETOPT?

die "USAGE: coverage_calculator.pl </PATH/TO/file.bam> </PATH/TO/DESTINATION/DIRECTORY> <col|whit> <m82|penn>\n" unless (@ARGV == 4);
my ($bam_file, $out_dir, $machine, $species) = @ARGV;

my $sam = Bio::DB::Sam->new(
	-bam  => $bam_file,
	-fasta => "/Volumes/SolexaRAID/Solexa_runs_Data/00.Downloaded_references_and_others/S_lycopersicum_chromosomes.2.40.fa"
); #my $sam

my @chromosomes = $sam->seq_ids;
foreach my $chr (@chromosomes) {
	open (COVERAGE_OUT,	 ">$out_dir/$species.$chr.coverage.$machine") or die "Cannot open >$out_dir/$species.$chr.coverage.$machine";
	$chr =~ tr/sl/SL/; # caps in chr names aren't consisent
	my $cov_pos = 1;
	my $chr_length = $sam->length($chr);
	print "Getting $species coverage data for $chr (length = $chr_length)\n";

	my $chunk_size = 999999;
	while ($cov_pos<=$chr_length) {
		if ($cov_pos+$chunk_size<$chr_length) {
			my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr,-start=>$cov_pos,-end=>$cov_pos+$chunk_size);
			my @coveragedata = $coverage->coverage;
			foreach my $cov_val (@coveragedata) {
				print COVERAGE_OUT "$chr,$cov_pos,$cov_val\n";
				$cov_pos++;
			} #foreach @coveragedata
		} else {
			my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr,-start=>$cov_pos,-end=>$chr_length);
			my @coveragedata = $coverage->coverage;
			foreach my $cov_val (@coveragedata) {
				print COVERAGE_OUT "$chr,$cov_pos,$cov_val\n";
				$cov_pos++;
			} #foreach @coveragedata
		} #ifelse
	} #while ($cov_pos<=$chr_length)
	close COVERAGE_OUT;
}
exit;

