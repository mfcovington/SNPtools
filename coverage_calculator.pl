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
	my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr);
	my @coveragedata       = $coverage->coverage;
	my $cov_pos = 0;
	foreach my $cov_val (@coveragedata) {
		$cov_pos++;
		print COVERAGE_OUT "$chr,$cov_pos,$cov_val\n";
		last if $cov_pos > 10;
	}
	last;
}
#print COVERAGE_OUT map { "$_ \n" } @coveragedata;
close COVERAGE_OUT;
exit;

