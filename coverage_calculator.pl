#!/usr/bin/perl
# coverage_calculator.v2.pl
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

# die "USAGE: coverage_calculator.pl </PATH/TO/file.bam> </PATH/TO/DESTINATION/DIRECTORY> <col|whit> <m82|penn>\n" unless (@ARGV == 4);
my ($ref_fa, $bam_file, $out_dir, $machine, $col, $whit, $m82, $pen, $species);
my $usage = "\nUSAGE: coverage_calculator.v2.pl --ref </PATH/TO/reference.fa> --bam </PATH/TO/file.bam> --out </PATH/TO/DESTINATION/DIRECTORY/> --col OR --whit --m82 OR --penn\n\n";

my $options = GetOptions (
	"ref=s"	=>	\$ref_fa,
	"bam=s"	=>	\$bam_file,
	"out=s" =>	\$out_dir,
	"col"	=> \$col,
	"whit"	=> \$whit,
	"m82"	=> \$m82,
	"pen"	=> \$pen
);

die $usage if defined $col && defined $whit;
die $usage if defined $m82 && defined $pen;

if ($col) {
	$machine = "col";
}elsif ($whit) {
	$machine = "whit";
}else{
	die "Something is wrong...\n";
}

if ($m82) {
	$species = "M82";
}elsif ($pen) {
	$species = "PEN";
}else{
	die "Something is wrong...\n";
}



my $sam = Bio::DB::Sam->new(
	-bam  => $bam_file,
	-fasta => $ref_fa
); #my $sam

my @chromosomes = $sam->seq_ids;
open (LOG, ">$out_dir/$species.$machine.log") or die "Cannot open >$out_dir/$species.$machine.log";
foreach my $chr (@chromosomes[7..12]) {
	$chr =~ tr/sl/SL/; # caps in chr names aren't consisent
	open (COVERAGE_OUT,	 ">$out_dir/$species.$chr.coverage.$machine") or die "Cannot open >$out_dir/$species.$chr.coverage.$machine";
	my $cov_pos = 1;
	my $chr_length = $sam->length($chr);
	print LOG "Getting $species coverage data for $chr (length = $chr_length)\n";

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
close LOG;
exit;

