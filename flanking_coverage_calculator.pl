#!/usr/bin/perl
# flanking_coverage_calculator_2.pl
# Mike Covington
# created: 2011-12-11
#
# Description: 
#
use strict; use warnings;
use File::Basename;


die "\n\tUSAGE: flanking_coverage_calculator_2.pl <SNP.CSV> <COL.COV> <WHIT.COV>\n\n" unless (@ARGV == 3);
my ($snp_in, $COL_in, $WHIT_in) = @ARGV;
#open files
open (SNP_IN, $snp_in) or die "Can't open $snp_in";
open (COL_COV_IN, $COL_in) or die "Can't open $COL_in";
open (WHIT_COV_IN, $WHIT_in) or die "Can't open $WHIT_in";
my ($filename, $directories, $suffix) = fileparse($snp_in, ".csv");
my $snp_out = $directories . $filename . ".nogap.gap.csv";
open (SNP_OUT, ">$snp_out") or die "Can't open $snp_out";

#read in first 17 lines
my (@COL_cov_chr, @COL_cov_pos, @COL_cov, @WHIT_cov_chr, @WHIT_cov_pos, @WHIT_cov, $COL_line, $WHIT_line);
for (my $count = 0; $count < 17; $count++) {
	$COL_line = <COL_COV_IN>;
	$WHIT_line = <WHIT_COV_IN>;
	chomp ($COL_line, $WHIT_line);
 	($COL_cov_chr[$count], $COL_cov_pos[$count], $COL_cov[$count]) = split(/,/, $COL_line);
 	($WHIT_cov_chr[$count], $WHIT_cov_pos[$count], $WHIT_cov[$count]) = split(/,/, $WHIT_line);
}

my $header = <SNP_IN>;
chomp $header;
print SNP_OUT join(",",$header, "nogap_cov_pos-8","nogap_cov_pos","nogap_cov_pos+8", "gap_cov_pos-8","gap_cov_pos","gap_cov_pos+8"), "\n";
while (my $snp_line = <SNP_IN>) {
	chomp $snp_line;
	my ($snp_chr, $snp_pos_unsplit, $snp_remainder) = split(/,/, $snp_line, 3); #read in CSV line
	my ($snp_pos, $snp_pos_index) = split(/\./, $snp_pos_unsplit);
	while ($snp_pos > $COL_cov_pos[8]) {
		shift @COL_cov_chr;
		shift @COL_cov_pos;
		shift @COL_cov;
		shift @WHIT_cov_chr;
		shift @WHIT_cov_pos;
		shift @WHIT_cov;
		$COL_line = <COL_COV_IN>;
		$WHIT_line = <WHIT_COV_IN>;
		chomp ($COL_line, $WHIT_line);
 		my @COL_line_elements = split(/,/, $COL_line);
 		my @WHIT_line_elements = split(/,/, $WHIT_line);
		push @COL_cov_chr, $COL_line_elements[0];
		push @COL_cov_pos, $COL_line_elements[1];
		push @COL_cov, $COL_line_elements[2];
		push @WHIT_cov_chr, $WHIT_line_elements[0];
		push @WHIT_cov_pos, $WHIT_line_elements[1];
		push @WHIT_cov, $WHIT_line_elements[2];
	}
	if ($snp_chr eq $COL_cov_chr[8] && $snp_pos == $COL_cov_pos[8]) {
		print SNP_OUT join (",", $snp_chr, $snp_pos_unsplit, $snp_remainder, $COL_cov[0], $COL_cov[8], $COL_cov[16], $WHIT_cov[0], $WHIT_cov[8], $WHIT_cov[16]), "\n";
	}else{
		die "Something strange going on here...";
	}
}

close SNP_IN;
close COL_COV_IN;
close WHIT_COV_IN;
close SNP_OUT;


