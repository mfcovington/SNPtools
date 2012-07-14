#!/usr/bin/env perl
# flanking_coverage_calculator.pl
# Mike Covington
# created: 2011-12-11
#
# Description: 
#
use strict;
use warnings;
use autodie;
use File::Basename;

# new version of coverage data is tab-delimited and only has data for regions with coverage

my $usage = <<USAGE_END;

USAGE:
flanking_coverage_calculator.pl
  --snp_file   </PATH/TO/snp_file.csv>
  --cov_prefix </PATH/TO/cov_file_prefix_only -- w/o .cov_(no)gap>
  --help

USAGE_END

my ( $snp_in, $cov_prefix, $help );
my $out_dir = "./";
my $options = GetOptions(
    "snp_file=s"   => \$snp_in,
    "cov_prefix=s" => \$cov_prefix,
    "help"         => \$help,
);

die $usage if $help;
die $usage unless defined $snp_in && defined $cov_prefix;

my ( $cov_gaps, $cov_nogaps ) = ( $cov_prefix, $cov_prefix );
$cov_gaps .= ".cov_gaps";
$cov_nogaps .= ".cov_nogaps";
my ($filename, $directories, $suffix) = fileparse($snp_in, ".csv");
my $snp_out = $directories . $filename . ".nogap.gap.csv";
#open files
open $snp_in_fh,     "<", $snp_in;
open $snp_out_fh,    ">", $snp_out;
open $cov_nogaps_fh, "<", $cov_nogaps;
open $cov_gaps_fh,   "<", $cov_gaps;

#read in first 17 lines
my (@COL_cov_chr, @COL_cov_pos, @COL_cov, @WHIT_cov_chr, @WHIT_cov_pos, @WHIT_cov, $COL_line, $WHIT_line);
for (my $count = 0; $count < 17; $count++) {
	$COL_line = <$cov_nogaps_fh>;
	$WHIT_line = <$cov_gaps_fh>;
	chomp ($COL_line, $WHIT_line);
 	($COL_cov_chr[$count], $COL_cov_pos[$count], $COL_cov[$count]) = split(/,/, $COL_line);
 	($WHIT_cov_chr[$count], $WHIT_cov_pos[$count], $WHIT_cov[$count]) = split(/,/, $WHIT_line);
}

my $header = <$snp_in_fh>;
chomp $header;
print $snp_out_fh join(",",$header, "nogap_cov_pos-8","nogap_cov_pos","nogap_cov_pos+8", "gap_cov_pos-8","gap_cov_pos","gap_cov_pos+8"), "\n";
while (my $snp_line = <$snp_in_fh>) {
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
		$COL_line = <$cov_nogaps_fh>;
		$WHIT_line = <$cov_gaps_fh>;
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
		print $snp_out_fh join (",", $snp_chr, $snp_pos_unsplit, $snp_remainder, $COL_cov[0], $COL_cov[8], $COL_cov[16], $WHIT_cov[0], $WHIT_cov[8], $WHIT_cov[16]), "\n";
	}else{
		die "Something strange going on here...";
	}
}

close $snp_in_fh;
close $cov_nogaps_fh;
close $cov_gaps_fh;
close $snp_out_fh;


