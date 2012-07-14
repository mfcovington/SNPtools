#!/usr/bin/perl
# FILE_NAME.pl
# Mike Covington
# created: 2011-12-07
#
# Description: 
#
#use strict; use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::DB::Sam;
use List::Util qw(sum);
use File::Basename;


die "\n\tUSAGE: XXXXXXX.pl </PATH/TO/file.bam> <col|whit> </PATH/TO/file.csv or /PATH/TO/file.nogap>\n\n" unless (@ARGV >= 3);
my ($bam_file, $machine, @files) = @ARGV;
die "$machine != 'col' or 'whit'" unless ($machine eq "col" || $machine eq "whit");

my $sam = Bio::DB::Sam->new(
	-bam  => $bam_file,
	-fasta => "/Volumes/SolexaRAID/Solexa_runs_Data/00.Downloaded_references_and_others/S_lycopersicum_chromosomes.2.40.fa"
); #my $sam

# my @chromosomes = $sam->seq_ids;
# foreach my $chr (@chromosomes) {
# 	open (COVERAGE_OUT,	 ">$out_dir/$species.$chr.coverage.$machine") or die "Cannot open >$out_dir/$species.$chr.coverage.$machine";
# 	$chr =~ tr/sl/SL/; # caps in chr names aren't consisent
# 	my $cov_pos = 1;
# 	my $chr_length = $sam->length($chr);
# 	print "Getting $species coverage data for $chr (length = $chr_length)\n";
# 


# 
# my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>'SL2.40ch00'); #,-start=>1128581,-end=>1128610);
# my @coveragedata       = $coverage->coverage;


foreach my $file_in (@files) { ## LOOP
	# specify and open files
	my ($filename, $directories) = fileparse($file_in, ".csv");
	my $file_out;

	if ($machine eq "col") {
		$file_out = $directories . $filename . ".nogap";
	}elsif ($machine eq "whit") {
		$file_out = $directories . $filename . ".gap.csv";
	}else{die "What?"}
	open (FILE_IN, $file_in) or die "Could not open $file_in\n";
	<FILE_IN>;
	my $row = <FILE_IN>;
	my @chr = split(/,/, $row);
	close FILE_IN;

	open (FILE_IN, $file_in) or die "Could not open $file_in\n";
	open (FILE_OUT, ">$file_out") or die "Could not open $file_out\n";

	print "Calculating flanking coverage for $chr[0]\n";
	my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr[0]); #,-start=>1128581,-end=>1128610);
	my @coveragedata = $coverage->coverage;

	print "Writing flanking coverage to $file_out\n";

	my $header = <FILE_IN>;
	chomp $header;
	print FILE_OUT join(",",$header, "gap_cov_pos-8","gap_cov_pos","gap_cov_pos+8"), "\n";
	foreach my $line (<FILE_IN>) {
		chomp $line;
		my @delim_line = split(/,/,$line);
		my $position = $delim_line[1];
	# 	die "$position isn't a number!" if $position =~ m/\D/; # die if contains non-digit character
	# 	die "$delim_line[1] isn't a number!" if $delim_line[1] =~ m/\D/; # die if contains non-digit character
		if ($position =~ m/\D/) {
			if ($position =~ m/\./) {
				my @split_pos = split(/\./,$position);
				$position = $split_pos[0];
			} else {die "$position isn't a number!"}
		}
		print FILE_OUT join(",",$line, $coveragedata[$position-8-1], $coveragedata[$position-1], $coveragedata[$position+8-1]), "\n";
	#	print "NEWDEL\t", $_, "\t", $coveragedata[$_-8-1], "\t", $coveragedata[$_-1], "\t", $coveragedata[$_+8-1], "\n";
	}
	
	
	
	close FILE_IN;
	close FILE_OUT;
}

exit;