#!/usr/bin/perl 
use strict;
use warnings;
use Bio::DB::Sam;
use Data::Dumper;
use List::Util qw(max sum);

my $ratio_threshold = 2; ##MFC

#THIS SCRIPT REMOVES SNPS WITH LOW COVERAGE AND SUPORTED ONLY BY ON SEGMENT OF THE READS
# my %files = (
# 		'chr0' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.FILTERED2.csv"},
# 		'chr1' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.2.2.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.2.2.csv.nogap.gap.FILTERED2.csv"},
# 		'chr0' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.FILTERED2.csv"},
# 		'chr0' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.FILTERED2.csv"},
# 		'chr0' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.FILTERED2.csv"},
# 		'chr0' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.FILTERED2.csv"},
# 		'chr0' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.FILTERED2.csv"},
# 		'chr0' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.FILTERED2.csv"},
# 		'chr0' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.FILTERED2.csv"},
# 		'chr0' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.FILTERED2.csv"},
# 		'chr0' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.FILTERED2.csv"},
# 		'chr0' =>	{'snp_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.csv",
# 					'out_file' => "01.2.SNP_table.PEN.final.Picard_and_GATK.1.1.csv.nogap.gap.FILTERED2.csv"},
# 
# 		'het' =>	{'snp_file' => "SNP_table.PEN.final.Picard_and_GATK.hets.1.1.csv",
# 					'out_file' => "02.2.SNP_table.PEN.final.Picard_and_GATK.hets.1.1.filtered.csv"},
# 			);



# foreach my $sp (keys %files){
my @files = @ARGV;
foreach my $file (@files) {
#OPEN SNP FILE AND GET HEADER AND HASH WTH SP => ARRAY OF POSITIONS OF THE COLUMNS WITH THE NUM OF READS FOR TAH SP IN THE FILE
# 	my $SNP_file = $files{$sp}->{'snp_file'};
# 	my $out_file = $files{$sp}->{'out_file'};
	my $SNP_file = $file;
	my $out_file = "$SNP_file.FILTERED.csv";
	

	print "\nNumber of SNPs to look at:\n";
	my $wc_l = `wc -l $SNP_file `;
	print "$wc_l\n";
	my @temp = split " ",$wc_l;
	my $unfilt_reads_number = $temp[0] - 1;
	open (SNPs, $SNP_file) or die "cannot open $SNP_file\n";
	my $head = <SNPs>;
	chomp $head;
	open(OUT, ">$out_file") or die "Can't open >$out_file for output: $!";
	print OUT "$head\n";
	open(OUT_KEEP, ">$out_file.quint_keep.csv") or die "Can't open >$out_file.quint_keep.csv for output: $!";
	open(OUT_KICK, ">$out_file.quint_kick.csv") or die "Can't open >$out_file.quint_kick.csv for output: $!";
	print OUT_KEEP "internal.terminal", $head, "\n";
	print OUT_KICK "low_del.threshold", $head, "\n";



	
	my ($counter,$counter_passed,$counter_keep,$counter_kick) = (0,0);
	my %bases_n;
	
	while (<SNPs>){
		$counter++;
#		last if($counter>40);
		print $counter," SNPs analyzed\n" if($counter%250000==0);
		#print "$_";
		my @elements = split ",",$_;
#		print "@elements\n";
		my ($scaff_name,$snp_pos,$reference,$snp_base) = @elements[0..2,8];

#		next if ($scaff_name ne "scaffold00395");
#		next if ($elements[1]!=2260033);

		#REMOVE THE DOT DECIMAL FROM THE INSERTIONS: 1001.1 => 1001
		my $short_snp_pos = $snp_pos;
		$short_snp_pos =~ s/\..*//;
	
		#############################
	#	#FILTER 1
	#	#REMOVE SNPS SUPPORTE BY LESS THAN 4 READS
		if(($elements[3] + $elements[4] + $elements[5] + $elements[6])<4){
			if($elements[7]<4){
#				print "BAD!!!\n";
				next;
			}
		}
# 
# 		#############################
# 	#	#FILTER 2
# 	#	#SKIP SNPs FOR WHICH ALL THE READS HAVE THE SNP IN THE SAME FIFTH (COUNTING FROM BOTH DIRECTIONS)
# 		if($elements[15]!=1 && $elements[23]!=1){
# 			print OUT join ",", @elements;
# 			$counter_passed++;
# #			print "close to juntion\t@elements" if($SNP_is_close_to_junction);
# 		}

		#############################
	#	#FILTER 2 REDONE BY MFC
	#	#SKIP SNPs FOR WHICH ALL THE READS HAVE THE SNP IN THE SAME "TERMINAL" FIFTH (COUNTING FROM BOTH DIRECTIONS) AND DON'T PASS FLANKING COVERAGE TEST
		if($elements[15]!=1 && $elements[23]!=1){ ## Keep if the SNPs are not all in one fifth of the reads
			print OUT join ",", @elements;
			$counter_passed++;
		} elsif (($elements[15]==1 && $elements[10]+$elements[11]+$elements[12] > 0) || ($elements[23]==1 && $elements[18]+$elements[19]+$elements[20] > 0)) { ## Keep if the SNPs are all in one fifth of the reads, but it is an internal quintile
			print OUT_KEEP "internal,",join(",", @elements);
			print OUT join ",", @elements;
			$counter_passed++;
			$counter_keep++;
		} elsif (($elements[27]/$elements[26]>$ratio_threshold && $elements[30]/$elements[29]<$ratio_threshold) || ($elements[25]/$elements[26]>$ratio_threshold && $elements[28]/$elements[29]<$ratio_threshold)) { ## Kick if the SNPs are in a terminal quintile and don't pass the flanking coverage test
			print OUT_KICK "threshold,",join(",", @elements);
			$counter_kick++;
		} elsif (($elements[8] eq "del") && ($elements[7] < 4)) {
			print OUT_KICK "low_del,",join(",", @elements);
			$counter_kick++;
		} else {
			print OUT_KEEP "terminal,",join(",", @elements);
			print OUT join ",", @elements;
			$counter_passed++;
			$counter_keep++;
		}

	}

close OUT;
close OUT_KEEP;
close OUT_KICK;
print "Filtered down to $counter_passed (",$unfilt_reads_number-$counter_passed ," reads filtered (",($unfilt_reads_number-$counter_passed)*100/$unfilt_reads_number," %)) \n---Done!;\n";
print "Kept: $counter_keep\nKicked: $counter_kick\n";
print "Saved: $counter_keep/",$unfilt_reads_number-$counter_passed+$counter_keep, " (", $counter_keep/($unfilt_reads_number-$counter_passed+$counter_keep)*100, " %)\n";
}

#
##subroutines

# sub round {
#     my($number) = shift;
#     return int($number + .5);
# }
# 
# sub print_hash{
# 	my %hash = @_;	
# 	my @keys = keys(%hash);
# 	if ($keys[0] =~ m/[^0-9.]/){ 
# 		foreach my $key (sort (keys(%hash))) {
# 		   print "$key => $hash{$key}\n";
# 		}
# 	}else{
# 		foreach my $key (sort {$a <=> $b} (keys(%hash))) {
# 		   print "$key => $hash{$key}\n";
# 		}
# 	}		
# }


exit;

