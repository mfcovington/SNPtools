#!/usr/bin/perl
# coverage_calculator.v2.pl
# Mike Covington
# created: 2011-12-05
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use Bio::DB::Sam;
use IO::File;

my $usage = <<USAGE_END;

USAGE:
coverage_calculator.v2.pl
  --ref </PATH/TO/reference.fa>
  --bam </PATH/TO/file.bam>
  --id <sample identifier for file.bam>
  --out </PATH/TO/DESTINATION/DIRECTORY/>
  --col OR --whit
  --help

USAGE_END

my ( $ref_fa, $bam_file, $machine, $col, $whit, $id, $help );
my $out_dir = "./";
my $options = GetOptions(
    "ref=s" => \$ref_fa,
    "bam=s" => \$bam_file,
    "id=s"  => \$id,
    "out=s" => \$out_dir,
    "col"   => \$col,
    "whit"  => \$whit,
    "help"  => \$help,
);

die $usage if $help;
die $usage if defined $col && defined $whit;
die $usage unless defined $ref_fa && defined $bam_file && defined $id;

if ($col) {
    $machine = "col";
}
elsif ($whit) {
    $machine = "whit";
}
else {
    die "Something is wrong...\n";
}

my $sam = Bio::DB::Sam->new(
    -bam   => $bam_file,
    -fasta => $ref_fa
);
my @chromosomes = $sam->seq_ids;

open my $log_fh, ">", "$out_dir/$id.$machine.log";
$log_fh->autoflush(1);

foreach my $chr (@chromosomes) {

    # $chr =~ tr/sl/SL/; # caps in chr names aren't consisent
    open my $cov_out_fh, ">", "$out_dir/$id.$chr.coverage.$machine";
    $cov_out_fh->autoflush(1);

    my $cov_pos    = 1;
    my $chr_length = $sam->length($chr);
    say $log_fh "Getting $id coverage data for $chr (length = $chr_length)";

    my $chunk_size = 999999;
    while ( $cov_pos <= $chr_length ) {

        if ( $cov_pos + $chunk_size < $chr_length ) {
            my ($coverage) = $sam->features(
                -type   => 'coverage',
                -seq_id => $chr,
                -start  => $cov_pos,
                -end    => $cov_pos + $chunk_size
            );
            my @coveragedata = $coverage->coverage;

            foreach my $cov_val (@coveragedata) {
                say $cov_out_fh "$chr,$cov_pos,$cov_val";
                $cov_pos++;
            }

        }
        else {
            my ($coverage) = $sam->features(
                -type   => 'coverage',
                -seq_id => $chr,
                -start  => $cov_pos,
                -end    => $chr_length
            );
            my @coveragedata = $coverage->coverage;

            foreach my $cov_val (@coveragedata) {
                say $cov_out_fh "$chr,$cov_pos,$cov_val";
                $cov_pos++;
            }

        }
    }
    close $cov_out_fh;
}
close $log_fh;
exit;

