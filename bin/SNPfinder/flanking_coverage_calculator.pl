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
use feature 'say';
use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);

# TODO: fix paths
use SNPtools::Coverage::DB::Main;

# TODO: Get proper (flanking) coverage for inserts. (.XX seems to be causing a problem)

my $usage = <<USAGE_END;

USAGE:
flanking_coverage_calculator.pl
  --snp_file   </PATH/TO/snp_file.csv>
  --sample_id
  --chromosome
  --flank_dist <flanking distance to use [8]>
  --cov_db_dir
  --help

USAGE_END

#get options
my ( $snp_in, $id, $chr, $cov_db_dir, $help );
my $flank_dist = 8;    # default
my $options = GetOptions(
    "snp_file=s"   => \$snp_in,
    "sample_id=s"  => \$id,
    "chromosome=s" => \$chr,
    "flank_dist=i" => \$flank_dist,
    "cov_db_dir=s" => \$cov_db_dir,
    "help"         => \$help,
);
die $usage if $help;
die $usage
  unless
    defined $snp_in &&
    defined $id     &&
    defined $chr;

#generate filenames
my ( $filename, $directories, $suffix ) = fileparse( $snp_in, ".csv" );
my $snp_out = $directories . $filename . ".nogap.gap.csv";

#open files
open my $snp_in_fh,     "<", $snp_in;
open my $snp_out_fh,    ">", $snp_out;

#build coverage hashes
my $dbi      = 'SQLite';
my $db       = "$cov_db_dir/coverage.db";
my $schema   = SNPtools::Coverage::DB::Main->connect("dbi:$dbi:$db");
my %cov_hash = build_cov_hash( ( $id, $chr ) );

#write header
my $header = <$snp_in_fh>;
chomp $header;
say $snp_out_fh join( ",",
    $header,
    "nogap_pos-$flank_dist", "nogap_pos", "nogap_pos+$flank_dist",
    "gap_pos-$flank_dist",   "gap_pos",   "gap_pos+$flank_dist" );

#write coverage
while ( my $snp_line = <$snp_in_fh> ) {
    chomp $snp_line;
    my ( $snp_chr, $snp_pos_unsplit, $snp_remainder ) = split( /,/, $snp_line, 3 );
    my ( $snp_pos, $snp_pos_index ) = split( /\./, $snp_pos_unsplit );

    # nogaps doesn't see deletions, must compensate
    my ( $del_cov ) = ${ [ split /,/, $snp_remainder ] }[5];

    my $lt_pos = $snp_pos - $flank_dist;
    my $rt_pos = $snp_pos + $flank_dist;

    say $snp_out_fh join( ",",
        $snp_chr,
        $snp_pos_unsplit,
        $snp_remainder,
        $cov_hash{$lt_pos}{nogap}             // 0,
        $cov_hash{$snp_pos}{nogap} + $del_cov // 0,
        $cov_hash{$rt_pos}{nogap}             // 0,
        $cov_hash{$lt_pos}{gap}               // 0,
        $cov_hash{$snp_pos}{gap}              // 0,
        $cov_hash{$rt_pos}{gap}               // 0 );
}

#close shop
close $snp_in_fh;
close $snp_out_fh;

sub build_cov_hash {
    my ( $sample_id, $chr ) = @_;

    my $rs = $schema->resultset('Coverage')->search(
        { 'sample_id' => $sample_id, 'chromosome' => $chr, },
        { select => [qw/ position gap_cov nogap_cov /] }
    );

    return
      map { $_->position => { gap => $_->gap_cov, nogap => $_->nogap_cov } }
      $rs->all;
}

exit;