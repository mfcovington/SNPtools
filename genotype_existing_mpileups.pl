#!/usr/bin/env perl
# genotype_mpileups.pl
# Mike Covington
# created: 2011-12-05
#
# Description:
#
use strict;
use warnings;
use Getopt::Long;

use genotyping_commander;

my $usage = <<USAGE_END;

USAGE:
snp_finder.pl
  --id           Sample identifier
  --par1         Parent 1 ID
  --par2         Parent 2 ID
  --out_dir      Output directory [current]
  --threads      Number of threads [1]
  --verbose
  --help

USAGE_END

my ( $id, $par1, $par2, $out_dir, $threads, $verbose, $help );
my $options = GetOptions(
    "id=s"      => \$id,
    "par1=s"    => \$par1,
    "par2=s"    => \$par2,
    "out_dir=s" => \$out_dir,
    "threads=i" => \$threads,
    "verbose"   => \$verbose,
    "help"      => \$help,
);

die $usage unless $options;
die $usage if $help;
die $usage
  unless defined $id

my $geno = genotyping_commander->new(
    id      => $id,
    par1    => $par1,
    par2    => $par2,
    out_dir => $out_dir,
    threads => $threads,
    verbose => $verbose,
);

$geno->genotype;

exit;
