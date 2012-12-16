#!/usr/bin/env perl
# genoplot_by_id.pl
# Mike Covington
# created: 2012-12-13
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use genoplot_commander;


my $usage = <<USAGE_END;

USAGE:
$0
  --id          Sample identifier
  --par1        Parent 1 ID
  --par2        Parent 2 ID
  --col_par1    Color for Parent 1 [Magenta]
  --col_par2    Color for Parent 2 [Green]
  --col_het     Color for Heterzygous [Black]
  --bam         Sample alignment file (.bam)
  --seq_list    OPTIONAL: Comma-delimted list of sequence IDs to analyze
                (By default, this list is generated from the bam file header.)
  --out_dir     Output directory [current]
  --threads     Number of threads [1]
  --no_nr       Use if noise reduction has not been performed
  --verbose
  --help

USAGE_END

my (
    $id,    $par1,     $par2,     $col_par1,   $col_par2,
    $col_het,  $bam_file, $seq_list, $out_dir, $threads,
    $no_nr, $verbose,  $help
);
my $options = GetOptions(
    "id=s"       => \$id,
    "par1=s"     => \$par1,
    "par2=s"     => \$par2,
    "col_par1=s"    => \$col_par1,
    "col_par2=s"    => \$col_par2,
    "col_het=s"     => \$col_het,
    "bam=s"      => \$bam_file,
    "seq_list=s" => \$seq_list,
    "out_dir=s"  => \$out_dir,
    "threads=i"  => \$threads,
    "no_nr"      => \$no_nr,
    "verbose"    => \$verbose,
    "help"       => \$help,
);

die $usage unless $options;
die $usage if $help;
die $usage
  unless defined $id
  && defined $par1
  && defined $par2
  && defined $bam_file;

my $genoplot = genoplot_commander->new(
    id       => $id,
    par1     => $par1,
    par2     => $par2,
    col_par1    => $col_par1,
    col_par2    => $col_par2,
    col_het     => $col_het,
    bam      => $bam_file,
    seq_list => $seq_list,
    out_dir  => $out_dir,
    threads  => $threads,
    verbose  => $verbose,
);

$genoplot->before_noise_reduction(1) if $no_nr;
$genoplot->genoplot_by_id;

exit;
