package SNPtools::Plot;
use namespace::autoclean;
use Moose;
extends 'SNPtools';
use MooseX::UndefTolerant;
use File::Basename;
use File::Path 'make_path';
use Parallel::ForkManager;
use autodie;
use Statistics::R;
use feature 'say';
use List::Util 'max';
use FindBin qw($Bin);
# use Data::Printer;

# TODO:
# - chromosome nicknames/abbreviations
# - fork to genoplot multiple IDs in parallel
# - when getting sequence names from custom seq_list/region,
#   validate against bam header before `return @seq_names;`
#    - do this for all relevant modules
# - fix x-axis labels when plotting regions
# - add option to change backgorund color

sub BUILD {
    my $self = shift;

    $self->_validity_tests;
}

sub genoplot_by_chr {
    my $self = shift;

    my $id          = $self->id;
    my $chromosome  = $self->ch;
    my $plot_format = $self->plot_format;
    my $plot_width  = $self->plot_width;
    my $plot_height = $self->plot_height;

    # $self->_make_tmp_dir;
    my $R = Statistics::R->new();

    $R->run_from_file("$Bin/Plot/genoplot_by_chr.build_df.R");
    $R->run_from_file("$Bin/Plot/genoplot_by_chr.build_plot.R");
    $R->run_from_file("$Bin/Plot/genoplot_by_chr.add_summary.R")
        if $self->plot_summary;
    $R->run(
        qq`ggsave(
          filename = paste($id, $chromosome, "$plot_format", sep = "."),
          plot = geno.plot,
          width = $plot_width,
          height = $plot_height)`
    );
}


sub genoplot_by_id {
    my $self = shift;

    my $id          = $self->id;
    my $par1        = $self->par1;
    my $par2        = $self->par2;
    my $col_par1    = $self->col_par1;
    my $col_par2    = $self->col_par2;
    my $col_het     = $self->col_het;
    my $plot_format = $self->plot_format;
    my $plot_width  = $self->plot_width;
    my $plot_height = $self->plot_height;
    my $plot_dir    = $self->_plot_dir;
    make_path($plot_dir);

    my @chromosomes = $self->get_seq_names;
    my @filenames;
    foreach my $chr (@chromosomes) {
        $self->_chromosome($chr);
        push @filenames, $self->_genotyped_dir . $self->_get_genofile;
    }

    my $R = Statistics::R->new();

    $R->set( 'filenames', \@filenames );
    $R->set( 'id',        $id );
    $R->set( 'par1',      $par1 );
    $R->set( 'par2',      $par2 );
    $R->set( 'col_par1',  $col_par1 );
    $R->set( 'col_par2',  $col_par2 );
    $R->set( 'col_het',   $col_het );

    if ( defined $self->region ) {
        my @seq_lengths = $self->get_seq_lengths unless $self->_region_end;
        my $max_length  = max @seq_lengths;
        my $start       = $self->_region_start || $self->_region_start(1);
        my $end         = $self->_region_end || $self->_region_end($max_length);
        $R->set( 'start', $start );
        $R->set( 'end',   $end );
        my $region = $self->region;
        say "  Building data frame for $id plot ($region)." if $self->verbose;
        $R->run_from_file("$Bin/Plot/genoplot_by_id.region.build_df.R");
        say "  Generating plot for $id ($region)." if $self->verbose;
        $R->run_from_file("$Bin/Plot/genoplot_by_id.build_plot.R");
    }
    else {
        say "  Building data frame for $id plot." if $self->verbose;
        $R->run_from_file("$Bin/Plot/genoplot_by_id.build_df.R");
        say "  Generating plot for $id." if $self->verbose;
        $R->run_from_file("$Bin/Plot/genoplot_by_id.build_plot.R");
    }

    my $plot_path = $self->_plot_path;

    $R->run(qq`setwd("$plot_dir")`);
    my $plot_name = "$plot_path.$plot_format" ;
    $plot_name =~ s|/{2,}|/|g;
    say "  Saving: $plot_name" if $self->verbose;
    $R->run(
        qq`ggsave(
          filename = "$plot_name",
          plot     = geno.plot,
          width    = $plot_width,
          height   = $plot_height)`
    );
}

sub get_seq_names {
    my $self = shift;

    my @seq_names;
    if ( defined $self->region ) {
        my $region = $self->region;
        my ( $chr, $start, $end ) = $region =~ /(.*):(\d*)-(\d*)/;
        my $region_error =
          "  ERROR: There is a problem with the region designation: $region\n";
        die $region_error unless $start || $end;
        die $region_error if $start && $end && $start > $end;
        $self->_region_start($start)  if $start;
        $self->_region_end($end) if $end;
        say "  Getting region of $chr from $start-$end" if $self->verbose;
        @seq_names = $chr;
    }
    elsif ( defined $self->seq_list ) {
        @seq_names = split /,/, $self->seq_list;
    }
    else {
        say "  Getting sequence names from bam file" if $self->verbose;
        my @header = $self->_get_header;
        @seq_names = map { $_ =~ m/\t SN: (.*) \t LN:/x } @header;
    }
    return @seq_names;
}

sub _get_genofile {
    my $self = shift;

    my $genofile = join ".", $self->id, $self->_chromosome, "genotyped";
    $genofile .= ".nr" unless $self->before_noise_reduction;
    return $genofile;
}

# around [qw(extract_mpileup genotype noise_reduction)] => sub {
#     my $orig = shift;
#     my $self = shift;

#     my @chromosomes = $self->get_seq_names;
#     my $pm = new Parallel::ForkManager($self->threads);
#     foreach my $chr (@chromosomes) {
#         $pm->start and next;
#         $self->_chromosome($chr);
#         $self->$orig(@_);
#         $pm->finish;
#     }
#     $pm->wait_all_children;
# };

has 'before_noise_reduction' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 0,
    lazy    => 1,
);

has 'plot_format' => (
    is      => 'ro',
    isa     => 'Str',
    default => 'png',
    lazy    => 1,
);

has 'plot_width' => (
    is      => 'ro',
    isa     => 'Num',
    default => 10,
    lazy    => 1,
);

has 'plot_height' => (
    is      => 'ro',
    isa     => 'Num',
    default => 8,
    lazy    => 1,
);

has 'col_par1' => (
    is      => 'ro',
    isa     => 'Str',
    default => 'orange',
    lazy    => 1,
);

has 'col_par2' => (
    is      => 'ro',
    isa     => 'Str',
    default => 'sky blue',
    lazy    => 1,
);

has 'col_het' => (
    is      => 'ro',
    isa     => 'Str',
    default => 'black',
    lazy    => 1,
);

has 'region' => (
    is  => 'rw',
    isa => 'Str',
);

has '_region_start' => (
    is  => 'rw',
    isa => 'Int',
);

has '_region_end' => (
    is  => 'rw',
    isa => 'Int',
);

sub _plot_path {
    my $self = shift;

    my $path = join "/", $self->_plot_dir, $self->id;
    $path .=
      "." . $self->_chromosome . "." . $self->_region_start . "-" . $self->_region_end
      if defined $self->region;
    $path .= ".before_nr" if $self->before_noise_reduction;
    return $path;
}

sub _make_dir {
    my $self = shift;

    my ( $filename, $dir_name ) = fileparse( $self->out_file );
    make_path( $dir_name );
}


sub _validity_tests {
    my $self = shift;

    $self->_validity_tests_samtools;
    # $self->_valid_fasta;
    $self->_valid_bam;
    # $self->_valid_bam_index;
}

sub _validity_tests_samtools {
    my $self = shift;

    $self->_valid_samtools_path;
    $self->_valid_samtools_version;
}

sub _get_header {
    my $self = shift;

    $self->_validity_tests();
    my $get_header_cmd = "samtools view -H " . $self->bam;
    my @header = `$get_header_cmd`;
    return @header;
}

sub _valid_bam {
    my $self = shift;

    if ( -e $self->bam and $self->bam =~ m/ \.bam$ /ix ) {
        say "  Found valid bam file: " . $self->bam if $self->verbose;
        return 1;
    }
    else {
        die "  Can't find valid bam file: " . $self->bam;
    }
}

sub _valid_bam_index {
    my $self = shift;

    my ( $bam_prefix, $bam_dir ) = fileparse( $self->bam, ".bam" );
    if ( -e "$bam_dir/$bam_prefix.bai" or -e "$bam_dir/$bam_prefix.bam.bai" ) {
        say "  Found valid index for " . $self->bam if $self->verbose;
        return 1;
    }
    else {
        say "  Can't find valid index for " . $self->bam;
        $self->bam_index;
    }
}

# sub _valid_fasta {
#     my $self = shift;

#     if ( -e $self->fasta and $self->fasta =~ m/ \.fasta$ | \.fa$ /ix ) {
#         say "  Found valid fasta file: " . $self->fasta if $self->verbose;
#         return 1;
#     }
#     else {
#         die "  Can't find valid fasta file: " . $self->fasta;
#     }
# }

sub _valid_samtools_path {
    my $self = shift;

    my $sam_status = `which samtools`;
    die "  Samtools not found in PATH" unless $sam_status =~ m/samtools/i;
    say "  Samtools looks good!" if $self->verbose;
}

sub _valid_samtools_version {
    my $self = shift;

    my @usage = `samtools 2>&1`; #change to samtools
    my $version;
    for (@usage) {
        $_ =~ m/Version: ([\d\.].*) /i;
        $version = $1 and last if $1;
    }
    my @version_parsed = $version =~ m/ (\d*) \. (\d*) \. (\d*) /x;
    say "  Samtools is version $version" if $self->verbose;
    die "  Need samtools version 0.1.XX+"
      unless ( $version_parsed[0] >= 1
        or $version_parsed[1] >= 2
        or $version_parsed[1] == 1 and $version_parsed[2] >= 18 );
}

__PACKAGE__->meta->make_immutable;
