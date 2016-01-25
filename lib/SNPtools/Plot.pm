package SNPtools::Plot;
use namespace::autoclean;
use Moose;
extends 'SNPtools';
use MooseX::UndefTolerant;
use File::Path 'make_path';
use Parallel::ForkManager;
use autodie;
use Statistics::R;
use feature 'say';
use List::Util 'max';
use FindBin qw($Bin);

# TODO:
# - chromosome nicknames/abbreviations
# - fork to genoplot multiple IDs in parallel
# - when getting sequence names from custom seq_list/region,
#   validate against bam header before `return @seq_names;`
#    - do this for all relevant modules
# - fix x-axis labels when plotting regions
# - add option to change backgorund color

my $bin_dir = "$Bin/../../bin";


# Public Attributes

has 'before_noise_reduction' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 0,
    lazy    => 1,
);

has 'col_het' => (
    is      => 'ro',
    isa     => 'Str',
    default => 'black',
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

has 'plot_format' => (
    is      => 'ro',
    isa     => 'Str',
    default => 'png',
    lazy    => 1,
);

has 'plot_height' => (
    is      => 'ro',
    isa     => 'Num',
    default => 8,
    lazy    => 1,
);

has 'plot_width' => (
    is      => 'ro',
    isa     => 'Num',
    default => 10,
    lazy    => 1,
);

has 'chr_pat' => (
    is  => 'rw',
    isa => 'Str',
);

has 'chr_sub' => (
    is  => 'rw',
    isa => 'Str',
);

has 'region' => (
    is  => 'rw',
    isa => 'Str',
);


# Private Attributes

has '_region_end' => (
    is  => 'rw',
    isa => 'Int',
);

has '_region_start' => (
    is  => 'rw',
    isa => 'Int',
);


# Public Methods

sub genoplot_by_chr {
    my $self = shift;

    my $id          = $self->id;
    my $chromosome  = $self->ch;
    my $plot_format = $self->plot_format;
    my $plot_width  = $self->plot_width;
    my $plot_height = $self->plot_height;

    # $self->_make_tmp_dir;
    my $R = Statistics::R->new();

    $R->run_from_file("$bin_dir/Plot/genoplot_by_chr.build_df.R");
    $R->run_from_file("$bin_dir/Plot/genoplot_by_chr.build_plot.R");
    $R->run_from_file("$bin_dir/Plot/genoplot_by_chr.add_summary.R")
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
    my $chr_pat     = $self->chr_pat;
    my $chr_sub     = $self->chr_sub;
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
        $R->run_from_file("$bin_dir/Plot/genoplot_by_id.region.build_df.R");
    }
    else {
        say "  Building data frame for $id plot." if $self->verbose;
        $R->run_from_file("$bin_dir/Plot/genoplot_by_id.build_df.R");
    }

    if ( defined $chr_pat && defined $chr_sub ) {
        say "  Renaming chromosomes." if $self->verbose;
        my $chr_ids = 'geno_df$chr';
        my $chromosome_renaming_command
            = "$chr_ids <- sub('$chr_pat', '$chr_sub', $chr_ids)";
        $R->run(qq`$chromosome_renaming_command`);
    }

    say "  Generating plot for $id." if $self->verbose;
    $R->run_from_file("$bin_dir/Plot/genoplot_by_id.build_plot.R");

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


# Private Methods

sub _get_genofile {
    my $self = shift;

    my $genofile = join ".", $self->id, $self->_chromosome, "genotyped";
    $genofile .= ".nr" unless $self->before_noise_reduction;
    return $genofile;
}

sub _plot_path {
    my $self = shift;

    my $path = join "/", $self->_plot_dir, $self->id;
    $path .=
      "." . $self->_chromosome . "." . $self->_region_start . "-" . $self->_region_end
      if defined $self->region;
    $path .= ".before_nr" if $self->before_noise_reduction;
    return $path;
}

__PACKAGE__->meta->make_immutable;
