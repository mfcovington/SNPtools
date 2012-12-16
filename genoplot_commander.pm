package genoplot_commander;
use Moose;
use MooseX::UndefTolerant;
use File::Basename;
use File::Path 'make_path';
use Parallel::ForkManager;
use autodie;
use Statistics::R;
use feature 'say';
# use Data::Printer;

# TODO:
# - chromosome nicknames/abbreviations
# - fork to genoplot multiple IDs in parallel

sub genoplot_by_chr {
    my $self = shift;

    my $id          = $self->id;
    my $chromosome  = $self->ch;
    my $plot_format = $self->plot_format;
    my $plot_width  = $self->plot_width;
    my $plot_height = $self->plot_height;

    # $self->_make_tmp_dir;
    my $R = Statistics::R->new();

    $R->run_from_file("genoplot_by_chr.build_df.R");
    $R->run_from_file("genoplot_by_chr.build_plot.R");
    $R->run_from_file("genoplot_by_chr.add_summary.R") if $self->plot_summary;
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
    my $plot_path   = $self->_plot_path;
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
    # $R->set( 'chromosomes', @chromosomes );
    $R->run_from_file("genoplot_by_id.build_df.R");
    $R->run_from_file("genoplot_by_id.build_plot.R");
    # $R->run_from_file("genoplot_by_chr.add_summary.R") if $self->plot_summary;
    $R->run(qq`setwd("$plot_dir")`);
    $R->run(
        qq`ggsave(
          filename = paste("$plot_path", "$plot_format", sep = "."),
          plot     = geno.plot,
          width    = $plot_width,
          height   = $plot_height)`
    );
}

sub bam_index {
    my $self = shift;

    $self->_validity_tests_samtools;
    $self->_valid_bam;
    say "  Building index for " . $self->bam if $self->verbose;
    my $samtools_cmd = "samtools index " . $self->bam;
    system( $samtools_cmd );
}

sub get_seq_names {
    my $self = shift;

    my @seq_names;
    if ( defined $self->seq_list ) {
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

has 'filename' => (
    is         => 'rw',
    isa        => 'Str',
);

has 'plot_format' => (
    is      => 'ro',
    isa     => 'Str',
    default => 'png',
);

has 'plot_width' => (
    is      => 'ro',
    isa     => 'Num',
    default => '8',
);

has 'plot_height' => (
    is      => 'ro',
    isa     => 'Num',
    default => '10',
);

has 'id' => (
    is  => 'rw',
    isa => 'Str',
);

has 'par1' => (
    is  => 'ro',
    isa => 'Str',
);

has 'par2' => (
    is  => 'ro',
    isa => 'Str',
);

has 'col_par1' => (
    is      => 'ro',
    isa     => 'Str',
    default => 'magenta',
    lazy    => 1,
);

has 'col_par2' => (
    is      => 'ro',
    isa     => 'Str',
    default => 'green',
    lazy    => 1,
);

has 'col_het' => (
    is      => 'ro',
    isa     => 'Str',
    default => 'black',
    lazy    => 1,
);

has 'bam' => (
    is  => 'ro',
    isa => 'Str',
);

has 'fasta' => (
    is  => 'ro',
    isa => 'Str',
);

has 'seq_list' => (
    is  => 'rw',
    isa => 'Str',
);

has '_chromosome' => (
    is  => 'rw',
    isa => 'Str',
);

has 'out_file' => (
    is  => 'rw',
    isa => 'Str',
);

has 'out_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => "./",
);

has 'threads' => (
    is      => 'rw',
    isa     => 'Int',
    default => 1,
);

has 'verbose' => (
    is      => 'ro',
    isa     => 'Bool',
    default => 0,
);

has '_plot_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => sub {
        my $self = shift;

        return $self->out_dir . "/genoplot/";
    },
    lazy => 1,
);

has '_genotyped_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => sub {
        my $self = shift;

        return $self->out_dir . "/genotyped/";
    },
    lazy => 1,
);

sub _plot_path {
    my $self = shift;

    my $path = join "/", $self->_plot_dir, $self->id;
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




no Moose;
__PACKAGE__->meta->make_immutable;