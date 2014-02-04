package SNPtools;
use namespace::autoclean;
use Moose;
use MooseX::UndefTolerant;
# use feature 'say';
# use Modern::Perl;
# use File::Basename;
# use File::Path 'make_path';
# use Parallel::ForkManager;
# use autodie;
# use Data::Printer;
# use FindBin qw($Bin);


# Public Attributes

has 'bam' => (
    is  => 'rw',
    isa => 'Str',
);

has 'fasta' => (
    is  => 'ro',
    isa => 'Str',
);

has 'id' => (
    is      => 'ro',
    isa     => 'Str',
    default => "unidentified",
    lazy    => 1,
);

has 'out_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => "./",
    lazy    => 1,
);

has 'out_file' => (
    is  => 'rw',
    isa => 'Str',
);

has 'par1' => (
    is      => 'ro',
    isa     => 'Str',
);

has 'par2' => (
    is      => 'ro',
    isa     => 'Str',
);

has 'seq_list' => (
    is  => 'rw',
    isa => 'Str',
);

has 'threads' => (
    is      => 'rw',
    isa     => 'Int',
    default => 1,
    lazy    => 1,
);

has 'verbose' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 0,
    lazy    => 1,
);


# Private Attributes

has '_chromosome' => (
    is  => 'rw',
    isa => 'Str',
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

has '_plot_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => sub {
        my $self = shift;

        return $self->out_dir . "/genoplot/";
    },
    lazy => 1,
);

has '_mpileup_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => sub {
        my $self = shift;

        return $self->out_dir . "/mpileup/";
    },
    lazy => 1,
);

has '_snp_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => sub {
        my $self = shift;

        return $self->out_dir . "/snp_master/";
    },
    lazy => 1,
);

__PACKAGE__->meta->make_immutable;
