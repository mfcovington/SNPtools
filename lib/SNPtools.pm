package SNPtools;
use namespace::autoclean;
use Moose;
use MooseX::UndefTolerant;
use feature 'say';
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

# Public Methods

sub bam_index {
    my $self = shift;

    $self->_validity_tests_samtools;
    $self->_valid_bam;
    say "  Building index for " . $self->bam if $self->verbose;
    my $samtools_cmd = "samtools index " . $self->bam;
    system( $samtools_cmd );
}

sub get_seq_lengths {
    my $self = shift;

    say "  Getting sequence lengths from bam file" if $self->verbose;
    my @header = $self->_get_header;
    my @seq_lengths = map { $_ =~ m/\t SN: .* \t LN: (.*)/x } @header;
    return @seq_lengths;
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

__PACKAGE__->meta->make_immutable;
