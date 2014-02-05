package SNPtools;
use namespace::autoclean;
use Moose;
use MooseX::UndefTolerant;
use File::Path 'make_path';
use feature 'say';
use File::Basename;
use autodie;


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

# Private Methods

sub _get_header {
    my $self = shift;

    $self->_validity_tests();
    my $get_header_cmd = "samtools view -H " . $self->bam;
    my @header = `$get_header_cmd`;
    return @header;
}

sub _make_dir {
    my $self = shift;

    my ( $filename, $dir_name ) = fileparse( $self->out_file );
    make_path( $dir_name );
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

sub _valid_fasta {
    my $self = shift;

    if ( -e $self->fasta and $self->fasta =~ m/ \.fasta$ | \.fa$ /ix ) {
        say "  Found valid fasta file: " . $self->fasta if $self->verbose;
        return 1;
    }
    else {
        die "  Can't find valid fasta file: " . $self->fasta;
    }
}

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

sub _validity_tests_samtools {
    my $self = shift;

    $self->_valid_samtools_path;
    $self->_valid_samtools_version;
}

__PACKAGE__->meta->make_immutable;
