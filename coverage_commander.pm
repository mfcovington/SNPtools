package coverage_commander;
use Moose;
use Modern::Perl;
use File::Basename;
use File::Path 'make_path';
use autodie;
use Data::Printer;

#TODO: check for presence of valid region!!!!
#TODO: require certain arguments to be defined
#TODO: generate log files??
#will it cause a problem if i look up a region that has no coverage?  will it return empty string, undef or 0?....looks like empty or undef

sub samtools_cmd_gaps {
    my $self = shift;

    my $samtools_cmd = "samtools mpileup" . $self->_region . $self->bam . " | cut -f1-2,4 > " . $self->out_file . ".gaps";
    return $samtools_cmd;
}

sub samtools_cmd_nogaps {
    my $self = shift;

    my $samtools_cmd = "samtools depth" . $self->_region . $self->bam . " > " . $self->out_file . ".nogaps";
    return $samtools_cmd;
}

sub get_seq_names {
    my $self = shift;

    say "  Getting sequence names from bam file" if $self->verbose;
    my @header = $self->_get_header;
    my @seq_names = map { $_ =~ m/\t SN: (.*) \t LN:/x } @header;
    return @seq_names;
}

sub get_seq_lengths {
    my $self = shift;

    say "  Getting sequence names from bam file" if $self->verbose;
    my @header = $self->_get_header;
    my @seq_lengths = map { $_ =~ m/\t SN: .* \t LN: (.*)/x } @header;
    return @seq_lengths;
}

sub get_coverage {
    my $self = shift;

    $self->_validity_tests();
    $self->_make_dir();

    if ( $self->gap ) {
        say "  Running: " . $self->samtools_cmd_gaps() if $self->verbose();
        system( $self->samtools_cmd_gaps );
    }
    if ( $self->nogap ) {
        say "  Running: " . $self->samtools_cmd_nogaps() if $self->verbose();
        system( $self->samtools_cmd_nogaps );
    }
}

has 'bam' => (
    is  => 'ro',
    isa => 'Str',
);

has 'chromosome' => (
    is  => 'rw',
    isa => 'Str',
);

has 'pos_start' => (
    is  => 'rw',
    isa => 'Int',
);

has 'pos_end' => (
    is  => 'rw',
    isa => 'Int',
);

has 'out_file' => (
    is  => 'rw',
    isa => 'Str',
);

has 'gap' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 1,
);

has 'nogap' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 1,
);

has 'verbose' => (
    is      => 'ro',
    isa     => 'Bool',
    default => 0,
);

sub _region {
    my $self = shift;

    my $region;
    given ( $self->chromosome ) {
        $region = " -r " . $self->chromosome . ":" . $self->pos_start . "-" . $self->pos_end . " "
            when defined
            and defined $self->pos_start
            and defined $self->pos_end
            and $self->pos_start < $self->pos_end;
        $region = " -r " . $self->chromosome . " " when defined;
        default { $region = " " }
    }

    return $region;
}

sub _make_dir {
    my $self = shift;

    my ( $filename, $dir_name ) = fileparse( $self->out_file );
    make_path( $dir_name );
}

sub _validity_tests {
    my $self = shift;

    $self->_valid_bam;
    $self->_valid_samtools_path;
    $self->_valid_samtools_version;
}

sub samtools_cmd_header {
    my $self = shift;

    my $samtools_cmd = "`samtools view -H " . $self->bam . "`";
    return $samtools_cmd;
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

    if ( -e $self->bam and $self->bam =~ m/.bam$/i ) {
        say "  Found valid bam file: " . $self->bam if $self->verbose;
        return 1;
    }
    else {
        die "  Can't find valid bam file: " . $self->bam;
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

no Moose;
__PACKAGE__->meta->make_immutable;

