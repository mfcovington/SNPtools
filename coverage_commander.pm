package coverage_commander;
use Moose;
use Modern::Perl;
use File::Basename;

#TODO: check for presence of valid region!!!!

has 'bam' => (
    is  => 'ro',
    isa => 'Str',
);

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

sub _validity_tests {
    my $self = shift;

    $self->_valid_bam;
    $self->_valid_samtools_path;
    $self->_valid_samtools_version;
}

has 'chromosome' => (
    is  => 'ro',
    isa => 'Str',
);

has 'pos_start' => (
    is  => 'ro',
    isa => 'Int',
);

has 'pos_end' => (
    is  => 'ro',
    isa => 'Int',
);

has 'include_intron_gaps' => (
    is      => 'ro',
    isa     => 'Bool',
    default => 0,
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
        $region = " -r $self->chromosome " when defined;
        default { $region = " " }
    }

    return $region;
}

sub samtools_cmd {
    my $self = shift;

    my $samtools_prefix = "samtools ";
    my $samtools_suffix = "";
    if ( $self->include_intron_gaps ) {
        $samtools_prefix .= "mpileup";
        $samtools_suffix  = " | cut -f1-2,4"; #watch memory usage. if this causes a significant increase in memory reqs, change it. 
    }
    else {
        $samtools_prefix .= "depth";
    }
    my $samtools_cmd = $samtools_prefix . $self->_region . $self->bam . $samtools_suffix;

    return $samtools_cmd;
}

sub get_coverage {
    my $self = shift;

    $self->_validity_tests();
    say "  Running: " . $self->samtools_cmd() if $self->verbose();
    system( $self->samtools_cmd );
}

no Moose;
__PACKAGE__->meta->make_immutable;

