package coverage_commander;
use Moose;
use Modern::Perl;
use File::Basename;

#TODO: check for presence of bam file and region!!!!
#TODO: check samtools version!!!!

has 'bam' => (
    is  => 'ro',
    isa => 'Str',
);

sub valid_bam {
    my $self = shift;

    if ( -e $self->bam and $self->bam =~ m/.bam$/i ) {
        say "  Found valid bam file: " . $self->bam if $self->verbose;
        return 1;
    }
    else {
        die "  Can't find valid bam file: " . $self->bam;
    }
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

sub region {
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
    my $samtools_cmd = $samtools_prefix . $self->region . $self->bam . $samtools_suffix;

    return $samtools_cmd;
}

sub get_coverage {
    my $self = shift;

    if ( $self->valid_bam() ) {
        say "  Running: " . $self->samtools_cmd() if $self->verbose();
        system( $self->samtools_cmd );
    }
}

no Moose;
__PACKAGE__->meta->make_immutable;

