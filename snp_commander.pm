package snp_commander;
use Moose;
use Modern::Perl;
use File::Basename;
use File::Path 'make_path';
use autodie;
# use Data::Printer;

#TODO: move common subroutines elsewhere (validation, mkdir, bam header-related, etc.)
#TODO: require certain arguments to be defined
#TODO: generate log files??
#TODO: verbose + very verbose

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

sub identify_snps {
    my $self = shift;

    $self->_validity_tests();
    $self->_make_dir();

    my $identify_snps_cmd_cmd = 
    "~/git.repos/snp_identification/01.1.SNP_calling_homos.pl \\
    --chromosome " . $self->chromosome . " \\
    --o " . $self->out_file . " \\
    --n_reads " . $self->cov_min . " \\
    --ref_freq " . $self->snp_min . " \\
    --indel_freq " . $self->indel_min . " \\
    --fasta_ref " . $self->fasta . " \\
    --bam_file " . $self->bam;

    say "  Running:\n  " . $identify_snps_cmd_cmd if $self->verbose();
    system( $identify_snps_cmd_cmd );
}

has 'bam' => (
    is  => 'ro',
    isa => 'Str',
);

has 'fasta' => (
    is  => 'ro',
    isa => 'Str',
);

has 'chromosome' => (
    is  => 'rw',
    isa => 'Str',
);

has 'cov_min' => (
    is      => 'rw',
    isa     => 'Int',
    default => 4,
);

has 'snp_min' => (
    is      => 'rw',
    isa     => 'Num',
    default => 0.33,
);

has 'indel_min' => (
    is      => 'rw',
    isa     => 'Num',
    default => 0.66,
);

has 'out_file' => (
    is  => 'rw',
    isa => 'Str',
);

has 'verbose' => (
    is      => 'ro',
    isa     => 'Bool',
    default => 0,
);

sub _make_dir {
    my $self = shift;

    my ( $filename, $dir_name ) = fileparse( $self->out_file );
    make_path( $dir_name );
}

sub _validity_tests {
    my $self = shift;

    $self->_valid_bam;
    $self->_valid_fasta;
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

no Moose;
__PACKAGE__->meta->make_immutable;

