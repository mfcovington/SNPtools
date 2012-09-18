package genotyping_commander;
use Moose;
use MooseX::UndefTolerant;
use Modern::Perl;
use File::Basename;
use File::Path 'make_path';
use Parallel::ForkManager;
use autodie;
# use Data::Printer;

sub samtools_cmd_mpileup {
    my $self = shift;

    my $samtools_cmd =
        "samtools mpileup -l "
      . $self->snp_dir . "/"
      . join( '.', "polyDB", $self->chromosome ) . " -f "
      . $self->fasta . " > "
      . $self->mpileup_dir . "/"
      . join( '.', $self->id, $self->chromosome, "mpileup" )

    return $samtools_cmd;
}

sub extract_mpileup {
    my $self = shift;

    say "  Running: " . $self->samtools_cmd_mpileup() if $self->verbose();
    system( $self->samtools_cmd_mpileup );
}

around 'extract_mpileup' => sub {
    my $orig = shift;
    my $self = shift;

    my @chromosomes = $self->get_seq_names;
    my $pm = new Parallel::ForkManager($self->threads);
    foreach my $chr (@chromosomes) {
        $pm->start and next;

        $self->chromosome($chr);

        $self->$orig(@_);

        $pm->finish;
    }
    $pm->wait_all_children;
};


sub genotype {
    my $self = shift;

    my $genotyping_cmd =
      "/Volumes/OuterColomaData/Mike/genotyping_pileups.custom_names.pl \\
    --pileup "  . $self->mpileup_dir . "/" . join( '.', $self->id, $self->chromosome, "mpileup" ) . " \\
    --snp "     . $self->snp_dir     . "/" . join( '.', "polyDB", $self->chromosome ) . " \\
    --out_dir " . $self->genotyped_dir . " \\
    --par1 "    . $self->par1 . " \\
    --par2 "    . $self->par2;

    say "  Running:\n  " . $genotyping_cmd if $self->verbose();
    system( $genotyping_cmd );
}

around 'genotype' => sub {
    my $orig = shift;
    my $self = shift;

    my @chromosomes = $self->get_seq_names;
    my $pm = new Parallel::ForkManager($self->threads);
    foreach my $chr (@chromosomes) {
        $pm->start and next;

        $self->chromosome($chr);

        $self->$orig(@_);

        $pm->finish;
    }
    $pm->wait_all_children;
};

has 'id' => (
    is  => 'ro',
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

has 'out_file' => (
    is  => 'rw',
    isa => 'Str',
);

has 'out_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => "./",
    lazy => 1,
);

has 'genotyped_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => {
        my $self = shift;

        return $self->out_dir . "/genotyped/";
    },
    lazy => 1,
);

has 'mpileup_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => {
        my $self = shift;

        return $self->out_dir . "/mpileup/";
    },
    lazy => 1,
);

has 'snp_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => {
        my $self = shift;

        return $self->out_dir . "/snp_master/";
    },
    lazy => 1,
);

has 'threads' => (
    is      => 'rw',
    isa     => 'Int',
    default => 1,
    lazy => 1,
);

has 'verbose' => (
    is      => 'ro',
    isa     => 'Bool',
    default => 0,
    lazy => 1,
);

sub _make_dir {
    my $self = shift;

    my ( $filename, $dir_name ) = fileparse( $self->out_file );
    make_path( $dir_name );
}

sub _validity_tests {
    my $self = shift;

    $self->_validity_tests_samtools;
    $self->_valid_fasta;
    $self->_valid_bam;
    $self->_valid_bam_index;
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


