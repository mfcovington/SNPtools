package SNPtools::SNPfinder;
use namespace::autoclean;
use Moose;
extends 'SNPtools';
use MooseX::UndefTolerant;
use Modern::Perl;
use File::Basename;
use File::Path 'make_path';
use Parallel::ForkManager;
use autodie;
# use Data::Printer;
use FindBin qw($Bin);

#TODO: move common subroutines elsewhere (validation, mkdir, bam header-related, etc.)
#TODO: require certain arguments to be defined
#TODO: generate log files??
#TODO: verbose + very verbose
#TODO: make _out_dir_snp, etc. subs
#TODO: Integrate 01.1.SNP_calling_homos.pl functionality

sub BUILD {
    my $self = shift;

    $self->_validity_tests;
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

    my $identify_snps_cmd =
      "$Bin/SNPfinder/01.1.SNP_calling_homos.pl \\
    --chromosome " . $self->_chromosome . " \\
    --o " . $self->out_file . " \\
    --n_reads " . $self->cov_min . " \\
    --ref_freq " . $self->snp_min . " \\
    --indel_freq " . $self->indel_min . " \\
    --fasta_ref " . $self->fasta . " \\
    --bam_file " . $self->bam;

    say "  Running:\n  " . $identify_snps_cmd if $self->verbose();
    system($identify_snps_cmd );
}

around 'identify_snps' => sub {
    my $orig = shift;
    my $self = shift;

    my @chromosomes = $self->get_seq_names;
    my $pm = new Parallel::ForkManager($self->threads);
    foreach my $chr (@chromosomes) {
        $pm->start and next;

        $self->_chromosome($chr);
        # my $cov_out = $self->out_dir . "/" . $self->id . ".snps." . $self->_chromosome;
        # $self->out_file($cov_out);
        $self->out_file( $self->out_dir . "/snps/" . $self->id . "." . $self->_chromosome . ".snps" );
        # $self->identify_snps;

        $self->$orig(@_);

        $pm->finish;
    }
    $pm->wait_all_children;
};

# TODO: Incorporate script into this module (or coverage module)
sub flanking_cov {
    my $self = shift;

    $self->_validity_tests();
    # $self->_make_dir();

    my $out_dir    = $self->out_dir;
    my $sample_id  = $self->id;
    my $chromosome = $self->_chromosome;
    my $snp_file   = "$out_dir/snps/$sample_id.$chromosome.snps.csv";
    my $cov_dir    = "$out_dir/coverage";

    my $flanking_cov_cmd = <<EOF;
$Bin/flanking_coverage_calculator.pl \\
    --snp_file   $snp_file   \\
    --sample_id  $sample_id  \\
    --chromosome $chromosome \\
    --cov_db_dir $cov_dir
EOF

#ADD FLANK DIST

    say "  Running:\n  " . $flanking_cov_cmd if $self->verbose();
    system($flanking_cov_cmd );
}

around 'flanking_cov' => sub {
    my $orig = shift;
    my $self = shift;

    my @chromosomes = $self->get_seq_names;
    my $pm = new Parallel::ForkManager($self->threads);
    foreach my $chr (@chromosomes) {
        $pm->start and next;

        $self->_chromosome($chr);
        # my $cov_out = $self->out_dir . "/" . $self->id . ".snps." . $self->_chromosome;
        # $self->out_file($cov_out);
        # $self->out_file( $self->out_dir . "/snps/" . $self->id . "." . $self->_chromosome . ".snps" );
        # $self->identify_snps;

        $self->$orig(@_);

        $pm->finish;
    }
    $pm->wait_all_children;
};

has '_chromosome' => (
    is  => 'rw',
    isa => 'Str',
);

has 'cov_min' => (
    is      => 'rw',
    isa     => 'Int',
    default => 4,
    lazy    => 1,
);

has 'snp_min' => (
    is      => 'rw',
    isa     => 'Num',
    default => 0.33,
    lazy    => 1,
);

has 'indel_min' => (
    is      => 'rw',
    isa     => 'Num',
    default => 0.66,
    lazy    => 1,
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

__PACKAGE__->meta->make_immutable;
