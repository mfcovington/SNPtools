package SNPtools::SNPfinder;
use namespace::autoclean;
use Moose;
extends 'SNPtools';
use MooseX::UndefTolerant;
use feature 'say';
use Parallel::ForkManager;
use autodie;
use FindBin qw($Bin);

#TODO: move common subroutines elsewhere (validation, mkdir, bam header-related, etc.)
#TODO: require certain arguments to be defined
#TODO: generate log files??
#TODO: verbose + very verbose
#TODO: make _out_dir_snp, etc. subs
#TODO: Integrate identify-polymorphisms.pl functionality

sub BUILD {
    my $self = shift;

    $self->_validity_tests;
}

my $bin_dir = "$Bin/../../bin";


# Public Attributes

has 'cov_min' => (
    is      => 'rw',
    isa     => 'Int',
    default => 4,
    lazy    => 1,
);

has 'indel_min' => (
    is      => 'rw',
    isa     => 'Num',
    default => 0.66,
    lazy    => 1,
);

has 'snp_min' => (
    is      => 'rw',
    isa     => 'Num',
    default => 0.33,
    lazy    => 1,
);


# Public Methods

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
$bin_dir/SNPfinder/flanking_coverage_calculator.pl \\
    --snp_file   $snp_file   \\
    --sample_id  $sample_id  \\
    --chromosome $chromosome \\
    --cov_db_dir $cov_dir
EOF

#ADD FLANK DIST

    say "  Running:\n  " . $flanking_cov_cmd if $self->verbose();
    system($flanking_cov_cmd );
}

sub identify_snps {
    my $self = shift;

    $self->_validity_tests();
    $self->_make_dir();

    my $identify_snps_cmd =
      "$bin_dir/SNPfinder/identify-polymorphisms.pl \\
    --chromosome " . $self->_chromosome . " \\
    --outputfile " . $self->out_file . " \\
    --min_cov " . $self->cov_min . " \\
    --min_snp_ratio " . $self->snp_min . " \\
    --min_ins_ratio " . $self->indel_min . " \\
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


# Private Methods

sub _validity_tests {
    my $self = shift;

    $self->_validity_tests_samtools;
    $self->_valid_fasta;
    $self->_valid_bam;
    $self->_valid_bam_index;
}

__PACKAGE__->meta->make_immutable;
