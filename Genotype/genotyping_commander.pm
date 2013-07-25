package genotyping_commander;
use Moose;
use MooseX::UndefTolerant;
use Modern::Perl;
use File::Basename;
use File::Path 'make_path';
use Parallel::ForkManager;
use Statistics::R;
use autodie;
# use Data::Printer;

#TODO:
# add relevant validity tests for extract_mpileup
# make so that validity tests are done once and remembered
# allow override of samtools version check
# TO DO: incorporate option to ignore indels (do for snp ID, too?) (see line 60)
# Update "  Need samtools version 0.1.XX+" in sub _valid_samtools_version
# Add method that returns full usage statement

sub extract_mpileup {
    my $self = shift;

    $self->_make_dir( $self->_mpileup_dir );

    my $samtools_cmd =
        "samtools mpileup \\
      -l ${ \$self->_snp_path } \\
      -f ${ \$self->fasta } \\
         ${ \$self->bam } \\
      >  ${ \$self->_pileup_path }";

    if ( ! -e $self->_snp_path ) {
        say "  SNP file not found: ${ \$self->_snp_path }" if $self->verbose();
        return;
    }
    else {
        say "  Running: $samtools_cmd" if $self->verbose();
        system( $samtools_cmd );
    }
}

sub genotype {
    my $self = shift;

    $self->_make_dir( $self->_genotyped_dir );

    my $genotyping_cmd =
      "./genotyping_pileups.pl \\
    --mpileup  ${ \$self->_pileup_path } \\
    --snp      ${ \$self->_snp_path } \\
    --par1_id  ${ \$self->par1 } \\
    --par2_id  ${ \$self->par2 } \\
    --out_file ${ \$self->_genotyped_path }";

    # TO DO: incorporate option to ignore indels (do for snp ID, too?):
    # $genotyping_cmd .= " --no_indels" if $no_indels;

    if ( ! -e $self->_snp_path ) {
        say "  SNP file not found: ${ \$self->_snp_path }" if $self->verbose();
        return;
    }
    elsif ( ! -e $self->_pileup_path ) {
        say "  Pileup file not found: ${ \$self->_pileup_path }" if $self->verbose();
        return;
    }
    else {
        say "  Running:\n  " . $genotyping_cmd if $self->verbose();
        system( $genotyping_cmd );
    }
}

sub noise_reduction {
    my $self = shift;

    my $R = Statistics::R->new();
    my $par1_genotyped = $self->_genotyped_dir . "/" . join( '.', $self->par1, $self->_chromosome, "genotyped" );
    my $par2_genotyped = $self->_genotyped_dir . "/" . join( '.', $self->par2, $self->_chromosome, "genotyped" );

    if ( ! -e $par1_genotyped ) {
        say "  Parent 1 genotype file not found: $par1_genotyped" if $self->verbose();
        return;
    }
    elsif ( ! -e $par2_genotyped ) {
        say "  Parent 2 genotype file not found: $par2_genotyped" if $self->verbose();
        return;
    }
    else {
        my $min_ratio     = $self->nr_ratio;
        my $polymorphisms = $self->_snp_path;
        $self->before_noise_reduction(0);
        my $polymorphisms_nr = $self->_snp_path;

        my $cmd_id_pos_pass_ratio = <<EOF;
PAR1 <- read.table("$par1_genotyped")
PAR2 <- read.table("$par2_genotyped")
PAR1_ratio <- PAR1[ , 3 ]/PAR1[ , 5 ]
PAR2_ratio <- PAR2[ , 4 ]/PAR2[ , 5 ]
pos_nr_PAR1 <- PAR1[ PAR1_ratio >= $min_ratio , 2 ]
pos_nr_PAR2 <- PAR2[ PAR2_ratio >= $min_ratio , 2 ]
pos_nr <- intersect( pos_nr_PAR1, pos_nr_PAR2 )
EOF
        $R->run($cmd_id_pos_pass_ratio);

        my $cmd_filter_and_write_nr_SNPs = <<EOF;
SNP <- read.table( "$polymorphisms", head = T )
SNP_nr <- SNP[ is.element( SNP\$pos, pos_nr) , ]
write.table( SNP_nr, file = "$polymorphisms_nr", quote = F, sep = "\t", row.names = F )
EOF
        $R->run($cmd_filter_and_write_nr_SNPs);
    }
};

around [qw(extract_mpileup genotype noise_reduction)] => sub {
    my $orig = shift;
    my $self = shift;

    my @chromosomes = $self->get_seq_names;
    my $pm = new Parallel::ForkManager($self->threads);
    foreach my $chr (@chromosomes) {
        $pm->start and next;
        $self->_chromosome($chr);
        $self->$orig(@_);
        $pm->finish;
    }
    $pm->wait_all_children;
};

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

has 'id' => (
    is  => 'rw',
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
    is  => 'rw',
    isa => 'Str',
);

has 'fasta' => (
    is  => 'ro',
    isa => 'Str',
);

has 'seq_list' => (
    is  => 'rw',
    isa => 'Str',
);

has '_chromosome' => (
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

has '_genotyped_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => sub {
        my $self = shift;

        return $self->out_dir . "/genotyped/";
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

has 'threads' => (
    is      => 'rw',
    isa     => 'Int',
    default => 1,
    lazy => 1,
);

has 'nr_ratio' => (
    is      => 'rw',
    isa     => 'Num',
    default => 0.7,
    lazy    => 1,
);

has 'verbose' => (
    is      => 'ro',
    isa     => 'Bool',
    default => 0,
    lazy => 1,
);

has 'before_noise_reduction' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 0,
    lazy    => 1,
);

sub _pileup_path {
    my $self = shift;

    return $self->_mpileup_dir . "/"
      . join( '.', $self->id, $self->_chromosome, $self->_mpileup_suffix );
}

sub _snp_path {
    my $self = shift;

    my $path = $self->_snp_dir . "/polyDB." . $self->_chromosome;
    $path .= ".nr" unless $self->before_noise_reduction;
    return $path;
}

sub _genotyped_path {
    my $self = shift;

    return $self->_genotyped_dir . "/"
      . join( '.', $self->id, $self->_chromosome, $self->_genotyped_suffix );
}

sub _mpileup_suffix {
    my $self = shift;

    my $suffix = "mpileup";
    $suffix .= ".nr" unless $self->before_noise_reduction;
    return $suffix;
}

sub _genotyped_suffix {
    my $self = shift;

    my $suffix = "genotyped";
    $suffix .= ".nr" unless $self->before_noise_reduction;
    return $suffix;
}

sub _make_dir {
    my $self = shift;
    my $dir_name = shift;

    ( my $filename, $dir_name ) = fileparse( $self->out_file ) unless defined $dir_name;
    make_path( $dir_name ) unless -e $dir_name;
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


