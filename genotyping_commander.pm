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

sub samtools_cmd_mpileup {
    my $self = shift;

    my $samtools_cmd =
        "samtools mpileup -l "
      . $self->snp_dir . "/"
      . join( '.', "polyDB", $self->chromosome ) . " -f "
      . $self->fasta . " "
      . $self->bam . " > "
      . $self->mpileup_dir . "/"
      . join( '.', $self->id, $self->chromosome, $self->_mpileup_suffix );

    return $samtools_cmd;
}

sub extract_mpileup {
    my $self = shift;

    $self->_make_dir( $self->mpileup_dir );

    say "  Running: " . $self->samtools_cmd_mpileup() if $self->verbose();
    system( $self->samtools_cmd_mpileup );
}

sub genotype {
    my $self = shift;

    $self->_make_dir( $self->genotyped_dir );

    my $pileup_path = $self->mpileup_dir . "/" . join( '.', $self->id, $self->chromosome, $self->_mpileup_suffix );
    my $snp_path    = $self->snp_dir     . "/" . join( '.', "polyDB", $self->chromosome );

    my $genotyping_cmd =
      "./genotyping_pileups.pl \\
    --pileup $pileup_path \\
    --snp $snp_path \\
    --out_dir " . $self->genotyped_dir;

    if ( !-e $pileup_path ) {
        say "  Pileup file not found: $pileup_path" if $self->verbose();
        return;
    }
    elsif ( !-e $snp_path ) {
        say "  SNP file not found: $snp_path" if $self->verbose();
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
    my $par1_genotyped = $self->genotyped_dir . "/" . join( '.', $self->par1, $self->chromosome, "genotyped" );
    my $par2_genotyped = $self->genotyped_dir . "/" . join( '.', $self->par2, $self->chromosome, "genotyped" );

    if ( !-e $par1_genotyped ) {
        say "  Parent 1 genotype file not found: $par1_genotyped" if $self->verbose();
        return;
    }
    elsif ( !-e $par2_genotyped ) {
        say "  Parent 2 genotype file not found: $par2_genotyped" if $self->verbose();
        return;
    }
    else {
        $R->run(qq`PAR1 <- read.table("$par1_genotyped")`);
        $R->run(qq`PAR2 <- read.table("$par2_genotyped")`);
        $R->run(q`PAR1_ratio <- PAR1[ , 3 ]/PAR1[ , 5 ]`);
        $R->run(q`PAR2_ratio <- PAR2[ , 4 ]/PAR2[ , 5 ]`);
        my $min_ratio = 0.7;
        $R->run(qq`pos_nr <- PAR1[ PAR1_ratio > $min_ratio & PAR2_ratio > $min_ratio , 2 ]`);
        my $polymorphisms = $self->snp_dir . "/" . join( '.', "polyDB", $self->chromosome );
        my $polymorphisms_nr = $polymorphisms . ".nr";
        $R->run(qq`SNP <- read.table( "$polymorphisms", head = T )`);
        $R->run(q`SNP_nr <- SNP[ is.element( SNP$pos, pos_nr) , ]`);
        $R->run(qq`write.table( SNP_nr, file = "$polymorphisms_nr", quote = F, sep = "\t", row.names = F )`);
    }
};

around [qw(extract_mpileup genotype noise_reduction)] => sub {
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

    say "  Getting sequence names from bam file" if $self->verbose;
    my @header = $self->_get_header;
    my @seq_names = map { $_ =~ m/\t SN: (.*) \t LN:/x } @header;
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
    default => sub {
        my $self = shift;

        return $self->out_dir . "/genotyped/";
    },
    lazy => 1,
);

has 'mpileup_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => sub {
        my $self = shift;

        return $self->out_dir . "/mpileup/";
    },
    lazy => 1,
);

has 'snp_dir' => (
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

sub _mpileup_suffix {
    my $self = shift;

    my $suffix = "mpileup";
    $suffix .= ".nr" unless $self->before_noise_reduction;
    return $suffix;
}

sub _make_dir {
    my $self = shift;
    my $dir_name = shift;

    ( my $filename, $dir_name ) = fileparse( $self->out_file ) unless defined $dir_name;
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


