package SNPtools::Coverage;
use Moose;
# use MooseX::UndefTolerant;
use Modern::Perl;
use File::Basename;
use File::Path 'make_path';
use Parallel::ForkManager;
use autodie;
# use Data::Printer;
use Capture::Tiny 'capture_stderr';
use SNPtools::Coverage::DB::Main;
use POSIX;
use DBI;

#TODO: check for presence of valid region!!!!
#TODO: require certain arguments to be defined
#TODO: generate log files??
#will it cause a problem if i look up a region that has no coverage?  will it return empty string, undef or 0?....looks like empty or undef
#TODO: Do I need to make defaults lazy and uncomment UndefTolerant?

sub samtools_cmd_gaps {
    my $self = shift;

    my $samtools_cmd = "samtools mpileup" . $self->_region . $self->bam . " | cut -f1-2,4 > " . $self->out_file . ".cov_gaps";
    return $samtools_cmd;
}

sub samtools_cmd_nogaps {
    my $self = shift;

    my $samtools_cmd = "samtools depth" . $self->_region . $self->bam . " > " . $self->out_file . ".cov_nogaps";
    return $samtools_cmd;
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

    say "  Getting sequence lengths from bam file" if $self->verbose;
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

sub get_coverage_all {
    my $self = shift;

    if ( $self->gap ) {
        say "  Running: " . $self->samtools_cmd_gaps() if $self->verbose();
        system( $self->samtools_cmd_gaps );
    }
    if ( $self->nogap ) {
        say "  Running: " . $self->samtools_cmd_nogaps() if $self->verbose();
        system( $self->samtools_cmd_nogaps );
    }
}

around 'get_coverage_all' => sub {
    my $orig = shift;
    my $self = shift;

    $self->_validity_tests();

    my @chromosomes = $self->get_seq_names;
    my $pm = new Parallel::ForkManager($self->threads);
    foreach my $chr (@chromosomes) {
        $pm->start and next;

        $self->_chromosome($chr);
        $self->out_file( $self->out_dir . "/coverage/" . $self->id . "." . $self->_chromosome . ".coverage" );
        $self->_make_dir();

        $self->$orig(@_);

        $pm->finish;
    }
    $pm->wait_all_children;
};

has 'cov_pos' => (
    is  => 'rw',
    isa => 'HashRef'
);

has 'cov_data' => (
    is  => 'ro',
    isa => 'ArrayRef'
);

has 'flank_dist' => (
    is      => 'rw',
    isa     => 'Int',
    default => 8,
);

has 'id' => (
    is      => 'ro',
    isa     => 'Str',
    default => "unidentified",
);

has 'bam' => (
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

has 'pos_start' => (
    is  => 'rw',
    isa => 'Int',
);

has 'pos_end' => (
    is  => 'rw',
    isa => 'Int',
);

has 'out_dir' => (
    is  => 'rw',
    isa => 'Str',
    default => "./",
);

has 'out_file' => (
    is  => 'rw',
    isa => 'Str',
);

has 'db' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 1,
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

has 'threads' => (
    is      => 'rw',
    isa     => 'Int',
    default => 1,
);

has 'verbose' => (
    is      => 'ro',
    isa     => 'Bool',
    default => 0,
);


sub add_positions {
    my $self = shift;

    my $chr = $self->_chromosome;
    my $flank_dist = $self->flank_dist;
    my $out_dir = $self->out_dir;
    my $sample_id = $self->id;

    my %cov_pos; # = $self->cov_pos;

# TODO: custom path
    open my $snps_fh, "<", "$out_dir/snps/$sample_id.$chr.snps.csv"; #/genotyping/snp_master/polyDB.$chr.nr";
    <$snps_fh>;
    while (<$snps_fh>) {
        my $snp_pos = [ split /,/ ]->[1];
        $cov_pos{$chr}{$snp_pos}                 = 1;
        $cov_pos{$chr}{ $snp_pos - $flank_dist } = 1;
        $cov_pos{$chr}{ $snp_pos + $flank_dist } = 1;
    }
    close $snps_fh;
    # print scalar keys $cov_pos{$chr}, "\n";
    $self->cov_pos( \%cov_pos );
}

sub get_coverage_db {
    my $self = shift;

    $self->add_positions;
    $self->populate_CoverageDB_by_chr;
}

around [ 'get_coverage_db', 'reciprocal_coverage' ] => sub {
    my $orig = shift;
    my $self = shift;

    $self->_validity_tests();

    my @chromosomes = $self->get_seq_names;
    my $pm = new Parallel::ForkManager( floor $self->threads / 2 );
    foreach my $chr (@chromosomes) {
        $pm->start and next;

        $self->_chromosome($chr);

        $self->$orig(@_);

        $pm->finish;
    }
    $pm->wait_all_children;
};

sub populate_CoverageDB_by_chr {
    my $self = shift;

    my $cov_dir     = $self->out_dir . "/coverage";
    my $chromosome  = $self->_chromosome;
    my $cov_pos_ref = $self->cov_pos;
    my $bam_file    = $self->bam;
    my $sample_id   = $self->id;

# TODO: Create DB table using DBIx::Class instead of DBI
    my $dbi    = 'SQLite';
    my $db     = "$cov_dir/coverage.db";
    my $schema = SNPtools::Coverage::DB::Main->connect("dbi:$dbi:$db");

    create_db( $db, $dbi, $cov_dir, $self->verbose ) unless -f $db;

    say "  Getting coverage for $chromosome" if $self->verbose;
    my $count = 1;
    my @cov_data;

    system("samtools index $bam_file") if ! -e "$bam_file.bai";

    my $sam_gap_cmd = "samtools mpileup -r $chromosome $bam_file | cut -f1-2,4";
    my $sam_nogap_cmd = "samtools depth -r $chromosome $bam_file";

    my $gap_fh;
    capture_stderr {    # suppress mpileup output sent to stderr
        open $gap_fh,   "-|", $sam_gap_cmd;
    };
    open my $nogap_fh, "-|", $sam_nogap_cmd;
    while ( my $gap_line = <$gap_fh> ) {
        my $nogap_line = <$nogap_fh>;
        chomp( $gap_line, $nogap_line );
        my ( $chr, $pos, $gap_cov ) = split /\t/, $gap_line;
        my $nogap_cov = [ split /\t/, $nogap_line ]->[2];
        if ( exists $$cov_pos_ref{$chr}{$pos} ) {
            $count++;
            push @cov_data, [ $sample_id, $chr, $pos, $gap_cov, $nogap_cov ];
        }

        populate_and_reset( \$count, \@cov_data, \$schema ) if $count % 100000 == 0;
    }
    close $gap_fh;
    close $nogap_fh;

    populate_and_reset( \$count, \@cov_data, \$schema ) if scalar @cov_data;
}

# TODO: Convert create_db into OO method
sub create_db {
    my ( $db, $dbi, $db_dir, $verbose ) = @_;
    say "  Creating coverage database: $db" if $verbose;
    make_path( $db_dir );
    my $dbh = DBI->connect("dbi:$dbi:$db");
    $dbh->do(<<'END_SQL');
CREATE TABLE coverage (
sample_id  TEXT    NOT NULL,
chromosome TEXT    NOT NULL,
position   INTEGER NOT NULL,
gap_cov    INTEGER NOT NULL,
nogap_cov  INTEGER,
PRIMARY KEY ( sample_id, chromosome, position )
);
END_SQL
    $dbh->disconnect();
}

sub populate_and_reset {
    my ( $count_ref, $cov_data_ref, $schema_ref ) = @_;
    $$count_ref = 1;
    $$schema_ref->populate(
        'Coverage',
        [
            [qw/sample_id chromosome position gap_cov nogap_cov/],
            @$cov_data_ref
        ]
    );
    @$cov_data_ref = ();
}

sub reciprocal_coverage {

    my $self = shift;
    my $chromosome = $self->_chromosome;
    my $out_dir = $self->out_dir;
    my $cov_dir = "$out_dir/coverage";

    # ----------
    # Run in parallel by chromosome
    # ----------
    # get SNP positions for PAR1
    # get SNP positions for PAR2
    # ----------

    my $par1 = 'R500';
    my $par2 = 'IMB211';
    my $par1_bam = 'sample-files/bam/R500.10kb.bam';
    my $par2_bam = 'sample-files/bam/IMB211.10kb.bam';
    my %bams = (
        $par1 => $par1_bam,
        $par2 => $par2_bam,
    );

    my %reciprocal_par = (
        $par1 => $par2,
        $par2 => $par1,
    );

    my %snp_pos;

    for my $id ( $par1, $par2 ) {
        my $snp_file = "$out_dir/snps/$id.$chromosome.snps.nogap.gap.FILTERED.csv";
        open my $snp_fh, "<", $snp_file;
        <$snp_fh>;
        while (<$snp_fh>) {
            my $pos = ( split /,/ )[1];
            $snp_pos{$id}{$pos} = 1;
        }
        close $snp_fh;
    }

    # ----------
    # get cov positions for PAR1
    # get cov positions for PAR2
    # ----------
    # compare lists and get:
    # ----------
    # SNP positions in PAR1 w/o PAR2 cov
    # SNP positions in PAR2 w/o PAR1 cov
    # ----------

    # my $cov_pos_ref = $self->cov_pos;
    # my $bam_file    = $self->bam;
    # my $sample_id   = $self->id;

    my $dbi    = 'SQLite';
    my $db     = "$cov_dir/coverage.db";
    my $schema = SNPtools::Coverage::DB::Main->connect("dbi:$dbi:$db");

    # say "  Getting reciprocal coverage for $chromosome" if $self->verbose;
    # my $count = 1;
    # my @cov_data;

    for my $id ( $par1, $par2 ) {

        my $other_par = $reciprocal_par{$id};

        my $all_ref = get_pos_from_cov_db( $id, $chromosome, \$schema );

        for my $pos ( @$all_ref ) {
            delete $snp_pos{$other_par}{$pos->position};
        }

    }

    # calculate coverage for PAR1/PAR2 at 'missing' positions
    # populate DB

}

sub get_pos_from_cov_db {

    my ( $sample_id, $chr, $schema_ref ) = @_;

    my $rs = $$schema_ref->resultset('Coverage')->search(
        {
            'sample_id'  => $sample_id,
            'chromosome' => $chr,
        },
        { select => [qw/ position /] }
    );

    my @all = $rs->all;

    return \@all;
}

sub _region {
    my $self = shift;

    my $region;
    given ( $self->_chromosome ) {
        $region = " -r " . $self->_chromosome . ":" . $self->pos_start . "-" . $self->pos_end . " "
            when defined
            and defined $self->pos_start
            and defined $self->pos_end
            and $self->pos_start < $self->pos_end;
        $region = " -r " . $self->_chromosome . " " when defined;
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

    $self->_validity_tests_samtools;
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

    # $self->_validity_tests();
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

