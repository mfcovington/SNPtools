package SNPtools;
use namespace::autoclean;
use Moose;
use MooseX::UndefTolerant;
# use feature 'say';
# use Modern::Perl;
# use File::Basename;
# use File::Path 'make_path';
# use Parallel::ForkManager;
# use autodie;
# use Data::Printer;
# use FindBin qw($Bin);

has 'id' => (
    is      => 'ro',
    isa     => 'Str',
    default => "unidentified",
    lazy    => 1,
);

has 'out_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => "./",
    lazy    => 1,
);

__PACKAGE__->meta->make_immutable;
