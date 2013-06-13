#!/usr/bin/env perl
use strict;
use warnings;
use feature 'say';
use CoverageDB::Main;

my $schema = CoverageDB::Main->connect('dbi:SQLite:db/coverage.db');


get_cov( ( 'test', 'SL2.40ch02', [49914000..49914020] ) );

sub get_cov {
    # my $pos_ref = shift;
    my ( $sample_id, $chr, $pos ) = @_;
    say "get_cov($sample_id, $chr, $pos):";
    my $rs =
      $schema->resultset('Coverage')->search(
        {
            'sample_id'  => $sample_id,
            'chromosome' => $chr,
            'position'   => $pos,
        },
        { select => [qw/ position gap_cov nogap_cov /] }
      );
    # while ( my $cd = $rs->next ) {
    #     print $cd->title . "\n";
    # }
    # my $cov = $rs->first;
    # say $cov->gap_cov;
    # say $cov->nogap_cov;
    say "# of hits: " . $rs->count;
    while ( my $cov = $rs->next ) {
        say $cov->position;
        say $cov->gap_cov;
        say $cov->nogap_cov;
    }

    use Data::Printer;
    p $rs;
}

