package CoverageDB::Main::Result::Coverage;
use base qw/DBIx::Class::Core/;

__PACKAGE__->table('coverage');
__PACKAGE__->add_columns(qw/sample_id chromosome position gap_cov nogap_cov/);
__PACKAGE__->set_primary_key( 'sample_id', 'chromosome', 'position' );

1;