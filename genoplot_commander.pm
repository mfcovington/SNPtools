package genoplot_commander;
use Moose;
use MooseX::UndefTolerant;
# use Modern::Perl;
# use File::Basename;
use File::Path 'make_path';
use Parallel::ForkManager;
use autodie;
use Statistics::R;
# use Data::Printer;

sub genoplot_by_chr {
    my $self = shift;

    my $id          = $self->id;
    my $chromosome  = $self->ch;
    my $plot_format = $self->plot_format;
    my $plot_width  = $self->plot_width;
    my $plot_height = $self->plot_height;

    # $self->_make_tmp_dir;
    my $R = Statistics::R->new();

    $R->run_from_file("genoplot_by_chr.build_df.R");
    $R->run_from_file("genoplot_by_chr.build_plot.R");
    $R->run_from_file("genoplot_by_chr.add_summary.R") if $self->plot_summary;
    $R->run(
        qq`ggsave(
          filename = paste($id, $chromosome, "$plot_format", sep = "."),
          plot = geno.plot,
          width = $plot_width,
          height = $plot_height)`
    );
}





no Moose;
__PACKAGE__->meta->make_immutable;