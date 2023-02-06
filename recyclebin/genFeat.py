import click
from gaussianutility.utils import feat_gen

@click.command(name = 'genFeat',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('file_name', type=str) # it must include a file format

def main(file_name):
    """Generate features of Gaussian input/output or XYZ file for ML application.
    The order of the features is:
       # Al, # H, # O, # Si, # generated water, path feature, connectivity feature

    The last one is a target value = formation Gibbs free energy per # monomer
    (@ 100 'C), if applicable (normally terminated Gaussian output file)

    """

    return_val = feat_gen(file_name)
    click.echo(return_val)

