import click
from gaussianutility.utils import freeze_low

@click.command(name = 'freezeLow',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('file_name', type=str) # it must include a file format

def main(file_name):
    """Freeze Carterian coordinates of atoms in the low ONIOM layer.    
    Currently, ".com" and ".gjf" formats are supported.
    """
    freeze_low(file_name)
