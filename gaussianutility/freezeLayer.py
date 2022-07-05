import click
from gaussianutility.utils import freeze_layer

@click.command(name = 'freezeLayer',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('file_name', type=str) # it must include a file format
@click.option('-i', '--layer-idx', metavar='<s>', required=False, help='Index of layer to freeze: H, M, L or their combination')

def main(file_name, layer_idx):
    """
    Freeze all the atoms in the specified ONIOM layer(s) for geometry optimization
    by setting freeze-code as -1
    High, middle, and low ONIOM layers are called H, M, and L, respectively
    Default is freezing the low layer only
    """    

    return_val = freeze_layer(file_name, layer_idx)
    click.echo(return_val)

