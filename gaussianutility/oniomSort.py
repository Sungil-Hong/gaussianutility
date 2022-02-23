import click
from gaussianutility.utils import ONIOM_sort

@click.command(name = 'oniomSort',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('file_name', type=str) # it must include a file format
@click.option('-i','--sort', metavar='<s>', required=False,
             help='sort index')


def main(file_name, sort_idx):
    """Sort ONIOM input file by the layer, from high to low.
    Additional index can be 'x', 'y', 'z' (by coordinate, ascending order)
    or 'Atom' (by atomic number, descending order)
    Sorting by 'Atom' is defualt
    """
    
    if not sort_idx:
        sort_idx = 'Atom'
        
    ONIOM_sort(file_name, sort_idx)
