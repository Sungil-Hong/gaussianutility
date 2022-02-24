import click
from gaussianutility.utils import ONIOM_sort

@click.command(name = 'oniomSort',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('file_name', type=str) # it must include a file format
@click.option('-s', '--sort-idx', metavar='<s>', required=False, help='index for sorting')

def main(file_name, sort_idx):
    """Sort ONIOM input file by the layer, from high to low.
    Additional index can be 'x', 'y', 'z' (by coordinate, ascending order)
    or 'Atom' (by atomic number, descending order)
    Sorting by 'Atom' is defualt
    """
        
    ONIOM_sort(file_name, sort_idx)
