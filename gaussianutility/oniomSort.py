import click
from gaussianutility.utils import ONIOM_sort

@click.command(name = 'oniomSort',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('file_name', type=str) # it must include a file format
@click.option('-s', '--sort-idx', metavar='<s>', required=False, help='index for sorting')
@click.option('-f', '--freeze-idx', metavar='<s>', required=False, help='index for freezing atom position for geometry optimization')

def main(file_name, sort_idx, freeze_idx):
    """Sort ONIOM input file by the layer, from high to low.
    Additional sorting index can be 'x', 'y', 'z' (by coordinate, ascending order)
    or 'Atom' (by atomic number, descending order)
    Sorting by 'Atom' is defualt
    
    Additional freezing index can be 'H', 'M', 'L', or their combination
    For example, 'ML' freezes coordinates of atoms in the middle and low ONIOM layers
    """
        
    ONIOM_sort(file_name, sort_idx, freeze_idx)
