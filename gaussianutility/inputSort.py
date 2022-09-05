import click
from gaussianutility.utils import input_sort

@click.command(name = 'inputSort',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('file_name', type=str) # it must include a file format
@click.option('-s', '--sort-idx', metavar='<s>', required=False, help='index for sorting')
@click.option('-f', '--freeze-idx', metavar='<s>', required=False, help='index for freezing atom position for geometry optimization')

def main(file_name, sort_idx, freeze_idx):
    """Sort input file by index; 'x', 'y', or 'z' coordinate, or 'Atom' symbol
    Sorting by 'Atom' is defualt
    
    Additional freezing index can be 'H', 'M', 'L', or their combination
    For example, 'ML' freezes coordinates of atoms in the middle and low ONIOM layers
    This is only for ONIOM type input structure
    """
        
    input_sort(file_name, sort_idx, freeze_idx)
