import click
from gaussianutility.utils import com2VASP

@click.command(name = 'com2vasp',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('com_file', type=str)
@click.option('-t', '--title', required=False, help='Title of the input file')
@click.option('-l', '--lattice_param', required=False, nargs = 3, help='List of lattice parameters, a, b, and c')

def main(com_file, title, lattice_param):
    """
    Convert .com file to .vasp file.
    Currently, only support cubin structure, where a,b,c represent lattice parameters.
    """

    com2VASP(com_file, title, lattice_param)
    
