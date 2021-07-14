import click
from gaussianutility.utils import extract_geom

@click.command(name = 'xgeom',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('file_name', type=str) # it must include a file format
@click.option('-o','--out-file', metavar='<s>', required=False,
             help='output file name')

def main(file_name, out_file):
    """Extract the last/optimized geometry from Gaussian output file
    and generate Gaussian input fle (.com or .gjf) or xyz file.
    
    The argument is the file name of a Gaussian output file, including ".out".
    
    Default output name and format is: input file name + "geom.com".
    
    To change, use -o option.
    
    Currently, ".com", ".gjf", and ".xyz" formats are supported.
    """
    if not out_file:
        out_file = file_name.rsplit(".",1)[0] + ".geom.com"
        
    extract_geom(file_name, out_file)

