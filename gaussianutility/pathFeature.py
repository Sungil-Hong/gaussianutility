import click
from gaussianutility.utils import path_feature
#from canela import __version__

@click.command(name = 'pathFeature',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
#@click.version_option(__version__)
@click.argument('file_name', type=str) # it must include a file format

def main(file_name):
    """
    Calculate a path feature, which is defined as a weighted average of two shortest paths of every Al-Al pair.
    This feature captures the different positioning or distribution of Al in structures.
    More weight is given to the shorter path.
    
    The input file format must be one among ".com", ".gjf", ".xyz", and ".out" (Gaussian output fule).
    
    Concept credit to Michael Cowan (mcowan92@gmail.com)
    """
    print(path_feature(file_name))
    
