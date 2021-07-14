#! usr/bin/env/python

import numpy as np
import ase.io
import networkx as nx
import featpro.utils as utils
from itertools import islice

def k_shortest_paths(G, source, target, k, weight=None):
    return list(
        islice(nx.shortest_simple_paths(G, source, target, weight=weight), k)
    )
    # Sourced from https://networkx.org/documentation/stable/index.html

def path_feature(file_name):
    """
    Calculate a path feature, which is defined as a weighted average of two shortest paths of every Al-Al pair.
    This feature captures the different positioning or distribution of Al in structures.
    More weight is given to the shorter path.
    
    The input file format must be one among ".com", ".gjf", ".xyz", and ".out" (Gaussian output fule).
    
    Concept credit to Michael Cowan (mcowan92@gmail.com)
    """
    informat = file_name.rsplit(".",1)[-1]
    if informat == 'com' or informat == 'gjf' or informat == 'xyz':
        atoms = ase.io.read(file_name)
    elif informat == 'out':
        atoms = ase.io.read(file_name, format = 'gaussian-out', index=-1)
        
    atomSymbols = np.array(atoms.get_chemical_symbols())
    atomList = np.unique(atomSymbols, return_counts=True)
    
    if 'Al' not in atomSymbols: # no aluminum
        return 'N/A'
    elif 'Al' in atomSymbols and atomList[1][0] == 1: # single aluminum
        return 'N/A'
    else:
        bonds = utils.get_bonds(atoms) # more than 1 aluminum

        Al_indices = [atom.index for atom in atoms if atom.symbol == 'Al']
        Al_first = Al_indices[0]
        Al_others = Al_indices[1:]

        G = nx.Graph()

        G.add_nodes_from(np.unique(bonds))
        G.add_edges_from(bonds)

        pathLength = 0
        for Al_target in Al_others:
            paths = k_shortest_paths(G, Al_first, Al_target, 2)
            pathLength += (len(paths[0])-1)**1.5 + (len(paths[1])-1)

        pathLength /= len(Al_others)

        return pathLength

    
