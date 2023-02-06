import networkx as nx
import featpro.utils as utils
from featpro.pos import connectivity
from itertools import islice

#=======================================================================
def k_shortest_paths(G, source, target, k, weight=None):
    return list(
        islice(nx.shortest_simple_paths(G, source, target, weight=weight), k)
    )
# Sourced from https://networkx.org/documentation/stable/index.html       


#=======================================================================
def path_feature(atoms):
    """
    Calculate a path feature, which is defined as a weighted average of two shortest paths of every Al-Al pair.
    This feature captures the different positioning or distribution of Al in structures.
    More weight is given to the shorter path.
    
    The input is ASE.atoms object
    
    Concept credit to Michael Cowan (mcowan92@gmail.com)
    """
    
    atoms = atoms.copy()
    
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

#=======================================================================
def feat_gen(file_name):
    """Print features and target for ML aplication in the order of
    {# of Al, # of H, # of O, # of Si, # of water produced, path feature, connectivity feature,
    formation Gibbs free E}
    from Gaussian output file (.out).
    
    If Gaussian input (.com or .gjf) or XYZ file is used, the Gibbs free energy will be missing.
    If the structure does not contain more than 1 Al, path feature returns N/A.
    """

    atoms = gen_atoms(file_name)
    informat = file_name.rsplit(".",1)[-1]
    return_var = np.array([file_name.rsplit(".",1)[:-1]], dtype=object)
    
    # Read number of atoms in the order of [Al, H, O, Si]
    atomSymbols = np.array(atoms.get_chemical_symbols())
    atomList = np.unique(atomSymbols, return_counts=True)
    
    if 'Al' not in atomList[0]:
        atomCount = np.concatenate(([0],atomList[1]))
    else:
        atomCount = atomList[1]

    # Number of produced water produced by condensation of monomers
    # This feature is intended to distinguish structue types
    no_water = (atomCount[0] + atomCount[-1])*4 - atomCount[-2]

    return_var = np.append(return_var, atomCount) # # of atoms
    return_var = np.append(return_var, no_water)  # # of produced water from condensation
    
    # Path feature
    pathval = path_feature(atoms)
    return_var = np.append(return_var, pathval) # path feature
    
    # Connectivity feature
    connectval = connectivity(atoms)
    return_var = np.append(return_var, connectval) # connectivity feature
    
    if informat == 'out':
        # Set target = formation energy
        # Read the Gibbs free energy from the Gaussian output file
            
        G_Si = -592.968165   # G of Si(OH)4 in hartree @ 100'C
        G_Al = -546.445      # G of Al(OH)3(H2O) in hartree @ 100'C
        G_water = -76.421375 # G of water in hartree @ 100'C
        
        lines = open(file_name,'r')
        for line in lines:
            if "Sum of electronic and thermal Free Energies" in line:
                gibbsE = float(line.split()[-1])

        # Normalize with the number of monomers and convert the unit from hartree to kJ/mol
        gibbsE -= G_Al*atomCount[0] + G_Si*atomCount[-1] - G_water*no_water
        gibbsE /= (atomCount[0] + atomCount[-1])/2625.4996394799 
        
    else: gibbsE = 'N/A'

    return_var = np.append(return_var, gibbsE) # formation Gibbs free E
    return_var_str = ', '.join(map(str,return_var))

    return return_var_str
    
    