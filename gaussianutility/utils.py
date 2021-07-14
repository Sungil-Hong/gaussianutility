#! usr/bin/env python

import numpy as np
import pandas as pd
import ase.io
import ase
from ase import atoms
import networkx as nx
from periodictable import elements
import featpro.utils as utils
from featpro.pos import connectivity
from itertools import islice

#=======================================================================
def extract_geom(file_name, out_file = 0):

    informat = file_name.rsplit(".",1)[-1]
    if informat != 'out':
        raise TypeError('The input file format must be .out')

    if not out_file:
        out_file = file_name.rsplit(".",1)[0] + ".geom.com"
        
    outformat = out_file.rsplit(".",1)[-1]
    if outformat !='com' and outformat !='gjf' and outformat !='xyz':
        raise TypeError('The output format must be .com, .gjf, or .xyz')
        
    # open the input file and read a route section

    lines = open(file_name).readlines()

    for idx, line in enumerate(lines):
        if "#" in line:
            idx_route = idx
            break

    route = []
    for line in lines[idx_route:]:
        route.append(line.strip('\n'))
        if "----------------------------" in line:
            break

    route = route[:-1]
    for idx in range(len(route)):
        if route[idx].startswith(' '):
            route[idx] = route[idx][1:]

    routeStr = ""
    routeStr = routeStr.join(route)
    if "geom=connectivity" in routeStr:
        routeStr = routeStr.replace("geom=connectivity","")

    oniom = 0
    if ("oniom" or "ONIOM" or "Oniom") in routeStr:
        oniom += 1

    # Read charge and multiplicity data

    for idx, line in enumerate(lines):
        if "Charge =" in line:
            idx_charge = idx
            break   

    if oniom == 0:
        charge, multiplicity = lines[idx_charge].split()[2], lines[idx_charge].split()[5]
        charAndMult = "{} {}".format(charge, multiplicity)

    if oniom == 1:
        chargeL, multiplicityL = lines[idx_charge].split()[2], lines[idx_charge].split()[5]
        chargeM, multiplicityM = lines[idx_charge+1].split()[2], lines[idx_charge+1].split()[5]
        chargeH, multiplicityH = lines[idx_charge+2].split()[2], lines[idx_charge+2].split()[5]

        charAndMult = "{} {} {} {} {} {}".format(chargeL, multiplicityL, chargeM, multiplicityM, chargeH, multiplicityH)    

    # Read oniom layer data

    iniGeom = []

    if oniom == 1:
        for line in lines[idx_charge+3:]:
            if len(line) <= 3 : break
            iniGeom.append(line.split())

        iniGeom = list(filter(None, iniGeom))
         
        shape_index = 0
        length = len(iniGeom[0])
        for i in iniGeom:
            if len(i) != length:
                shape_index += 1
                break
    
        if shape_index == 0:
            df_iniGeom = pd.DataFrame(iniGeom, columns = ['element','X','Y','Z','oniom level'])
        else:
            df_iniGeom = pd.DataFrame(iniGeom, columns = ['element','atomic type','X','Y','Z','oniom level','connected atom','level of connected atom','val1','val2'])
    
        oniom_layer = pd.DataFrame(df_iniGeom, columns = ['oniom level'])

    # Read final geometry, convert format, and write a new file

    for idx, line in enumerate(reversed(lines)):
        if "Coordinates (Angstroms)" in line:
            break

    idx_geo = len(lines)-idx-1
    geom = []

    for line in lines[idx_geo+3:]:
        geom.append(line.split())
        if "----------------------------" in line:
            break

    df_geom = pd.DataFrame(geom[:-1], columns = ['index','atomic number','atomic type','X','Y','Z'])
    df_geom = df_geom.drop(columns=['index', 'atomic type'])

    elem = []
    for el in elements:
        elem.append([el.number, el.symbol])

    df_elem = pd.DataFrame(elem, columns = ['atomic number', 'symbol'])

    atomNo = df_geom['atomic number']
    atomSb = []
    for i in atomNo:
        atomSb.append(df_elem.loc[df_elem['atomic number'] == int(i)]['symbol'].item())

    df_geom = df_geom.drop(columns='atomic number')
    df_geom.insert(0, 'symbol', atomSb)

    if oniom == 1:
        df_geom['oniom level'] = oniom_layer
        
    no_elem = str(len(df_geom))

    if outformat == "com" or outformat == "gjf":
        output = open(out_file, 'w')
        output.write(routeStr + '\n\n')
        output.write(out_file + '\n\n')
        output.write(charAndMult + '\n')
        output.write(df_geom.to_string(index=False, header=False))
        output.write('\n\n')
        output.close()
    elif outformat == "xyz":
        output = open(out_file, 'w')
        output.write(no_elem + '\n')
        output.write(out_file + '\n')
        output.write(df_geom.to_string(index=False, header=False))
        output.write('\n')
        output.close()

    
#=======================================================================    
def gen_atoms(file_name):
    """Read Gaussian input (.com or .gjf), XYZ file, or Gaussian output (.out)
    and return ASE.atoms object
    """
    
    informat = file_name.rsplit(".",1)[-1]
    
    if informat == 'com' or informat == 'gjf' or informat == 'xyz':
        atoms = ase.io.read(file_name)
    elif informat == 'out':
        try:
            atoms = ase.io.read(file_name, format = 'gaussian-out', index=-1)
        except StopIteration: 
            # Sometimes, ase.io.read does not read a gaussian output file correctly.
            # Generate an atom object manually.
            lines = open(file_name).readlines()
            for idx, line in enumerate(reversed(lines)):
                if "Coordinates (Angstroms)" in line:
                    break  
            idx_geo = len(lines)-idx-1
            geom = []

            for line in lines[idx_geo+3:]:
                geom.append(line.split())
                if "----------------------------" in line:
                    break

            df_geom = pd.DataFrame(geom[:-1], columns = ['index','atomic number','atomic type','X','Y','Z'])

            elem = []
            for el in elements:
                elem.append([el.number, el.symbol])

            df_elem = pd.DataFrame(elem, columns = ['atomic number', 'symbol'])

            atomNo = df_geom['atomic number']
            atomSb = []
            for i in atomNo:
                atomSb.append(df_elem.loc[df_elem['atomic number'] == int(i)]['symbol'].item())

            df_geom = df_geom.drop(columns=['index', 'atomic number', 'atomic type'])
            atomPs = pd.DataFrame.to_numpy(df_geom)
            atoms = ase.Atoms(symbols = atomSb, positions=atomPs)
            
    else: raise TypeError('The input file format must be .out, .com, .gjf, or .xyz')
        
    return atoms


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

    return_var = np.append(return_var, atomCount)
    return_var = np.append(return_var, no_water)
    # (# of Al, # of H, # of O, # of Si, # of water produced)
    
    # Path feature
    pathval = path_feature(atoms)
    return_var = np.append(return_var, pathval)
    # (# of Al, # of H, # of O, # of Si, # of water produced, path feature)
    
    # Connectivity feature
    connectval = connectivity(atoms)
    return_var = np.append(return_var, connectval)
    # (# of Al, # of H, # of O, # of Si, # of water produced, path feature, connectivity feature)
    
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

    return_var = np.append(return_var, gibbsE)
    # (# of Al, # of H, # of O, # of Si, # of water produced, path feature, connectivity feature, formation Gibbs free E)

    return_var_str = ', '.join(map(str,return_var))

    return return_var_str
    
