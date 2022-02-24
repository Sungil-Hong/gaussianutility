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

        for idx, line in enumerate(lines[idx_charge:]):
            if not "Charge" in line:
                idx_charge_end = idx_charge+idx
                break

        charAndMult = []

        for line in lines[idx_charge: idx_charge_end]:
                charAndMult.append(str(line.split()[2]))
                charAndMult.append(str(line.split()[5]))

        charAndMult = " ".join(charAndMult)

   # Read oniom layer data

    if oniom == 1:

        iniGeom = []
        for line in lines[idx_charge_end:]:
            if len(line) <= 2: break
            iniGeom.append(line.split())

        iniGeom = list(filter(None, iniGeom))

        oniom_layer = []
        for line in iniGeom:
            oniom_layer.append(line[5])

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
    
#=======================================================================
def freeze_low(file_name):
    """Freeze Cartesian coordinates of atoms in the low ONIOM layer.
    Only .com and .gjf file is supported.
    """
    informat = file_name.rsplit(".",1)[-1]
    if not informat == 'com' or informat == 'gjf':
        raise TypeError('The input file format must be Gaussian input (com or gjf)')
        
    lines = open(file_name).readlines()

    for idx, line in enumerate(lines):
        if "#" in line:
            idx_route = idx
            route = line
            break

    if len(lines[idx_route+1]) >=2:
        route = route.replace('\n', ' ')
        route += lines[idx_route+1]
        idx_route += 1

    if not ("oniom" or "ONIOM" or "Oniom") in route:
        raise AttributeError("Must have ONIOM formulation")

    if not "modredundant" in route:
        lines[idx_route] = lines[idx_route].replace('opt', 'opt=modredundant')

    connect = 0
    if "geom=connectivity" in route:
        connect += 1

    geom = []

    for idx2, line in enumerate(lines[idx_route+5:]):
        geom.append(line.split())
        geomEndIdx = idx_route+idx2+5
        if len(line) <= 2: break

    # if connect == 1:
    #     for id3, line in enumerate(lines[])
    geom = list(filter(None, geom))

    lineLen = 0
    for line in geom:
        if len(line) > lineLen: lineLen = len(line)

    clmnName = []
    for i in range(lineLen):
        name = "Column {}".format(i)
        clmnName.append(name)

    clmnName[5] = 'ONIOM layer'
    df_geom = pd.DataFrame(geom, columns= clmnName)
    oniom_layer = pd.DataFrame(df_geom, columns = ['ONIOM layer'])

    Low_level = np.array(oniom_layer.index[oniom_layer['ONIOM layer']=='L'])+1

    if connect == 0:
        editfile = open(file_name, 'w')
        for line in lines[:geomEndIdx]:
            editfile.write(line)
        editfile.write('\n')
        for elem in Low_level:
            editfile.write('X ' + str(elem) + ' F\n')
        editfile.write('\n')
        editfile.close()

    if connect == 1:
        for idx, line in enumerate(lines[geomEndIdx+1:]):
            if len(line) <= 2:
                connectEndIdx = geomEndIdx+idx+1
                break

        editfile = open(file_name, 'w')
        for line in lines[:connectEndIdx+1]:
            editfile.write(line)
        for elem in Low_level:
            editfile.write('X ' + str(elem) + ' F\n')
        editfile.write('\n')
        editfile.close()



#=======================================================================
def ONIOM_sort(file_name, sort_idx = 0):
    
    informat = file_name.rsplit(".",1)[-1]
    if not informat == 'com' or informat == 'gjf':
        raise TypeError('The input file format must be Gaussian input (com or gjf)')
 
    inputfile = open(file_name,'r')
    
    if not sort_idx:
        sort_idx = 'Atom'
    if not sort_idx in ['Atom','x','y','z']:
        raise ValueError(" Sort index must be \'x\', \'y\', \'z\', or \'Atom\' ")

    lines = open(file_name).readlines()
    for idx, line in enumerate(lines):
        if "#" in line:
            idx_route = idx
            route = line
            break
            
    if len(lines[idx_route+1]) >=2:
        route = route.replace('\n', ' ')
        route += lines[idx_route+1]
        idx_route += 1
        
    title_and_spin = lines[idx_route+2:idx_route+5]

    if not ("oniom" or "ONIOM" or "Oniom") in route:
        raise AttributeError("Must have ONIOM formulation")        
        
    if "geom=connectivity" in route:
        route = route.replace("geom=connectivity","")
    
    geom=[]
    for idx2, line in enumerate(lines[idx_route+5:]):
        geom.append(line.split())
        geomEndIdx = idx_route+idx2+5
        if len(line) <= 2: break
     
    lineLen = 0
    for line in geom:
        if len(line) > lineLen: lineLen = len(line)

    if lineLen == 5:
        df_geom = pd.DataFrame(geom, columns = ['Atom','x','y','z','ONIOM_layer'])
    elif lineLen > 5:
        clmnName = []
        for i in range(lineLen):
            name = "C{}".format(i)
            clmnName.append(name)

        df_geom = pd.DataFrame(geom, columns = clmnName)
        df_geom = df_geom.iloc[ : ,0:6]
        df_geom = df_geom.rename({'C0':'Atom', 'C1':'Index', 'C2':'x', 'C3':'y', 'C4':'z', 'C5':'ONIOM_layer'}, axis='columns')

    df_layer_sort = pd.DataFrame({'ONIOM_layer': ['H','M','L']})
    layer_sort_mapping = df_layer_sort.reset_index().set_index('ONIOM_layer')

    df_geom['ONIOM_layer_num'] = df_geom['ONIOM_layer'].map(layer_sort_mapping['index'])

    if sort_idx == 'Atom':
        elem = []
        for el in elements:
            elem.append(el.symbol)

        df_atom_sort = pd.DataFrame(elem[1:], columns = ['symbol'])
        atom_sort_mapping = df_atom_sort.reset_index().set_index('symbol')

        df_geom['Atomic_num'] = df_geom['Atom'].map(atom_sort_mapping['index'])
        df_geom = df_geom.sort_values(by=['ONIOM_layer_num', 'Atomic_num'], ascending = [1,0])
        df_geom = df_geom.drop(columns=['ONIOM_layer_num', 'Atomic_num'])
    else:
        df_geom = df_geom.sort_values(by=['ONIOM_layer_num', sort_idx])
        df_geom = df_geom.drop(columns=['ONIOM_layer_num'])

    outfile = open(file_name, 'w')
    outfile.write(route)
    outfile.write('\n')
    outfile.write(''.join(title_and_spin))
    outfile.write(df_geom.to_string(index=False, header=False))
    outfile.write('\n\n')
    outfile.close()
