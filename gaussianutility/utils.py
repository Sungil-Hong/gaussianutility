#! usr/bin/env python

import numpy as np
import pandas as pd
import ase.io
import ase
from ase import atoms
from periodictable import elements

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
def read_input(file_name):
    """Return input file format, route section, title (including spin/multiplicity), and geometry.
    'route' and 'title_and_spin' are string.
    If input is in ONIOM format, geom is pandas dataframe; otherwise, list
    """
    informat = file_name.rsplit(".",1)[-1]

    inputfile = open(file_name,'r')
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

    geom=[]
    for idx2, line in enumerate(lines[idx_route+5:]):
        geom.append(line.split())
        geomEndIdx = idx_route+idx2+5
        if len(line) <= 2: break
    geom = list(filter(None, geom))
    
    rest = []
    for idx3, line in enumerate(lines[geomEndIdx:]):
        rest.append(line)
 
    if ("oniom" or "ONIOM" or "Oniom") in route: 
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
    
    else:
        lineLen = 0
        for line in geom:
            if len(line) > lineLen: lineLen = len(line)
        if lineLen != 4:
            raise ValueError("Check input file")

        df_geom = pd.DataFrame(geom, columns = ['Atom','x','y','z'])

    return informat, route, title_and_spin, df_geom, rest

#=======================================================================
def input_sort(file_name, sort_idx = 0):
    
    if not sort_idx:
        sort_idx = 'Atom'
    if not sort_idx in ['Atom','x','y','z']:
        raise ValueError(" Sort index must be \'x\', \'y\', \'z\', or \'Atom\' ")  
        
    informat, route, title_and_spin, df_geom, rest = read_input(file_name)
    
    if not informat == 'com' or informat == 'gjf':
        raise TypeError("The input file format must be Gaussian input (com or gjf)") 
        
    oniom_idx = 0
    if ("oniom" or "ONIOM" or "Oniom") in route:
        oniom_idx = 1
     
    if "geom=connectivity" in route:
        route = route.replace("geom=connectivity","")
    
    if oniom_idx == 1:
        df_layer_sort = pd.DataFrame({'ONIOM_layer': ['L','M','H',]})
        layer_sort_mapping = df_layer_sort.reset_index().set_index('ONIOM_layer')
        df_geom['ONIOM_layer_num'] = df_geom['ONIOM_layer'].map(layer_sort_mapping['index'])

    if sort_idx == 'Atom':
        elem = []
        for el in elements:
            elem.append(el.symbol)

        df_atom_sort = pd.DataFrame(elem[1:], columns = ['symbol'])
        atom_sort_mapping = df_atom_sort.reset_index().set_index('symbol')
        df_geom['Atomic_num'] = df_geom['Atom'].map(atom_sort_mapping['index'])
        
        if oniom_idx == 0:
            df_geom = df_geom.sort_values(by=['Atomic_num'], ascending = 0)
            df_geom = df_geom.drop(columns=['Atomic_num'])
            
        elif oniom_idx == 1:
            df_geom = df_geom.sort_values(by=['ONIOM_layer_num', 'Atomic_num'], ascending = [1,0])
            df_geom = df_geom.drop(columns=['ONIOM_layer_num', 'Atomic_num'])
            
    else:
        if oniom_idx == 0:
            df_geom = df_geom.sort_values(by=sort_idx)
            
        elif oniom_idx == 1:
            df_geom = df_geom.sort_values(by=['ONIOM_layer_num', sort_idx])
            df_geom = df_geom.drop(columns=['ONIOM_layer_num'])
    
    
        df_geom = df_geom[['Atom', 'Index', 'x', 'y', 'z', 'ONIOM_layer']]
        
    
    outfile = open(file_name, 'w')
    outfile.write(route)
    outfile.write('\n')
    outfile.write(''.join(title_and_spin))
    outfile.write(df_geom.to_string(index=False, header=False))
    outfile.write('\n\n')
    outfile.close()

#=======================================================================
def freeze_layer(file_name, layer_idx = 0):
    if not layer_idx:
        layer_idx = 'L'
        
    layers = list(layer_idx)

    if ('H' not in layers) and  ('M' not in layers) and ('L' not in layers):
        raise ValueError("Layer indext should be composed of \'H\', \'M\', \'L\'")
 
    informat, route, title_and_spin, df_geom, rest = read_input(file_name)

    if not informat == 'com' or informat == 'gjf':
        raise TypeError('The input file format must be Gaussian input (com or gjf)')

    if not ("oniom" or "ONIOM" or "Oniom") in route:
        raise AttributeError("Must have ONIOM formulation")

    for layer in layers:
        df_geom.loc[df_geom["ONIOM_layer"] == layer, "Index"] = -1

    outfile = open(file_name, 'w')
    outfile.write(route)
    outfile.write('\n')
    outfile.write(''.join(title_and_spin))
    outfile.write(df_geom.to_string(index=False, header=False))
    outfile.write('\n')
    outfile.write(''.join(rest))
    outfile.close()

