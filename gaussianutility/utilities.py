#!/usr/bin/env python3

import numpy as np
import pandas as pd
from periodictable import elements

"""

This file contains functions that can be used in the other scripts.
The utilities herein are not used by itself but inside another script.

"""

def readinput(file_name):
    """
    This reads an Gaussian inputfile and decomposes it into multiple elements 
    as a return for the use in other scripts.
    Return route section, title, charge and multiplicity, and and geometry.
    Geometry is a pandas dataframe and all the others are strings
    """
    # Check input file
    name, input_format = file_name.rsplit(".", 1)
    if input_format not in ('com', 'gjf'):
        raise TypeError('The input file format must be .com or .gjf')

    # Find the elements of the Gaussian input file
    with open(file_name, 'r') as inputfile:
        lines = inputfile.readlines()

    # Route section
    for idx, line in enumerate(lines):
        if "#" in line:
            idx_route = idx
            route = line
            break

    if len(lines[idx_route+1]) >=2:
        route = route.replace('\n', ' ')
        route += lines[idx_route+1]
        idx_route += 1

    connectIdx = "geom=connectivity" in route
            
    # Title and charge/multiplicity
    title = lines[idx_route+2]
    charge_mult = lines[idx_route+4]

    # Geometry
    geom=[]
    for idx, line in enumerate(lines[idx_route+5:]):
        geom.append(line.split())
        geomEndIdx = idx_route+idx+5
        if len(line) <= 2: break
        
    geom = list(filter(None, geom))
    
    connectivity = []
    if connectIdx:
        for idx, line in enumerate(lines[geomEndIdx:]):
            if line.startswith("1") or line.startswith(" 1"):
                break
        
        connectIdx = geomEndIdx + idx
        for line in lines[connectIdx:]:
            connectivity.append(line)
            if line.startswith(str(len(geom))) or line.startswith(" "+str(len(geom))):
                break

    if any(keyword in route.lower() for keyword in ["oniom"]):
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
            
            if np.sum(np.mod(np.array(df_geom.iloc[:,1], dtype=float),1)) == 0:
                df_geom = df_geom.iloc[:,0:6]
                df_geom = df_geom.rename({'C0':'Atom', 'C1':'Index', 'C2':'x', 'C3':'y', 'C4':'z', 'C5':'ONIOM_layer'}, axis='columns')
            else:
                df_geom = df_geom.iloc[:,0:5]
                df_geom = df_geom.rename({'C0':'Atom', 'C1':'x', 'C2':'y', 'C3':'z', 'C4':'ONIOM_layer'}, axis='columns')
                
    else:
        lineLen = 0
        for line in geom:
            if len(line) > lineLen: lineLen = len(line)
        if lineLen == 4:
            df_geom = pd.DataFrame(geom, columns = ['Atom','x','y','z'])
        elif lineLen > 4:
            clmnName = []
            for i in range(lineLen):
                name = "C{}".format(i)
                clmnName.append(name)

            df_geom = pd.DataFrame(geom, columns = clmnName)
            
            if np.sum(np.mod(np.array(df_geom.iloc[:,1], dtype=float),1)) == 0:
                df_geom = df_geom.iloc[:,0:5]
                df_geom = df_geom.rename({'C0':'Atom', 'C1':'Index', 'C2':'x', 'C3':'y', 'C4':'z'}, axis='columns')
            else:
                df_geom = df_geom.iloc[:,0:4]
                df_geom = df_geom.rename({'C0':'Atom', 'C1':'x', 'C2':'y', 'C3':'z'}, axis='columns')            
                
    df_geom['Atom'] = np.array([atomSb.split('-')[0] for atomSb in list(df_geom['Atom'])])
    
    return route, title, charge_mult, df_geom, connectivity


def readoutput(file_name, stepIdx = -1):
    """
    This reads an Gaussian outputfile and decomposes it into multiple elements 
    as a return for the use in other scripts.
    Return route section, title, spin and multiplicity, and and geometry.
    Geometry is a pandas dataframe and all the others are strings
    """
    # Check input file
    name, input_format = file_name.rsplit(".", 1)
    if input_format != 'out':
        raise TypeError('The input file format must be .out')

    # Read file
    with open(file_name, 'r') as inputfile:
        lines = inputfile.readlines()

    # Read computational chemistry method and remove unnecessary keywords 
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
    
    if 'qst' in routeStr:
        routeStr = routeStr.replace('=qst3', '').replace('qst3', '')
        routeStr = routeStr.replace('=qst2', '').replace('qst2', '')

    if "geom=connectivity" in routeStr:
        routeStr = routeStr.replace("geom=connectivity","")

    oniom = any(word in routeStr.lower() for word in ["oniom"])

    # Read charge and multiplicity information
    for idx, line in enumerate(lines):
        if "Charge =" in line:
            idx_charge = idx
            idx_charge_end = idx+1
            break

    if oniom:
        for idx, line in enumerate(lines[idx_charge:]):
            if not "Charge" in line:
                idx_charge_end = idx_charge+idx
                break

        charge_mult = []
        for line in lines[idx_charge: idx_charge_end]:
            charge_mult.append(str(line.split()[2]))
            charge_mult.append(str(line.split()[5]))

        charge_mult = " ".join(charge_mult)

    else:
        charge, multiplicity = lines[idx_charge].split()[2], lines[idx_charge].split()[5]
        charge_mult = "{} {}".format(charge, multiplicity)

    # Read oniom layer data
    iniGeom = []
    for line in lines[idx_charge_end:]:
        if len(line) <= 2: break
        iniGeom.append(line.split())

    iniGeom = list(filter(None, iniGeom))

    if oniom:
        oniomIdx = 0
        if len(iniGeom[0]) == 5:
            oniom_layer = [line[4] for line in iniGeom]
        elif len(iniGeom[0]) > 5:
            oniom_layer = [line[5] for line in iniGeom]

    # Read atom index data (0: optimized; -1:frozen)
    try:
        if isinstance(int(iniGeom[0][1]),int):
            indexFlag = True
    except ValueError:
        indexFlag = False
        
    if indexFlag:
        indices = [line[1] for line in iniGeom]

    # Read final geometry
    if stepIdx == "L":
        energies = []
        if oniom:
            for idx, line in enumerate(lines):
                if "ONIOM: extrapolated energy" in line:
                    energies.append(float(line.split()[-1]))
                    
            if energies[-1] == np.min(energies): stepIdx = -1
            else: stepIdx = np.argmin(np.array(energies))
            
        else:
            for idx, line in enumerate(lines):
                if "SCF Done:  " in line:
                    energies.append(float(line.split()[4]))
                    
            if energies[-1] == np.min(energies): stepIdx = -1
            else: stepIdx = np.argmin(np.array(energies))

    else: stepIdx = int(stepIdx)
    
    optLineNo = []
    standardOri = False
    for idx, line in enumerate(lines):
        if "Standard orientation" in line: 
            standardOri = True
            break
            
    for idx, line in enumerate(lines):        
        if standardOri:
            if "Standard orientation" in line:
                optLineNo.append(idx+1)
        else:
            if "Input orientation" in line:
                optLineNo.append(idx+1)

    if any(word in routeStr for word in ["freq","Freq","FREQ"]):
        optLineNo.pop(-1)

    geom = []
    for line in lines[optLineNo[stepIdx]+4:]:
        geom.append(line.split())
        if "----------------------------" in line:
            break

    # Convert geometry information to Pandas dataframe
    df_geom = pd.DataFrame(geom[:-1], columns = ['index','atomic number','atomic type','x','y','z'])
    df_geom = df_geom.drop(columns=['index', 'atomic type'])

    if oniom:
        df_geom['ONIOM_layer'] = oniom_layer
        
    if indexFlag:
        df_geom.insert(1, 'Index', indices)
        
    elem = np.array([[el.number, el.symbol] for el in elements])
    elem = {int(elem[i,0]): elem[i,1] for i in range(len(elem))}
    atomSb = [elem[atomNo] for atomNo in np.array(df_geom['atomic number'], dtype=int)]

    df_geom = df_geom.drop(columns='atomic number')
    df_geom.insert(0, 'Atom', atomSb)

    return routeStr, charge_mult, df_geom
    
