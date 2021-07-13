import click
import pandas as pd
import sys
import numpy as np
from periodictable import elements
#from canela import __version__

@click.command(name='exgeom',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
#@click.version_option(__version__)
@click.argument('file_name', type=str) # it must include a file format
file_name_input = file_name.copy()
out_name_default = file_name_input.rsplit(".",1)[0] + ".geom.com"
@click.argument('out_name', type=str, default=out_name_default)

#def exgeom(file_name, file_format):
def exgeom(file_name):
    """Extract the last/optimized geometry from Gaussian output file
    and generate Gaussian input fle (.com or .gjf) or xyz file.
    
    The first argument is the file name of a Gaussian output file including ".out"
    The second argument is output file name including desired format (optional)
    Default output name and format = input file name + "geom.com"
    Currently ".com" and ".xyz" are supported
    """

    outformat = out_name.rsplit(".",1)[-1]
    
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
        
    no_elem = len(df_geom)

    if outformat == "com":
        output = open(out_name, 'w')
        output.write(routeStr + '\n\n')
        output.write(out_name + '\n\n')
        output.write(charAndMult + '\n')
        output.write(df_geom.to_string(index=False, header=False))
        output.write('\n\n')
        output.close()
    elif outformat == "xyz":
        output = open(out_name, 'w')
        output.write(no_elem + '\n')
        output.write(out_name + '\n')
        output.write(df_geom.to_string(index=False, header=False))
        output.write('\n')
        output.close()

## Extract the final geometry from the Gaussian output file as an input file
## Arguments = file names you want to extract the geometry
## The default name of the extracted geometry file is "file name of output".geom.com
## If you want to specify the file name, see the following example (valid only with a single file)
## Examples:
#### xgeom methane.out
#### xgeom methane.out methanol.out ethane.out
#### xgeom *.out
#### xgeom methane.out - methane.geometry.com


