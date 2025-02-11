#!/usr/bin/env python3

import pandas as pd
import numpy as np
from periodictable import elements
import argparse
from argparse import RawTextHelpFormatter
from gaussianutility.utilities import readinput

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Convert Gaussian input file to .vasp file\n\n"
                     "Return file_name.vasp: VASP geometry file (POSCAR)",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument("file_name", help='Gaussian input file (.com or .gjf)')
    parser.add_argument('-c', '--ctype', \
    help='Coordination system type; d for direct and c for cartesian', default='d')
    
    return parser.parse_args()

def com_2_vasp(file_name, ctype):
    if ctype not in ('d', 'c'):
        raise ValueError('The coordination type must be d for direct or c for cartesian')

    # Read geometry from a Gaussian input file
    _, _, _, df_geom, _ = readinput(file_name)

    if 'Index' in df_geom.columns:
        indices = df_geom['Index']
        df_geom = df_geom[['Atom', 'x', 'y', 'z']]
    else:
        df_geom = df_geom[['Atom', 'x', 'y', 'z']]

    if not 'Tv' in df_geom['Atom'].values:
        raise ValueError('The provided Gaussian input file does not contain lattice information')

    lattice = df_geom[-3:].drop(['Atom'], axis=1)
    df_geom = df_geom.drop(df_geom[df_geom['Atom'] == 'Tv'].index)
    atomlist = np.array(df_geom['Atom'])
        
    if ctype == 'd':
        cart_coords = np.array(df_geom[['x','y','z']], dtype=float)
        lattice_inv = np.linalg.inv(np.array(lattice, dtype=float))
        direct_coords = np.dot(cart_coords, lattice_inv)

        df_geom['x'] = direct_coords[:,0]
        df_geom['y'] = direct_coords[:,1]
        df_geom['z'] = direct_coords[:,2]

    if 'indices' in list(locals()):
        selective = []
        for idx in indices:
            if idx == '0':
                selective.append(['T', 'T', 'T'])
            elif idx == '-1':
                selective.append(['F', 'F', 'F'])

        selective = np.array(selective)
        df_geom['selective_x'] = selective[:,0]
        df_geom['selective_y'] = selective[:,1]
        df_geom['selective_z'] = selective[:,2]

    uniqatomlist = []
    uniqatomcount = []

    for elem in atomlist:
        if elem not in uniqatomlist:
            uniqatomlist.append(elem)
            uniqatomcount.append(np.count_nonzero(atomlist == elem))

    uniqatomlist = np.array(uniqatomlist)
    uniqatomcount= np.array(uniqatomcount, dtype=str)

    # Write VASP input file
    name = ".".join(file_name.rsplit(".",1)[0:-1])
    out_file = name + ".vasp"
    df_geom = df_geom.drop(['Atom'], axis=1)

    with open(out_file, 'w') as output:
        output.write(name + '\n')
        output.write('1.0\n')
        output.write(lattice.to_csv(index=False, header=False, sep='\t'))
        output.write('\t'+'\t'.join(uniqatomlist)+'\n')
        output.write('\t'+'\t'.join(uniqatomcount)+'\n')
        
        if 'indices' in list(locals()):
            output.write('Selective dynamics\n')

        if ctype == 'd':
            output.write('Direct\n')
        elif ctype == 'c':
            output.write('Cartesian\n')
            
        output.write(df_geom.to_csv(index=False, header=False, sep='\t'))
        output.write('\n')
        

def main():
    args = parse_args()
    com_2_vasp(args.file_name, args.ctype)

if __name__ == "__main__":
    main()


