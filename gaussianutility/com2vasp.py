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
    return parser.parse_args()

def com_2_vasp(file_name):
    # Read geometry from a Gaussian input file
    _, _, _, df_geom, _ = readinput(file_name)
    df_geom = df_geom[['Atom', 'x', 'y', 'z']]
    df_geom['Atom'] = np.array([i.rsplit("-",1)[0] for i in list(df_geom['Atom'])])
    if not 'Tv' in df_geom['Atom'].values:
        raise ValueError('The Gaussian input file does not contain lattice information')

    # Convert Gaussian input to VASP input
    elem = [el.symbol for el in elements]
    df_atom_sort = pd.DataFrame(elem[1:], columns = ['symbol'])
    atom_sort_mapping = df_atom_sort.reset_index().set_index('symbol')
    df_geom['Atomic_num'] = df_geom['Atom'].map(atom_sort_mapping['index'])

    df_geom = df_geom.sort_values(by=['Atomic_num'], ascending = 0)
    df_geom = df_geom.drop(columns=['Atomic_num'])

    # Read lattice information
    lattice = df_geom[-3:].drop(['Atom'], axis=1) 

    df_geom = df_geom.drop(df_geom[df_geom['Atom'] == 'Tv'].index)
    df_geom = df_geom.drop(df_geom[df_geom['Atom'] == '?s'].index)

    atomlist = np.array(df_geom['Atom'])
    coord = df_geom[['x','y','z']]

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
    with open(out_file, 'w') as output:
        output.write(name + '\n')
        output.write('1.0\n')
        output.write(lattice.to_csv(index=False, header=False, sep='\t'))
        output.write('   '+' '.join(uniqatomlist)+'\n')
        output.write('   '+' '.join(uniqatomcount)+'\n')
        output.write('Cartesian\n')
        output.write(coord.to_csv(index=False, header=False, sep='\t'))
        output.write('\n')

def main():
    args = parse_args()
    com_2_vasp(args.file_name)

if __name__ == "__main__":
    main()


