#!/usr/bin/env python

import numpy as np
import itertools
import mendeleev as md
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Convert vasp input file to Gaussian input file (.com)\n"
                     "The computational method in the Gaussian calculation must be modified manually if needed\n\n"
                     "Return file_name.com: Gaussian input file",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument("file_name", help='VASP input file (.vasp)')
    return parser.parse_args()

def vasp_2_com(file_name):
    # Read geometry information and make Pandas dataframe
    lines = open(file_name).readlines()

    elem_type = lines[5].split()
    elem_num = np.array(list(map(int, lines[6].split())))
    elect_num = np.array([md.element(e).electrons for e in elem_type])

    multiplicity = 1 if sum(elect_num * elem_num) % 2 == 0 else 2

    elem_list = np.array(list(itertools.chain.from_iterable([ [e] * n for e, n in zip(elem_type, elem_num)])))
    
    geom = np.array([line.split() for line in lines[8:] if len(line) >= 2], dtype='float_')

    lattice_scailing_factor = float(lines[1].split()[0])
    lattice_vector = [line.split() for line in lines[2:5]]
    lattice_vector = np.array(lattice_vector, dtype = float)

    geomMod = lines[7].split()[0]
    if geomMod.lower() == 'cartesian':
        geom *= lattice_scailing_factor
    elif geomMod.lower() == 'direct':
        geom = np.dot(lattice_vector, geom.T).T * lattice_scailing_factor
        
    df_geom = pd.DataFrame(geom, columns = ['x','y','z'])
    df_geom.insert(loc=0, column='element', value=elem_list)

    # Write .com file
    out_file = file_name.rsplit(".",1)[0] + ".com"
    with open(out_file, 'w') as output:
        output.write('# pbepbe/3-21g/auto\n\n')
        output.write(f'{file_name.rsplit(".", 1)[0]}\n\n')
        output.write(f'0 {multiplicity}\n')
        output.write(df_geom.to_string(index=False, header=False) + '\n')
        output.write('Tv ' + lines[2])
        output.write('Tv ' + lines[3])
        output.write('Tv ' + lines[4])
        output.write('\n')

def main():
    args = parse_args()
    file_name = args.file_name

    if not file_name.endswith(".vasp"):
        raise ValueError('The input structure must be a .vasp file')

    vasp_2_com(file_name)

if __name__ == "__main__":
    main()

