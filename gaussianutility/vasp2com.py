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

    parser.add_argument("file_name", help='VASP input file, e.g., POSCAR')
    return parser.parse_args()

def vasp_2_com(file_name):
    # Read geometry information and make Pandas dataframe
    lines = open(file_name).readlines()
    
    title = lines[0]
    
    lattice_scailing_factor = float(lines[1].split()[0])
    
    lattice_vector = [line.split() for line in lines[2:5]]
    lattice_vector = np.array(lattice_vector, dtype = float)
    
    elem_type = lines[5].split()
    
    elem_num = np.array(list(map(int, lines[6].split())))
    elect_num = np.array([md.element(e).electrons for e in elem_type])
    multiplicity = 1 if sum(elect_num * elem_num) % 2 == 0 else 2
    
    elem_list = np.array(list(itertools.chain.from_iterable([ [e] * n for e, n in zip(elem_type, elem_num)])))
    
    if lines[7].split()[0].lower() == 'selective':
        geom_lines = []
        for line in lines[9:]:
            if len(line.split()) <= 2:
                break
                
            geom_lines.append(line.split())
    
        geom_lines = np.array(geom_lines)
        geom = np.array(geom_lines[:,:3], dtype=float)
    
        if lines[8].split()[0].lower() == 'cartesian':
            geom *= lattice_scailing_factor
        elif lines[8].split()[0].lower() == 'direct':
            geom = np.dot(geom, lattice_vector) * lattice_scailing_factor
        else:
            raise TypeError('Cannot read geometry type parameter - It should be either cartesian or direct.')
            
        flags = geom_lines[:,3:]
        flags_for_gaussian = []
        for flag in flags:
            if np.unique(flag) == 'T':
                flags_for_gaussian.append('0')
            elif np.unique(flag) == 'F':
                flags_for_gaussian.append('-1')
            else:
                raise TypeError('Supported selective dynamics flags are only T,T,T and F,F,F.')
    
        df_geom = pd.DataFrame(geom, columns = ['x', 'y', 'z'])
        df_geom.insert(loc=0, column='flag', value=flags_for_gaussian)
        df_geom.insert(loc=0, column='element', value=elem_list)
        
    else:
        geom = np.array([line.split() for line in lines[8:] if len(line.split()) >= 2], dtype='float_')
        if lines[7].split()[0].lower() == 'cartesian':
            geom *= lattice_scailing_factor
        elif lines[7].split()[0].lower() == 'direct':
            geom = np.dot(geom, lattice_vector) * lattice_scailing_factor
        else:
            raise TypeError('Cannot read geometry type parameter - It should be either cartesian or direct.')
    
        df_geom = pd.DataFrame(geom, columns = ['x', 'y', 'z'])
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

    vasp_2_com(file_name)

if __name__ == "__main__":
    main()

