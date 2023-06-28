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
                     "Return file_name.com: Gaussian input file"
                     , formatter_class=RawTextHelpFormatter)

    parser.add_argument("file_name", help='VASP input file (.vasp)')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    file_name = args.file_name
        
    # Check input file
    name, input_format = file_name.rsplit(".",1)[0], file_name.rsplit(".",1)[-1]
    if input_format != "vasp":
        raise ValueError('The input structure must be .vasp file')

    # Read geometry information and make Pandas dataframe
    lines = open(file_name).readlines()

    elem_type = lines[5].split()
    elem_num = np.array(list(map(int, lines[6].split())))
    elect_num = []

    for e in elem_type:
        elect_num.append(md.element(e).electrons)
        
    elect_num = np.array(elect_num)

    if sum(elect_num * elem_num)%2 == 0:
        multiplicity = 1
    else: multiplicity = 2

    elem_list = []
    for i in range(len(elem_type)):
        elem_list.append((elem_type[i],)*elem_num[i])
        
    elem_list = np.array(list(itertools.chain(*elem_list)))

    geom = []
    for line in lines[8:]:
        if len(line) < 2: break
        geom.append(line.split())

    geom = np.array(geom, dtype='float_')

    lattice_scailing_factor = float(lines[1].split()[0])
    lattice_vector = [line.split() for line in lines[2:5]]
    lattice_vector = np.array(lattice_vector, dtype = float)

    geomMod = lines[7].split()[0]
    if geomMod == 'Cartesian' or geomMod == 'cartesian':
        geom = geom*lattice_scailing_factor
    elif geomMod == 'Direct' or geomMod == 'direct':
        geom = np.dot(lattice_vector, geom.T).T
        geom = geom*lattice_scailing_factor
        
    df_geom = pd.DataFrame(geom, columns = ['x','y','z'])
    df_geom.insert(loc=0, column='element', value=elem_list)

    # Write .com file
    out_file = file_name.rsplit(".",1)[0] + ".com"
    output = open(out_file, 'w')
    output.write('# pbepbe/3-21g/auto\n')
    output.write('\n')
    output.write(name)
    output.write('\n\n')
    output.write('0 '+str(multiplicity))
    output.write('\n')
    output.write(df_geom.to_string(index=False, header=False))
    output.write('\n')
    output.write('Tv ' + lines[2])
    output.write('Tv ' + lines[3])
    output.write('Tv ' + lines[4])
    output.write('\n')
    output.close()


