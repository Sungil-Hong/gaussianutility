#!/usr/bin/env python3

import numpy as np
import pandas as pd
from periodictable import elements
import argparse
from argparse import RawTextHelpFormatter
from gaussianutility.utilities import readoutput

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Extract geometry from Gaussian output file (.out)"+\
                     "to xyz structure file\n"
                     "The ONIOM partition and connectivity information will be lost\n\n"
                     "Return file_name.xyz: xyz structure file"
                     , formatter_class=RawTextHelpFormatter)

    parser.add_argument('file_name', help='Gaussian output file (.out)')
    parser.add_argument('-i', '--index', nargs='?', const=1, \
    help='Optimization step index;\n0 for input, -1 for the last geometry,\n'+\
          'and L for the lowest energy geometry', default=-1)
    args = parser.parse_args()
    return args
    
def main():
    args = parse_args()
    file_name = args.file_name
    stepIdx = args.index

    # Read Gaussian output file
    _, _, df_geom = readoutput(file_name, stepIdx)
    df_geom = df_geom[['Atom', 'x', 'y', 'z']]
    df_geom['Atom'] = np.array([i.rsplit("-",1)[0] for i in list(df_geom['Atom'])])
    elem_len = str(len(df_geom))

    # Write xyz file
    name = ".".join(file_name.rsplit(".",1)[0:-1])
    out_file = name + ".xyz"
    output = open(out_file, 'w')
    output.write(elem_len + '\n')
    output.write(out_file + '\n')
    output.write(df_geom.to_string(index=False, header=False))
    output.write('\n')
    output.close()

