#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter
from gaussianutility.utilities import readinput

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Convert Gaussian input file (.com or .gjf) to xyz file\n\n"
                     "Return file_name.xyz: xyz structure file",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument("file_name", help='Gaussian input file (.com or .gjf)')
    return parser.parse_args()
    
def com_2_xyz(file_name):
    # Read geometry from a Gaussian input file
    _, _, _, df_geom, _ = readinput(file_name)
    df_geom = df_geom[['Atom', 'x', 'y', 'z']]
    df_geom['Atom'] = np.array([i.rsplit("-",1)[0] for i in list(df_geom['Atom'])])
    elem_len = str(len(df_geom))

    # Write xyz file
    name = ".".join(file_name.rsplit(".",1)[0:-1])
    out_file = name + ".xyz"
    with open(out_file, 'w') as output:
        output.write(elem_len + '\n')
        output.write(out_file + '\n')
        output.write(df_geom.to_string(index=False, header=False) + '\n')

def main():
    args = parse_args() 
    com_2_xyz(args.file_name)

if __name__ == "__main__":
    main()

