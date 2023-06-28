#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter
from gaussianutility.utilities import readinput

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Freeze a given ONIOM layer(s) in a Gaussian input structure during\n"
                     "optimization calculation by adding indices (-1) next to atom symbols\n\n"
                     "Return file_name.com: Modified input structure"
                     , formatter_class=RawTextHelpFormatter)

    parser.add_argument("file_name", help='Gaussian input file (.com or .gjf)')
    parser.add_argument('-i', '--index', nargs='?', const=1, \
    help="'H', 'M', or 'L' for freezing a high, medium, or low ONIOM layer, or their combination\n"+\
         "'L' is default", default='L')
    args = parser.parse_args()
    return args
    
def main():
    args = parse_args()
    file_name = args.file_name
    freezeIdx = args.index

    # Read a Gaussian input file
    route, title, charge_mult, df_geom, connectivity = readinput(file_name)

    if ("oniom" or "ONIOM" or "Oniom") not in route:
        raise TypeError("Input structure must be in ONIOM scheme")
        
    if "Index" in list(df_geom.columns):
        df_geom["Index"] = np.zeros(len(df_geom), dtype = int)
        
    if "Index" not in list(df_geom.columns):
        initIdx = np.zeros(len(df_geom), dtype = int)
        df_geom.insert(1, "Index", initIdx)
        
    # Read index for freezing
    freezeIdx = list(freezeIdx)    
    for idx in freezeIdx:
        if idx not in ['H', 'M', 'L']:
            raise ValueError("Index for layer to freeze must be H, M, L or their combination")
            
    # Add indicies to freeze atoms
    for layer in freezeIdx:
        df_geom.loc[df_geom["ONIOM_layer"] == layer, "Index"] = -1
        
    # Write Gaussian input file
    outfile = open(file_name, 'w')
    outfile.write(route + '\n')
    outfile.write(title + '\n')
    outfile.write(charge_mult)
    outfile.write(df_geom.to_string(index=False, header=False))
    outfile.write('\n\n')
    outfile.writelines(connectivity)
    outfile.write('\n')
    outfile.close()

