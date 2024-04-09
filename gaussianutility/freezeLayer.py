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
                     "Return file_name.com: Modified input structure",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument("file_name", help='Gaussian input file (.com or .gjf)')
    parser.add_argument('-i', '--index', nargs='?', const=1, \
    help="'H', 'M', or 'L' for freezing a high, medium, or low ONIOM layer, or their combination\n"+\
         "'U' unfreezes all the atoms. 'L' is default", default='L')
    args = parser.parse_args()
    return args
    
def freeze_layer(file_name, index):
    # Read a Gaussian input file
    route, title, charge_mult, df_geom, connectivity = readinput(file_name)

    if "oniom" not in route.lower():
        raise TypeError("Input structure must be in ONIOM scheme")
        
    if "Index" in df_geom.columns:
        df_geom["Index"] = 0
    else:
        df_geom.insert(1, "Index", 0)
        
    # Read index for freezing
    freezeIdx = list(freezeIdx)    
    if 'U' in freezeIdx:
        if len(freezeIdx) > 1:
            raise ValueError("U cannot be combined with other layer indices.")

    else:
        for layer in freezeIdx:
            if layer not in ['H', 'M', 'L']:
                raise ValueError("Layer index must be H, M, L, their combination, or U")
                        
            df_geom.loc[df_geom["ONIOM_layer"] == layer, "Index"] = -1

    # Write Gaussian input file
    with open(file_name, 'w') as output:
        output.write(f"{route}\n")
        output.write(f"{title}\n")
        output.write(charge_mult)
        output.write(df_geom.to_csv(index=False, header=False, sep='\t'))
        output.write('\n\n')
        output.writelines(connectivity)
        output.write('\n')

def main():
    args = parse_args()
    freeze_layer(args.file_name, args.index)

if __name__ == "__main__":
    main()
    

