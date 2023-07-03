#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter
from gaussianutility.utilities import readoutput

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Extract geometry from Gaussian output file (.out)\n"
                     "to Gaussian input file (.com)\n"
                     "This support ONIOM-type calculations\n"
                     "Connectivity information will be lost\n\n"
                     "Return file_name_geom.com: Gaussian input file"
                     , formatter_class=RawTextHelpFormatter)

    parser.add_argument('file_name', help='Gaussian output file (.out)')
    parser.add_argument('-i', '--index', nargs='?', const=1, \
    help='Optimization step index;\n0 for input, -1 for the last geometry, '+\
         'and L for the lowest energy geometry', default=-1)
    args = parser.parse_args()
    return args
    
def main():
    args = parse_args()
    file_name = args.file_name
    stepIdx = args.index

    # Read Gaussian output file
    routeStr, charge_mult, df_geom = readoutput(file_name, stepIdx)

    # Write .com file
    name = ".".join(file_name.rsplit(".",1)[0:-1])
    out_file = name + "_geom.com"
    output = open(out_file, 'w')
    output.write(routeStr + '\n\n')
    output.write(out_file + '\n\n')
    output.write(charge_mult + '\n')
    output.write(df_geom.to_string(index=False, header=False))
    output.write('\n')
    output.close()


