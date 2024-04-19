#!/usr/bin/env python3

import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import mendeleev as md

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Convert xyz file to Gaussian input file (.com)\n\n"
                     "Return file_name.com: Gaussian input file",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument("file_name", help='XYZ structure file (.xyz)')
    return parser.parse_args()

def xyz_2_com(file_name):
    # Read geometry from XYZ file
    lines = open(file_name).readlines()

    # Find the line where the XYZ coordinates start
    for idx, line in enumerate(lines[:]):
        if len(line.split()) == 4:
            try:
                _ = [float(val.replace('n','')) for val in line.split()[1:]]
                break
            except ValueError:
                continue

    # Decide multiplicity assuming overall neutral charge
    elemt_list = [line.split()[0] for line in lines[idx:]]
    elemt_type, elemt_count = np.unique(elemt_list, return_counts=True)
    elemt_electrons = [md.element(elemt).electrons for elemt in elemt_type]
    total_electrons = sum(np.array(elemt_electrons) * elemt_count)
    multiplicity = 1 if total_electrons % 2 == 0 else 2

    # Write Gaussian input file
    name = ".".join(file_name.rsplit(".",1)[0:-1])
    out_file = name + ".com"
    with open(out_file, 'w') as output:
        output.write('# hf/3-21g\n\n')
        output.write(f'{out_file}\n\n')
        output.write(f'0 {multiplicity}\n')
        for line in lines[idx:]:
            output.write(line)
        output.write('\n')

def main():
    args = parse_args()
    file_name = args.file_name

    if not file_name.endswith(".xyz"):
        raise ValueError('The input structure must be a .xyz file')
    
    xyz_2_com(file_name)

if __name__ == "__main__":
    main()


