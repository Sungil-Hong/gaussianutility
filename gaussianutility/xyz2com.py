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
    with open(file_name, 'r') as f:
        lines = f.readlines()

    # Atoms start at line 2
    atoms_start = 2
    atoms_lines = [line for line in lines[atoms_start:] if len(line.split()) >= 4]

    # Check if this is an extended XYZ file by determining whether the first atom
    # line contains more than four tokens (atom, x, y, z)
    is_extended = len(atoms_lines[0].split()) > 4 

    # Decide multiplicity using only the element type from each atom line
    elemt_list = [line.split()[0] for line in atoms_lines] 
    elemt_type, elemt_count = np.unique(elemt_list, return_counts=True)
    elemt_electrons = [md.element(elemt).electrons for elemt in elemt_type]
    total_electrons = sum(np.array(elemt_electrons) * elemt_count)
    multiplicity = 1 if total_electrons % 2 == 0 else 2

    # Write Gaussian input file using only the first four tokens: element, x, y, z
    name = ".".join(file_name.rsplit(".",1)[0:-1])
    out_file = name + ".com"
    with open(out_file, 'w') as output:
        output.write('# hf/3-21g\n\n')
        output.write(f'{out_file}\n\n')
        output.write(f'0 {multiplicity}\n')
        for line in atoms_lines:  # MODIFIED: iterate over atoms_lines
            tokens = line.split()
            if len(tokens) >= 4:
                element = tokens[0]
                x, y, z = tokens[1:4]
                output.write(f' {element:2s}    {x:>15s} {y:>15s} {z:>15s}\n')
        output.write('\n')

def main():
    args = parse_args()
    file_name = args.file_name

    if not file_name.endswith(".xyz"):
        raise ValueError('The input structure must be a .xyz file')
    
    xyz_2_com(file_name)

if __name__ == "__main__":
    main()


