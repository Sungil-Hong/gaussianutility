#!/usr/bin/env python3

import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import mendeleev as md

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Convert xyz file to Gaussian input file (.com)\n\n"
                     "Return file_name.com: Gaussian input file"
                     , formatter_class=RawTextHelpFormatter)

    parser.add_argument("file_name", help='XYZ structure file (.xyz)')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    file_name = args.file_name

    # Read geometry from a Gaussian input file
    lines = open(file_name).readlines()
    
    for idx, line in enumerate(lines[:]):
        if len(line.split()) == 4:
            isxyz = True
            for val in line.split()[1:]:
                 isxyz *= val.replace('.', '', 1).isdigit()

            if isxyz: break

    # Decide multiplicity assuming overall neutral charge
    elemts = []
    for line in lines[idx:]:
        elemts.append(line.split()[0])
        
    elemts, counts = np.unique(elemts, return_counts=True)
        
    elem_elect = []

    for elem in elemts:
        elem_elect.append(md.element(elem).electrons)

    total_elect = sum(np.array(elem_elect) * counts)

    if total_elect%2 == 0:
        multiplicity = 1
    else: multiplicity = 2

    # Write Gaussian input file
    name = ".".join(file_name.rsplit(".",1)[0:-1])
    out_file = name + ".com"
    output = open(out_file, 'w')
    output.write('# hf/3-21g\n\n')
    output.write(out_file + '\n\n')
    output.write('0 ' + str(multiplicity) + '\n')
    for line in lines[idx:]:
        output.write(line)
    output.write('\n')
    output.close()


