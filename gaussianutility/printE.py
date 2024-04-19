#!/usr/bin/env python3

import numpy as np
import argparse
from argparse import RawTextHelpFormatter

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Print electronic energy, electronic energy with zero-point correction,\n"
                     "enthalpy, and Gibbs free energy from Gaussian output file(s) after frequency calculation\n"
                     "with some other structural information\n\n"
                     "Print the following values on terminal:\n"
                     "  File_name, Stoichiometry(*charge *multiplicity), E, E+ZPE, H, G, *number of imaginary frequencies\n"
                     "    *charge: printed with a sign (+ or -) if not charge-neutral\n"
                     "    *multiplicity: printed if not singlet\n"
                     "    *number of imaginary frequencies: printied if exist"
                     , formatter_class=RawTextHelpFormatter)

    parser.add_argument('file_name', nargs='+', help='Gaussian output file(s) (.out)')
    args = parser.parse_args()
    return args
    
def print_E(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    stoich, E, EZPE, H, G, imagf = '', '', '', '', '', ''

    for line in lines:
        if "Stoichiometry" in line: 
            stoich = line.split()[1]
            break
    
    for idx, line in reversed(list(enumerate(lines))):
        if "SCF Done" in line: 
            E = line.split()[4]
            break
        
    for line in lines[len(lines)-idx:]:
        if "zero-point Energies" in line: EZPE = line.split()[6]
        if "thermal Enthalpies" in line: H = line.split()[6]
        if "thermal Free Energies" in line: G = line.split()[7]
        if "imaginary frequencies ignored" in line: imagf = line.split()[0]
        else: imagf = '0'
    
    return {
        'file_name': file_path.rsplit("/",1)[-1].rsplit(".",1)[0],
        'stoich': stoich,
        'E': E,
        'EZPE': EZPE,
        'H': H,
        'G': G,
        'imagf': imagf
    }

def main():
    args = parse_args()
    file_names = args.file_name

    for file_name in file_names:
        result = print_E(file_name)
        print(f"{result['file_name']}  {result['stoich']}  {result['E']}  {result['EZPE']}  "
              f"{result['H']}  {result['G']}  {result['imagf']}")

if __name__ == '__main__':
    main()
        