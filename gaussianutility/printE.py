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
    
def main():
    args = parse_args()
    file_name = args.file_name

    #file_name = sys.argv[1:]
    for file in file_name:
        inputfile = open(file)
        lines = inputfile.readlines()
        inputfile.close()

        E, EZPE, H, G, imagf = '', '', '', '', ''

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
        
        print(file.rsplit("/",1)[-1].rsplit(".",1)[0] + "  " + str(stoich) + "  " + str(E) + "  " + str(EZPE)\
        + "  " + str(H) + "  " + str(G)+ "  "  + str(imagf))

