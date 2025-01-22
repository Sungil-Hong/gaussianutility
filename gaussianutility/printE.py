#!/usr/bin/env python3

import numpy as np
import argparse
from argparse import RawTextHelpFormatter

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Print electronic energy and thermochemistry results (if applicable)from Gaussian output file(s)\n"
                     "Print the following values on terminal:\n"
                     "  File_name, Stoichiometry(*charge *multiplicity)\n"
                     "     Calculation type, E, E+ZPE, H, G, number of imaginary frequencies\n\n"
                     "*charge: print with a sign (+ or -) if not neutral\n"
                     "*multiplicity: print if not in singlet spin state\n"
                     , formatter_class=RawTextHelpFormatter)

    parser.add_argument('file_name', nargs='+', help='Gaussian output file(s) (.out or .log)')
    args = parser.parse_args()
    return args
    
def read_E(lines):
    stoich, E, EZPE, H, G, imagf = '', '', '', '', '', ''

    for line in lines:
        if "Stoichiometry" in line: 
            stoich = line.split()[1]
            break
    
    for idx, line in reversed(list(enumerate(lines))):
        if "SCF Done" in line: 
            E = line.split()[4]
            break
        
    job_type = None
    for line in lines[idx:]:
        if "Thermochemistry" in line:
            job_type = 'frequency'
            break
            
        if "Stationary point found" in line:
            job_type = 'optimization'
            break
        
    if job_type == None:
        job_type = 'single-point'

    if job_type == 'frequency':
        for line in lines[idx:]:
            if "zero-point Energies" in line: EZPE = line.split()[6]
            if "thermal Enthalpies" in line: H = line.split()[6]
            if "thermal Free Energies" in line: G = line.split()[7]
            if "imaginary frequencies ignored" in line: imagf = line.split()[0]
            else: imagf = '0'
        
        return {
            'stoich': stoich,
            'job type': job_type,
            'E': E,
            'EZPE': EZPE,
            'H': H,
            'G': G,
            'imagf': imagf
        }

    else:
        return {
            'stoich': stoich,
            'job type': job_type,
            'E': E
        }

def main():
    args = parse_args()
    file_names = args.file_name

    for file_name in file_names:
        with open(file_name, 'r') as file:
            lines = file.readlines()

        try:
            lines[-1]
        except IndexError:
            print("!!!Caution: " + file_name + " is an empty file.")
            continue

        if "Normal" not in lines[-1].split():
            calc_E = [read_E(lines)]
            print("Results of " + file_name + ", " + calc_E[0]['stoich'])
            print("!!!Caution: The calculation does not seem to be normally ternimated!!!")
            print(f"   1 calc: {calc_E[0]['job type']}  {calc_E[0]['E']}")

        else:        
        ### Divide multiple calculation steps ###
            linked_jobs_idx = [0]
            for idx, line in enumerate(lines):
                if line.startswith(" Normal"):
                    linked_jobs_idx.append(idx)
        
            ### Read data of each job ###
            calc_lines = [lines[linked_jobs_idx[i]:linked_jobs_idx[i+1]+2] for i in range(len(linked_jobs_idx)-1)]
            calc_E = [read_E(calc) for calc in calc_lines]
            print("Results of " + file_name + ", " + calc_E[0]['stoich'])
        
            for i in range(len(calc_E)):
                result = calc_E[i]
                if result['job type'] == 'frequency':
                    print(f"   {i+1} calc: {result['job type']}  {result['E']}  {result['EZPE']}  {result['H']}  {result['G']}  {result['imagf']}")
                else:
                    print(f"   {i+1} calc: {result['job type']}  {result['E']}")
    
                
    if "Normal" not in lines[-1].split():
        print("!!!Caution: The calculation does not seem to be normally ternimated!!!")


if __name__ == '__main__':
    main()
        
