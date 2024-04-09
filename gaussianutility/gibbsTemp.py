#!/usr/bin/env python3

import numpy as np
import math as m
import sys
import argparse
from argparse import RawTextHelpFormatter

# Thermodynamic constants
kB = 1.380649e-23 # J/K
h = 6.62607015e-34 # Js
Rgas = 8.3144626 # J/mol/K
Na = 6.02214076e23

def parse_args():
    parser = argparse.ArgumentParser(
        description= """
        Calculate Gibbs free energy (in Hartrees) at different tempeature(s) from Gaussian output file (.out)
        The Gaussian job must be normally terminated with frequency calculation
        Print temperatures and Gibbs free energies on terminal.
        """,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument('file_name', help='Gaussian output file (.out)')
    parser.add_argument('T1', type=float, help='Temperature to calculate G (in K)')
    parser.add_argument('T2', type=float, nargs='?', help='Upper bound of temperature range (in K)')
    parser.add_argument('step_number', type=int, nargs='?', help="Number of steps of temperature between T1 and T2 to calculate G", default=20)
    args = parser.parse_args()
    return args
    
# Calculate contributions from translational, rotational, vibrational, and eletronic motions
# The equations are sourced from "Thermochemistry in Gaussian" written by Joseph W. Ochterski (2000)
# https://gaussian.com/wp-content/uploads/dl/thermo.pdf
def Contributions(temp, mass, press, vibTemp, multiplicity, rho_r, theta_r):
    ## Translational part
    qt = m.pow(2*m.pi*mass*kB*temp/h**2, 3/2) * kB*temp/press
    St = Rgas*(m.log(qt)+1+3/2)
    Et = 3/2*Rgas*temp

    ## Rotational part
    if len(theta_r) > 1:    
        qr = m.pi**(1/2)/rho_r*m.sqrt((temp**3/np.prod(theta_r)))
    elif len(theta_r) == 1:
        qr = m.pi**(1/2)/rho_r*m.sqrt((temp/(theta_r[0])))
        
    if 'qr' in locals():
        Sr = Rgas*(m.log(qr)+3/2)
        Er = 3/2*Rgas*temp
    else:
        Sr = 0
        Er = 0

    ## Vibrational part
    summ1 = 0 # for entropy contribution
    summ2 = 0 # for energy contribution
    for val in vibTemp:
        summ1 += val/temp/(m.exp(val/temp) - 1) - m.log(1-m.exp(-val/temp))
        summ2 += val*(1/2 + 1/(m.exp(val/temp)-1))

    Sv, Ev = Rgas*np.array([summ1, summ2])

    ## Electronic part
    Se = Rgas * m.log(multiplicity)
    Ee = 0
    Stot = St + Sr + Sv + Se
    Etot = Et + Er + Ev + Ee
    
    return Stot, Etot # return total entropy and energy corrections
    
def gibbs_temp(file_name, T1, T2, step_number):
    #Set temperature based on the provided arguments
    if T2:
        temperature = np.linspace(T1, T2, step_number)
    else:
        temperature = np.array([T1])

    # Exctract thermochemistry result from the output file
    with open(file_name, 'r') as inFile:
        lines = inFile.readlines()

    multiplicity = None
    thermochem = []
    collecting_thermochem = False

    for line in lines:
        if "Multiplicity" in line:
            multiplicity = int(line.split()[5]) 
        elif "Thermochemistry" in line:
            collecting_thermochem = True
        elif "Total Bot" in line:
            break
        elif collecting_thermochem:
            thermochem.append(line.strip())
    
    if not thermochem:
        raise ValueError("Required information not found in the output file.")

    # Extract required informations from the thermochemistry results
    for idx, line in enumerate(thermochem):
        if "Pressure" in line:
            press = float(line.split()[4])*101325 # Pa
        elif "Molecular mass" in line:
            mass = float(line.split()[2])*1.6605391e-27 # kg
        elif "Rotational symmetry number" in line:\
            rho_r = float(line.split()[3]) 
        elif "Rotational temperatures" in line: 
            if len(line.split()) < 5: 
                theta_r = np.array([float(line.split()[3])])
            else:
                theta_r = np.float64(np.array(line.split()[3:6]))
        elif "Vibrational temperature" in line:
            vibTempBeginIdx = idx
        elif "Zero-point correction" in line: 
            zpCorr = float(line.split()[2])*2625.4996394799e3 # J/mol
            vibTempEndIdx = idx
        elif "Sum of electronic and zero-point Energies" in line:
            ElectE = float(line.split()[6])*2625.4996394799e3 - zpCorr # J/mol

    vibTemp = []
    for line in thermochem[vibTempBeginIdx:vibTempEndIdx]:
        for elem in line.split():
            try: vibTemp.append(float(elem))
            except: ValueError
    
    vibTemp = np.array(vibTemp)

    # Arrays of entropy and energy corrections in the given temperature (range)
    StotArr=[]
    EtotArr=[]
    if 'rho_r' not in locals():
        rho_r = []
        theta_r = []

    for temp in temperature:
        Stot, Etot = Contributions(temp, mass, press, vibTemp, multiplicity, rho_r, theta_r)
        StotArr.append(Stot)
        EtotArr.append(Etot)

    StotArr = np.array(StotArr)
    EtotArr = np.array(EtotArr)
    
    Hcorr = EtotArr + kB*temperature*Na # Thermal corrections to Enthalpy
    Gcorr = Hcorr - temperature * StotArr # Thermal corrections to Free energy
    Gibbs = (ElectE + Gcorr) / 2625.4996394799e3 # Final Gibbs free energy in Hartree

    temperature = np.ndarray.tolist(np.round(temperature,decimals=6))
    Gibbs = np.ndarray.tolist(np.round(Gibbs,decimals=6))

    print("Temperature [K]: " + str(temperature).strip('[]'))
    print("Gibbs free energy [Hartree]: " + str(Gibbs).strip('[]'))

def main():
    args = parse_args()
    gibbs_temp(args.file_name, args.T1, args.T2, args.step_number)

if __name__ == "__main__":
    main()
    
