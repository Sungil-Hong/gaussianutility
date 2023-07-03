#!/usr/bin/env python3

import numpy as np
import math as m
import sys

# Thermodynamic constants
kB = 1.380649e-23 # J/K
h = 6.62607015e-34 # Js
Rgas = 8.3144626 # J/mol/K
Na = 6.02214076e23

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
    
    return Stot, Etot # return total entropy correction and energy correction
    
def main():
    # Print error/help messages that looks like argparse functionality
    if len(sys.argv) == 1:
        print("usage: gibbsTemp [-h] file_name temp1 [temp2 [step_number]]")
        print("gibbsTemp: error: the following arguments are required: file_name temp1")
        sys.exit()

    inFile = sys.argv[1]
    if inFile in ['-h', '--help']:
        print("usage: gibbsTemp [-h] file_name temp1 [temp2 [step_number]]")
        print("\n")
        print("Calculated Gibbs free energy at different tempeature(s) from Gaussian output file (.out)")
        print("The Gaussian job must be normally terminated with frequency calculation\n")
        print("Print the following values on terminal:")
        print("     Single G value, if one temp is provided")
        print("     List of 20 G values between temp1 and temp2 if both are provided")
        print("     List of G values as many as step_number if temp1, temp2, and step_number are provided\n")
        print("positional arguments:")
        print("  file_name              Gaussian output file (.out)")
        print("  temp1                  (starting) temperature for G calculation")
        print("  [temp2]                End temperature for G calculation")
        print("  [step_number]          Number of steps for G calculation bewteen temp1 and temp2")
        print("                         Default is 20\n")
        print("options:")
        print("  -h, --help             show this help message and exit")
        sys.exit()

    #Set temperature based on the provided arguments
    l = len(sys.argv)
    if l == 2: print("Provide temperature argument(s)"); sys.exit(1)
    elif l == 3: temperature = np.array([float(sys.argv[2])])
    elif l == 4: temperature = np.linspace(float(sys.argv[2]), float(sys.argv[3]), 20)
    elif l == 5: temperature = np.linspace(float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]))
    else: print("Too many arguments"); sys.exit(1)

    # Exctract thermochemistry result from the output file
    inFileRead = open(inFile,'r')
    for idx, line in enumerate(inFileRead):
        if "Thermochemistry" in line: beginIdx = idx
        if "Total Bot" in line: endIdx = idx; break
    inFileRead.close()

    thermochem = []   
    inFileRead = open(inFile,'r')
    for idx, line in enumerate(inFileRead):
        if "Multiplicity" in line: multiplicity = int(line.split()[5]) 
        if idx >= beginIdx-1 and idx < endIdx-3:
            thermochem.append(line)

    # Extract required informations from the thermochemistry output
    for idx, line in enumerate(thermochem):
        if "Pressure" in line:
            press = float(line.split()[4])*101325 # Pa
        if "Molecular mass" in line:
            mass = float(line.split()[2])*1.6605391e-27 # kg
        if "Rotational symmetry number" in line:\
            rho_r = float(line.split()[3]) 
        if "Rotational temperatures" in line: 
            if len(line.split()) < 5: 
                theta_r = np.array([float(line.split()[3])])
            else:
                theta_r = np.float64(np.array(line.split()[3:6]))
        if "Vibrational temperature" in line:
            vibTempBeginIdx = idx
        if "Zero-point correction" in line: 
            zpCorr = float(line.split()[2])*2625.4996394799e3 # J/mol
            vibTempEndIdx = idx
        if "Sum of electronic and zero-point Energies" in line:
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


