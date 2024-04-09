#!/usr/bin/env python3

import math as m
import numpy as np
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter
import matplotlib.pyplot as plt

##### Change these lines to modify X range or make corrections #####
##### For UV-Vis #####
UV_wl_range = np.linspace(100,500,801) # nm
UV_alpha = 1 #UV-Vis wavelength scaling factor; default
#UV_alpha = 1.12 #for M062X; ref: https://doi.org/10.26434/chemrxiv-2022-pfvhc-v3
UV_delta = 0.4/m.sqrt(2) #UV-Vis bandwith parameter, eV; default
#UV_delta = 0.21  #for M062X; ref: https://doi.org/10.26434/chemrxiv-2022-pfvhc-v3

##### For IR and Raman #####
IR_wv_range = np.linspace(0,4000,8001) # cm-1
IR_wn_factor = 1 # default
#IR_wn_factor = 0.947 #For M062X; ref: J. Phys. Chem. A, 121(11), 2265 (2017)#
IR_HWHM = 4 # IR peak half-width at half-max, cm-1 # default

# Define functions to extract spectrum data from output file(s)
def uv_vis(file_name):
    readfile = open(file_name).readlines()
    wavelengths = []
    strengths = []
    
    for idx, line in enumerate(reversed(readfile)):
        if "Excitation energies and oscillator strengths" in line:
            break
            
    ex_idx = len(readfile)-idx-1
    
    for line in readfile[ex_idx:]:
        if 'Excited State' in line:
            wavelengths.append(float(line.split()[6]))
            strengths.append(float(line.split()[8].split('=')[1]))
            
    if len(wavelengths) == 0:
        raise TypeError('The Gaussian job does not look like excited state calculations')

    wavelengths = np.array(wavelengths)
    strengths = np.array(strengths)
    
    return wavelengths, strengths


def ir(file_name):
    readfile = open(file_name, 'r')
    frequencies = []
    intensities = []
    
    for line in readfile:
        if 'Frequencies --' in line:
            frequencies.append(line.split()[2:])
        if 'IR Inten' in line:
            intensities.append(line.split()[3:])
            
    if len(intensities) == 0:
        raise TypeError('The Gaussian job does not look like conntaining normal IR information')
            
    frequencies = np.array(frequencies, dtype='float').flatten() #cm-1
    intensities = np.array(intensities, dtype='float').flatten() #km/mole

    return frequencies, intensities

def raman(file_name):
    readfile = open(file_name, 'r')
    frequencies = []
    intensities = []
    
    for line in readfile:
        if 'Frequencies --' in line:
            frequencies.append(line.split()[2:])
        if 'Raman Activ' in line:
            intensities.append(line.split()[3:])
            
    if len(intensities) == 0:
        raise TypeError('The Gaussian job does not look like conntaining normal IR information')
            
    frequencies = np.array(frequencies, dtype='float').flatten() #cm-1
    intensities = np.array(intensities, dtype='float').flatten() 

    return frequencies, intensities

# Constants
h = 6.6261e-34 # J/s
c = 299792458  # m/s
J2eV = 6.02214076e23/96.48530749925793/1000

# Converstion between energy and wavelength
def wvl2E(wvl): # wvl in nm
    E = h*c/(wvl/1e9)*J2eV
    return E # E in eV

def E2wvl(E): #E in eV
    wvl = h*c*1e9/(E/J2eV)
    return wvl # wvl in nm

# Gaussian distribution for UV-Vis
def uvGauss(x, f, wv):
    #Reference: https://gaussian.com/uvvisplot/
#    return 13.062973*f*3099.6*m.exp(-((1/x-1/wv)/0.00032262)**2)
    return 13.062973*f*E2wvl(m.sqrt(2)*UV_delta)*m.exp(-((wvl2E(x)-wvl2E(wv)/UV_alpha)/UV_delta/m.sqrt(2))**2)

# Cauchy distribution for normal and raman IR
def cauchy(x, mu):
    gamma = IR_HWHM
    return 1/m.pi/gamma/(1+((x-mu)/gamma)**2)

def parse_args():
    parser = argparse.ArgumentParser(
        description = """
        Produce Normal-IR, Raman, or UV-Vis spectrum from Gaussian output file (.out)
        For Raman spectrum, the Gaussian calculation should be done with freq=raman keyword
        For UV-Vis spectrum, excited state calculations (like TD-DFT and EOM-CCSD) should be done
        Linear combination of multiple spectra for mixtures is available by adding ratio arguments
        This case, the sum of the ratios must be 1

        Examples of command line usage are:
            spectrum ir file1.out
            spectrum uv file1.out file2.out -r 0.8 0.2

        Return file_name_type.png: Spectrum image
               file_name_type.csv: Sepctrum data
        """,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument('type', 
                        choices=["uv", "ir", "raman"],
                        help="Type of spectrum; 'uv' for UV-Vis, 'ir' for normal IR, and 'raman' for Raman")
    parser.add_argument('file_name', nargs='+', 
                        help="Gaussian output files (.out)")
    parser.add_argument('-r', '--ratio', nargs='?', 
                        help="Ratios of input structures (required when multiple structures are provided)")
    args = parser.parse_args()

    return args

def main():
    # Read provided arguments
    args = parse_args()
    type_name = args.type
    file_names = args.file_name
    ratios = [float(val) for val in args.ratio] if args.ratio else [1]

    if len(file_names) != len(ratios):
        raise ValueError("Ratio(s) of one or more species is not provided")
    
    if sum(ratios) > 1.0001 or sum(ratios) < 0.9999:
        raise ValueError('Sum of the ratios of the species should be 1')

    names = []
    for file_name in file_names:
        name, informat = file_name.rsplit(".", 1)
        names.append(name)
        if informat not in ("out"):
            raise TypeError("The input file format must be .out")

    # Plot spectra
    # UV-Vis
    if type_name == 'uv':
        spectrum_save = '_'.join(names) + '_uv_vis.png'
        data_save = '_'.join(names) + '_uv_vis.csv'
        
        f = plt.figure()
        f.set_figwidth(12)
        f.set_figheight(8)
        
        xrange = UV_wl_range
        curve_per_file = []
        
        for file_name in file_names:
            wavelengths, strengths = uv_vis(file_name)
            
            curve_per_x = []
            for idx, wv in enumerate(wavelengths):
                wv_array = np.ones(len(xrange)) * wv
                st_array = np.ones(len(xrange)) * strengths[idx]
                
                curve = np.array(list(map(uvGauss, xrange, st_array, wv_array)))
                curve_per_x.append(curve)
            
            final_curve = np.sum(curve_per_x, axis=0)
            curve_per_file.append(final_curve)
            
            if len(file_names) == 1:
                plt.plot(xrange, final_curve, 'b')
            else:
                blue = 0.8
                alpha = 0.6
                red = np.random.uniform(low=0, high=1)
                green = np.random.uniform(low=0, high=1)
                plt.plot(xrange, final_curve, color = (red, green, blue), alpha=alpha)
        
        if len(file_names) == 1:
            f.legend(names, loc='upper right', bbox_to_anchor=(0.9, 0.88))
            
        else:
            summed_curve = np.dot(np.array(curve_per_file).T, np.array([ratios]).T)
            plt.plot(xrange, summed_curve, 'b')
            f.legend(names + ['Mixture'], loc='upper right', bbox_to_anchor=(0.9, 0.88))
            
        plt.xlabel('Wavelength, $\lambda$  $(nm)$', fontsize=20)
        plt.ylabel('Absorbance $(L/mol/cm)$', fontsize=20)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.savefig(spectrum_save)
        
        if len(file_names) == 1:
            uv_vis_df = pd.DataFrame(
                np.hstack((np.array([xrange]).T,np.array(curve_per_file).T)),
                columns = ['wavelength'] + names)
        else:
            uv_vis_df = pd.DataFrame(
                np.hstack((np.array([xrange]).T,np.array(curve_per_file).T, summed_curve)),
                columns = ['wavelength']+ names + ['Mixture'])
        uv_vis_df.to_csv(data_save, index=False)

    # IR: both normal IR and raman
    if type_name in ['ir', 'raman']:
        if type_name == 'ir':
            spectrum_save = '_'.join(names) + '_ir.png'
            data_save = '_'.join(names) + '_ir.csv'
        elif type_name == 'raman':
            spectrum_save = '_'.join(names) + '_raman.png'
            data_save = '_'.join(names) + '_raman.csv'
        
        f = plt.figure()
        f.set_figwidth(12)
        f.set_figheight(8)
        
        xrange = IR_wv_range
        curve_per_file = []
        
        for file_name in file_names:
            if type_name == 'ir':
                frequencies, intensities = ir(file_name)
                frequencies *= IR_wn_factor # Correction factor
                intensities *= 1e5/2.24e3
            elif type_name == 'raman':
                frequencies, intensities = raman(file_name)
                
            curve_per_x = []
            for idx, mu in enumerate(frequencies):
                mu_array = np.ones(len(xrange)) * mu
                curve = np.array(list(map(cauchy, xrange, mu_array))) * intensities[idx]
                curve_per_x.append(curve)
            
            final_curve = np.sum(curve_per_x, axis=0)
            curve_per_file.append(final_curve)
            
            if len(file_names) == 1:
                plt.plot(xrange, final_curve, 'b')
            else:
                blue = 0.8
                alpha = 0.6
                red = np.random.uniform(low=0, high=1)
                green = np.random.uniform(low=0, high=1)
                plt.plot(xrange, final_curve, color = (red, green, blue), alpha=alpha)
        
        if len(file_names) == 1:
            f.legend(name, loc='upper right', bbox_to_anchor=(0.9, 0.88))
            
        else:
            summed_curve = np.dot(np.array(curve_per_file).T, np.array([ratios]).T)
            plt.plot(xrange, summed_curve, 'b')
            f.legend(names + ['Mixture'], loc='upper right', bbox_to_anchor=(0.9, 0.88))
            
        plt.xlabel('Wavenumber, $\\nu$  $(cm^{-1})$', fontsize=20)
        plt.ylabel('Absorbance (a.u.)', fontsize=20)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.savefig(spectrum_save)
        
        if len(file_names) == 1:
            ir_df = pd.DataFrame(
                np.hstack((np.array([xrange]).T,np.array(curve_per_file).T)),
                columns = ['wavenumber'] + names)
        else:
            ir_df = pd.DataFrame(
                np.hstack((np.array([xrange]).T,np.array(curve_per_file).T, summed_curve)),
                columns = ['wavenumber'] + names + ['Mixture'])
        ir_df.to_csv(data_save, index=False)

if __name__ == "__main__":
    main()
    
