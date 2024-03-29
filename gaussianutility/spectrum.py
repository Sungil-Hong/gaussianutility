#!/usr/bin/env python3

import sys
import math as m
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

##### Change these lines to modify X range or make corrections #####

UV_wl_range = np.linspace(100,500,801) # nm
UV_alpha = 1 #UV-Vis wavelength scaling factor; default
#UV_alpha = 1.12 #for M062X; ref: https://doi.org/10.26434/chemrxiv-2022-pfvhc-v3
UV_delta = 0.4/m.sqrt(2) #UV-Vis bandwith parameter, eV; default
#UV_delta = 0.21  #for M062X; ref: https://doi.org/10.26434/chemrxiv-2022-pfvhc-v3

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

h = 6.6261e-34 # J/s
c = 299792458  # m/s
J2eV = 6.02214076e23/96.48530749925793/1000

def wvl2E(wvl): # wvl in nm
    E = h*c/(wvl/1e9)*J2eV
    return E # E in eV

def E2wvl(E): #E in eV
    wvl = h*c*1e9/(E/J2eV)
    return wvl # wvl in nm

def uvGauss(x, f, wv):
    #Gaussian distribution for UV-Vis
    #Reference: https://gaussian.com/uvvisplot/
#    return 13.062973*f*3099.6*m.exp(-((1/x-1/wv)/0.00032262)**2)
    return 13.062973*f*E2wvl(m.sqrt(2)*UV_delta)*m.exp(-((wvl2E(x)-wvl2E(wv)/UV_alpha)/UV_delta/m.sqrt(2))**2)

def cauchy(x, mu):
    # Cauchy distribution for normal and raman IR
    gamma = IR_HWHM
    return 1/m.pi/gamma/(1+((x-mu)/gamma)**2)

def main():
    # Print error/help messages that looks like argparse functionality
    argvLen = len(sys.argv)
    if argvLen == 1:
        print("usage: spectrum [-h] type file_name1 [ratio1 [file_name2 ratio2 ...]]")
        print("spectrum: error: the following arguments are required: type file_name1")
        sys.exit()
        
    type_name = sys.argv[1]
    if type_name in ['-h', '--help']:
        print("usage: freezeLayer [-h] type file_name1 [ratio1 [file_name2 ratio2 ...]]")
        print("\n")
        print("Produce Normal-IR, Raman, or UV-Vis spectrum from Gaussian output file (.out)")
        print("For Raman spectrum, the Gaussian calculation should be done with freq=raman keyword")
        print("For UV-Vis spectrum, excited state calculations (like TD-DFT and EOM-CCSD) should be done\n")
        print("Return file_name_type.png: Spectrum image")
        print("       file_name_type.csv: Sepctrum data\n")
        print("positional arguments:")
        print("  type                   Type of spectrum:")
        print("                         'uv' for UV-Vis, 'ir' for normal IR, and 'raman' for Raman")
        print("  file_name              Gaussian output file (.out)")
        print("  [ratio]                Required to be provided for plotting a spectrum of a mixture")
        print("                         The number of the components is not restricted")
        print("                         However, the summation of the ratios must be 1")
        print("                         Ex) spectrum uv isomer1.out 0.6 isomer2.out 0.2 isomer3.out 0.2\n")
        print("options:")
        print("  -h, --help             show this help message and exit")
        sys.exit()

    # Read provided arguments
    if type_name not in ('uv', 'ir', 'raman'):
        raise TypeError('The type of the spectrum must be uv, ir, or raman')

    if argvLen == 3:
        file_name = sys.argv[2]
        name, input_format = file_name.rsplit(".",1)[0], file_name.rsplit(".",1)[-1]
        file_name = [file_name]
        name = [name]
        
        if input_format != 'out':
            raise TypeError('The input file format must be .out')
        
    elif argvLen > 3 and argvLen%2 == 1:
        raise ValueError('Ratio of one or more species is not provided')
        
    elif argvLen > 3 and argvLen%2 == 0:
        file_name = []
        ratio = []
        
        for i in range(argvLen-2):
            if i % 2 == 0:
                file_name.append(sys.argv[i+2])
                ratio.append(float(sys.argv[i+3]))
        
        name = []    
        for fn in file_name:
            name.append(fn.rsplit(".",1)[0])
            input_format = fn.rsplit(".",1)[-1]
            if input_format != 'out':
                raise TypeError('The input file format must be .out')
                
        if sum(ratio) > 1.0001 or sum(ratio) < 0.9999:
            raise ValueError('Sum of the ratios of the species should be 1')
            
    else: raise ValueError('Not enough arguments are provided')
     

    # Plot spectra
    # UV-Vis
    if type_name == 'uv':
        spectrum_save = '_'.join(name) + '_uv_vis.png'
        data_save = '_'.join(name) + '_uv_vis.csv'
        wavelengths = []
        strengths = []
        
        blue = 0.8
        alpha = 0.6

        f = plt.figure()
        f.set_figwidth(12)
        f.set_figheight(8)
        
        xrange = UV_wl_range
        curve_per_file = []
        
        for fn in file_name:
            wavelengths, strengths = uv_vis(fn)
            
            curve_per_x = []
            for idx, wv in enumerate(wavelengths):
                wv_array = np.ones(len(xrange)) * wv
                st_array = np.ones(len(xrange)) * strengths[idx]
                
                curve = np.array(list(map(uvGauss,xrange,st_array,wv_array)))
                curve_per_x.append(curve)
            
            final_curve = np.sum(curve_per_x, axis=0)
            curve_per_file.append(final_curve)
            
            if argvLen == 3:
                plt.plot(xrange, final_curve, 'b')
            elif argvLen > 3:
                red = np.random.uniform(low=0, high=1)
                green = np.random.uniform(low=0, high=1)
                plt.plot(xrange, final_curve, color = (red, green, blue), alpha=alpha)
        
        if argvLen == 3:
            f.legend(name, loc='upper right', bbox_to_anchor=(0.9, 0.88))
            
        elif argvLen > 3:
            summed_curve = np.dot(np.array(curve_per_file).T, np.array([ratio]).T)
            plt.plot(xrange, summed_curve, 'b')
            f.legend(name+['Mixture'], loc='upper right', bbox_to_anchor=(0.9, 0.88))
            
        plt.xlabel('Wavelength, $\lambda$  $(nm)$', fontsize=20)
        plt.ylabel('Absorbance $(L/mol/cm)$', fontsize=20)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.savefig(spectrum_save)
        
        if argvLen == 3:
            uv_vis_df = pd.DataFrame(np.hstack((np.array([xrange]).T,np.array(curve_per_file).T)),\
                        columns = ['wavelength']+name)
        elif argvLen > 3:
            uv_vis_df = pd.DataFrame(np.hstack((np.array([xrange]).T,np.array(curve_per_file).T,\
                        summed_curve)), columns = ['wavelength']+name+['Mixture'])
        uv_vis_df.to_csv(data_save, index=False)

    # IR: both normal IR and raman
    if type_name in ['ir', 'raman']:
        if type_name == 'ir':
            spectrum_save = '_'.join(name) + '_ir.png'
            data_save = '_'.join(name) + '_ir.csv'
            
        elif type_name == 'raman':
            spectrum_save = '_'.join(name) + '_raman.png'
            data_save = '_'.join(name) + '_raman.csv'

        frequencies = []
        intensities = []
        
        blue = 0.8
        alpha = 0.6

        f = plt.figure()
        f.set_figwidth(12)
        f.set_figheight(8)
        
        xrange = IR_wv_range
        curve_per_file = []
        
        for fn in file_name:
            if type_name == 'ir':
                frequencies, intensities = ir(fn)
                frequencies *= IR_wn_factor # Correction factor
                intensities *= 1e5/2.24e3
            elif type_name == 'raman':
                frequencies, intensities = raman(fn)
                
            curve_per_x = []
            for idx, mu in enumerate(frequencies):
                mu_array = np.ones(len(xrange))*mu
                curve = np.array(list(map(cauchy,xrange,mu_array))) * intensities[idx]
                curve_per_x.append(curve)
            
            final_curve = np.sum(curve_per_x, axis=0)
            curve_per_file.append(final_curve)
            
            if argvLen == 3:
                plt.plot(xrange, final_curve, 'b')
            elif argvLen > 3:
                red = np.random.uniform(low=0, high=1)
                green = np.random.uniform(low=0, high=1)
                plt.plot(xrange, final_curve, color = (red, green, blue), alpha=alpha)
        
        if argvLen == 3:
            f.legend(name, loc='upper right', bbox_to_anchor=(0.9, 0.88))
            
        elif argvLen > 3:
            summed_curve = np.dot(np.array(curve_per_file).T, np.array([ratio]).T)
            plt.plot(xrange, summed_curve, 'b')
            f.legend(name+['Mixture'], loc='upper right', bbox_to_anchor=(0.9, 0.88))
            
        plt.xlabel('Wavenumber, $\\nu$  $(cm^{-1})$', fontsize=20)
        plt.ylabel('Absorbance (a.u.)', fontsize=20)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.savefig(spectrum_save)
        
        if argvLen == 3:
            ir_df = pd.DataFrame(np.hstack((np.array([xrange]).T,np.array(curve_per_file).T)),\
                        columns = ['wavenumber']+name)
        elif argvLen > 3:
            ir_df = pd.DataFrame(np.hstack((np.array([xrange]).T,np.array(curve_per_file).T,\
                        summed_curve)), columns = ['wavenumber']+name+['Mixture'])
        ir_df.to_csv(data_save, index=False)

        
