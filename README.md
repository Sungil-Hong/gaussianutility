# Gaussian Utility

The Gaussian Utility is a Python package that provides useful utilities for working with Gaussian 09/16 software. 

Gaussian is a computer program used by chemists, chemical engineers, biochemists, physicists, and other scientists for performing quantum chemistry calculations.

## Purpose and Functionalities

The Gaussian Utility package offers several functionalities to enhance the usage of Gaussian software.
Especially, this package can address ONIOM-type Gaussian input/output files, which is not widely available in ASE.

Most of the scripts that work on a Gaussian output file are written based Density Functional Theory calculations.
Hence, the compatibility with different types of calculations, e.g., coupled cluster, Hartree-Fock, semi-empirical, and Møller–Plesset perturbation methodes has not been confirmed yet.

These functionalities include:

1. Convert input structure file type:
   - `com2vasp`: Converts Gaussian input files (.com) to VASP input files (.vasp).
   - `com2xyz`: Converts Gaussian input files (.com) to XYZ coordinate files (.xyz).
   - `vasp2com`: Converts VASP input files (.vasp) to Gaussian input files (.com).
   
2. Extract structure from output file:
   - `out2com`: Generate Gaussian input files (.com) from Gaussian output files (.out).
   - `out2xyz`: Generate XYZ coordinate files (.xyz) from Gaussian output files (.out).

3. Spectrum utilities:
   - `spectrum`: Produces Normal-IR, Raman, or UV-Vis spectra from Gaussian output files (.out).

4. Other utilities:
   - `gibbsTemp`: Calculates the Gibbs free energies at given temperatures from a Gaussian output file with frequency calculation.
   - `printE`: Prints the energy from a Gaussian output file.
   - `freezeLayer`: Freezes a specific ONIOM layer of atoms in a Gaussian input file with ONIOM formulation.
   - `sortInput`: Sorts the atoms in a Gaussian input file based on their atomic numbers, x, y, or z coordiate, or ONIOM layer.

## Installation

Gaussian Utility is a pip installable package:
   ```
   pip install gaussianutility
   ```

Alternatively, it can be installed from GitHub following these steps:

1. Clone the GitHub repository:
   ```
   git clone https://github.com/Sungil-Hong/gaussianutility.git
   ```

2. Navigate to the cloned repository:
   ```
   cd gaussianutility
   ```

3. Install the package using `setuptools`:
   ```
   pip install .
   ```

## Dependencies

The Gaussian Utility package has the following dependencies:

- numpy>=1.17.2
- pandas>=1.1.0
- periodictable>=1.6.1
- ase
- argparse>=1.1
- matplotlib>=3.5.3
- mendeleev>=0.12.0

These dependencies will be automatically installed when you install the package using `setuptools`.

## Usage

The Gaussian Utility package provides command-line utilities that can be accessed by running the corresponding command.
Here are some examples of how to use the utilities:

1. Convert a Gaussian input file to a VASP input file:
   ```
   com2vasp input.com
   ```

2. Convert a Gaussian input file to an XYZ coordinate file:
   ```
   com2xyz input.com
   ```

3. Convert a VASP input file to a Gaussian input file:
   ```
   vasp2com output.vasp
   ```

4. Extract optimized geometry to a Gaussian input file from a Gaussian output file:
   ```
   out2com output.out
   ```

5. Extract optimized geometry to an XYZ coordinate file:
   ```
   out2xyz output.out
   ```

6. Produce a spectrum from a Gaussian output file:
   ```
   spectrum type file_name1 [ratio1 [file_name2 ratio2 ...]]
   ```

7. Freeze a specific ONIOM layer of atoms in a Gaussian input file:
   ```
   freezeLayer [-i [layer index]] file_name
   ```

8. Calculate the Gibbs free energy at a given temperature:
   ```
   gibbsTemp file_name temp1 [temp2 [step_number]]
   ```

9. Print the energy of a Gaussian output file:
   ```
   printE file_name
   ```

10. Sort the atoms in a Gaussian input file based on their atomic numbers:
    ```
    sortInput [-s [sort index]] [-o [order index]] file_name
    ```

## Authors

The Gaussian Utility package was developed by Sungil Hong. For any inquiries or issues, you can contact Sungil Hong via email at suh33@pitt.edu.

## License

The Gaussian Utility package is distributed under the MIT License. Please refer to the [LICENSE](LICENSE) file for more information.

## Additional Information

For more information about the Gaussian Utility package, you can visit the [GitHub repository](https://github.com/Sungil-Hong/gaussianutility). The repository contains the source code, documentation, and examples of usage for the package.
