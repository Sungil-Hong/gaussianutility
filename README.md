# Gaussian Utility

[![License Badge](https://img.shields.io/github/license/Sungil-Hong/gaussianutility)](https://github.com/loevlie/GPT4Readability/blob/main/LICENSE)
[![Issues Badge](https://img.shields.io/github/issues/Sungil-Hong/gaussianutility)](https://github.com/loevlie/GPT4Readability/issues)
[![Pull Requests Badge](https://img.shields.io/github/issues-pr/Sungil-Hong/gaussianutility)](https://github.com/loevlie/GPT4Readability/pulls)
[![Contributors Badge](https://img.shields.io/github/contributors/Sungil-Hong/gaussianutility)](https://github.com/loevlie/GPT4Readability/graphs/contributors)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dwyl/esta/issues)

![Logo](GaussianUtility_Logo.png)

The Gaussian Utility is a Python package that provides useful utilities for working with Gaussian 09/16 software. 

Gaussian is a computer program used by chemists, chemical engineers, biochemists, physicists, and other scientists for performing quantum chemistry calculations.

The Gaussian Utility package offers several functionalities to enhance the usage of Gaussian software.
Especially, this package can address ONIOM-type Gaussian input/output files, which is not widely available in ASE.

Most of the scripts that work on a Gaussian output file are written based Density Functional Theory calculations.
Hence, the compatibility with different types of calculations, e.g., coupled cluster, Hartree-Fock, semi-empirical, and Møller–Plesset perturbation methodes has not been confirmed yet.

The package includes the following scripts:

1. Convert input structure file type:
   - `com2vasp`: Converts Gaussian input file (.com) to VASP input file (.vasp).
   - `com2xyz`: Converts Gaussian input file (.com) to XYZ coordinate file (.xyz).
   - `vasp2com`: Converts VASP input file (.vasp) to Gaussian input file (.com).
   - `xyz2com`: Converts XYZ coordinate file (.xyz) to Gaussian input file (.com).
   
2. Extract structure from output file:
   - `out2com`: Generate Gaussian input file (.com) from Gaussian output file (.out).
   - `out2xyz`: Generate XYZ coordinate file (.xyz) from Gaussian output file (.out).

3. Spectrum utilities:
   - `spectrum`: Produces Normal-IR, Raman, or UV-Vis spectra from Gaussian output file (.out).

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
- matplotlib>=3.5.1
- mendeleev>=0.9.0

These dependencies will be automatically installed when you install the package using `setuptools`.

## Usage

The Gaussian Utility package provides command-line utilities that can be accessed by running the corresponding command.
Here are some examples of how to use the utilities:

1. Convert a structure file format:
   ```
   com2vasp structure.com
   com2xyz structure.com
   vasp2com structure.vasp
   xyz2com structure.xyz
   ```

2. Extract optimized geometry from a Gaussian output file:
   ```
   out2com output.out
   out2xyz output.out
   ```

3. Produce a spectrum from a Gaussian output file:
   ```
   spectrum type file_name1 [file_name2 ... [-r ratio1 ratio2 ...]]
   ```

4. Freeze a specific ONIOM layer of atoms in a Gaussian input file:
   ```
   freezeLayer [-i [layer index]] file_name
   ```

5. Calculate the Gibbs free energy at a given temperature:
   ```
   gibbsTemp file_name temp1 [temp2 [step_number]]
   ```

6. Print the energy of a Gaussian output file:
   ```
   printE file_name
   ```

7. Sort the atoms in a Gaussian input file based on their atomic numbers:
    ```
    sortInput [-s [sort index]] [-o [order index]] file_name
    ```
## Test
In the test folder, you'll find files for testing each command with the naming convention "command_structure.extension". For example, to test the "com2vasp" command, follow these steps:
1. Navigate to the test directory:

```
cd test
```

2. Run the desired command with the corresponding test file. For example, to test the "com2vasp" command with the file "com2vasp_ZSM5.com", use the following command:
```
com2vasp com2vasp_ZSM5.com
```

3. Verify that the command executes as expected and produces the correct output.

## Authors

The Gaussian Utility package was developed by Sungil Hong. For any inquiries or issues, you can contact Sungil Hong via email at s.hong@pitt.edu.

## License

The Gaussian Utility package is distributed under the MIT License. Please refer to the [LICENSE](LICENSE) file for more information.

## Additional Information

For more information about the Gaussian Utility package, you can visit the [GitHub repository](https://github.com/Sungil-Hong/gaussianutility). The repository contains the source code, documentation, and examples of usage for the package.

## Acknowledgment

The Gaussian Utility package has developed over multiple computational chemistry research projects conducted in the Computer-Aided Nano and Energy Lab (CANELa), led by Prof. Mpourmpakis at the University of Pittsburgh.
The projects have been financially supported by U.S. Department of Energy Nuclear Energy University Program (DOE-NEUP) (18-15496), and National Science Foundation (NSF) (1920623),
and computationally supported by the Center for Research Computing (CRC) at the University of Pittsburgh for performing DFT calculations.

