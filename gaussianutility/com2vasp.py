#!/usr/bin/env python3

import pandas as pd
import numpy as np
from periodictable import elements
import argparse
import re
from argparse import RawTextHelpFormatter
from gaussianutility.utilities import readinput

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Convert Gaussian input file to .vasp file\n\n"
                     "Return file_name.vasp: VASP geometry file (POSCAR)",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument("file_name", help='Gaussian input file (.com or .gjf)')
    parser.add_argument('-c', '--ctype', \
    help='Coordination system type; d for direct and c for cartesian', default='d')
    
    return parser.parse_args()

def com_2_vasp(file_name, ctype):
    if ctype not in ('d', 'c'):
        raise ValueError('The coordination type must be d for direct or c for cartesian')

    # Read geometry from a Gaussian input file
    _, _, _, df_geom, _ = readinput(file_name)

    if 'Index' in df_geom.columns:
        indices = df_geom['Index']
        df_geom = df_geom[['Atom', 'x', 'y', 'z']]
    else:
        df_geom = df_geom[['Atom', 'x', 'y', 'z']]

    if not 'Tv' in df_geom['Atom'].values:
        raise ValueError('The provided Gaussian input file does not contain lattice information')

    lattice = df_geom[-3:].drop(['Atom'], axis=1)
    df_geom = df_geom.drop(df_geom[df_geom['Atom'] == 'Tv'].index)
    atoms = np.array(df_geom['Atom'])

    if ctype == 'd':
        cart_coords = np.array(df_geom[['x','y','z']], dtype=float)
        lattice_inv = np.linalg.inv(np.array(lattice, dtype=float))
        direct_coords = np.dot(cart_coords, lattice_inv)

        df_geom['x'] = direct_coords[:,0]
        df_geom['y'] = direct_coords[:,1]
        df_geom['z'] = direct_coords[:,2]

    # Fix the atomic position if the given index is -1 using selective dynamics scheme
    if 'indices' in list(locals()):
        selective = [] 
        for idx in indices:
            if idx == '0':
                selective.append(['T', 'T', 'T'])
            elif idx == '-1':
                selective.append(['F', 'F', 'F'])

        selective = np.array(selective)
        df_geom['selective_x'] = selective[:,0]
        df_geom['selective_y'] = selective[:,1]
        df_geom['selective_z'] = selective[:,2]

    # Find unique atoms for defining species list and atom numbers
    uniqatomlist = []
    uniqatomcount = []

    for elem in atoms:
        if elem not in uniqatomlist:
            uniqatomlist.append(elem)
            uniqatomcount.append(np.count_nonzero(atoms == elem))

    uniqatomlist = np.array(uniqatomlist)
    uniqatomcount= np.array(uniqatomcount, dtype=str)

    # Define list of sequences of atoms and atom counts for sorting if the sequence of a single element appears more than once
    # For example, if the list of atoms comes as C, O, H, C, it should be sorted to be C, O, H
    atom_type = []
    for i in range(len(atoms)-1):
        if atoms[i] not in atom_type:
            atom_type.append(atoms[i])
            
        elif atoms[i+1] != atoms[i]:
            atom_type.append(atoms[i+1])

    atom_count = np.zeros(len(atom_type))
    atoms_buffer = atoms.copy()

    for idx in range(len(atom_type)):
        for idx2, atom in enumerate(atoms_buffer):
            if atom == atom_type[idx]:
                atom_count[idx] += 1
            else:
                atoms_buffer = atoms_buffer[int(idx2):]
                break

    atom_count = np.array(atom_count, dtype='int')

    def cumulative_sum_loop(numbers):
        cumulative = []
        current_sum = 0
        for num in numbers:
            current_sum += num
            cumulative.append(current_sum)
        return cumulative

    atom_count_cumul = cumulative_sum_loop(atom_count)
    # First and last index of each chunk of elements
    starts = [0] + atom_count_cumul[:-1] 
    ends = atom_count_cumul              

    # Segment the dataframe to each chunk of element
    segments = []
    for start, end in zip(starts, ends):
        segment = df_geom[start:end]  # Slicing the DataFrame
        segments.append(segment) # Storing each segment

    # Assign unique labels to duplicate atom types
    # For instance, C1, O, H, C2
    counts = {}
    atom_type_with_counts = []

    for atom in atom_type:
        counts[atom] = counts.get(atom, 0) + 1
        if atom_type.count(atom) == 1:
            # If the atom type is unique, keep it as is
            atom_label = atom
        else:
            # If there are duplicates, append a count to the atom type
            atom_label = f"{atom}{counts[atom]}"
        
        atom_type_with_counts.append(atom_label)

    # Pair each segment with its atom type label
    segment_pairs = list(zip(atom_type_with_counts, segments))

    # Extract base atom types in order of first occurrence
    base_atoms_in_order = []
    for atom in atom_type:
        if atom not in base_atoms_in_order:
            base_atoms_in_order.append(atom)
            
    # Create a mapping of base atom types to their sorting order
    base_atom_order = {atom: index for index, atom in enumerate(base_atoms_in_order)}

    # Define a sorting key function
    def sort_key(item):
        atom_label, _ = item
        # Extract base atom type and count
        if atom_label in base_atom_order:
            # Unique atom types without counts (e.g., 'O')
            base_atom = atom_label
            count = 0
        else:
            # Atom types with counts (e.g., 'C1', 'H2')
            match = re.match(r"([A-Za-z]+)(\d+)", atom_label)
            if match:
                base_atom = match.group(1)
                count = int(match.group(2))
            else:
                base_atom = atom_label
                count = 0
        return (base_atom_order[base_atom], count)

    # Sort the segments based on the custom order
    sorted_segment_pairs = sorted(segment_pairs, key=sort_key)

    # Extract the sorted segments
    sorted_segments = [segment for _, segment in sorted_segment_pairs]

    # Reassemble the DataFrame
    final_df = pd.concat(sorted_segments, ignore_index=True)

    # Write VASP input file
    final_df = final_df.drop(['Atom'], axis=1)
    name = ".".join(file_name.rsplit(".",1)[0:-1])
    out_file = name + ".vasp"

    with open(out_file, 'w') as output:
        output.write(name + '\n')
        output.write('1.0\n')
        output.write(lattice.to_csv(index=False, header=False, sep='\t'))
        output.write('\t'+'\t'.join(uniqatomlist)+'\n')
        output.write('\t'+'\t'.join(uniqatomcount)+'\n')

        if 'indices' in list(locals()):
            output.write('Selective dynamics\n')

        if ctype == 'd':
            output.write('Direct\n')
        elif ctype == 'c':
            output.write('Cartesian\n')
        
        output.write(final_df.to_csv(index=False, header=False, sep='\t'))
        output.write('\n')

        

def main():
    args = parse_args()
    com_2_vasp(args.file_name, args.ctype)

if __name__ == "__main__":
    main()


